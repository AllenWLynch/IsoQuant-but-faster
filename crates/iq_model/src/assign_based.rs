//! U14: assignment-based models (`assignment_based_constructor.py:32-256`, design §7.5).
//!
//! One pass over the reads partitions them into known mono-exon, known spliced (non-FL) and
//! novel mono-exon candidates; then models are emitted in 4.1 order: known mono-exon
//! (`construct_monoexon_isoforms`), known non-FL (`construct_nonfl_isoforms`, unless
//! `fl_only`), novel mono-exon (`construct_monoexon_novel`, with `--report-novel-unspliced`).
//!
//! Order: every map 4.1 iterates here is insertion-ordered by read order (mirrored with
//! `IndexMap`); the sets it uses (`detected_known_isoforms`, `known_isoforms_in_graph_ids`,
//! `novel_monoexon`) are only tested for membership, so no hash order reaches the output.
//! The defaultdict reads that insert zero entries (`read_assignment_counts[read_id]` at
//! `:83`, the support counters) have no effect on content and are not mirrored.
//!
//! Differences from the 4.1 code shape (identical results):
//! - mono-exon coverage `[0]*t_len` bitmap (`:65-74`) → union of the reads' clipped spans;
//! - `cluster_monoexons` (`:197-207`) repeated `max` over the dict → one sort by
//!   `(count desc, insertion asc)` (counts never change) + a `BTreeMap` for range removal;
//! - `is_internal_monoexonic_read` → binary search over terminal exons sorted by the
//!   coordinate compared against `apa_delta`;
//! - D-M1: no ids are drawn (4.1 draws two per surviving cluster, even for rejected models).

use crate::ctx::{Ctx, DetectedKnown};
use crate::fixture_util::next_novel_gene;
use iq_types::indexmap::IndexMap;
use iq_types::model::vertex_code;
use iq_types::{
    BlockStore, EventType, IsoRef, IsoformIndex, Iv, ModelKind, ModelStore, Pos, ReadIdx, ReadView, Strand,
    TranscriptModel, TxIdx,
};
use std::collections::{BTreeMap, HashMap};

/// `construct_assignment_based_isoforms(read_assignment_storage)` over all block reads.
pub fn construct_assignment_based(ms: &mut ModelStore, detected_known: &mut DetectedKnown, ctx: &Ctx<'_>) {
    construct_assignment_based_reads(ms, detected_known, ctx, 0..ctx.store.len() as ReadIdx);
}

#[derive(Default)]
struct Spliced {
    reads: Vec<ReadIdx>,
    left: u32,
    right: u32,
}

struct Mono {
    reads: Vec<ReadIdx>,
    /// Clipped read spans `[s, e)` in transcript-relative coordinates (non-empty only).
    spans: Vec<(i64, i64)>,
    t_len: i64,
    polya: u32,
}

#[inline]
fn has_event(r: &ReadView<'_>, ty: EventType) -> bool {
    r.first_match().is_some_and(|m| m.events().iter().any(|e| e.ty == ty))
}

/// `construct_assignment_based_isoforms` over the reads `reads` (the 4.1 argument).
pub fn construct_assignment_based_reads(
    ms: &mut ModelStore,
    detected_known: &mut DetectedKnown,
    ctx: &Ctx<'_>,
    reads: impl IntoIterator<Item = ReadIdx>,
) {
    let p = ctx.params;
    let apa_delta = p.matching.apa_delta as i64;
    let mapq_cutoff = p.filter.simple_alignments_mapq_cutoff;
    let index = ctx.index;
    let iso_of: HashMap<TxIdx, usize> = index
        .iso
        .iter()
        .enumerate()
        .filter_map(|(i, r)| match r.key {
            IsoRef::Ref(t) => Some((t, i)),
            IsoRef::Model(_) => None,
        })
        .collect();

    let mut spliced: IndexMap<TxIdx, Spliced> = IndexMap::new();
    let mut mono: IndexMap<TxIdx, Mono> = IndexMap::new();
    let mut novel_mono: Vec<ReadIdx> = Vec::new();

    for ri in reads {
        let r = ctx.store.view(ri as usize);
        let n_exons = r.n_corrected_exons();
        if n_exons <= 2 && (r.multimapper() || r.mapq() < mapq_cutoff) {
            continue;
        }
        // `if not read_assignment` (:46) is dead: ReadAssignment has no __bool__/__len__.
        let atype = r.atype();
        if n_exons == 1 && r.polya_found() && !r.multimapper() && (atype.is_unassigned() || atype.is_inconsistent()) {
            novel_mono.push(ri);
        }
        if !atype.is_unique() {
            continue;
        }
        // A unique assignment always names an annotated isoform of this block's index; 4.1
        // would crash (KeyError) otherwise, so such reads are skipped.
        let Some(IsoRef::Ref(tx)) = r.first_match().and_then(|m| m.tx()) else { continue };
        if detected_known.contains(&tx) {
            continue;
        }
        let Some(&iso) = iso_of.get(&tx) else { continue };
        let exons = index.iso_exons(iso);
        let strand = index.iso[iso].strand;
        if has_event(&r, EventType::MonoExonMatch) {
            debug_assert!(index.iso_introns(iso).is_empty());
            let ts = exons[0].s as i64;
            let t_len = exons[exons.len() - 1].e as i64 - ts + 1;
            let e = mono.entry(tx).or_insert_with(|| Mono { reads: Vec::new(), spans: Vec::new(), t_len, polya: 0 });
            e.reads.push(ri);
            let s = (r.start() as i64 - ts).max(0);
            let en = (r.end() as i64 - ts + 1).min(t_len);
            if s < en {
                e.spans.push((s, en));
            }
            let polya_ev = if strand == Strand::Plus { EventType::CorrectPolyaSiteRight } else { EventType::CorrectPolyaSiteLeft };
            if has_event(&r, polya_ev) {
                e.polya += 1;
            }
        } else if exons.len() > 1 {
            let e = spliced.entry(tx).or_default();
            e.reads.push(ri);
            if p.requires_polya_for_construction && strand == Strand::Minus {
                if has_event(&r, EventType::CorrectPolyaSiteLeft) {
                    e.left += 1;
                }
            } else if (exons[0].s as i64 - r.start() as i64).abs() <= apa_delta {
                e.left += 1;
            }
            if p.requires_polya_for_construction && strand == Strand::Plus {
                if has_event(&r, EventType::CorrectPolyaSiteRight) {
                    e.right += 1;
                }
            } else if (exons[exons.len() - 1].e as i64 - r.end() as i64).abs() <= apa_delta {
                e.right += 1;
            }
        }
    }

    // construct_monoexon_isoforms (:209-228)
    for (&tx, m) in &mut mono {
        let count = m.reads.len() as u32;
        m.spans.sort_unstable();
        let mut covered = 0i64;
        let mut cur: Option<(i64, i64)> = None;
        for &(s, e) in &m.spans {
            match cur {
                Some((cs, ce)) if s <= ce => cur = Some((cs, ce.max(e))),
                _ => {
                    if let Some((cs, ce)) = cur {
                        covered += ce - cs;
                    }
                    cur = Some((s, e));
                }
            }
        }
        if let Some((cs, ce)) = cur {
            covered += ce - cs;
        }
        let coverage = covered as f64 / m.t_len as f64;
        if count < p.model.min_known_count
            || coverage < p.min_mono_exon_coverage
            || (p.model.require_monoexonic_polya && m.polya == 0)
        {
            continue;
        }
        if detected_known.contains(&tx) {
            continue;
        }
        let key = add_known_model(ms, index, iso_of[&tx], tx);
        detected_known.insert(tx);
        for &ri in &m.reads {
            ms.save_assigned_read(ri, key);
        }
    }

    // construct_nonfl_isoforms (:230-256)
    if !p.model.fl_only {
        let in_graph: std::collections::HashSet<TxIdx> = ctx.paths.known_isoforms_in_graph.values().copied().collect();
        for (&tx, s) in &spliced {
            if detected_known.contains(&tx) || !in_graph.contains(&tx) {
                continue;
            }
            if (s.reads.len() as u32) < p.model.min_known_count || s.left < 1 || s.right < 1 {
                continue;
            }
            let key = add_known_model(ms, index, iso_of[&tx], tx);
            detected_known.insert(tx);
            ms.nonfl_known_reads.insert(tx, s.reads.clone());
            for &ri in &s.reads {
                ms.save_assigned_read(ri, key);
            }
        }
    }

    if p.report_novel_unspliced {
        construct_monoexon_novel(ms, ctx, &novel_mono);
    }
}

/// `TranscriptModel.from_reference_transcript` + append (`:221-222`, `:248-249`).
fn add_known_model(ms: &mut ModelStore, index: &IsoformIndex, iso: usize, tx: TxIdx) -> u32 {
    let key = ms.new_key();
    let rec = &index.iso[iso];
    ms.models.push(TranscriptModel {
        key,
        strand: rec.strand,
        gene: rec.gene,
        kind: ModelKind::Known(tx),
        exons: index.iso_exons(iso).to_vec(),
        intron_path: Box::new([]),
        polya_confirmed: false,
        tss_confirmed: false,
        source_known: None,
    });
    key
}

/// Terminal exons (`collect_terminal_exons_from_graph`, `:108-121`): `(polya, polyt)` as
/// `(intron.end, polyA pos)` and `(polyT pos, intron.start)`, in graph order.
fn terminal_exons(g: &iq_types::IntronGraph) -> (Vec<Iv>, Vec<Iv>) {
    let mut polya = Vec::new();
    let mut polyt = Vec::new();
    for (k, vs) in &g.outgoing {
        for v in vs {
            if v.0 == vertex_code::POLYA {
                polya.push(Iv::new(k.1 as Pos, v.1 as Pos));
            }
        }
    }
    for (k, vs) in &g.incoming {
        for v in vs {
            if v.0 == vertex_code::POLYT {
                polyt.push(Iv::new(v.1 as Pos, k.0 as Pos));
            }
        }
    }
    (polya, polyt)
}

/// `is_internal_monoexonic_read` over terminal exons pre-sorted by the coordinate compared
/// with `apa_delta` (`e` for forward, `s` otherwise).
fn is_internal(sorted: &[Iv], start: i64, end: i64, apa_delta: i64, delta: i64, forward: bool) -> bool {
    if forward {
        // |e.e - end| <= apa_delta && start >= e.s - delta
        let lo = sorted.partition_point(|x| (x.e as i64) < end - apa_delta);
        sorted[lo..].iter().take_while(|x| x.e as i64 <= end + apa_delta).any(|x| start >= x.s as i64 - delta)
    } else {
        // |e.s - start| <= apa_delta && end <= e.e + delta
        let lo = sorted.partition_point(|x| (x.s as i64) < start - apa_delta);
        sorted[lo..].iter().take_while(|x| x.s as i64 <= start + apa_delta).any(|x| end <= x.e as i64 + delta)
    }
}

/// `cluster_monoexons` (`:197-207`): greedy, largest group first (ties: first inserted),
/// absorbing every group within `±apa_delta` (ascending position). Returns
/// `(top position, reads)` in cluster creation order.
fn cluster_monoexons(groups: IndexMap<Pos, Vec<ReadIdx>>, apa_delta: i64) -> Vec<(Pos, Vec<ReadIdx>)> {
    let mut order: Vec<usize> = (0..groups.len()).collect();
    order.sort_by_key(|&i| (std::cmp::Reverse(groups[i].len()), i));
    let mut remaining: BTreeMap<i64, usize> = groups.keys().enumerate().map(|(i, &p)| (p as i64, i)).collect();
    let mut groups: Vec<(Pos, Vec<ReadIdx>)> = groups.into_iter().collect();
    let mut out = Vec::new();
    for i in order {
        let top = groups[i].0 as i64;
        if !remaining.contains_key(&top) {
            continue;
        }
        let hit: Vec<(i64, usize)> = remaining.range(top - apa_delta..=top + apa_delta).map(|(&p, &j)| (p, j)).collect();
        let mut reads = Vec::new();
        for (p, j) in hit {
            remaining.remove(&p);
            reads.append(&mut groups[j].1);
        }
        out.push((top as Pos, reads));
    }
    out
}

/// `construct_monoexon_novel` (`:137-154`) + `generate_monoexon_from_clustered` (`:156-195`).
fn construct_monoexon_novel(ms: &mut ModelStore, ctx: &Ctx<'_>, novel: &[ReadIdx]) {
    let p = ctx.params;
    let apa_delta = p.matching.apa_delta as i64;
    let delta = p.matching.delta as i64;
    let (mut polya_ex, mut polyt_ex) = terminal_exons(ctx.graph);
    polya_ex.sort_unstable_by_key(|x| x.e);
    polyt_ex.sort_unstable_by_key(|x| x.s);
    let store: &BlockStore = ctx.store;
    let mut polya_reads: IndexMap<Pos, Vec<ReadIdx>> = IndexMap::new();
    let mut polyt_reads: IndexMap<Pos, Vec<ReadIdx>> = IndexMap::new();
    for &ri in novel {
        let r = store.view(ri as usize);
        let pa = r.polya();
        let (s, e) = (r.start() as i64, r.end() as i64);
        if pa.ext_a != -1 && !is_internal(&polya_ex, s, e, apa_delta, delta, true) {
            polya_reads.entry(pa.ext_a).or_default().push(ri);
        }
        if pa.ext_t != -1 && !is_internal(&polyt_ex, s, e, apa_delta, delta, false) {
            polyt_reads.entry(pa.ext_t).or_default().push(ri);
        }
    }
    for (groups, forward) in [(polya_reads, true), (polyt_reads, false)] {
        for (three, reads) in cluster_monoexons(groups, apa_delta) {
            if (reads.len() as u32) < p.model.min_novel_count {
                continue;
            }
            let coords = if forward {
                let five = reads.iter().map(|&r| store.view(r as usize).start()).min().unwrap();
                Iv::new(five, three)
            } else {
                let five = reads.iter().map(|&r| store.view(r as usize).end()).max().unwrap();
                Iv::new(three, five)
            };
            // intersection_len(exon, coords) > interval_len(coords) / 2
            let len2 = coords.e as i64 - coords.s as i64 + 1;
            let overlaps = |x: &Iv| {
                let inter = (x.e.min(coords.e) as i64 - x.s.max(coords.s) as i64 + 1).max(0);
                2 * inter > len2
            };
            if ms.models.iter().any(|m| m.exons.iter().any(overlaps)) {
                continue;
            }
            let key = ms.new_key();
            let gene = next_novel_gene(ms);
            ms.models.push(TranscriptModel {
                key,
                strand: if forward { Strand::Plus } else { Strand::Minus },
                gene,
                kind: ModelKind::Nnic,
                exons: vec![coords],
                intron_path: Box::new([]),
                polya_confirmed: false,
                tss_confirmed: false,
                source_known: None,
            });
            for &ri in &reads {
                ms.save_assigned_read(ri, key);
            }
        }
    }
}

pub(crate) fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    reg.add::<fixture::U14Construct>();
}

pub(crate) mod fixture {
    use super::*;
    use crate::ctx::ReadAssigner;
    use crate::fixture_util::{compare_stores, StageInput};
    use iq_types::fixture::schema::{StageStoreIn, StoreOut};
    use iq_types::fixture::{FixResult, FixtureCtx, UnitFixture, Verdict};
    use iq_types::{Assignment, CombinedProfile, Predictions};

    /// U14 never assigns.
    struct NoAssigner;
    impl ReadAssigner for NoAssigner {
        fn assign(&self, _: &IsoformIndex, _: &CombinedProfile, _: bool, _: u64) -> Option<Assignment> {
            unreachable!("U14 does not call the assigner")
        }
    }

    pub struct U14Construct;
    impl UnitFixture for U14Construct {
        const UNIT: &'static str = "U14";
        const FN: &'static str = "construct_assignment_based_isoforms";
        type In = StageStoreIn;
        type Out = StoreOut;
        fn run(input: &StageStoreIn, ctx: &FixtureCtx) -> FixResult<StoreOut> {
            let si = StageInput::load(input, ctx)?;
            let mut ms = si.store.clone();
            let mut dk = si.detected_known.clone();
            let preds = Predictions::default();
            let c = Ctx {
                params: &si.params,
                index: &si.index,
                store: &si.block.store,
                graph: &si.graph,
                paths: &si.paths,
                predictions: &preds,
                assigner: &NoAssigner,
            };
            match &input.reads {
                Some(r) => construct_assignment_based_reads(&mut ms, &mut dk, &c, r.iter().copied()),
                None => construct_assignment_based(&mut ms, &mut dk, &c),
            }
            Ok(StoreOut { store: si.store_out(&ms, &dk) })
        }
        fn compare(got: &StoreOut, want: &StoreOut, ctx: &FixtureCtx) -> Verdict {
            compare_stores(&got.store, &want.store, ctx)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use iq_types::{
        Assignment, AssignmentType, GeneKey, IntronGraph, IsoRec, IsoformMatch, MatchClass, MatchEvent, Params, PathStorage,
        PolyAInfo, Predictions, ReadRecord, Vertex, NO_POLYA,
    };

    struct NoAssigner;
    impl crate::ctx::ReadAssigner for NoAssigner {
        fn assign(&self, _: &IsoformIndex, _: &iq_types::CombinedProfile, _: bool, _: u64) -> Option<Assignment> {
            unreachable!()
        }
    }

    /// tx0: mono-exon [100,199] '+'; tx1: spliced [(1000,1100),(1200,1300)] '-'.
    fn index() -> IsoformIndex {
        let mut ix = IsoformIndex::default();
        let mut add = |t: u32, strand, exons: &[Iv]| {
            let introns: Vec<Iv> = exons.windows(2).map(|w| Iv::new(w[0].e + 1, w[1].s - 1)).collect();
            ix.iso.push(IsoRec {
                key: IsoRef::Ref(t),
                gene: GeneKey::Ref(0),
                strand,
                span: Iv::new(exons[0].s, exons[exons.len() - 1].e),
                exon_off: ix.iso_exons.len() as u32,
                n_exon: exons.len() as u32,
                intron_off: ix.iso_introns.len() as u32,
                n_intron: introns.len() as u32,
            });
            ix.iso_exons.extend_from_slice(exons);
            ix.iso_introns.extend(introns);
        };
        add(0, Strand::Plus, &[Iv::new(100, 199)]);
        add(1, Strand::Minus, &[Iv::new(1000, 1100), Iv::new(1200, 1300)]);
        ix
    }

    fn read(idx: u32, exons: &[Iv], atype: AssignmentType, tx: Option<u32>, events: &[EventType], polya: PolyAInfo) -> ReadRecord {
        let matches = match tx {
            Some(t) => [IsoformMatch::new(
                MatchClass::FullSpliceMatch,
                Some(GeneKey::Ref(0)),
                Some(IsoRef::Ref(t)),
                events.iter().map(|&e| MatchEvent::new(e)),
                Strand::Plus,
                0.0,
            )]
            .into_iter()
            .collect(),
            None => Default::default(),
        };
        ReadRecord {
            idx,
            name: format!("r{idx}"),
            exons: exons.to_vec(),
            corrected_exons: exons.to_vec(),
            boundary_deletions: vec![],
            polya,
            polya_found: polya.ext_a != -1 || polya.ext_t != -1,
            strand: Strand::Plus,
            mapped_strand: Strand::Plus,
            mapq: 60,
            multimapper: false,
            group: 0,
            assignment: Assignment::new(atype, matches),
            genomic_region: Iv::new(1, 100000),
        }
    }

    fn run(recs: &[ReadRecord], params: &Params, graph: &IntronGraph, paths: &PathStorage, dk: &mut DetectedKnown, ms: &mut ModelStore) {
        let ix = index();
        let store = BlockStore::from_records(recs);
        let preds = Predictions::default();
        let c = Ctx { params, index: &ix, store: &store, graph, paths, predictions: &preds, assigner: &NoAssigner };
        construct_assignment_based(ms, dk, &c);
    }

    fn params() -> Params {
        let mut p = Params::ont();
        p.report_novel_unspliced = true;
        p
    }

    #[test]
    fn known_mono_exon_coverage_gate() {
        use EventType::*;
        let p = params();
        let u = AssignmentType::Unique;
        // covers [100,150] and [140,180] -> 81/100 >= 0.75
        let recs = [
            read(0, &[Iv::new(90, 150)], u, Some(0), &[MonoExonMatch], NO_POLYA),
            read(1, &[Iv::new(140, 180)], u, Some(0), &[MonoExonMatch, CorrectPolyaSiteRight], NO_POLYA),
        ];
        let (mut dk, mut ms) = (DetectedKnown::default(), ModelStore::new(2));
        run(&recs, &p, &IntronGraph::default(), &PathStorage::default(), &mut dk, &mut ms);
        assert_eq!(ms.models.len(), 1);
        assert_eq!(ms.models[0].kind, ModelKind::Known(0));
        assert_eq!(ms.transcript_read_ids[&0], vec![0, 1]);
        assert!(dk.contains(&0));
        assert!(ms.nonfl_known_reads.is_empty());

        // [100,150] alone: 51/100 < 0.75
        let (mut dk, mut ms) = (DetectedKnown::default(), ModelStore::new(1));
        run(&recs[..1], &p, &IntronGraph::default(), &PathStorage::default(), &mut dk, &mut ms);
        assert!(ms.models.is_empty());

        // already detected: skipped at read level
        let (mut dk, mut ms) = (DetectedKnown::from_iter([0]), ModelStore::new(2));
        run(&recs, &p, &IntronGraph::default(), &PathStorage::default(), &mut dk, &mut ms);
        assert!(ms.models.is_empty());

        // require_monoexonic_polya: '+' needs correct_polya_site_right
        let mut p2 = p.clone();
        p2.model.require_monoexonic_polya = true;
        let (mut dk, mut ms) = (DetectedKnown::default(), ModelStore::new(2));
        run(&recs, &p2, &IntronGraph::default(), &PathStorage::default(), &mut dk, &mut ms);
        assert_eq!(ms.models.len(), 1);
        let recs0 = [recs[0].clone(), read(1, &[Iv::new(140, 180)], u, Some(0), &[MonoExonMatch], NO_POLYA)];
        let (mut dk, mut ms) = (DetectedKnown::default(), ModelStore::new(2));
        run(&recs0, &p2, &IntronGraph::default(), &PathStorage::default(), &mut dk, &mut ms);
        assert!(ms.models.is_empty());
    }

    #[test]
    fn known_nonfl_needs_graph_and_both_ends() {
        let p = params();
        let u = AssignmentType::Unique;
        let mut paths = PathStorage::default();
        paths.known_isoforms_in_graph.insert(vec![Vertex(1101, 1199)].into_boxed_slice(), 1);
        let full = [Iv::new(1010, 1100), Iv::new(1200, 1290)];
        let recs = [read(0, &full, u, Some(1), &[EventType::Fsm], NO_POLYA), read(1, &full, u, Some(1), &[EventType::Fsm], NO_POLYA)];
        let (mut dk, mut ms) = (DetectedKnown::default(), ModelStore::new(2));
        run(&recs, &p, &IntronGraph::default(), &paths, &mut dk, &mut ms);
        assert_eq!(ms.models.len(), 1);
        assert_eq!(ms.nonfl_known_reads[&1], vec![0, 1]);
        assert_eq!(ms.construction_assignment, vec![Some(0), Some(0)]);
        // not in graph
        let (mut dk, mut ms) = (DetectedKnown::default(), ModelStore::new(2));
        run(&recs, &p, &IntronGraph::default(), &PathStorage::default(), &mut dk, &mut ms);
        assert!(ms.models.is_empty());
        // right end too far (> apa_delta): no right support
        let short = [Iv::new(1010, 1100), Iv::new(1200, 1240)];
        let recs = [read(0, &short, u, Some(1), &[], NO_POLYA)];
        let (mut dk, mut ms) = (DetectedKnown::default(), ModelStore::new(1));
        run(&recs, &p, &IntronGraph::default(), &paths, &mut dk, &mut ms);
        assert!(ms.models.is_empty());
        // requires_polya: '-' strand needs correct_polya_site_left for left support
        let mut p2 = p.clone();
        p2.requires_polya_for_construction = true;
        let recs = [read(0, &full, u, Some(1), &[], NO_POLYA)];
        let (mut dk, mut ms) = (DetectedKnown::default(), ModelStore::new(1));
        run(&recs, &p2, &IntronGraph::default(), &paths, &mut dk, &mut ms);
        assert!(ms.models.is_empty());
        let recs = [read(0, &full, u, Some(1), &[EventType::CorrectPolyaSiteLeft], NO_POLYA)];
        let (mut dk, mut ms) = (DetectedKnown::default(), ModelStore::new(1));
        run(&recs, &p2, &IntronGraph::default(), &paths, &mut dk, &mut ms);
        assert_eq!(ms.models.len(), 1);
    }

    fn pa(a: Pos, t: Pos) -> PolyAInfo {
        PolyAInfo { ext_a: a, ext_t: t, int_a: -1, int_t: -1 }
    }

    #[test]
    fn novel_mono_exon() {
        let mut p = params();
        p.model.min_novel_count = 2;
        let n = AssignmentType::Noninformative;
        let recs = [
            read(0, &[Iv::new(5000, 5500)], n, None, &[], pa(5500, -1)),
            read(1, &[Iv::new(4900, 5510)], n, None, &[], pa(5510, -1)),
            read(2, &[Iv::new(4950, 5490)], AssignmentType::Inconsistent, None, &[], pa(5490, -1)),
            // polyT cluster of two reads at 8000/8020
            read(3, &[Iv::new(8000, 8400)], n, None, &[], pa(-1, 8000)),
            read(4, &[Iv::new(8020, 8500)], n, None, &[], pa(-1, 8020)),
            // singleton polyT: below min_novel_count
            read(5, &[Iv::new(20000, 20400)], n, None, &[], pa(-1, 20000)),
            // overlaps the first novel model by > half: rejected
            read(6, &[Iv::new(4900, 5300)], n, None, &[], pa(-1, 4900)),
            read(7, &[Iv::new(4910, 5300)], n, None, &[], pa(-1, 4910)),
        ];
        let (mut dk, mut ms) = (DetectedKnown::default(), ModelStore::new(recs.len()));
        run(&recs, &p, &IntronGraph::default(), &PathStorage::default(), &mut dk, &mut ms);
        assert_eq!(ms.models.len(), 2, "{:?}", ms.models);
        assert_eq!(ms.models[0].exons, vec![Iv::new(4900, 5500)]);
        assert_eq!(ms.models[0].strand, Strand::Plus);
        assert_eq!(ms.models[0].kind, ModelKind::Nnic);
        // cluster read order: by position ascending within top±apa_delta
        assert_eq!(ms.transcript_read_ids[&0], vec![2, 0, 1]);
        assert_eq!(ms.models[1].exons, vec![Iv::new(8000, 8500)]);
        assert_eq!(ms.models[1].strand, Strand::Minus);
        assert_eq!(ms.transcript_read_ids[&1], vec![3, 4]);
        assert_ne!(ms.models[0].gene, ms.models[1].gene);

        // off without --report-novel-unspliced
        p.report_novel_unspliced = false;
        let (mut dk, mut ms) = (DetectedKnown::default(), ModelStore::new(recs.len()));
        run(&recs, &p, &IntronGraph::default(), &PathStorage::default(), &mut dk, &mut ms);
        assert!(ms.models.is_empty());

        // internal read: a polyA terminal exon (5000 after intron (4000,4999)) ending near the read end
        p.report_novel_unspliced = true;
        let mut g = IntronGraph::default();
        g.outgoing.insert(Vertex(4000, 4800), vec![Vertex(vertex_code::POLYA, 5505)]);
        let (mut dk, mut ms) = (DetectedKnown::default(), ModelStore::new(recs.len()));
        run(&recs, &p, &g, &PathStorage::default(), &mut dk, &mut ms);
        // polyA reads are internal (drop), so the 4900 polyT cluster no longer overlaps anything
        let ex: Vec<_> = ms.models.iter().map(|m| (m.strand, m.exons[0])).collect();
        assert_eq!(ex, vec![(Strand::Minus, Iv::new(8000, 8500)), (Strand::Minus, Iv::new(4900, 5300))]);
    }

    /// Literal port of 4.1 `cluster_monoexons` (repeated max over the dict).
    fn cluster_naive(mut g: IndexMap<Pos, Vec<ReadIdx>>, apa: i64) -> Vec<(Pos, Vec<ReadIdx>)> {
        let mut out: Vec<(Pos, Vec<ReadIdx>)> = Vec::new();
        while !g.is_empty() {
            let mut best = 0;
            for i in 1..g.len() {
                if g[i].len() > g[best].len() {
                    best = i;
                }
            }
            let top = *g.get_index(best).unwrap().0;
            let mut acc = Vec::new();
            for pos in (top as i64 - apa)..=(top as i64 + apa) {
                if let Some(v) = g.shift_remove(&(pos as Pos)) {
                    acc.extend(v);
                }
            }
            out.push((top, acc));
        }
        out
    }

    fn lcg(s: &mut u64) -> u64 {
        *s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        *s >> 33
    }

    #[test]
    fn cluster_matches_naive() {
        let mut s = 7u64;
        for _ in 0..300 {
            let mut g: IndexMap<Pos, Vec<ReadIdx>> = IndexMap::new();
            let n = lcg(&mut s) % 40;
            for r in 0..n {
                let pos = 1000 + (lcg(&mut s) % 400) as Pos;
                g.entry(pos).or_default().push(r as ReadIdx);
            }
            let apa = (lcg(&mut s) % 60) as i64;
            assert_eq!(cluster_monoexons(g.clone(), apa), cluster_naive(g, apa));
        }
    }

    #[test]
    fn is_internal_matches_naive() {
        let mut s = 11u64;
        for _ in 0..2000 {
            let n = lcg(&mut s) % 8;
            let mut ex: Vec<Iv> = (0..n)
                .map(|_| {
                    let a = (lcg(&mut s) % 500) as Pos;
                    Iv::new(a, a + (lcg(&mut s) % 300) as Pos)
                })
                .collect();
            let st = (lcg(&mut s) % 600) as i64;
            let en = st + (lcg(&mut s) % 300) as i64;
            let (apa, d) = (50, 6);
            let fw = ex.iter().any(|e| (e.e as i64 - en).abs() <= apa && st >= e.s as i64 - d);
            let rv = ex.iter().any(|e| (e.s as i64 - st).abs() <= apa && en <= e.e as i64 + d);
            ex.sort_unstable_by_key(|x| x.e);
            assert_eq!(is_internal(&ex, st, en, apa, d, true), fw);
            ex.sort_unstable_by_key(|x| x.s);
            assert_eq!(is_internal(&ex, st, en, apa, d, false), rv);
        }
    }

    #[test]
    fn next_novel_gene_after_max() {
        let mut ms = ModelStore::new(0);
        assert_eq!(next_novel_gene(&ms), GeneKey::Novel(0));
        let m = TranscriptModel {
            key: 0,
            strand: Strand::Plus,
            gene: GeneKey::Novel(4),
            kind: ModelKind::Nnic,
            exons: vec![Iv::new(1, 2)],
            intron_path: Box::new([]),
            polya_confirmed: false,
            tss_confirmed: false,
            source_known: None,
        };
        ms.models.push(m);
        assert_eq!(next_novel_gene(&ms), GeneKey::Novel(5));
    }
}
