//! U19: `ModelReadCounter.build_model_read_assignments`
//! (`model_construction/model_read_counter.py:44-170`).
//!
//! One final assignment per block read against the block's finished models, emitted in read
//! order through a sink (nothing is retained). Categories, checked in this order:
//! 1. **construction member** (`construction_assignment[r]` names a live model) and
//! 2. **reference-unique member of a kept known isoform**: `_model_read_assignment_single`.
//!    A known model with a unique/umd source match on it reuses that match. Otherwise the
//!    read is re-profiled against a single-model index with `delta = max(delta,
//!    graph_clustering_distance)` (`Params::final_single_delta`). A consistent result whose
//!    first match is the model is kept, anything else is forced onto the model as
//!    `unique_minor_difference` / `full_splice_match` without events.
//! 3. **discarded** (reference-unique to a dropped known isoform): no transcript; a `genic`
//!    match on the isoform's block gene with `gene_atype = unique` if that gene still has a
//!    model, else no match.
//! 4. **leftover**: assigned against `from_models(all models)` (`quick_mode = false`),
//!    or `noninformative` when the block has no models.
//!
//! Structure (performance):
//! - 4.1 builds a single-model `GeneInfo` + `LongReadAssigner` per model, lazily. iqrs
//!   routes every read first (cheap), builds the index and `ProfileBuilder` of each model
//!   that some read needs once per call (in parallel), and shares them.
//! - Reads go in chunks of [`CHUNK`]; within a chunk the per-read evaluation is pure, so it
//!   runs on rayon when the chunk has at least `par_threshold` reads, and results are
//!   emitted in read order. The fixture runner checks parallel == serial.
//!
//! Model ids: the Rust model index names each model by its key (`IsoRef::Model`). A reused
//! source match names the reference transcript (`IsoRef::Ref`); it is renamed to the key
//! of the known model (4.1 uses the same string for both). The model index's provisional
//! ids only order the matches of ambiguous reads, which every consumer treats as a set.
//!
//! `ReadAssigner`: leftovers use `ctx.assigner`. The single-model calls need a different
//! junction delta than the run's assigner, which the trait cannot carry, so
//! [`final_assignments`] builds a [`LiveAssigner`] with the widened delta for them.

use crate::ctx::{Ctx, ReadAssigner};
use crate::fl::LiveAssigner;
use crate::reassign::reference_unique_known_isoform;
use iq_assign::profile::{ProfileBuilder, Scratch};
use iq_assign::AssignParams;
use iq_types::params::MatchingParams;
use iq_types::{
    Assignment, AssignmentType, FinalAssign, GeneKey, IsoRef, IsoformIndex, IsoformMatch, MatchClass, Matches,
    ModelIdx, ModelKind, ModelStore, ReadIdx, ReadView, Strand, TxIdx,
};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

/// Minimum chunk size (reads) for rayon.
pub const PAR_THRESHOLD: usize = 256;
/// Reads per chunk (bounds the results held before emission).
pub const CHUNK: usize = 1 << 14;

/// Calls `sink(read_idx, assignment)` for every block read, in read order.
pub fn final_assignments(ms: &ModelStore, ctx: &Ctx<'_>, sink: impl FnMut(u32, FinalAssign)) {
    let reads: Vec<ReadIdx> = (0..ctx.store.len() as ReadIdx).collect();
    let full = full_index(ms, ctx);
    let single = single_assigner(ctx);
    for_each_final(ms, ctx, &reads, &full, &single, PAR_THRESHOLD, |a| compact(&a), sink);
}

/// `from_models(all models)` with provisional ids (empty index when there are no models).
pub fn full_index(ms: &ModelStore, ctx: &Ctx<'_>) -> IsoformIndex {
    if ms.models.is_empty() {
        IsoformIndex::default()
    } else {
        iq_annot::index_from_models(ctx.index.chr, &ms.models, ctx.params.matching.delta)
    }
}

/// Matching params of the single-model path (`_graph_tolerance_args`).
pub fn single_matching(ctx: &Ctx<'_>) -> MatchingParams {
    MatchingParams { delta: ctx.params.final_single_delta(), ..ctx.params.matching.clone() }
}

/// The single-model assigner (`LongReadAssigner(gi, _graph_tolerance_args, quick_mode=False)`).
pub fn single_assigner(ctx: &Ctx<'_>) -> LiveAssigner {
    LiveAssigner { params: AssignParams::new(single_matching(ctx), ctx.params.splice.clone(), false) }
}

/// Compacts an assignment for the counters: the model of every match that
/// names one, in match order.
pub fn compact(a: &Assignment) -> FinalAssign {
    FinalAssign {
        atype: a.atype,
        gene_atype: a.gene_atype,
        tx: a
            .matches
            .iter()
            .filter_map(|m| match m.tx {
                Some(IsoRef::Model(k)) => Some(k),
                _ => None,
            })
            .collect(),
    }
}

/// Category of a read (decided from the store alone).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Route {
    /// `_model_read_assignment_single` on `models[pos]`; `reuse` = the source match is reused.
    Single { pos: usize, reuse: bool },
    /// Reference-unique to a dropped known isoform; its block gene if it has a surviving model.
    Discarded(Option<GeneKey>),
    /// Leftover against the full model index.
    Full,
    /// Leftover in a block without models.
    Noninformative,
}

/// Per-call lookup tables (4.1 `model_by_id`, `model_ids`, `surviving_genes`, `gene_id_map`).
struct Router<'a> {
    ms: &'a ModelStore,
    /// model key → storage position (dict semantics: last wins).
    by_key: HashMap<ModelIdx, usize>,
    /// known transcript → storage position of its model (last wins).
    by_known: HashMap<TxIdx, usize>,
    surviving_genes: HashSet<GeneKey>,
    /// block `gene_id_map` (last wins).
    ref_gene: HashMap<TxIdx, GeneKey>,
}

impl<'a> Router<'a> {
    fn new(ms: &'a ModelStore, block_index: &IsoformIndex) -> Self {
        let mut by_key = HashMap::new();
        let mut by_known = HashMap::new();
        for (i, m) in ms.models.iter().enumerate() {
            by_key.insert(m.key, i);
            if let ModelKind::Known(t) = m.kind {
                by_known.insert(t, i);
            }
        }
        let surviving_genes = ms.models.iter().map(|m| m.gene).collect();
        let mut ref_gene = HashMap::new();
        for rec in &block_index.iso {
            if let IsoRef::Ref(t) = rec.key {
                ref_gene.insert(t, rec.gene);
            }
        }
        Router { ms, by_key, by_known, surviving_genes, ref_gene }
    }

    fn single(&self, v: &ReadView<'_>, pos: usize) -> Route {
        let reuse = match self.ms.models[pos].kind {
            ModelKind::Known(t) if v.atype().is_unique() => v.matches().any(|m| m.tx() == Some(IsoRef::Ref(t))),
            _ => false,
        };
        Route::Single { pos, reuse }
    }

    fn route(&self, r: ReadIdx, v: &ReadView<'_>) -> Route {
        if let Some(k) = self.ms.construction_assignment.get(r as usize).copied().flatten() {
            if let Some(&pos) = self.by_key.get(&k) {
                return self.single(v, pos);
            }
        }
        match reference_unique_known_isoform(v) {
            Some(IsoRef::Ref(t)) => match self.by_known.get(&t) {
                Some(&pos) => self.single(v, pos),
                None => Route::Discarded(self.ref_gene.get(&t).copied().filter(|g| self.surviving_genes.contains(g))),
            },
            // a reference-step assignment never names a model; 4.1 would find no such id
            Some(IsoRef::Model(_)) => Route::Discarded(None),
            None if self.ms.models.is_empty() => Route::Noninformative,
            None => Route::Full,
        }
    }
}

fn single_match(class: MatchClass, gene: Option<GeneKey>, tx: Option<IsoRef>) -> Matches {
    let mut m = Matches::new();
    m.push(IsoformMatch { gene, tx, tx_strand: Strand::Unknown, class, penalty: 0.0, events: Default::default() });
    m
}

/// Core of U19 over `reads` (the 4.1 argument, in order). `full` = `from_models(ms.models)`
/// (any id naming; unused without models), `single` assigns against the single-model
/// indices (built here with `single_matching` profiles). Each assignment is mapped by `map`
/// (in parallel) and handed to `emit` in read order. Parallel above `par_threshold` reads per
/// chunk (`0` = always, `usize::MAX` = never); the output does not depend on it.
#[allow(clippy::too_many_arguments)]
pub fn for_each_final<T: Send>(
    ms: &ModelStore,
    ctx: &Ctx<'_>,
    reads: &[ReadIdx],
    full: &IsoformIndex,
    single: &dyn ReadAssigner,
    par_threshold: usize,
    map: impl Fn(Assignment) -> T + Sync,
    mut emit: impl FnMut(u32, T),
) {
    let router = Router::new(ms, ctx.index);
    let store = ctx.store;

    // single-model indices needed by some read (4.1: built lazily, cached per model)
    let mut needed = vec![false; ms.models.len()];
    for &r in reads {
        if let Route::Single { pos, reuse: false } = router.route(r, &store.view(r as usize)) {
            needed[pos] = true;
        }
    }
    let chr = ctx.index.chr;
    let delta = ctx.params.matching.delta;
    let build = |(i, &need): (usize, &bool)| {
        need.then(|| iq_annot::index_from_models(chr, std::slice::from_ref(&ms.models[i]), delta))
    };
    let n_needed = needed.iter().filter(|&&b| b).count();
    let single_idx: Vec<Option<IsoformIndex>> = if n_needed >= 8 && par_threshold != usize::MAX {
        needed.par_iter().enumerate().map(build).collect()
    } else {
        needed.iter().enumerate().map(build).collect()
    };
    let single_params = single_matching(ctx);
    let single_builders: Vec<Option<ProfileBuilder<'_>>> =
        single_idx.iter().map(|ix| ix.as_ref().map(|ix| ProfileBuilder::new(ix, &single_params))).collect();
    let full_builder = (!ms.models.is_empty()).then(|| ProfileBuilder::new(full, &ctx.params.matching));

    let eval = |r: ReadIdx, scratch: &mut Scratch| -> Assignment {
        let v = store.view(r as usize);
        match router.route(r, &v) {
            Route::Single { pos, reuse } => {
                let model = &ms.models[pos];
                let key = IsoRef::Model(model.key);
                if reuse {
                    let ModelKind::Known(t) = model.kind else { unreachable!() };
                    let src = v.matches().find(|m| m.tx() == Some(IsoRef::Ref(t))).expect("routed reuse");
                    let mut m = src.to_match();
                    m.tx = Some(key);
                    let mut ms_ = Matches::new();
                    ms_.push(m);
                    return Assignment::new(v.atype(), ms_);
                }
                let (Some(ix), Some(b)) = (&single_idx[pos], &single_builders[pos]) else { unreachable!() };
                let exons = v.corrected_exons();
                let prof = b.build(&exons, v.polya(), scratch);
                if let Some(a) = single.assign(ix, &prof, false, r as u64) {
                    if a.atype.is_consistent() && a.matches.first().is_some_and(|m| m.tx == Some(key)) {
                        return a;
                    }
                }
                Assignment::new(
                    AssignmentType::UniqueMinorDifference,
                    single_match(MatchClass::FullSpliceMatch, Some(model.gene), Some(key)),
                )
            }
            Route::Discarded(Some(g)) => Assignment {
                atype: AssignmentType::Discarded,
                gene_atype: AssignmentType::Unique,
                matches: single_match(MatchClass::Genic, Some(g), None),
            },
            Route::Discarded(None) => Assignment::empty(AssignmentType::Discarded),
            Route::Full => {
                let b = full_builder.as_ref().expect("models exist");
                let exons = v.corrected_exons();
                let prof = b.build(&exons, v.polya(), scratch);
                // 4.1 assign_to_isoform never returns None here (it would crash on it)
                ctx.assigner
                    .assign(full, &prof, false, r as u64)
                    .unwrap_or_else(|| Assignment::empty(AssignmentType::Noninformative))
            }
            Route::Noninformative => Assignment::empty(AssignmentType::Noninformative),
        }
    };

    let mut scratch = Scratch::new();
    let mut out: Vec<T> = Vec::new();
    for chunk in reads.chunks(CHUNK) {
        if chunk.len() >= par_threshold {
            chunk.par_iter().map_init(Scratch::new, |s, &r| map(eval(r, s))).collect_into_vec(&mut out);
        } else {
            out.clear();
            out.extend(chunk.iter().map(|&r| map(eval(r, &mut scratch))));
        }
        for (&r, t) in chunk.iter().zip(out.drain(..)) {
            emit(r, t);
        }
    }
}

/// Registers the U19 fixture runner.
pub(crate) fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    reg.add::<fixture::U19Final>();
}

pub(crate) mod fixture {
    //! U19 runner (plumbing: [`crate::fixture_util`]). Nested assignments run live
    //! ([`LiveAssigner`], U6 exact). The full index uses the real model ids (4.1 sorts the
    //! matches of ambiguous reads by id string). The runner runs serial and forced-parallel
    //! and fails the case if they differ. Output: the full assignment per read (stronger
    //! than the `(atype, gene_atype, [(tx, gene)])` contract).
    use super::*;
    use crate::fixture_util::StageInput;
    use iq_types::fixture::schema::{AssignmentFixture, StageStoreIn, U19AssignmentsOut, U19ReadAssignment};
    use iq_types::fixture::{FixError, FixResult, FixtureCtx, UnitFixture};
    use iq_types::Predictions;

    pub struct U19Final;

    impl UnitFixture for U19Final {
        const UNIT: &'static str = "U19";
        const FN: &'static str = "build_model_read_assignments";
        // U6 deviations reach U19 through the nested assignments
        type In = StageStoreIn;
        type Out = U19AssignmentsOut;

        fn run(input: &StageStoreIn, ctx: &FixtureCtx) -> FixResult<U19AssignmentsOut> {
            let si = StageInput::load(input, ctx)?;
            let p = &si.params;
            let assigner = LiveAssigner { params: AssignParams::new(p.matching.clone(), p.splice.clone(), false) };
            let preds = Predictions::default();
            let c = Ctx {
                params: p,
                index: &si.index,
                store: &si.block.store,
                graph: &si.graph,
                paths: &si.paths,
                predictions: &preds,
                assigner: &assigner,
            };
            let reads: Vec<ReadIdx> = match &input.reads {
                Some(r) => r.clone(),
                None => (0..si.block.store.len() as ReadIdx).collect(),
            };
            let ms = &si.store;
            let names: Vec<&str> = ms.models.iter().map(|m| si.ids.models.resolve(m.key)).collect();
            let full = if ms.models.is_empty() {
                IsoformIndex::default()
            } else {
                iq_annot::index_from_models_named(si.index.chr, &ms.models, p.matching.delta, &names)
            };
            let single = single_assigner(&c);
            let run = |thr: usize| {
                let mut v = Vec::with_capacity(reads.len());
                for_each_final(ms, &c, &reads, &full, &single, thr, |a| a, |r, a| v.push((r, a)));
                v
            };
            let serial = run(usize::MAX);
            if serial != run(0) {
                return Err(FixError("U19: serial and parallel final assignments differ".into()));
            }
            // production entry (provisional ids): same compact result up to match order
            if reads.len() == si.block.store.len() && reads.iter().enumerate().all(|(i, &r)| i == r as usize) {
                let sorted = |mut f: FinalAssign| {
                    f.tx.sort_unstable();
                    f
                };
                let mut prod = Vec::with_capacity(reads.len());
                final_assignments(ms, &c, |r, fa| prod.push((r, sorted(fa))));
                let want: Vec<(u32, FinalAssign)> = serial.iter().map(|(r, a)| (*r, sorted(compact(a)))).collect();
                if prod != want {
                    return Err(FixError(
                        "U19: final_assignments (provisional ids) differs from the named-index run".into(),
                    ));
                }
            }
            let assignments = serial
                .iter()
                .map(|(r, a)| U19ReadAssignment { read_idx: *r, assignment: AssignmentFixture::from_core(a, &si.ids) })
                .collect();
            Ok(U19AssignmentsOut { assignments })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use iq_types::{
        BlockStore, CombinedProfile, IntronGraph, IsoRec, Iv, Params, PathStorage, Predictions, ReadRecord,
        TranscriptModel,
    };
    use std::sync::Mutex;

    /// Returns the scripted assignment for `read_key`, records which index sizes it saw.
    struct Scripted {
        by_read: HashMap<u64, Assignment>,
        calls: Mutex<Vec<(u64, usize)>>,
    }

    impl ReadAssigner for Scripted {
        fn assign(
            &self,
            index: &IsoformIndex,
            _p: &CombinedProfile,
            quick_mode: bool,
            read_key: u64,
        ) -> Option<Assignment> {
            assert!(!quick_mode);
            self.calls.lock().unwrap().push((read_key, index.n_iso()));
            self.by_read.get(&read_key).cloned()
        }
    }

    const EXONS: [Iv; 2] = [Iv::new(100, 200), Iv::new(300, 400)];

    fn model(key: ModelIdx, kind: ModelKind, gene: GeneKey) -> TranscriptModel {
        TranscriptModel {
            key,
            strand: Strand::Plus,
            gene,
            kind,
            exons: EXONS.to_vec(),
            intron_path: Default::default(),
            polya_confirmed: false,
            tss_confirmed: false,
            source_known: None,
        }
    }

    fn uniq(tx: IsoRef, gene: GeneKey) -> Assignment {
        Assignment::new(AssignmentType::Unique, single_match(MatchClass::FullSpliceMatch, Some(gene), Some(tx)))
    }

    fn record(idx: u32, a: Assignment) -> ReadRecord {
        ReadRecord {
            idx,
            name: format!("r{idx}"),
            exons: EXONS.to_vec(),
            corrected_exons: EXONS.to_vec(),
            boundary_deletions: vec![],
            polya: Default::default(),
            polya_found: false,
            strand: Strand::Plus,
            mapped_strand: Strand::Plus,
            mapq: 60,
            multimapper: false,
            group: 0,
            assignment: a,
            genomic_region: Iv::new(0, 1000),
        }
    }

    fn block_index(txs: &[(TxIdx, GeneKey)]) -> IsoformIndex {
        let iso = txs
            .iter()
            .map(|&(t, g)| IsoRec {
                key: IsoRef::Ref(t),
                gene: g,
                strand: Strand::Plus,
                span: Iv::new(100, 400),
                exon_off: 0,
                n_exon: 0,
                intron_off: 0,
                n_intron: 0,
            })
            .collect();
        IsoformIndex { iso, ..Default::default() }
    }

    /// Runs `for_each_final` (serial and parallel; asserts equal) with `single` and `full`
    /// scripted assigners.
    fn run(
        ms: &ModelStore,
        recs: &[ReadRecord],
        bidx: &IsoformIndex,
        single: &Scripted,
        full: &Scripted,
    ) -> Vec<Assignment> {
        let params = Params::pacbio();
        let store = BlockStore::from_records(recs);
        let (graph, paths, preds) = (IntronGraph::default(), PathStorage::default(), Predictions::default());
        let ctx = Ctx {
            params: &params,
            index: bidx,
            store: &store,
            graph: &graph,
            paths: &paths,
            predictions: &preds,
            assigner: full,
        };
        let fidx = full_index(ms, &ctx);
        let reads: Vec<ReadIdx> = (0..recs.len() as ReadIdx).collect();
        let go = |thr| {
            let mut v = Vec::new();
            for_each_final(ms, &ctx, &reads, &fidx, single, thr, |a| a, |r, a| v.push((r, a)));
            v
        };
        let s = go(usize::MAX);
        assert_eq!(s, go(0));
        assert!(s.iter().enumerate().all(|(i, (r, _))| i == *r as usize));
        s.into_iter().map(|(_, a)| a).collect()
    }

    fn scripted(by_read: impl IntoIterator<Item = (u64, Assignment)>) -> Scripted {
        Scripted { by_read: by_read.into_iter().collect(), calls: Mutex::new(vec![]) }
    }

    const G_REF: GeneKey = GeneKey::Ref(7);
    const G_NOV: GeneKey = GeneKey::Novel(0);

    #[test]
    fn categories() {
        // models: 0 = novel NIC (gene G_NOV), 1 = known tx 10 (gene G_REF)
        let mut ms = ModelStore::new(7);
        ms.models.push(model(0, ModelKind::Nic, G_NOV));
        ms.models.push(model(1, ModelKind::Known(10), G_REF));
        ms.next_key = 2;
        ms.construction_assignment[0] = Some(0); // construction, re-profiled: consistent → kept
        ms.construction_assignment[1] = Some(0); // construction, re-profiled: inconsistent → forced
        ms.construction_assignment[2] = Some(5); // dead model → falls through to known tx 10 (reuse)
        let bidx = block_index(&[(10, G_REF), (11, G_REF), (12, GeneKey::Ref(8))]);
        let recs = vec![
            record(0, Assignment::empty(AssignmentType::Inconsistent)),
            record(1, Assignment::empty(AssignmentType::Noninformative)),
            record(2, uniq(IsoRef::Ref(10), G_REF)),
            record(3, uniq(IsoRef::Ref(11), G_REF)), // dropped known, gene survives
            record(4, uniq(IsoRef::Ref(12), GeneKey::Ref(8))), // dropped known, gene gone
            record(5, Assignment::empty(AssignmentType::Inconsistent)), // leftover → full
            record(6, Assignment::empty(AssignmentType::Ambiguous)), // leftover, assigner None
        ];
        let kept = uniq(IsoRef::Model(0), G_NOV);
        let single = scripted([(0, kept.clone()), (1, Assignment::empty(AssignmentType::Inconsistent))]);
        let amb = Assignment::new(
            AssignmentType::Ambiguous,
            [IsoRef::Model(0), IsoRef::Model(1)]
                .iter()
                .map(|&t| single_match(MatchClass::FullSpliceMatch, Some(G_NOV), Some(t))[0].clone())
                .collect(),
        );
        let full = scripted([(5, amb.clone())]);
        let got = run(&ms, &recs, &bidx, &single, &full);

        assert_eq!(got[0], kept);
        assert_eq!(got[1].atype, AssignmentType::UniqueMinorDifference);
        assert_eq!(got[1].gene_atype, AssignmentType::UniqueMinorDifference); // derived as in ReadAssignment.__init__
        assert_eq!(got[1].matches[0].tx, Some(IsoRef::Model(0)));
        assert_eq!(got[1].matches[0].gene, Some(G_NOV));
        assert_eq!(got[1].matches[0].class, MatchClass::FullSpliceMatch);
        assert!(got[1].matches[0].events.is_empty());
        // reuse: tx renamed to the known model's key, gene kept
        assert_eq!(got[2], uniq(IsoRef::Model(1), G_REF));
        assert_eq!(got[3].atype, AssignmentType::Discarded);
        assert_eq!(got[3].gene_atype, AssignmentType::Unique);
        assert_eq!(
            (got[3].matches[0].tx, got[3].matches[0].gene, got[3].matches[0].class),
            (None, Some(G_REF), MatchClass::Genic)
        );
        assert_eq!(got[4], Assignment::empty(AssignmentType::Discarded));
        assert_eq!(got[5], amb);
        assert_eq!(got[6], Assignment::empty(AssignmentType::Noninformative));
        assert_eq!(compact(&got[5]).tx.as_slice(), &[0, 1]);
        assert!(compact(&got[3]).tx.is_empty());

        // single-model calls see a one-model index, full calls the whole model index
        let sc = single.calls.lock().unwrap().clone();
        assert!(sc.iter().all(|&(r, n)| n == 1 && r <= 1));
        let fc = full.calls.lock().unwrap().clone();
        assert!(fc.iter().all(|&(r, n)| n == 2 && r >= 5));
    }

    #[test]
    fn consistent_on_other_model_is_forced() {
        // single assigner returns a consistent match whose first tx is not the model
        let mut ms = ModelStore::new(1);
        ms.models.push(model(0, ModelKind::Nnic, G_NOV));
        ms.next_key = 1;
        ms.construction_assignment[0] = Some(0);
        let single = scripted([(0, uniq(IsoRef::Model(3), G_NOV))]);
        let got = run(
            &ms,
            &[record(0, Assignment::empty(AssignmentType::Noninformative))],
            &block_index(&[]),
            &single,
            &scripted([]),
        );
        assert_eq!(got[0].atype, AssignmentType::UniqueMinorDifference);
        assert_eq!(got[0].matches[0].tx, Some(IsoRef::Model(0)));
    }

    #[test]
    fn known_construction_without_source_match_reprofiles() {
        // construction-bound to known model tx 10, but the source match is on tx 11 → re-profile
        let mut ms = ModelStore::new(1);
        ms.models.push(model(0, ModelKind::Known(10), G_REF));
        ms.next_key = 1;
        ms.construction_assignment[0] = Some(0);
        let kept = Assignment::new(
            AssignmentType::UniqueMinorDifference,
            single_match(MatchClass::FullSpliceMatch, Some(G_REF), Some(IsoRef::Model(0))),
        );
        let single = scripted([(0, kept.clone())]);
        let got =
            run(&ms, &[record(0, uniq(IsoRef::Ref(11), G_REF))], &block_index(&[(11, G_REF)]), &single, &scripted([]));
        assert_eq!(got[0], kept);
        assert_eq!(single.calls.lock().unwrap().len(), 2); // serial + parallel run
    }

    #[test]
    fn no_models() {
        let ms = ModelStore::new(2);
        let recs =
            [record(0, Assignment::empty(AssignmentType::Inconsistent)), record(1, uniq(IsoRef::Ref(10), G_REF))];
        let got = run(&ms, &recs, &block_index(&[(10, G_REF)]), &scripted([]), &scripted([]));
        assert_eq!(got[0], Assignment::empty(AssignmentType::Noninformative));
        // gene has no surviving model
        assert_eq!(got[1], Assignment::empty(AssignmentType::Discarded));
    }
}
