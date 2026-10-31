//! U18: `TranscriptToGeneJoiner` (`model_construction/transcript_to_gene_joiner.py:17-125`).
//!
//! Joins novel genes into overlapping annotated genes, or into each other, by rewriting the
//! models' gene keys. Scoring is verbatim: same strand only; shared introns →
//! `1 + shared/union`; otherwise `overlap / min(len)` if ≥ 0.5; never two annotated genes.
//!
//! ## Layout (equivalent to 4.1, cheaper)
//! - Genes are kept in 4.1 `gene_to_transcripts` insertion order: annotated genes by first
//!   isoform in index order, then genes of non-known models in storage order. A pair's
//!   position in the 4.1 `scores` dict is the `(pos_a, pos_b)` (`pos_a < pos_b`) of its genes
//!   in that order, so `max(scores, key=get)` = the first maximal pair in that key order.
//! - Only positive scores are stored (`BTreeMap` by position). A zero score never wins
//!   (`<= 0 → break`). The number of live 4.1 pairs, needed only for the loop condition, is
//!   a closed formula: a pair dies exactly when one of its genes is merged away (always a
//!   non-annotated gene), so live pairs = all pairs of live genes that are not both annotated.
//! - Pairs whose genes' hulls (region ∪ intron extents) are disjoint score 0 exactly and are
//!   never evaluated.
//! - Transcript sets are not materialised: a merged gene points at its survivor, and each
//!   model's final gene is the root of its starting gene.
//!
//! ## Merge direction
//! 4.1 sorts the pair by id string; if the first is annotated it survives, else the second
//! survives. With one annotated gene the annotated one always survives. Between two novel
//! genes the survivor is the one with the larger 4.1 id string (`novel_gene_<chr>_<N>`);
//! [`join_genes_by`] takes that order as a comparator. iqrs has no gene strings at this
//! stage (D-M1), so [`join_genes`] uses the block-local key (creation order, = numeric `N`
//! order). The choice only decides which key the merged group keeps; in all fixture cases
//! the resulting partition of models into genes is identical either way.
//!
//! ## Last pair (matches IsoQuant 4.1)
//! 4.1 loops `while len(self.scores) > 1`, so a last remaining pair is never merged.

use iq_types::{GeneKey, IsoRef, IsoformIndex, Iv, ModelKind, Strand, TranscriptModel};
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap};

/// `STRONG_OVERLAP_FRACTION`.
const STRONG_OVERLAP_FRACTION: f64 = 0.5;

struct Gene {
    key: GeneKey,
    /// `gene_id in gene_info.gene_strands`.
    annotated: bool,
    strand: Strand,
    region: Iv,
    /// Has `region`/`strand` (4.1 `gene_id in self.gene_regions`).
    placed: bool,
    /// Sorted, deduplicated.
    introns: Vec<Iv>,
    /// Region ∪ intron extents; disjoint hulls ⇒ score 0.
    hull: Iv,
    /// Survivor position this gene was merged into (`u32::MAX` = alive).
    merged_into: u32,
}

impl Gene {
    fn update_hull(&mut self) {
        // `introns` is sorted, so the first has the smallest start.
        let s = self.introns.first().map_or(self.region.s, |i| i.s.min(self.region.s));
        let e = self.introns.iter().fold(self.region.e, |e, i| e.max(i.e));
        self.hull = Iv::new(s, e);
    }
}

/// Size of the intersection of two sorted, deduplicated interval lists.
fn shared_count(a: &[Iv], b: &[Iv]) -> usize {
    let (mut i, mut j, mut n) = (0, 0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            Ordering::Less => i += 1,
            Ordering::Greater => j += 1,
            Ordering::Equal => {
                n += 1;
                i += 1;
                j += 1;
            }
        }
    }
    n
}

/// `count_score` (`:56-74`).
fn score(g1: &Gene, g2: &Gene) -> f64 {
    if g1.strand != g2.strand || !g1.hull.overlaps(&g2.hull) {
        return 0.0;
    }
    let shared = shared_count(&g1.introns, &g2.introns);
    if shared > 0 {
        let union = g1.introns.len() + g2.introns.len() - shared;
        return 1.0 + shared as f64 / union.max(1) as f64;
    }
    let (r1, r2) = (g1.region, g2.region);
    let inter = (r1.e.min(r2.e) as i64 - r1.s.max(r2.s) as i64 + 1).max(0);
    if inter <= 0 {
        return 0.0;
    }
    let frac = inter as f64 / r1.len().min(r2.len()) as f64;
    if frac >= STRONG_OVERLAP_FRACTION {
        frac
    } else {
        0.0
    }
}

/// `junctions_from_blocks`: `(b[i].e + 1, b[i+1].s - 1)`.
fn push_junctions(exons: &[Iv], out: &mut Vec<Iv>) {
    out.extend(exons.windows(2).map(|w| Iv::new(w[0].e + 1, w[1].s - 1)));
}

/// `fn join_genes(&mut [TranscriptModel], &IsoformIndex)`: novel genes are ordered by their
/// block-local key (see the module docs, "Merge direction").
pub fn join_genes(models: &mut [TranscriptModel], index: &IsoformIndex) {
    join_genes_by(models, index, |a, b| a.cmp(&b));
}

/// [`join_genes`] with the order of novel gene keys (`GeneKey::Novel(k)`) that stands in for
/// 4.1's id-string order: when two novel genes merge, the greater one survives.
pub fn join_genes_by(models: &mut [TranscriptModel], index: &IsoformIndex, novel_cmp: impl Fn(u32, u32) -> Ordering) {
    // --- __init__ ---------------------------------------------------------------------------
    let annot: HashMap<GeneKey, usize> = index.genes.iter().enumerate().map(|(i, g)| (g.key, i)).collect();
    let mut pos: HashMap<GeneKey, u32> = HashMap::new();
    let mut genes: Vec<Gene> = Vec::new();
    let new_gene = |key: GeneKey| -> Gene {
        match annot.get(&key) {
            Some(&i) => {
                let r = &index.genes[i];
                Gene {
                    key,
                    annotated: true,
                    strand: r.strand,
                    region: r.region,
                    placed: true,
                    introns: Vec::new(),
                    hull: r.region,
                    merged_into: u32::MAX,
                }
            }
            None => Gene {
                key,
                annotated: false,
                strand: Strand::Unknown,
                region: Iv::default(),
                placed: false,
                introns: Vec::new(),
                hull: Iv::default(),
                merged_into: u32::MAX,
            },
        }
    };
    // Annotated transcripts (`gene_id_map` order = index order).
    let mut tx_gene: HashMap<u32, u32> = HashMap::new();
    for (i, rec) in index.iso.iter().enumerate() {
        let p = *pos.entry(rec.gene).or_insert_with(|| {
            genes.push(new_gene(rec.gene));
            (genes.len() - 1) as u32
        });
        let g = &mut genes[p as usize];
        if !g.placed {
            // Gene of an isoform missing from `gene_strands` (4.1 would fail on it later).
            g.placed = true;
            g.strand = rec.strand;
            g.region = rec.span;
        }
        g.introns.extend_from_slice(index.iso_introns(i));
        if let IsoRef::Ref(t) = rec.key {
            tx_gene.insert(t, p);
        }
    }
    // Novel models.
    for m in models.iter() {
        if matches!(m.kind, ModelKind::Known(_)) {
            continue;
        }
        let p = *pos.entry(m.gene).or_insert_with(|| {
            genes.push(new_gene(m.gene));
            (genes.len() - 1) as u32
        });
        let g = &mut genes[p as usize];
        let span = m.span();
        if !g.placed {
            g.placed = true;
            g.region = span;
            g.strand = m.strand;
        } else {
            g.region = Iv::new(g.region.s.min(span.s), g.region.e.max(span.e));
            debug_assert_eq!(g.strand, m.strand, "joiner: strand mismatch within gene");
        }
        push_junctions(&m.exons, &mut g.introns);
    }
    for g in &mut genes {
        g.introns.sort_unstable();
        g.introns.dedup();
        g.update_hull();
    }

    // --- count_scores -------------------------------------------------------------------------
    let n = genes.len();
    let n_annot = genes.iter().filter(|g| g.annotated).count() as u64;
    let mut n_novel = n as u64 - n_annot;
    let mut scores: BTreeMap<(u32, u32), f64> = BTreeMap::new();
    // Sweep over hull starts: only hull-overlapping pairs can score > 0.
    let mut by_start: Vec<u32> = (0..n as u32).collect();
    by_start.sort_unstable_by_key(|&i| genes[i as usize].hull.s);
    for (k, &i) in by_start.iter().enumerate() {
        let gi = &genes[i as usize];
        for &j in &by_start[k + 1..] {
            let gj = &genes[j as usize];
            if gj.hull.s > gi.hull.e {
                break;
            }
            if gi.annotated && gj.annotated {
                continue;
            }
            let s = score(gi, gj);
            if s > 0.0 {
                scores.insert((i.min(j), i.max(j)), s);
            }
        }
    }

    // --- join_transcripts -----------------------------------------------------------------------
    // 4.1 `len(self.scores)`: pairs of live genes, not both annotated.
    let live_pairs = |n_novel: u64| n_novel * n_novel.saturating_sub(1) / 2 + n_novel * n_annot;
    // Matches IsoQuant 4.1: `while len(self.scores) > 1` never merges the last pair.
    while live_pairs(n_novel) >= 2 {
        // `max(scores, key=get)`: first maximal entry in key order.
        let mut best: Option<((u32, u32), f64)> = None;
        for (&k, &s) in &scores {
            if best.is_none_or(|(_, b)| s > b) {
                best = Some((k, s));
            }
        }
        let Some(((a, b), _)) = best else { break };
        let (ga, gb) = (&genes[a as usize], &genes[b as usize]);
        let a_survives = match (ga.annotated, gb.annotated) {
            (true, _) => true,
            (_, true) => false,
            _ => match (ga.key, gb.key) {
                (GeneKey::Novel(x), GeneKey::Novel(y)) => novel_cmp(x, y) == Ordering::Greater,
                // Unannotated `Ref` keys (gene missing from the index): by key.
                (x, y) => x > y,
            },
        };
        let (g1, g2) = if a_survives { (a, b) } else { (b, a) };
        merge(&mut genes, g1, g2);
        n_novel -= 1;
        // Rescore: drop pairs with g2, recompute pairs with g1 (4.1 keeps dict positions).
        scores.retain(|&(x, y), _| x != g2 && y != g2);
        let s1 = &genes[g1 as usize];
        for (x, gx) in genes.iter().enumerate() {
            let x = x as u32;
            if x == g1 || gx.merged_into != u32::MAX || (s1.annotated && gx.annotated) {
                continue;
            }
            let k = (g1.min(x), g1.max(x));
            let s = score(s1, gx);
            if s > 0.0 {
                scores.insert(k, s);
            } else {
                scores.remove(&k);
            }
        }
    }

    // --- rewrite gene ids -------------------------------------------------------------------
    let root = |mut p: u32| {
        while genes[p as usize].merged_into != u32::MAX {
            p = genes[p as usize].merged_into;
        }
        genes[p as usize].key
    };
    for m in models.iter_mut() {
        let start = match m.kind {
            ModelKind::Known(t) => tx_gene.get(&t).copied(),
            _ => pos.get(&m.gene).copied(),
        };
        if let Some(p) = start {
            m.gene = root(p);
        }
    }
}

/// `merge_genes(gene1, gene2)` (`:85-101`) minus the score update.
fn merge(genes: &mut [Gene], g1: u32, g2: u32) {
    let (lo, hi) = genes.split_at_mut(g1.max(g2) as usize);
    let (a, b) = if g1 < g2 { (&mut lo[g1 as usize], &mut hi[0]) } else { (&mut hi[0], &mut lo[g2 as usize]) };
    a.region = Iv::new(a.region.s.min(b.region.s), a.region.e.max(b.region.e));
    let mut introns = Vec::with_capacity(a.introns.len() + b.introns.len());
    introns.extend_from_slice(&a.introns);
    introns.extend_from_slice(&b.introns);
    introns.sort_unstable();
    introns.dedup();
    a.introns = introns;
    a.update_hull();
    b.merged_into = g1;
    b.introns = Vec::new();
}

pub(crate) fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    reg.add::<fixture::U18Join>();
}

pub(crate) mod fixture {
    //! U18 runner. No new genes or models are created, so output gene ids are the input
    //! 4.1 strings; novel-novel merge direction uses the 4.1 id strings.
    use super::*;
    use iq_types::fixture::convert::ModelExtras;
    use iq_types::fixture::schema::{ModelFixture, ModelsOut, U18JoinIn};
    use iq_types::fixture::{FixResult, FixtureCtx, UnitFixture};

    pub struct U18Join;

    impl UnitFixture for U18Join {
        const UNIT: &'static str = "U18";
        const FN: &'static str = "join_transcripts";
        type In = U18JoinIn;
        type Out = ModelsOut;

        fn run(input: &U18JoinIn, ctx: &FixtureCtx) -> FixResult<ModelsOut> {
            let index = ctx.index(&input.index)?;
            let mut ids = ctx.ids().clone();
            ids.reset_block();
            let (mut models, extras): (Vec<TranscriptModel>, Vec<ModelExtras>) =
                input.models.iter().map(|m| m.to_core(&mut ids)).unzip();
            join_genes_by(&mut models, &index, |a, b| ids.novel_genes.resolve(a).cmp(ids.novel_genes.resolve(b)));
            let models = models.iter().zip(&extras).map(|(m, e)| ModelFixture::from_core(m, &ids, e)).collect();
            Ok(ModelsOut { models })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use iq_types::{GeneRec, IsoRec};

    fn model(key: u32, gene: GeneKey, exons: &[(u32, u32)]) -> TranscriptModel {
        TranscriptModel {
            key,
            strand: Strand::Plus,
            gene,
            kind: ModelKind::Nnic,
            exons: exons.iter().map(|&(s, e)| Iv::new(s as _, e as _)).collect(),
            intron_path: Default::default(),
            polya_confirmed: false,
            tss_confirmed: false,
            source_known: None,
        }
    }

    fn index_with_gene(region: (u32, u32)) -> IsoformIndex {
        // A gene is joinable only if it has a transcript (`gene_to_transcripts`).
        let r = Iv::new(region.0 as _, region.1 as _);
        IsoformIndex {
            iso: vec![IsoRec {
                key: IsoRef::Ref(0),
                gene: GeneKey::Ref(0),
                strand: Strand::Plus,
                span: r,
                exon_off: 0,
                n_exon: 1,
                intron_off: 0,
                n_intron: 0,
            }],
            iso_exons: vec![r],
            genes: vec![GeneRec { key: GeneKey::Ref(0), strand: Strand::Plus, region: r }],
            ..Default::default()
        }
    }

    /// A single (annotated, novel) pair is never merged (4.1 `len(scores) > 1`).
    #[test]
    fn single_pair() {
        let index = index_with_gene((100, 200));
        let mut ms = vec![model(0, GeneKey::Novel(0), &[(120, 150), (170, 220)])];
        join_genes(&mut ms, &index);
        assert_eq!(ms[0].gene, GeneKey::Novel(0));
    }

    /// Two novel genes sharing an intron: the comparator-greater one survives; weak overlap
    /// (< 0.5) does not merge.
    #[test]
    fn novel_pair_direction_and_weak_overlap() {
        let index = IsoformIndex::default();
        let mut ms = vec![
            model(0, GeneKey::Novel(0), &[(100, 150), (200, 250)]),
            model(1, GeneKey::Novel(1), &[(120, 150), (200, 300)]),
            model(2, GeneKey::Novel(2), &[(280, 600)]),
        ];
        join_genes_by(&mut ms, &index, |a, b| b.cmp(&a));
        assert_eq!(ms[0].gene, GeneKey::Novel(0));
        assert_eq!(ms[1].gene, GeneKey::Novel(0));
        assert_eq!(ms[2].gene, GeneKey::Novel(2));
    }
}
