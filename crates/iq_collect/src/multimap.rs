//! U9: multimapper resolution (`assignment/multimap_resolver.py:32-210`), the grouping of
//! `prepare_multimapper_dict` (`processed_read_manager.py:20-58`) and the application of
//! the result in `FullAssignmentLoader.load_read_assignment` (`assignment_loader.py:50-70`).
//!
//! Flow at the barrier:
//! 1. [`duplicate_groups`] over every block store of the run (in chromosome, block order):
//!    groups of reads sharing a QNAME, in 4.1 list order.
//! 2. [`BasicRead::from_view`] per member, [`resolve`] per group.
//! 3. [`apply_resolution`] per block: patch types / multimapper flags and drop suspended
//!    reads (re-indexing the block like the 4.1 loader does).
//!
//! As in 4.1, `BasicReadAssignment.penalty_score` is always 0 (`min(0.0, penalty)`; see
//! [`BasicRead::from_view`]).

use iq_types::{AssignmentType, BlockStore, GeneKey, IsoRef, Iv, MultimapStrategy, Pos, ReadView};
use std::collections::HashMap;

/// `BasicReadAssignment` fields used by the resolver.
#[derive(Clone, Debug, PartialEq)]
pub struct BasicRead {
    /// (block, row) of the read in its `BlockStore` (caller-defined numbering).
    pub block: u32,
    pub row: u32,
    /// Chromosome (any consistent numbering; part of the duplicate test).
    pub chr: u32,
    /// `exons[0][0]`, `exons[-1][1]` (0 if no exons).
    pub start: Pos,
    pub end: Pos,
    pub genomic_region: Iv,
    pub multimapper: bool,
    pub polya_found: bool,
    pub atype: AssignmentType,
    pub gene_atype: AssignmentType,
    /// `min(0, matches[0].penalty)`, i.e. 0 in practice (as in 4.1).
    pub penalty: f64,
    /// Distinct genes of the matches, sorted (a set).
    pub genes: Vec<GeneKey>,
    /// Distinct transcripts of the matches, sorted (a set).
    pub isoforms: Vec<IsoRef>,
}

impl BasicRead {
    /// `BasicReadAssignment(read_assignment)` (`isoform_assignment.py:550-591`).
    pub fn from_view(v: &ReadView<'_>, block: u32, chr: u32) -> Self {
        let (start, end) = v.exons_span();
        let mut genes: Vec<GeneKey> = v.matches().filter_map(|m| m.gene()).collect();
        genes.sort_unstable();
        genes.dedup();
        let mut isoforms: Vec<IsoRef> = v.matches().filter_map(|m| m.tx()).collect();
        isoforms.sort_unstable();
        isoforms.dedup();
        let p0 = v.first_match().map_or(0.0, |m| m.penalty());
        // Matches IsoQuant 4.1: `BasicReadAssignment.penalty_score` is effectively always 0.
        let penalty = if v.n_matches() > 0 { p0.min(0.0) } else { 0.0 };
        BasicRead {
            block,
            row: v.row() as u32,
            chr,
            start,
            end,
            genomic_region: v.genomic_region(),
            multimapper: v.multimapper(),
            polya_found: v.polya_found(),
            atype: v.atype(),
            gene_atype: v.gene_atype(),
            penalty,
            genes,
            isoforms,
        }
    }

    /// `BasicReadAssignment.__eq__` (same read name is implied by the group).
    fn same(&self, o: &BasicRead) -> bool {
        self.chr == o.chr && self.start == o.start && self.end == o.end && self.isoforms == o.isoforms
    }
}

/// `MultimapResolver.resolve` for one group of same-name reads (in 4.1 list order).
/// Groups of size ≤ 1 are left unchanged.
pub fn resolve(group: &mut [BasicRead], strategy: MultimapStrategy) {
    if group.len() <= 1 {
        return;
    }
    match strategy {
        MultimapStrategy::IgnoreMultimapper => {
            for a in group.iter_mut() {
                a.atype = AssignmentType::Suspended;
            }
        }
        MultimapStrategy::Merge => {
            // 4.1 passes a set to find_duplicates, which indexes it (TypeError for > 1
            // informative assignment); the sorted index list is the evident intent.
            let keep: Vec<usize> =
                (0..group.len()).filter(|&i| group[i].atype != AssignmentType::Noninformative).collect();
            filter_assignments(group, keep);
        }
        MultimapStrategy::TakeBest => select_best_assignment(group),
    }
}

fn select_best_assignment(g: &mut [BasicRead]) {
    let (mut primary_unique, mut consistent, mut inconsistent, mut primary_inconsistent, mut noninformative) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for (i, a) in g.iter().enumerate() {
        if a.atype.is_inconsistent() {
            if !a.multimapper {
                primary_inconsistent.push(i);
            }
            inconsistent.push(i);
        } else if a.atype.is_consistent() {
            consistent.push(i);
            if !a.multimapper && a.atype != AssignmentType::Ambiguous {
                primary_unique.push(i);
            }
        } else {
            noninformative.push(i);
        }
    }
    if !primary_unique.is_empty() {
        filter_assignments(g, primary_unique);
    } else if !consistent.is_empty() {
        filter_assignments(g, consistent);
    } else if !primary_inconsistent.is_empty() {
        select_best_inconsistent(g, primary_inconsistent);
    } else if !inconsistent.is_empty() {
        select_best_inconsistent(g, inconsistent);
    } else if !noninformative.is_empty() {
        select_noninformative(g, &noninformative);
    }
    // else: 4.1 returns `assignment_list[0:1]` (unreachable: every type is in one class)
}

fn select_best_inconsistent(g: &mut [BasicRead], idx: Vec<usize>) {
    if idx.len() <= 1 {
        return filter_assignments(g, idx);
    }
    let best = idx.iter().map(|&i| g[i].penalty).fold(f64::INFINITY, f64::min);
    let keep = idx.into_iter().filter(|&i| g[i].penalty == best).collect();
    filter_assignments(g, keep)
}

fn select_noninformative(g: &mut [BasicRead], idx: &[usize]) {
    let ov = |a: &BasicRead| -> i64 {
        // intersection_len(genomic_region, (start, end))
        let s = a.genomic_region.s.max(a.start) as i64;
        let e = a.genomic_region.e.min(a.end) as i64;
        (e - s + 1).max(0)
    };
    let max_ov = idx.iter().map(|&i| ov(&g[i])).max().unwrap_or(0);
    let mut best = None;
    let mut min_start = i64::MAX;
    for &i in idx {
        if ov(&g[i]) == max_ov && (g[i].genomic_region.s as i64) < min_start {
            min_start = g[i].genomic_region.s as i64;
            best = Some(i);
        }
    }
    filter_assignments(g, best.into_iter().collect())
}

/// `find_duplicates`: keeps the first of each group of equal assignments.
fn find_duplicates(g: &[BasicRead], idx: Vec<usize>) -> Vec<usize> {
    if idx.len() <= 1 {
        return idx;
    }
    let mut discarded = vec![false; idx.len()];
    let mut out = Vec::with_capacity(idx.len());
    for i in 0..idx.len() {
        if discarded[i] {
            continue;
        }
        out.push(idx[i]);
        for j in i + 1..idx.len() {
            if !discarded[j] && g[idx[i]].same(&g[idx[j]]) {
                discarded[j] = true;
            }
        }
    }
    out
}

/// `filter_assignments`: keep `keep` (deduplicated), mark the rest suspended; several
/// genes / isoforms among the kept ones make them ambiguous multimappers.
fn filter_assignments(g: &mut [BasicRead], keep: Vec<usize>) {
    let keep = find_duplicates(g, keep);
    let mut genes: Vec<GeneKey> = Vec::new();
    let mut isos: Vec<IsoRef> = Vec::new();
    for &i in &keep {
        genes.extend_from_slice(&g[i].genes);
        isos.extend_from_slice(&g[i].isoforms);
    }
    genes.sort_unstable();
    genes.dedup();
    isos.sort_unstable();
    isos.dedup();
    let change_tx = isos.len() > 1;
    let change_gene = genes.len() > 1;
    for (i, a) in g.iter_mut().enumerate() {
        if keep.contains(&i) {
            let amb = if a.atype.is_inconsistent() {
                AssignmentType::InconsistentAmbiguous
            } else {
                AssignmentType::Ambiguous
            };
            if change_tx {
                a.atype = amb;
                a.multimapper = true;
            }
            if change_gene {
                a.gene_atype = amb;
                a.multimapper = true;
            }
        } else {
            a.atype = AssignmentType::Suspended;
            a.gene_atype = AssignmentType::Suspended;
        }
    }
}

// ------------------------------------------------------------------------------------------
// Barrier helpers
// ------------------------------------------------------------------------------------------

/// Groups of reads sharing a QNAME across `stores` (given in chromosome, then block order),
/// each group in 4.1 list order (store order, then row order). Hash collisions are split by
/// comparing names. Groups are returned ordered by their first member.
pub fn duplicate_groups(stores: &[&BlockStore]) -> Vec<Vec<(u32, u32)>> {
    let mut count: HashMap<u64, u32> = HashMap::new();
    for st in stores {
        for v in st.iter() {
            *count.entry(v.name_hash()).or_default() += 1;
        }
    }
    count.retain(|_, c| *c > 1);
    if count.is_empty() {
        return Vec::new();
    }
    // (hash, name) → group index
    let mut groups: Vec<Vec<(u32, u32)>> = Vec::new();
    let mut by_hash: HashMap<u64, Vec<usize>> = HashMap::new();
    for (b, st) in stores.iter().enumerate() {
        for v in st.iter() {
            let h = v.name_hash();
            if !count.contains_key(&h) {
                continue;
            }
            let cands = by_hash.entry(h).or_default();
            let name = v.name();
            let hit = cands.iter().copied().find(|&gi| {
                let (b0, r0) = groups[gi][0];
                stores[b0 as usize].view(r0 as usize).name() == name
            });
            match hit {
                Some(gi) => groups[gi].push((b as u32, v.row() as u32)),
                None => {
                    cands.push(groups.len());
                    groups.push(vec![(b as u32, v.row() as u32)]);
                }
            }
        }
    }
    groups.retain(|g| g.len() > 1);
    groups
}

/// Applies resolved reads of one block: overwrites types and multimapper flags, and drops
/// suspended reads (rows re-indexed densely, like `ReadAssignmentLoader.get_next`).
/// `resolved` must all belong to `store`. Returns the store unchanged when nothing was
/// resolved.
pub fn apply_resolution(mut store: BlockStore, resolved: &[BasicRead]) -> BlockStore {
    if resolved.is_empty() {
        return store;
    }
    let mut drop = vec![false; store.len()];
    for r in resolved {
        let i = r.row as usize;
        if r.atype == AssignmentType::Suspended {
            drop[i] = true;
        } else {
            store.set_types(i, r.atype, r.gene_atype);
            store.set_multimapper(i, r.multimapper);
        }
    }
    if !drop.iter().any(|&d| d) {
        return store;
    }
    let mut out = BlockStore::new(store.region);
    for (i, &d) in drop.iter().enumerate() {
        if !d {
            let mut rec = store.record(i);
            rec.idx = out.len() as u32;
            out.push(&rec);
        }
    }
    out
}

/// Per-chromosome result of [`resolve_chr_duplicates`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct DupOutcome {
    /// 4.1 `total_assignments` of the chromosome: reads kept after resolution.
    pub total: u64,
    /// 4.1 `polya_assignments`: kept reads with `polyA_found`.
    pub polya: u64,
    /// Same-name groups of ≥ 2 reads that were resolved.
    pub groups: u64,
    /// Reads suspended (and removed from their stores).
    pub suspended: u64,
}

/// `ProcessedReadsManagerNoSecondary.finalize(chr)` (`processed_read_manager.py:20-87,
/// 165-186`) plus the loader-side application (`assignment_loader.py:50-70`) for one
/// chromosome: `stores` = all its blocks in genomic (= 4.1 dump) order.
///
/// - Same-name groups are built by real name, members in (block, row) order.
///   `candidates` (sorted name hashes) prefilter the reads to look at: it must contain every
///   name that can occur more than once among `stores`, i.e. the pass-A duplicate
///   candidates **and** the reads collected by more than one block (4.1 double processing,
///   [`crate::BlockCollector::finish_shared`]; such a record is not a pass-A duplicate).
///   `None` considers every read.
/// - Groups of ≥ 2 are resolved with [`resolve`]; types and multimapper flags are patched.
/// - **Suspended reads are removed** from their store, and the remaining rows are re-indexed
///   densely. This matches what `ReadAssignmentLoader.get_next` hands to every construct
///   and quant stage, so downstream code needs no filter.
pub fn resolve_chr_duplicates(
    stores: &mut [BlockStore],
    candidates: Option<&[u64]>,
    strategy: MultimapStrategy,
) -> DupOutcome {
    let mut out = DupOutcome::default();
    let is_cand = |h: u64| candidates.is_none_or(|c| c.binary_search(&h).is_ok());
    // groups by name, in first-appearance order
    let mut groups: Vec<Vec<(u32, u32)>> = Vec::new();
    {
        let mut by_name: HashMap<&str, usize> = HashMap::new();
        for (b, st) in stores.iter().enumerate() {
            for v in st.iter() {
                if !is_cand(v.name_hash()) {
                    out.total += 1;
                    out.polya += v.polya_found() as u64;
                    continue;
                }
                let gi = *by_name.entry(v.name()).or_insert_with(|| {
                    groups.push(Vec::new());
                    groups.len() - 1
                });
                groups[gi].push((b as u32, v.row() as u32));
            }
        }
    }
    let mut patches: Vec<Vec<BasicRead>> = vec![Vec::new(); stores.len()];
    for g in &groups {
        if g.len() == 1 {
            let (b, r) = g[0];
            out.total += 1;
            out.polya += stores[b as usize].view(r as usize).polya_found() as u64;
            continue;
        }
        let mut reads: Vec<BasicRead> =
            g.iter().map(|&(b, r)| BasicRead::from_view(&stores[b as usize].view(r as usize), b, 0)).collect();
        resolve(&mut reads, strategy);
        out.groups += 1;
        for r in reads {
            if r.atype == AssignmentType::Suspended {
                out.suspended += 1;
            } else {
                out.total += 1;
                out.polya += r.polya_found as u64;
            }
            patches[r.block as usize].push(r);
        }
    }
    for (st, p) in stores.iter_mut().zip(&patches) {
        if !p.is_empty() {
            *st = apply_resolution(std::mem::take(st), p);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use AssignmentType as A;

    fn br(atype: A, start: Pos, isos: &[u32], genes: &[u32], region: (i32, i32)) -> BasicRead {
        BasicRead {
            block: 0,
            row: 0,
            chr: 0,
            start,
            end: start + 100,
            genomic_region: Iv::new(region.0, region.1),
            multimapper: false,
            polya_found: false,
            atype,
            gene_atype: atype,
            penalty: 0.0,
            genes: genes.iter().map(|&g| GeneKey::Ref(g)).collect(),
            isoforms: isos.iter().map(|&t| IsoRef::Ref(t)).collect(),
        }
    }

    #[test]
    fn equal_copies_collapse() {
        let mut g = vec![br(A::Unique, 10, &[1], &[1], (0, 50)), br(A::Unique, 10, &[1], &[1], (51, 200))];
        resolve(&mut g, MultimapStrategy::TakeBest);
        assert_eq!(g[0].atype, A::Unique);
        assert!(!g[0].multimapper);
        assert_eq!(g[1].atype, A::Suspended);
    }

    #[test]
    fn unequal_copies_ambiguous() {
        let mut g = vec![br(A::Unique, 10, &[1], &[1], (0, 50)), br(A::Unique, 10, &[2], &[1], (51, 200))];
        resolve(&mut g, MultimapStrategy::TakeBest);
        assert_eq!((g[0].atype, g[0].gene_atype, g[0].multimapper), (A::Ambiguous, A::Unique, true));
        assert_eq!(g[1].atype, A::Ambiguous);
    }

    #[test]
    fn inconsistent_penalty_and_noninformative() {
        let mut g = vec![br(A::Inconsistent, 10, &[1], &[1], (0, 50)), br(A::Inconsistent, 20, &[2], &[2], (0, 50))];
        g[0].penalty = 2.0;
        g[1].penalty = 1.0;
        resolve(&mut g, MultimapStrategy::TakeBest);
        assert_eq!(g[0].atype, A::Suspended);
        assert_eq!(g[1].atype, A::Inconsistent);
        // noninformative: largest overlap with region, then lowest region start
        let mut g = vec![br(A::Noninformative, 10, &[], &[], (60, 300)), br(A::Noninformative, 10, &[], &[], (0, 300))];
        resolve(&mut g, MultimapStrategy::TakeBest);
        assert_eq!(g[0].atype, A::Suspended);
        assert_eq!(g[1].atype, A::Noninformative);
    }
}
