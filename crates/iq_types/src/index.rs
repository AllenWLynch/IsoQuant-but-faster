//! `IsoformIndex` (the 4.1 `GeneInfo` equivalent) and read profiles (U5 output
//! = U6 input).
//!
//! Builders live elsewhere: `IsoformIndex::from_genes` / `from_models` in `iq_annot` (U1),
//! and `IsoformIndex::from_dump` (Python dump) in `fixture::convert`.

use crate::assign::{GeneKey, IsoRef};
use crate::prim::{ChrIdx, Iv, PolyAInfo, Pos, Strand};
use serde::{Deserialize, Serialize};

/// Dense `i8` profile matrix `[iso × feature]` with 4.1 `FeatureProfiles.set_profiles`
/// semantics (`gene_info.py:96-118`): `1` = feature in isoform, `-1` = absent but inside the
/// isoform span, `-2` = outside the isoform span.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct ProfMatrix {
    pub n_iso: u32,
    pub n_feat: u32,
    /// Row-major, `n_iso * n_feat`.
    pub data: Vec<i8>,
}

impl ProfMatrix {
    pub fn new(n_iso: usize, n_feat: usize) -> Self {
        ProfMatrix { n_iso: n_iso as u32, n_feat: n_feat as u32, data: vec![-1; n_iso * n_feat] }
    }
    #[inline]
    pub fn row(&self, iso: usize) -> &[i8] {
        let f = self.n_feat as usize;
        &self.data[iso * f..(iso + 1) * f]
    }
    #[inline]
    pub fn row_mut(&mut self, iso: usize) -> &mut [i8] {
        let f = self.n_feat as usize;
        &mut self.data[iso * f..(iso + 1) * f]
    }
}

/// Isoform spans sorted by start with a prefix-max of ends, for overlap/containment
/// queries in O(log I + k).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SortedSpans {
    /// Isoform indices (into `IsoformIndex::iso`) sorted by `(span.s, iso index)`.
    pub order: Vec<u32>,
    pub starts: Vec<Pos>,
    /// `prefix_max_end[k] = max(span.e of order[0..=k])`.
    pub prefix_max_end: Vec<Pos>,
}

impl SortedSpans {
    pub fn build(spans: &[Iv]) -> Self {
        let mut order: Vec<u32> = (0..spans.len() as u32).collect();
        order.sort_by_key(|&i| (spans[i as usize].s, i));
        let starts = order.iter().map(|&i| spans[i as usize].s).collect();
        let mut prefix_max_end = Vec::with_capacity(order.len());
        let mut m = Pos::MIN;
        for &i in &order {
            m = m.max(spans[i as usize].e);
            prefix_max_end.push(m);
        }
        SortedSpans { order, starts, prefix_max_end }
    }

    /// Indices (into `IsoformIndex::iso`) of spans overlapping `q`, in ascending isoform
    /// order (= `iso` order, so iteration order matches 4.1 dict order).
    pub fn overlapping(&self, spans: &[Iv], q: Iv) -> Vec<u32> {
        // candidates have start <= q.e
        let hi = self.starts.partition_point(|&s| s <= q.e);
        // first position whose prefix max end >= q.s
        let lo = self.prefix_max_end[..hi].partition_point(|&e| e < q.s);
        let mut out: Vec<u32> = self.order[lo..hi].iter().copied().filter(|&i| spans[i as usize].e >= q.s).collect();
        out.sort_unstable();
        out
    }
}

/// One isoform of an [`IsoformIndex`].
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct IsoRec {
    pub key: IsoRef,
    pub gene: GeneKey,
    pub strand: Strand,
    /// `(exons[0].s, exons[-1].e)`.
    pub span: Iv,
    pub exon_off: u32,
    pub n_exon: u32,
    pub intron_off: u32,
    pub n_intron: u32,
}

/// One gene of an [`IsoformIndex`] (`gene_strands`, gene regions).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct GeneRec {
    pub key: GeneKey,
    pub strand: Strand,
    pub region: Iv,
}

/// Profile kinds, index into the `[_; 3]` arrays of [`IsoformIndex`].
pub mod prof_kind {
    pub const INTRON: usize = 0;
    pub const EXON: usize = 1;
    pub const SPLIT_EXON: usize = 2;
}

/// Per-block isoform index (4.1 `GeneInfo` with `prepare_profiles=True`, restricted to the
/// fields collect and construct read). Shared read-only (`Arc`) by U5-U19.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct IsoformIndex {
    pub chr: ChrIdx,
    /// `GeneInfo.start/end`.
    pub start: Pos,
    pub end: Pos,
    pub delta: i32,
    /// Order = GeneInfo dict insertion order (genes by start, transcripts by start); for
    /// model indices, storage order.
    pub iso: Vec<IsoRec>,
    /// Rank of each isoform's id string in Python codepoint order; every
    /// `sorted(isoform_ids)` in 4.1 uses it. For model indices without real ids:
    /// creation order.
    pub id_rank: Vec<u32>,
    /// Flat per-isoform exons / introns; `IsoRec` holds offsets.
    pub iso_exons: Vec<Iv>,
    pub iso_introns: Vec<Iv>,
    /// Sorted feature unions (`intron_profiles.features`, `exon_profiles.features`,
    /// `split_exon_profiles.features`; `gene_info.py:646-647, 830-875`).
    pub introns: Vec<Iv>,
    pub exons: Vec<Iv>,
    pub split_exons: Vec<Iv>,
    pub intron_prof: ProfMatrix,
    pub exon_prof: ProfMatrix,
    pub split_prof: ProfMatrix,
    /// Per isoform, per kind (`prof_kind`): `profile_ranges` = first/last+1 index with value 1.
    pub prof_range: Vec<[(u32, u32); 3]>,
    /// Per isoform, per kind: index range of features overlapping the isoform span.
    pub span_feat: Vec<[(u32, u32); 3]>,
    pub genes: Vec<GeneRec>,
    pub span_tree: SortedSpans,
}

impl IsoformIndex {
    #[inline]
    pub fn n_iso(&self) -> usize {
        self.iso.len()
    }
    /// True if the index has no isoforms (4.1 `GeneInfo.empty()`).
    pub fn is_empty(&self) -> bool {
        self.iso.is_empty()
    }
    #[inline]
    pub fn iso_exons(&self, i: usize) -> &[Iv] {
        let r = &self.iso[i];
        &self.iso_exons[r.exon_off as usize..(r.exon_off + r.n_exon) as usize]
    }
    #[inline]
    pub fn iso_introns(&self, i: usize) -> &[Iv] {
        let r = &self.iso[i];
        &self.iso_introns[r.intron_off as usize..(r.intron_off + r.n_intron) as usize]
    }
    /// Features of a profile kind (`prof_kind`).
    pub fn features(&self, kind: usize) -> &[Iv] {
        match kind {
            prof_kind::INTRON => &self.introns,
            prof_kind::EXON => &self.exons,
            _ => &self.split_exons,
        }
    }
    /// Profile matrix of a kind (`prof_kind`).
    pub fn prof(&self, kind: usize) -> &ProfMatrix {
        match kind {
            prof_kind::INTRON => &self.intron_prof,
            prof_kind::EXON => &self.exon_prof,
            _ => &self.split_prof,
        }
    }
    /// Position of an isoform by key (linear; indices are small).
    pub fn find_iso(&self, key: IsoRef) -> Option<usize> {
        self.iso.iter().position(|r| r.key == key)
    }
    /// Isoform indices whose span overlaps `q`, in `iso` order.
    pub fn overlapping_isoforms(&self, q: Iv) -> Vec<u32> {
        let spans: Vec<Iv> = self.iso.iter().map(|r| r.span).collect();
        self.span_tree.overlapping(&spans, q)
    }

    /// Recomputes the derived fields `id_rank` (from `ids`, one per isoform), `span_feat`
    /// and `span_tree` from the primary fields. Used by every builder.
    pub fn finish_derived(&mut self, ids: &[&str]) {
        assert_eq!(ids.len(), self.iso.len());
        self.id_rank = id_ranks(ids);
        let spans: Vec<Iv> = self.iso.iter().map(|r| r.span).collect();
        self.span_tree = SortedSpans::build(&spans);
        self.span_feat = spans
            .iter()
            .map(|sp| {
                let f = |feats: &[Iv]| feature_window(feats, *sp);
                [f(&self.introns), f(&self.exons), f(&self.split_exons)]
            })
            .collect();
    }
}

/// Ranks of `ids` under Python string ordering (Unicode code point order; byte-wise UTF-8
/// comparison is equivalent). Equal ids get equal ranks.
pub fn id_ranks(ids: &[&str]) -> Vec<u32> {
    let mut order: Vec<usize> = (0..ids.len()).collect();
    order.sort_by(|&a, &b| ids[a].as_bytes().cmp(ids[b].as_bytes()));
    let mut rank = vec![0u32; ids.len()];
    let mut r = 0u32;
    for (k, &i) in order.iter().enumerate() {
        if k > 0 && ids[order[k - 1]] != ids[i] {
            r += 1;
        }
        rank[i] = r;
    }
    rank
}

/// `[lo, hi)` index window of sorted features overlapping `span` (features are sorted by
/// start; ends are assumed sorted too, which holds for 4.1 feature unions of disjoint-start
/// sorted lists only approximately — so the window is computed conservatively by scanning).
pub fn feature_window(feats: &[Iv], span: Iv) -> (u32, u32) {
    let mut lo = None;
    let mut hi = 0usize;
    for (k, f) in feats.iter().enumerate() {
        if f.overlaps(&span) {
            if lo.is_none() {
                lo = Some(k);
            }
            hi = k + 1;
        }
    }
    match lo {
        Some(l) => (l as u32, hi as u32),
        None => (0, 0),
    }
}

// ------------------------------------------------------------------------------------------
// Read profiles (U5 → U6, U7)
// ------------------------------------------------------------------------------------------

/// 4.1 `MappedReadProfile` (`assignment/long_read_profiles.py:20-26`), dense form.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct MappedReadProfile {
    /// One value per known feature of the index (`1`, `-1`, `0`, `-2`).
    pub gene_profile: Vec<i8>,
    /// One value per read feature.
    pub read_profile: Vec<i8>,
    pub read_features: Vec<Iv>,
    /// `gene_profile_range` = `(first, last + 1)` indices into `gene_profile`.
    pub gene_profile_range: (i32, i32),
}

/// 4.1 `CombinedReadProfiles` (`long_read_profiles.py:29-37`). `exon` is `None` in scope
/// (`count_exons = false`). CAGE is not ported.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct CombinedProfile {
    pub intron: MappedReadProfile,
    pub exon: Option<MappedReadProfile>,
    pub split_exon: MappedReadProfile,
    pub polya: PolyAInfo,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ranks_codepoint_order() {
        // 'transcript:ENS…' sorts after 'transcript12…' (':' = 0x3a > '1' = 0x31)
        let ids = ["transcript:ENSDART1", "transcript12.chr8.nic", "transcript9.chr8.nic", "transcript12.chr8.nic"];
        assert_eq!(id_ranks(&ids), vec![2, 0, 1, 0]);
    }

    #[test]
    fn spans_overlap() {
        let spans = [Iv::new(10, 20), Iv::new(5, 100), Iv::new(30, 40), Iv::new(15, 16)];
        let t = SortedSpans::build(&spans);
        assert_eq!(t.overlapping(&spans, Iv::new(18, 31)), vec![0, 1, 2]);
        assert_eq!(t.overlapping(&spans, Iv::new(101, 200)), Vec::<u32>::new());
        assert_eq!(t.overlapping(&spans, Iv::new(16, 16)), vec![0, 1, 3]);
    }

    #[test]
    fn window() {
        let f = [Iv::new(1, 5), Iv::new(10, 20), Iv::new(30, 40), Iv::new(50, 60)];
        assert_eq!(feature_window(&f, Iv::new(15, 35)), (1, 3));
        assert_eq!(feature_window(&f, Iv::new(41, 49)), (0, 0));
    }
}
