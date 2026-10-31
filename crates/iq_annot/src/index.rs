//! The run-wide annotation index (design §2.2): flat tables, built once, `Arc`-shared.
//!
//! Layout:
//! - `genes` sorted by `(chr, start, file order)`;
//! - `txs` grouped by gene (gene `g` owns `txs[g.tx_off .. g.tx_off + g.n_tx]`), within a gene
//!   ordered by `(start, file order)` = gffutils `children(order_by='start')`;
//! - `exons` flat per transcript, after the CDS-as-exon rule (§4.2);
//! - `other` flat per transcript (CDS/start_codon/stop_codon/UTR children, file order);
//! - all strings in one arena (`StrRef`).
//!
//! Attribute strings are kept raw (the GFF3 column 9); the gffutils view
//! ([`AnnotationIndex::gene_attributes`]) and the GTF attribute string 4.1 prints
//! ([`AnnotationIndex::gene_attrs`], `set_gene_attributes`) are derived on demand.

use crate::attrs;
use iq_types::smallvec::SmallVec;
use iq_types::{quirk, ChrIdx, DevId, GeneIdx, Iv, Iv0, Pos, Strand, TxIdx};
use std::collections::HashMap;

/// Gene selection mode (design §4.2, D-A1).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Default)]
pub enum GeneTypes {
    /// Every parentless feature with transcript children that have exon/CDS children
    /// (default; deviation D-A1).
    #[default]
    All,
    /// 4.1 behaviour: type exactly `gene`, transcripts of `DEFAULT_RNA_TYPES`. Used for
    /// fixtures and parity runs.
    Isoquant,
}

/// `(offset, len)` into the string arena.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash)]
pub(crate) struct StrRef {
    pub off: u32,
    pub len: u32,
}

/// One gene row.
#[derive(Clone, Debug)]
pub(crate) struct GeneRow {
    pub chr: ChrIdx,
    pub iv: Iv,
    pub strand: Strand,
    /// Index into `AnnotationIndex::types`.
    pub ty: u16,
    /// Index into `AnnotationIndex::sources`.
    pub source: u16,
    pub id: StrRef,
    /// Raw GFF3 attribute column.
    pub attrs: StrRef,
    pub tx_off: u32,
    pub n_tx: u32,
}

/// One transcript row (a transcript under two genes appears once per gene).
#[derive(Clone, Debug)]
pub(crate) struct TxRow {
    pub gene: GeneIdx,
    /// The transcript row's own GFF coordinates (not the exon span).
    pub iv: Iv,
    pub strand: Strand,
    pub ty: u16,
    pub source: u16,
    pub id: StrRef,
    pub attrs: StrRef,
    pub exon_off: u32,
    pub n_exon: u32,
    pub other_off: u32,
    pub n_other: u32,
    /// Position among the gene's transcripts in gffutils primary-key order (= unordered
    /// `children()` iteration; 4.1 dict order of `isoform_strands`, `gene_id_map`, ...).
    pub file_rank: u32,
}

/// Flat annotation index. Field layout is the owner's choice; the query API below is the
/// contract other units use.
#[derive(Debug, Default)]
pub struct AnnotationIndex {
    pub(crate) gene_types: GeneTypes,
    pub(crate) arena: String,
    pub(crate) chroms: Vec<String>,
    pub(crate) chrom_len: Vec<Option<u64>>,
    pub(crate) chrom_map: HashMap<String, ChrIdx>,
    pub(crate) types: Vec<String>,
    pub(crate) sources: Vec<String>,
    pub(crate) genes: Vec<GeneRow>,
    pub(crate) txs: Vec<TxRow>,
    pub(crate) exons: Vec<Iv>,
    /// `(iv, type index)`.
    pub(crate) other: Vec<(Iv, u16)>,
    /// Per chromosome: `[lo, hi)` range of `genes`.
    pub(crate) chr_genes: Vec<(u32, u32)>,
    /// Running max of gene ends within each chromosome's range (same indexing as `genes`).
    pub(crate) prefix_max_end: Vec<Pos>,
    /// `(chr, start, end, strand)` → `exon_id` (§8.1; last exon row in file order wins).
    pub(crate) exon_ids: HashMap<(ChrIdx, Pos, Pos, Strand), StrRef>,
    /// Per chromosome: sorted novel-id numbers taken by input ids (`id_policy.py:20-45`).
    pub(crate) forbidden: Vec<Vec<u64>>,
}

impl AnnotationIndex {
    #[inline]
    pub(crate) fn s(&self, r: StrRef) -> &str {
        &self.arena[r.off as usize..(r.off + r.len) as usize]
    }

    /// Gene selection mode this index was built with.
    pub fn gene_types(&self) -> GeneTypes {
        self.gene_types
    }
    pub fn n_chroms(&self) -> usize {
        self.chroms.len()
    }
    pub fn n_genes(&self) -> usize {
        self.genes.len()
    }
    pub fn n_txs(&self) -> usize {
        self.txs.len()
    }

    /// Chromosome index by name.
    pub fn chrom(&self, name: &str) -> Option<ChrIdx> {
        self.chrom_map.get(name).copied()
    }
    pub fn chrom_name(&self, chr: ChrIdx) -> &str {
        &self.chroms[chr as usize]
    }
    /// Length from `##sequence-region`, if the file declared one.
    pub fn chrom_len(&self, chr: ChrIdx) -> Option<u64> {
        self.chrom_len.get(chr as usize).copied().flatten()
    }
    /// Genes of a chromosome, as an index range (sorted by `(start, file order)`).
    pub fn chrom_genes(&self, chr: ChrIdx) -> std::ops::Range<GeneIdx> {
        match self.chr_genes.get(chr as usize) {
            Some(&(lo, hi)) => lo..hi,
            None => 0..0,
        }
    }

    /// Genes with `start <= iv.e && end >= iv.s` (1-based closed), sorted by
    /// `(start, file order)` (§4.3). O(log G + k).
    pub fn genes_overlapping(&self, chr: ChrIdx, iv: Iv) -> SmallVec<[GeneIdx; 4]> {
        let mut out = SmallVec::new();
        let Some(&(lo, hi)) = self.chr_genes.get(chr as usize) else { return out };
        let (lo, hi) = (lo as usize, hi as usize);
        let genes = &self.genes[lo..hi];
        let end = genes.partition_point(|g| g.iv.s <= iv.e);
        let start = self.prefix_max_end[lo..lo + end].partition_point(|&e| e < iv.s);
        for (k, g) in genes[start..end].iter().enumerate() {
            if g.iv.e >= iv.s {
                out.push((lo + start + k) as GeneIdx);
            }
        }
        out
    }

    /// Genes for a region as IsoQuant passes it (`AbstractAlignmentStorage.region`,
    /// 0-based closed): the gene query of `get_gene_info_for_region` /
    /// `get_genes_in_region` (`alignment_processor.py:583-605`).
    ///
    /// D-R1 lives here: iqrs converts the region to 1-based (`s+1, e+1`); 4.1 passed the
    /// 0-based tuple to gffutils unchanged (and gffutils treats `start == 0` as "no start
    /// bound", leaving the strict test `gene.start < end`).
    pub fn genes_for_region0(&self, chr: ChrIdx, r: Iv0) -> SmallVec<[GeneIdx; 4]> {
        if quirk(DevId::R1) {
            if r.s == 0 {
                // gffutils `region()`: falsy start → only `start < end` (interface.py:752-767)
                let range = self.chrom_genes(chr);
                let lo = range.start as usize;
                let n = self.genes[lo..range.end as usize].partition_point(|g| (g.iv.s as i64) < r.e);
                return (lo..lo + n).map(|g| g as GeneIdx).collect();
            }
            self.genes_overlapping(chr, Iv::new(r.s as Pos, r.e as Pos))
        } else {
            self.genes_overlapping(chr, Iv::new((r.s + 1) as Pos, (r.e + 1) as Pos))
        }
    }

    pub fn gene_id(&self, g: GeneIdx) -> &str {
        self.s(self.genes[g as usize].id)
    }
    pub fn gene_span(&self, g: GeneIdx) -> (ChrIdx, Iv, Strand) {
        let r = &self.genes[g as usize];
        (r.chr, r.iv, r.strand)
    }
    /// GFF feature type of the gene row (`gene`, `ncRNA_gene`, `pseudogene`, ...).
    pub fn gene_type(&self, g: GeneIdx) -> &str {
        &self.types[self.genes[g as usize].ty as usize]
    }
    pub fn gene_source(&self, g: GeneIdx) -> &str {
        &self.sources[self.genes[g as usize].source as usize]
    }
    /// Raw GFF3 attribute column of the gene row.
    pub fn gene_raw_attrs(&self, g: GeneIdx) -> &str {
        self.s(self.genes[g as usize].attrs)
    }
    /// gffutils `attributes` of the gene (key order = file order).
    pub fn gene_attributes(&self, g: GeneIdx) -> Vec<(String, Vec<String>)> {
        attrs::parse_gff3_attributes(self.gene_raw_attrs(g))
    }
    /// The GTF attribute string 4.1 prints for the gene (`set_gene_attributes`,
    /// `gene_info.py:576-583`).
    pub fn gene_attrs(&self, g: GeneIdx) -> String {
        attrs::format_attrs(&self.gene_attributes(g), attrs::GENE_SKIP)
    }
    /// Transcripts of a gene, as an index range (ordered by `(start, file order)`).
    pub fn gene_txs(&self, g: GeneIdx) -> std::ops::Range<TxIdx> {
        let r = &self.genes[g as usize];
        r.tx_off..r.tx_off + r.n_tx
    }

    pub fn tx_id(&self, t: TxIdx) -> &str {
        self.s(self.txs[t as usize].id)
    }
    pub fn tx_gene(&self, t: TxIdx) -> GeneIdx {
        self.txs[t as usize].gene
    }
    pub fn tx_strand(&self, t: TxIdx) -> Strand {
        self.txs[t as usize].strand
    }
    /// The transcript row's own coordinates.
    pub fn tx_row_span(&self, t: TxIdx) -> Iv {
        self.txs[t as usize].iv
    }
    pub fn tx_type(&self, t: TxIdx) -> &str {
        &self.types[self.txs[t as usize].ty as usize]
    }
    pub fn tx_source(&self, t: TxIdx) -> &str {
        &self.sources[self.txs[t as usize].source as usize]
    }
    pub fn tx_raw_attrs(&self, t: TxIdx) -> &str {
        self.s(self.txs[t as usize].attrs)
    }
    pub fn tx_attributes(&self, t: TxIdx) -> Vec<(String, Vec<String>)> {
        attrs::parse_gff3_attributes(self.tx_raw_attrs(t))
    }
    /// The GTF attribute string 4.1 prints for the transcript (`gene_info.py:584-589`).
    pub fn tx_attrs(&self, t: TxIdx) -> String {
        attrs::format_attrs(&self.tx_attributes(t), attrs::TX_SKIP)
    }
    /// The exon attribute string 4.1 prints for a known exon (`feature_attributes` exon
    /// keys, `gene_info.py:590-596`: transcript attributes minus the exon skip list).
    pub fn tx_exon_attrs(&self, t: TxIdx) -> String {
        attrs::format_attrs(&self.tx_attributes(t), attrs::EXON_SKIP)
    }
    /// Exons of a transcript (after the CDS-as-exon rule, §4.2). Empty for transcripts
    /// without exon/CDS children (4.1 keeps those in strand/gene maps only).
    pub fn tx_exons(&self, t: TxIdx) -> &[Iv] {
        let r = &self.txs[t as usize];
        &self.exons[r.exon_off as usize..(r.exon_off + r.n_exon) as usize]
    }
    /// `other_features` of a transcript: `(iv, type)` of its CDS/start_codon/stop_codon/UTR
    /// children in file order (`set_other_features`, `gene_info.py:554-560`).
    pub fn tx_other_features(&self, t: TxIdx) -> impl Iterator<Item = (Iv, &str)> + '_ {
        let r = &self.txs[t as usize];
        self.other[r.other_off as usize..(r.other_off + r.n_other) as usize]
            .iter()
            .map(|&(iv, ty)| (iv, self.types[ty as usize].as_str()))
    }
    /// Position of the transcript among its gene's transcripts in file order.
    pub fn tx_file_rank(&self, t: TxIdx) -> u32 {
        self.txs[t as usize].file_rank
    }

    /// Annotation `exon_id` of an exon (`FeatureIdStorage`, `id_policy.py:56-66`). For
    /// coordinates shared by several exon rows, the **last** row in file order wins (4.1 dict
    /// overwrite; design §2.2/§8.1 "first" is wrong).
    pub fn exon_id(&self, chr: ChrIdx, iv: Iv, strand: Strand) -> Option<&str> {
        self.exon_ids.get(&(chr, iv.s, iv.e, strand)).map(|&r| self.s(r))
    }
    pub fn n_exon_ids(&self) -> usize {
        self.exon_ids.len()
    }
    /// Novel-id numbers taken by input ids on `chr` (sorted; `ExcludingIdDistributor`).
    pub fn forbidden_ids(&self, chr: ChrIdx) -> &[u64] {
        self.forbidden.get(chr as usize).map(|v| v.as_slice()).unwrap_or(&[])
    }

    /// Approximate heap footprint in bytes.
    pub fn heap_bytes(&self) -> usize {
        use std::mem::size_of;
        self.arena.capacity()
            + self.chroms.iter().map(|s| s.capacity() + size_of::<String>()).sum::<usize>()
            + self.genes.capacity() * size_of::<GeneRow>()
            + self.txs.capacity() * size_of::<TxRow>()
            + self.exons.capacity() * size_of::<Iv>()
            + self.other.capacity() * size_of::<(Iv, u16)>()
            + self.prefix_max_end.capacity() * 4
            + self.exon_ids.capacity() * (size_of::<((ChrIdx, Pos, Pos, Strand), StrRef)>() + 1)
            + self.forbidden.iter().map(|v| v.capacity() * 8).sum::<usize>()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse::parse_gff3_reader;
    use iq_types::{with_quirks, QuirkSet, QUIRKS41};

    const GFF: &str = "1\ts\tgene\t100\t900\t.\t+\t.\tID=A\n1\ts\tgene\t50\t60\t.\t+\t.\tID=B\n\
                       1\ts\tgene\t100\t200\t.\t-\t.\tID=C\n1\ts\tgene\t5000\t6000\t.\t+\t.\tID=D\n";

    fn ids(a: &AnnotationIndex, v: &[GeneIdx]) -> Vec<String> {
        v.iter().map(|&g| a.gene_id(g).to_owned()).collect()
    }

    #[test]
    fn overlap_query() {
        let a = parse_gff3_reader(GFF.as_bytes(), GeneTypes::Isoquant).unwrap();
        // sorted by (start, file order)
        assert_eq!(ids(&a, &a.genes_overlapping(0, Iv::new(1, 10000))), ["B", "A", "C", "D"]);
        assert_eq!(ids(&a, &a.genes_overlapping(0, Iv::new(201, 4999))), ["A"]);
        assert_eq!(ids(&a, &a.genes_overlapping(0, Iv::new(60, 60))), ["B"]);
        assert!(a.genes_overlapping(0, Iv::new(901, 4999)).is_empty());
        assert!(a.genes_overlapping(1, Iv::new(1, 10)).is_empty());
    }

    /// D-R1: 4.1 queries gffutils with the 0-based tuple; iqrs converts to 1-based.
    #[test]
    fn d_r1_region_shift() {
        let a = parse_gff3_reader(GFF.as_bytes(), GeneTypes::Isoquant).unwrap();
        let q = |r: Iv0| ids(&a, &a.genes_for_region0(0, r));
        // 0-based [900, 950] = 1-based [901, 951]: A ends at 900
        assert!(with_quirks(QuirkSet::NONE, || q(Iv0 { s: 900, e: 950 })).is_empty());
        // 0-based [4998, 4998] = 1-based 4999: D starts at 5000
        assert!(with_quirks(QuirkSet::NONE, || q(Iv0 { s: 4998, e: 4998 })).is_empty());
        assert_eq!(with_quirks(QuirkSet::NONE, || q(Iv0 { s: 4999, e: 4999 })), ["D"]);
        assert_eq!(with_quirks(QuirkSet::NONE, || q(Iv0 { s: 0, e: 99 })), ["B", "A", "C"]);
        if QUIRKS41 {
            let r1 = QuirkSet::only(DevId::R1);
            assert_eq!(with_quirks(r1, || q(Iv0 { s: 900, e: 950 })), ["A"]);
            assert_eq!(with_quirks(r1, || q(Iv0 { s: 4998, e: 4999 })), Vec::<String>::new());
            assert_eq!(with_quirks(r1, || q(Iv0 { s: 4999, e: 4999 })), Vec::<String>::new());
            assert_eq!(with_quirks(r1, || q(Iv0 { s: 5000, e: 5000 })), ["D"]);
            // gffutils falsy start: only `start < end`
            assert_eq!(with_quirks(r1, || q(Iv0 { s: 0, e: 100 })), ["B"]);
        }
    }
}
