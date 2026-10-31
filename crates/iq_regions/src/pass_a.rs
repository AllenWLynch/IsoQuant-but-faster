//! Pass A: stream every record of a chromosome once, cut clusters at gaps not
//! bridged by a shared gene (`AlignmentCollector.process`), split each cluster at coverage
//! valleys (`forward_alignments` → `split_coverage_regions`), and annotate each block with
//! max coverage, the coverage cap, its genes and an estimate of the records it processes
//! (owned primary records plus the boundary-straddling reads it shares with a neighbour).
//!
//! Memory: per cluster two flat `u32` arrays over its 256-bp bins (coverage difference array
//! and the ownership midpoint histogram), the spans of reads longer than 32 kb, and the
//! duplicate-name Bloom filter. Nothing else is kept per read.

use iq_types::read::{flags, COVERAGE_BIN};
use iq_types::{Block, Cluster, GeneIdx, Iv0, RecordLite};

use crate::dups::DupFilter;
use crate::owner::{mid_key, needs_explicit, OwnerCounter};
use crate::split::{coverage_cutoff, skip_read_fraction, split_coverage_regions_into, CoverageBins};

/// Gene overlap query of one chromosome, as 4.1 issues it (`get_genes_in_region` /
/// `get_gene_info_for_region`, `alignment_processor.py:583-605`).
///
/// Production: `iq_annot::AnnotationIndex::genes_for_region0(chr, r)`, which passes the
/// 0-based tuple as 4.1 does.
/// Without a gene DB use [`NoGenes`] (4.1 returns the empty set).
pub trait GeneOracle {
    /// Appends the genes overlapping the 0-based closed storage region `r`, in
    /// `(start, file order)` order, to `out`.
    fn genes_for_region0(&self, r: Iv0, out: &mut Vec<GeneIdx>);
}

/// No gene database: every gap splits, every block is intergenic.
pub struct NoGenes;
impl GeneOracle for NoGenes {
    fn genes_for_region0(&self, _r: Iv0, _out: &mut Vec<GeneIdx>) {}
}

impl<G: GeneOracle + ?Sized> GeneOracle for &G {
    fn genes_for_region0(&self, r: Iv0, out: &mut Vec<GeneIdx>) {
        (**self).genes_for_region0(r, out)
    }
}

/// What pass A reads per record.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct PassARecord {
    pub pos0: i64,
    /// 0-based exclusive end (`reference_end`).
    pub end0: i64,
    pub flag: u16,
    /// `iq_types::name_hash(qname)`; only read for primary records (0 otherwise is fine).
    pub name_hash: u64,
}

impl PassARecord {
    /// Pass-A view of a BAM record, or `None` when pysam's `reference_end` is `None`
    /// (unmapped flag or no CIGAR). 4.1 would crash on such a record (`None - 1` in
    /// `alignment_is_not_adjacent`/`add_alignment`); iqrs skips it.
    pub fn from_lite(r: &RecordLite<'_>) -> Option<PassARecord> {
        if r.is_unmapped() || r.cigar.is_empty() {
            return None;
        }
        let primary = r.flag & (flags::SECONDARY | flags::SUPPLEMENTARY) == 0;
        Some(PassARecord {
            pos0: r.pos0,
            end0: r.end0(),
            flag: r.flag,
            name_hash: if primary { iq_types::name_hash(r.qname) } else { 0 },
        })
    }
    #[inline]
    pub fn is_primary(&self) -> bool {
        self.flag & (flags::SECONDARY | flags::SUPPLEMENTARY | flags::UNMAPPED) == 0
    }
}

/// Per-chromosome pass-A settings.
#[derive(Copy, Clone, Debug)]
pub struct PassAConfig {
    /// Coverage cap cutoff of this chromosome ([`coverage_cutoff`]); ≤ 0 disables the cap.
    pub coverage_cutoff: i64,
    /// Expected number of primary records (sizes the duplicate filter; 0 = unknown).
    pub expected_primary: usize,
}

impl PassAConfig {
    /// `chr_len` is the FASTA record length (4.1 `len(self.chr_record)`).
    pub fn for_chr(chr_name: &str, chr_len: u64, small_chr_max_coverage: i64, usual_gene_max_coverage: i64) -> Self {
        PassAConfig {
            coverage_cutoff: coverage_cutoff(chr_name, chr_len, small_chr_max_coverage, usual_gene_max_coverage),
            expected_primary: 0,
        }
    }
}

/// Record counts of a pass (4.1 `alignment_stat_counter`, plus skipped records).
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct PassAStats {
    pub records: u64,
    pub primary: u64,
    pub secondary: u64,
    pub supplementary: u64,
    /// Records without a reference end (unmapped/no CIGAR), not seen by the clustering.
    pub skipped: u64,
}

/// Output of pass A for one chromosome.
#[derive(Clone, Debug, Default)]
pub struct PassAResult {
    pub clusters: Vec<Cluster>,
    /// Sorted name hashes of primary records that may occur more than once (superset).
    pub dup_candidates: Vec<u64>,
    pub stats: PassAStats,
    /// Largest heap footprint of the per-cluster arrays (bytes), for reporting.
    pub peak_cluster_bytes: usize,
}

/// Streaming pass A over one chromosome. Feed records in BAM (coordinate) order with
/// [`PassA::push`], then call [`PassA::finish`].
pub struct PassA<G: GeneOracle> {
    genes: G,
    cfg: PassAConfig,
    /// `storage.region` (0-based closed), `None` when empty.
    region: Option<Iv0>,
    n: u64,
    /// Bin index of `diff[0]` / `own.hist[0]`.
    base: i64,
    /// Coverage difference array (wrapping), later prefix-summed in place.
    diff: Vec<u32>,
    own: OwnerCounter,
    /// Smallest / largest touched bin (4.1 coverage_dict min/max key).
    bin_lo: i64,
    bin_hi: i64,
    ga: Vec<GeneIdx>,
    gb: Vec<GeneIdx>,
    splits: Vec<(Iv0, u64)>,
    regs: Vec<Iv0>,
    owned: Vec<u64>,
    clusters: Vec<Cluster>,
    dups: DupFilter,
    stats: PassAStats,
    peak_cluster_bytes: usize,
}

impl<G: GeneOracle> PassA<G> {
    pub fn new(genes: G, cfg: PassAConfig) -> Self {
        let dups = if cfg.expected_primary > 0 {
            DupFilter::with_capacity(cfg.expected_primary)
        } else {
            DupFilter::default()
        };
        PassA {
            genes,
            cfg,
            region: None,
            n: 0,
            base: 0,
            diff: Vec::new(),
            own: OwnerCounter::default(),
            bin_lo: i64::MAX,
            bin_hi: i64::MIN,
            ga: Vec::new(),
            gb: Vec::new(),
            splits: Vec::new(),
            regs: Vec::new(),
            owned: Vec::new(),
            clusters: Vec::new(),
            dups,
            stats: PassAStats::default(),
            peak_cluster_bytes: 0,
        }
    }

    /// Feeds one BAM record (skips records without a reference end, see
    /// [`PassARecord::from_lite`]).
    #[inline]
    pub fn push_lite(&mut self, r: &RecordLite<'_>) {
        match PassARecord::from_lite(r) {
            Some(rec) => self.push(rec),
            None => {
                self.stats.records += 1;
                self.stats.skipped += 1;
            }
        }
    }

    /// Feeds one record (`AlignmentCollector.process` loop body, `:299-319`).
    pub fn push(&mut self, rec: PassARecord) {
        self.stats.records += 1;
        if rec.flag & flags::SECONDARY != 0 {
            self.stats.secondary += 1;
        } else if rec.flag & flags::SUPPLEMENTARY != 0 {
            self.stats.supplementary += 1;
        } else {
            self.stats.primary += 1;
        }
        let iv = Iv0 { s: rec.pos0, e: rec.end0 - 1 };
        if let Some(r) = self.region {
            // alignment_is_not_adjacent (:101-103): not overlaps(region, read)
            let adjacent = !(r.e < iv.s || r.s > iv.e);
            if !adjacent && !self.shares_gene(r, iv) {
                self.emit();
            }
        }
        self.add(iv, rec);
    }

    /// `len(get_genes_in_region(region) & get_genes_in_region(read)) != 0`. The read's
    /// (smaller) query runs first; an empty result decides without the region query.
    fn shares_gene(&mut self, region: Iv0, read: Iv0) -> bool {
        self.gb.clear();
        self.genes.genes_for_region0(read, &mut self.gb);
        if self.gb.is_empty() {
            return false;
        }
        self.ga.clear();
        self.genes.genes_for_region0(region, &mut self.ga);
        if self.ga.is_empty() {
            return false;
        }
        if self.gb.len() <= 8 {
            return self.gb.iter().any(|g| self.ga.contains(g));
        }
        self.ga.sort_unstable();
        self.gb.iter().any(|g| self.ga.binary_search(g).is_ok())
    }

    /// Makes `[lo, hi]` (bin keys) addressable in `diff` / `own.hist`.
    #[inline]
    fn ensure(&mut self, lo: i64, hi: i64) {
        if self.diff.is_empty() {
            self.base = lo;
        } else if lo < self.base {
            // records out of coordinate order (not produced by a sorted BAM): shift right
            let shift = (self.base - lo) as usize;
            self.diff.splice(0..0, std::iter::repeat_n(0, shift));
            self.own.hist.splice(0..0, std::iter::repeat_n(0, shift));
            self.base = lo;
        }
        let need = (hi - self.base + 1) as usize;
        if need > self.diff.len() {
            self.diff.resize(need, 0);
            self.own.hist.resize(need, 0);
        }
    }

    /// `AbstractAlignmentStorage.add_alignment` (`:126-136`) + ownership accounting.
    #[inline]
    fn add(&mut self, iv: Iv0, rec: PassARecord) {
        let bin_s = iv.s.div_euclid(COVERAGE_BIN);
        let bin_e = iv.e.div_euclid(COVERAGE_BIN);
        // midpoint keys lie in [bin_s - 1, bin_e]; the difference array needs bin_e + 1
        self.ensure(bin_s.min(bin_e) - 1, bin_s.max(bin_e) + 1);
        if bin_s <= bin_e {
            let a = (bin_s - self.base) as usize;
            let b = (bin_e + 1 - self.base) as usize;
            self.diff[a] = self.diff[a].wrapping_add(1);
            self.diff[b] = self.diff[b].wrapping_sub(1);
            self.bin_lo = self.bin_lo.min(bin_s);
            self.bin_hi = self.bin_hi.max(bin_e);
        }
        self.region = Some(match self.region {
            None => iv,
            Some(r) => Iv0 { s: r.s.min(iv.s), e: r.e.max(iv.e) },
        });
        self.n += 1;
        if rec.is_primary() {
            if needs_explicit(iv.s, iv.e) {
                self.own.long.push(iv);
            } else {
                let k = (mid_key(iv.s, iv.e) - self.base) as usize;
                self.own.hist[k] += 1;
            }
            self.dups.insert(rec.name_hash);
        }
    }

    /// `forward_alignments` (`:324-335`) + `process_alignments_in_region` bookkeeping.
    fn emit(&mut self) {
        let Some(region) = self.region.take() else { return };
        let bytes = self.diff.capacity() * 4 + self.own.hist.capacity() * 4 + self.own.long.capacity() * 16;
        self.peak_cluster_bytes = self.peak_cluster_bytes.max(bytes);

        // coverage = prefix sum of the difference array, in place
        let mut acc = 0u32;
        for d in self.diff.iter_mut() {
            acc = acc.wrapping_add(*d);
            *d = acc;
        }
        let counts: &[u32] = if self.bin_lo <= self.bin_hi {
            &self.diff[(self.bin_lo - self.base) as usize..=(self.bin_hi - self.base) as usize]
        } else {
            &[]
        };
        split_coverage_regions_into(region, &CoverageBins { first_bin: self.bin_lo, counts }, self.n, &mut self.splits);

        self.regs.clear();
        self.regs.extend(self.splits.iter().map(|x| x.0));
        self.owned.clear();
        self.owned.resize(self.regs.len(), 0);
        self.own.assign(self.base, &self.regs, &mut self.owned);
        // 4.1 double processing: a block also processes the reads straddling its boundaries.
        // Each such read touches the bin holding the boundary, so that bin's coverage bounds
        // their number (all records counted; an upper bound, used only as a cost estimate).
        if self.regs.len() > 1 {
            let (diff, base, lo, hi) = (&self.diff, self.base, self.bin_lo, self.bin_hi);
            let cov_at = |pos: i64| -> u64 {
                let b = pos.div_euclid(COVERAGE_BIN);
                if b < lo || b > hi {
                    0
                } else {
                    diff[(b - base) as usize] as u64
                }
            };
            for k in 0..self.regs.len() - 1 {
                let c = cov_at(self.regs[k].e).max(cov_at(self.regs[k + 1].s));
                self.owned[k] += c;
                self.owned[k + 1] += c;
            }
        }

        let mut blocks = Vec::with_capacity(self.splits.len());
        for (i, &(reg, max_cov)) in self.splits.iter().enumerate() {
            let mut gene_ids = Vec::new();
            self.genes.genes_for_region0(reg, &mut gene_ids);
            blocks.push(Block {
                region: reg,
                max_cov: max_cov.min(u32::MAX as u64) as u32,
                skip: skip_read_fraction(max_cov, self.cfg.coverage_cutoff),
                gene_ids,
                owned_reads: self.owned[i],
            });
        }
        self.clusters.push(Cluster { region, blocks, n_records: self.n });

        // storage.reset()
        self.n = 0;
        self.diff.clear();
        self.own.clear();
        self.bin_lo = i64::MAX;
        self.bin_hi = i64::MIN;
    }

    /// Clusters emitted so far (complete ones only).
    pub fn clusters(&self) -> &[Cluster] {
        &self.clusters
    }

    /// Flushes the last cluster (`:321-323`).
    pub fn finish(mut self) -> PassAResult {
        self.emit();
        let peak = self.peak_cluster_bytes + self.dups.heap_bytes();
        PassAResult {
            clusters: std::mem::take(&mut self.clusters),
            dup_candidates: self.dups.finish(),
            stats: self.stats,
            peak_cluster_bytes: peak,
        }
    }
}

/// `fn pass_a(records, genes, cfg) -> PassAResult` (U3).
pub fn pass_a<G: GeneOracle>(
    records: impl IntoIterator<Item = PassARecord>,
    genes: G,
    cfg: PassAConfig,
) -> PassAResult {
    let mut p = PassA::new(genes, cfg);
    for r in records {
        p.push(r);
    }
    p.finish()
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Genes(Vec<(i64, i64)>); // 1-based closed gene spans, GeneIdx = position
    impl GeneOracle for Genes {
        fn genes_for_region0(&self, r: Iv0, out: &mut Vec<GeneIdx>) {
            for (i, &(s, e)) in self.0.iter().enumerate() {
                let (qs, qe) = (r.s + 1, r.e + 1); // 1-based query
                if s <= qe && e >= qs {
                    out.push(i as GeneIdx);
                }
            }
        }
    }
    fn rec(s: i64, e: i64) -> PassARecord {
        PassARecord { pos0: s, end0: e + 1, flag: 0, name_hash: (s as u64) << 20 ^ e as u64 }
    }
    const CFG: PassAConfig = PassAConfig { coverage_cutoff: -1, expected_primary: 0 };

    #[test]
    fn gap_splits_without_genes() {
        let r = pass_a([rec(100, 200), rec(150, 300), rec(400, 500)], NoGenes, CFG);
        assert_eq!(r.clusters.len(), 2);
        assert_eq!(r.clusters[0].region, Iv0 { s: 100, e: 300 });
        assert_eq!(r.clusters[0].n_records, 2);
        assert_eq!(r.clusters[1].blocks[0].owned_reads, 1);
        assert!(r.dup_candidates.is_empty());
    }

    #[test]
    fn shared_gene_bridges_gap() {
        // gene 1-based [250, 450] touches both sides
        let g = Genes(vec![(250, 450)]);
        let r = pass_a([rec(100, 300), rec(400, 500), rec(10_000, 10_100)], &g, CFG);
        assert_eq!(r.clusters.len(), 2);
        assert_eq!(r.clusters[0].region, Iv0 { s: 100, e: 500 });
        assert_eq!(r.clusters[0].blocks[0].gene_ids, vec![0]);
        // a gene ending exactly at the 0-based region end + 1 (1-based e) still bridges
        let g = Genes(vec![(301, 401)]);
        assert_eq!(pass_a([rec(100, 300), rec(400, 500)], &g, CFG).clusters.len(), 1);
        let g = Genes(vec![(302, 401)]);
        assert_eq!(pass_a([rec(100, 300), rec(400, 500)], &g, CFG).clusters.len(), 2);
    }

    #[test]
    fn coverage_matches_python_dict() {
        // bins: [0..3] by first read, [2..5] by second; no split (few reads)
        let r = pass_a([rec(10, 3 * 256 + 5), rec(2 * 256, 5 * 256)], NoGenes, CFG);
        assert_eq!(r.clusters.len(), 1);
        assert_eq!(r.clusters[0].blocks[0].max_cov, 2);
    }

    #[test]
    fn duplicates_and_secondaries() {
        let mut a = rec(100, 200);
        a.name_hash = 42;
        let mut b = rec(150, 250);
        b.name_hash = 42;
        let mut c = rec(160, 260);
        c.name_hash = 7;
        c.flag = flags::SECONDARY;
        let mut d = rec(170, 270);
        d.name_hash = 7;
        let r = pass_a([a, b, c, d], NoGenes, CFG);
        assert_eq!(r.dup_candidates, vec![42]);
        assert_eq!(r.clusters[0].blocks[0].owned_reads, 3);
        assert_eq!(r.stats.secondary, 1);
    }

    #[test]
    fn unsorted_input_grows_left() {
        let r = pass_a([rec(5000, 6000), rec(100, 5500)], NoGenes, CFG);
        assert_eq!(r.clusters.len(), 1);
        assert_eq!(r.clusters[0].region, Iv0 { s: 100, e: 6000 });
        assert_eq!(r.clusters[0].blocks[0].max_cov, 2);
        assert_eq!(r.clusters[0].blocks[0].owned_reads, 2);
    }
}
