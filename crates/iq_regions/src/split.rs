//! `split_coverage_regions` (`alignment_processor.py:640-682`) and the coverage cap of
//! `process_alignments_in_region` (`:337-375`), ported verbatim.

use iq_types::read::{ABS_COV_VALLEY, COVERAGE_BIN, REL_COV_VALLEY, SMALL_CHR_IDS, SMALL_CHR_LEN};
use iq_types::Iv0;

/// `AlignmentCollector.MIN_REGION_LEN`.
pub const MIN_REGION_LEN: i64 = 32768;
/// `AlignmentCollector.MAX_REGION_LEN`.
pub const MAX_REGION_LEN: i64 = 1_048_576;
/// `AlignmentCollector.MIN_READS_TO_SPLIT`.
pub const MIN_READS_TO_SPLIT: u64 = 1024;

/// Dense view of a cluster's 256-bp coverage bins (4.1 `coverage_dict`).
///
/// `counts[i]` is the coverage of bin `first_bin + i`. `first_bin` and the last bin must be
/// the smallest and largest *touched* bins (4.1 `sorted(coverage_dict.keys())[0]` and
/// `[-1]`); bins in between that no record touched are 0, exactly like a defaultdict miss.
#[derive(Clone, Copy, Debug)]
pub struct CoverageBins<'a> {
    pub first_bin: i64,
    pub counts: &'a [u32],
}

impl CoverageBins<'_> {
    #[inline]
    fn get(&self, bin: i64) -> u64 {
        let i = bin - self.first_bin;
        if i >= 0 && (i as usize) < self.counts.len() {
            self.counts[i as usize] as u64
        } else {
            0
        }
    }
}

/// `split_coverage_regions(genomic_region, alignment_storage)`: splits a cluster at coverage
/// valleys. Returns `(region, max_coverage)` per block; regions are 4.1's mixed coordinates
/// (`max(bin*256 + 1, region.s)` .. `min(pos*256, region.e)`), recorded verbatim.
///
/// 4.1 crashes on an empty coverage dict (every record had zero reference length on a bin
/// boundary); iqrs returns the whole region with coverage 0.
pub fn split_coverage_regions(region: Iv0, bins: &CoverageBins<'_>, n_reads: u64) -> Vec<(Iv0, u64)> {
    let mut out = Vec::new();
    split_coverage_regions_into(region, bins, n_reads, &mut out);
    out
}

/// [`split_coverage_regions`] into a reused buffer (cleared first).
pub fn split_coverage_regions_into(region: Iv0, bins: &CoverageBins<'_>, n_reads: u64, out: &mut Vec<(Iv0, u64)>) {
    out.clear();
    if bins.counts.is_empty() {
        out.push((region, 0));
        return;
    }
    if region.e - region.s + 1 < MIN_REGION_LEN || n_reads < MIN_READS_TO_SPLIT {
        let max_cov = bins.counts.iter().copied().max().unwrap_or(0) as u64;
        out.push((region, max_cov));
        return;
    }
    let last = bins.first_bin + bins.counts.len() as i64 - 1;
    let min_bins = MIN_REGION_LEN / COVERAGE_BIN;
    let max_bins = MAX_REGION_LEN / COVERAGE_BIN;

    let mut current_start = bins.first_bin;
    let mut pos = current_start + 1;
    let mut max_cov = bins.get(current_start);
    let mut abs_cov_valley = ABS_COV_VALLEY as u64;
    let mut rel_cov_valley = REL_COV_VALLEY;
    let mut prev_upd_pos = pos;

    while pos <= last {
        loop {
            if pos > last {
                break;
            }
            if pos - current_start >= min_bins {
                // `coverage_dict[pos] > max(abs_cov_valley, max_cov * rel_cov_valley)`; Python
                // compares int with float exactly, and every value here is < 2^53.
                let c = bins.get(pos);
                let is_valley = c <= abs_cov_valley || (c as f64) <= (max_cov as f64) * rel_cov_valley;
                if is_valley {
                    break;
                }
            }
            if pos - prev_upd_pos > max_bins {
                // increase valley cut-offs every MAX_REGION_LEN bases to avoid extra long regions
                prev_upd_pos = pos;
                abs_cov_valley += ABS_COV_VALLEY as u64;
                rel_cov_valley += REL_COV_VALLEY;
            }
            max_cov = max_cov.max(bins.get(pos));
            pos += 1;
        }
        out.push((
            Iv0 { s: (current_start * COVERAGE_BIN + 1).max(region.s), e: (pos * COVERAGE_BIN).min(region.e) },
            max_cov,
        ));
        current_start = pos;
        prev_upd_pos = pos;
        abs_cov_valley = ABS_COV_VALLEY as u64;
        rel_cov_valley = REL_COV_VALLEY;
        max_cov = bins.get(current_start);
        pos = (current_start + 1).min(last + 1);
    }
}

/// 4.1 small-chromosome test (`process_alignments_in_region`, `:340-341`): `chr_len` is the
/// FASTA record length.
pub fn is_small_chr(chr_name: &str, chr_len: u64) -> bool {
    SMALL_CHR_IDS.contains(&chr_name) || chr_len < SMALL_CHR_LEN
}

/// Coverage cutoff for a chromosome (`--max_coverage_small_chr` / `--max_coverage_normal_chr`).
pub fn coverage_cutoff(chr_name: &str, chr_len: u64, small_chr_max_coverage: i64, usual_gene_max_coverage: i64) -> i64 {
    if is_small_chr(chr_name, chr_len) {
        small_chr_max_coverage
    } else {
        usual_gene_max_coverage
    }
}

/// `skip_read_fraction`: `ceil(max_coverage / cutoff)` if `max_coverage > cutoff > 0`, else 1.
/// (4.1 divides in floats; for operands < 2^31 the float quotient of a non-integral ratio
/// never rounds to an integer, so integer ceil is identical.)
pub fn skip_read_fraction(max_coverage: u64, cutoff: i64) -> u32 {
    if cutoff > 0 && max_coverage > cutoff as u64 {
        max_coverage.div_ceil(cutoff as u64).min(u32::MAX as u64) as u32
    } else {
        1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bins_from(first: i64, v: &[u32]) -> Vec<(Iv0, u64)> {
        let r = Iv0 { s: first * 256, e: (first + v.len() as i64) * 256 - 1 };
        split_coverage_regions(r, &CoverageBins { first_bin: first, counts: v }, 5000)
    }

    #[test]
    fn short_or_few_reads_not_split() {
        let v = vec![5u32; 100];
        let r = Iv0 { s: 0, e: 100 * 256 - 1 };
        assert_eq!(split_coverage_regions(r, &CoverageBins { first_bin: 0, counts: &v }, 5000), vec![(r, 5)]);
        let v = vec![5u32; 300];
        let r = Iv0 { s: 0, e: 300 * 256 - 1 };
        assert_eq!(split_coverage_regions(r, &CoverageBins { first_bin: 0, counts: &v }, 1023), vec![(r, 5)]);
    }

    #[test]
    fn valley_splits() {
        // 200 bins of 100, one bin of 1 (valley), 200 bins of 50
        let mut v = vec![100u32; 200];
        v.push(1);
        v.extend(std::iter::repeat_n(50, 200));
        let out = bins_from(10, &v);
        assert_eq!(out.len(), 2);
        assert_eq!(out[0], (Iv0 { s: 10 * 256 + 1, e: 210 * 256 }, 100)); // start = bin*256+1 (4.1 mixed coords)
                                                                          // the valley bin starts the next region; its own coverage seeds max_cov
        assert_eq!(out[1], (Iv0 { s: 210 * 256 + 1, e: 411 * 256 - 1 }, 50));
    }

    #[test]
    fn valley_inside_min_bins_ignored() {
        let mut v = vec![100u32; 50];
        v.push(0);
        v.extend(std::iter::repeat_n(100, 200));
        assert_eq!(bins_from(0, &v).len(), 1);
    }

    #[test]
    fn empty_coverage_does_not_panic() {
        let r = Iv0 { s: 256, e: 255 };
        assert_eq!(split_coverage_regions(r, &CoverageBins { first_bin: 0, counts: &[] }, 1), vec![(r, 0)]);
    }

    #[test]
    fn skip_fraction() {
        assert_eq!(skip_read_fraction(100, -1), 1);
        assert_eq!(skip_read_fraction(100, 100), 1);
        assert_eq!(skip_read_fraction(101, 100), 2);
        assert_eq!(skip_read_fraction(300, 100), 3);
        assert!(is_small_chr("chrM", 10_000_000));
        assert!(is_small_chr("8", 499_999));
        assert!(!is_small_chr("8", 500_000));
    }
}
