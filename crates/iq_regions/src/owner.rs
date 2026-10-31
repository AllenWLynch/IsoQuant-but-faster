//! Read ownership (design §5.2).
//!
//! 4.1 re-fetches every split region, so a record overlapping two blocks is processed in
//! both, and iqrs does the same (U8 decides with the fetch window, see
//! `iq_collect::ClusterOwner`). The single **owner block** (largest overlap with
//! `[pos0, end0 - 1]`, ties to the earlier block) remains as a partition of the reads for
//! cost estimates: pass A counts owned primary records per block without storing reads (see
//! [`OwnerCounter`]) and adds the boundary-straddling reads separately.

use iq_types::read::COVERAGE_BIN;
use iq_types::{Block, Iv0};

use crate::split::MIN_REGION_LEN;

/// Owner of `read` among `regions` (sorted, pairwise disjoint): the region with the largest
/// overlap with the 0-based closed read span; ties go to the earlier region.
///
/// A record that overlaps no region (only possible for degenerate spans, e.g. zero reference
/// length, or the first base of a cluster that starts exactly on a bin boundary, which 4.1's
/// `bin*256 + 1` region start leaves out) goes to the first region ending at or after its
/// start, clamped to the last region. Panics if `regions` is empty.
pub fn owner_of(read: Iv0, regions: &[Iv0]) -> usize {
    owner_by(read, regions.len(), |i| regions[i])
}

/// [`owner_of`] over blocks (design `owner_block(read_iv0, &[Block])`).
pub fn owner_block(read: Iv0, blocks: &[Block]) -> usize {
    owner_by(read, blocks.len(), |i| blocks[i].region)
}

#[inline]
fn owner_by(read: Iv0, n: usize, reg: impl Fn(usize) -> Iv0) -> usize {
    assert!(n > 0, "owner_block: no blocks");
    // first region ending at or after read.s
    let (mut lo, mut hi) = (0usize, n);
    while lo < hi {
        let mid = (lo + hi) / 2;
        if reg(mid).e < read.s {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    let first = lo;
    let mut best = first.min(n - 1);
    let mut best_ov = 0i64;
    let mut j = first;
    while j < n {
        let r = reg(j);
        if r.s > read.e {
            break;
        }
        let ov = read.e.min(r.e) - read.s.max(r.s) + 1;
        if ov > best_ov {
            best_ov = ov;
            best = j;
        }
        j += 1;
    }
    best
}

/// Blocks a record is fetched into under 4.1 semantics: every block whose region overlaps
/// the read span.
pub fn overlapping_blocks(read: Iv0, blocks: &[Block]) -> std::ops::Range<usize> {
    let lo = blocks.partition_point(|b| b.region.e < read.s);
    let hi = lo + blocks[lo..].partition_point(|b| b.region.s <= read.e);
    lo..hi
}

/// Reads with a span longer than this are resolved explicitly; shorter ones by midpoint.
///
/// Why the midpoint is exact for short reads: split_coverage_regions' blocks are contiguous
/// and every *inner* block spans `[c*256 + 1, c'*256]` with `c' - c >= MIN_REGION_LEN / 256`,
/// i.e. ≥ 32768 bp. A span of ≤ 32768 bp therefore crosses at most one boundary `p` (block
/// `k` ends at `p`, `k+1` starts at `p+1`). Then overlap_k = p - s + 1 ≥ overlap_{k+1} = e - p
/// iff `floor((s+e)/2) <= p`: the owner is the block containing the midpoint, ties included.
pub const LONG_SPAN: i64 = MIN_REGION_LEN;

/// Reads resolved explicitly by [`owner_of`]: long spans, and degenerate (empty) spans whose
/// "no overlap" fallback is not midpoint-shaped.
#[inline]
pub(crate) fn needs_explicit(s: i64, e: i64) -> bool {
    e - s + 1 > LONG_SPAN || e < s
}

/// Midpoint key: positions `[c*256 + 1, c*256 + 256]` map to `c`, so the inner block
/// `[c_k*256 + 1, c_{k+1}*256]` is exactly the key range `c_k .. c_{k+1}`.
#[inline]
pub(crate) fn mid_key(s: i64, e: i64) -> i64 {
    let m = (s + e).div_euclid(2);
    (m - 1).div_euclid(COVERAGE_BIN)
}

/// Counts owned records per block without per-read storage, except for the rare reads
/// longer than [`LONG_SPAN`]. Histogram indices are relative to a caller-chosen base key.
#[derive(Default, Debug)]
pub struct OwnerCounter {
    /// `hist[i]` = short reads whose midpoint key is `base + i`.
    pub(crate) hist: Vec<u32>,
    pub(crate) long: Vec<Iv0>,
}

impl OwnerCounter {
    pub fn clear(&mut self) {
        self.hist.clear();
        self.long.clear();
    }

    /// Adds the owned counts to `owned[k]` for blocks `regions` (sorted, contiguous as
    /// produced by `split_coverage_regions`), with `hist[0]` at key `base`.
    pub fn assign(&self, base: i64, regions: &[Iv0], owned: &mut [u64]) {
        let n = regions.len();
        debug_assert_eq!(owned.len(), n);
        let mut k = 0usize;
        // last midpoint key of block k: keys up to (e-1)/256 have midpoints ≤ e when e = c*256
        let last_key = |k: usize| (regions[k].e - 1).div_euclid(COVERAGE_BIN);
        for (i, &c) in self.hist.iter().enumerate() {
            if c == 0 {
                continue;
            }
            let key = base + i as i64;
            while k + 1 < n && key > last_key(k) {
                k += 1;
            }
            owned[k] += c as u64;
        }
        for &r in &self.long {
            owned[owner_of(r, regions)] += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::split::{split_coverage_regions, CoverageBins};
    use proptest::prelude::*;

    fn regs(v: &[(i64, i64)]) -> Vec<Iv0> {
        v.iter().map(|&(s, e)| Iv0 { s, e }).collect()
    }

    #[test]
    fn basic_owner() {
        let r = regs(&[(0, 99), (100, 199), (200, 299)]);
        assert_eq!(owner_of(Iv0 { s: 10, e: 20 }, &r), 0);
        assert_eq!(owner_of(Iv0 { s: 90, e: 120 }, &r), 1);
        assert_eq!(owner_of(Iv0 { s: 80, e: 119 }, &r), 0); // tie 20/20 → earlier
        assert_eq!(owner_of(Iv0 { s: 50, e: 250 }, &r), 1); // full middle block wins
        assert_eq!(owner_of(Iv0 { s: 500, e: 600 }, &r), 2); // no overlap → clamp
        assert_eq!(owner_of(Iv0 { s: 100, e: 99 }, &r), 1); // zero length
    }

    /// Builds a random cluster the way pass A does and checks, for every read, that its owner
    /// overlaps it (when any block does) and is unique/deterministic, and that the
    /// histogram counts equal brute-force owner counts.
    fn check_cluster(reads: &[(i64, i64)]) -> Result<usize, TestCaseError> {
        let s0 = reads.iter().map(|r| r.0).min().unwrap();
        let e0 = reads.iter().map(|r| r.1).max().unwrap();
        let region = Iv0 { s: s0, e: e0 };
        let b_lo = s0.div_euclid(256);
        let b_hi = e0.div_euclid(256);
        let mut cov = vec![0u32; (b_hi - b_lo + 1) as usize];
        let base = b_lo - 1;
        let mut oc = OwnerCounter { hist: vec![0; (b_hi - b_lo + 3) as usize], long: vec![] };
        for &(s, e) in reads {
            for b in s.div_euclid(256)..=e.div_euclid(256) {
                cov[(b - b_lo) as usize] += 1;
            }
            if needs_explicit(s, e) {
                oc.long.push(Iv0 { s, e });
            } else {
                oc.hist[(mid_key(s, e) - base) as usize] += 1;
            }
        }
        let splits =
            split_coverage_regions(region, &CoverageBins { first_bin: b_lo, counts: &cov }, reads.len() as u64);
        let rg: Vec<Iv0> = splits.iter().map(|x| x.0).collect();
        // contiguity of the regions that makes the midpoint rule exact
        for w in rg.windows(2) {
            prop_assert_eq!(w[0].e + 1, w[1].s);
        }
        let mut brute = vec![0u64; rg.len()];
        for &(s, e) in reads {
            let iv = Iv0 { s, e };
            let o = owner_of(iv, &rg);
            prop_assert_eq!(o, owner_of(iv, &rg));
            let ovs: Vec<i64> = rg.iter().map(|r| e.min(r.e) - s.max(r.s) + 1).collect();
            let mx = *ovs.iter().max().unwrap();
            if mx > 0 {
                prop_assert_eq!(ovs[o], mx);
                prop_assert_eq!(ovs.iter().position(|&x| x == mx).unwrap(), o);
            }
            brute[o] += 1;
        }
        let mut got = vec![0u64; rg.len()];
        oc.assign(base, &rg, &mut got);
        prop_assert_eq!(got, brute);
        Ok(rg.len())
    }

    fn cluster_strategy() -> impl Strategy<Value = Vec<(i64, i64)>> {
        // dense blocks of reads separated by low-coverage stretches, plus some long reads
        (1usize..5, any::<u64>()).prop_map(|(n_islands, seed)| {
            let mut x = seed | 1;
            let mut rnd = move |m: i64| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                (x % m as u64) as i64
            };
            let mut reads = Vec::new();
            let mut at = 1000 + rnd(512);
            for _ in 0..n_islands {
                let len = 40_000 + rnd(60_000);
                for _ in 0..(1100 / n_islands + 50) {
                    let s = at + rnd(len);
                    let l = match rnd(20) {
                        0 => 30_000 + rnd(80_000),
                        1 => rnd(3),
                        _ => 200 + rnd(5000),
                    };
                    reads.push((s, s + l - 1));
                }
                // sparse bridge so the cluster stays one region with valleys
                let gap = 2_000 + rnd(30_000);
                for k in 0..rnd(3) {
                    let s = at + len + k * 700;
                    reads.push((s, s + 5000 + gap));
                }
                at += len + gap;
            }
            // bin-aligned starts exercise the `bin*256 + 1` hole
            reads.push((at.div_euclid(256) * 256, at.div_euclid(256) * 256 + 300));
            reads
        })
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(200))]
        #[test]
        fn ownership_unique_overlapping_and_counts_exact(reads in cluster_strategy()) {
            check_cluster(&reads)?;
        }
    }

    /// The generator must actually produce split clusters, or the property is vacuous.
    #[test]
    fn generator_produces_multi_block_clusters() {
        use proptest::strategy::ValueTree;
        let mut runner = proptest::test_runner::TestRunner::deterministic();
        let mut multi = 0;
        for _ in 0..100 {
            let reads = cluster_strategy().new_tree(&mut runner).unwrap().current();
            if check_cluster(&reads).unwrap() > 1 {
                multi += 1;
            }
        }
        assert!(multi >= 30, "only {multi}/100 clusters split");
    }
}
