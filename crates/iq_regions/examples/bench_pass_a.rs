//! Pass-A micro-benchmark over synthetic coordinate-sorted records.
//!
//! `cargo run --release -p iq_regions --example bench_pass_a [n_reads]` (default 10M).
//! Records are generated on the fly (not stored), so the peak heap is pass A's own: 2,000
//! loci of 20-120 kb on a ~360 Mb span, plus one vtg-like locus with 30% of all reads,
//! 1-4 kb spans, 3% spliced spans of 20-80 kb, 10% secondary/supplementary, 0.5% repeated
//! names. Each primary record's 36-byte name is hashed (xxh3) like the real path.

use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use iq_regions::{GeneOracle, PassA, PassAConfig, PassARecord};
use iq_types::{GeneIdx, Iv0};

struct Counting;
static CUR: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);
unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        let p = unsafe { System.alloc(l) };
        if !p.is_null() {
            let c = CUR.fetch_add(l.size(), Ordering::Relaxed) + l.size();
            PEAK.fetch_max(c, Ordering::Relaxed);
        }
        p
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        unsafe { System.dealloc(p, l) };
        CUR.fetch_sub(l.size(), Ordering::Relaxed);
    }
}
#[global_allocator]
static A: Counting = Counting;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, m: u64) -> u64 {
        self.next() % m
    }
}

/// Sorted 1-based gene spans; binary search on start plus a max-gene-length window.
struct Genes {
    iv: Vec<(i64, i64)>,
    max_len: i64,
}
impl GeneOracle for Genes {
    fn genes_for_region0(&self, r: Iv0, out: &mut Vec<GeneIdx>) {
        let (qs, qe) = (r.s + 1, r.e + 1);
        let lo = self.iv.partition_point(|g| g.0 < qs - self.max_len);
        let hi = self.iv.partition_point(|g| g.0 <= qe);
        for (i, g) in self.iv[lo..hi].iter().enumerate() {
            if g.1 >= qs {
                out.push((lo + i) as GeneIdx);
            }
        }
    }
}

fn main() {
    let n: u64 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(10_000_000);
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    // loci: (start, len, n_reads)
    let n_loci = 2000u64;
    let giant = n * 3 / 10;
    let per = (n - giant) / n_loci;
    let mut loci = Vec::new();
    let mut at = 10_000i64;
    for k in 0..n_loci {
        let len = 20_000 + rng.below(100_000) as i64;
        let reads = if k == n_loci / 2 { giant } else { per / 2 + rng.below(per) };
        loci.push((at, len, reads));
        at += len + 100_000 + rng.below(100_000) as i64;
    }
    let genes = Genes { iv: loci.iter().map(|&(s, l, _)| (s + 1, s + l)).collect(), max_len: 120_000 };

    let base = CUR.load(Ordering::Relaxed);
    PEAK.store(base, Ordering::Relaxed);
    let t = Instant::now();
    let mut p = PassA::new(&genes, PassAConfig { coverage_cutoff: -1, expected_primary: 0 });
    let mut name = *b"m64012_190920_173625/000000000/ccs\0\0";
    let mut total = 0u64;
    let mut serial = 0u64;
    for &(s, len, reads) in &loci {
        // sorted starts: cumulative uniform increments over the locus
        let step = (len as u64 * 2 / reads.max(1)).max(1);
        let mut pos = s;
        for _ in 0..reads {
            pos += rng.below(step + 1) as i64;
            let span = if rng.below(33) == 0 { 20_000 + rng.below(60_000) } else { 1_000 + rng.below(3_000) } as i64;
            let flag = match rng.below(20) {
                0 => 0x100,
                1 => 0x800,
                _ => 0,
            };
            serial += 1;
            let id = if rng.below(200) == 0 { serial / 2 } else { serial };
            name[21..30].copy_from_slice(format_args_9(id).as_slice());
            let name_hash = if flag == 0 { iq_types::name_hash(&name) } else { 0 };
            p.push(PassARecord { pos0: pos, end0: pos + span, flag, name_hash });
            total += 1;
        }
    }
    let res = p.finish();
    let dt = t.elapsed().as_secs_f64();
    let peak = PEAK.load(Ordering::Relaxed) - base;
    let n_blocks: usize = res.clusters.iter().map(|c| c.blocks.len()).sum();
    let owned: u64 = res.clusters.iter().flat_map(|c| &c.blocks).map(|b| b.owned_reads).sum();
    println!("records        {total}");
    println!("clusters       {}  blocks {n_blocks}", res.clusters.len());
    println!("owned primary  {owned} (stats.primary {})", res.stats.primary);
    println!("dup candidates {}", res.dup_candidates.len());
    println!("time           {dt:.3} s");
    println!("throughput     {:.1} M records/s", total as f64 / dt / 1e6);
    println!("peak heap      {:.2} MB ({:.2} B/record)", peak as f64 / 1e6, peak as f64 / total as f64);
}

/// 9 decimal digits without allocating.
fn format_args_9(mut v: u64) -> [u8; 9] {
    let mut b = [b'0'; 9];
    for i in (0..9).rev() {
        b[i] = b'0' + (v % 10) as u8;
        v /= 10;
    }
    b
}
