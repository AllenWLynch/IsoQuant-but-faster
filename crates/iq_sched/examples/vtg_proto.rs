//! M1 vtg memory prototype. Measurement tool, not production code.
//!
//! ```text
//! vtg_proto passa <bam> <chr> [threads]   # pass A only: time, rec/s, peak RSS, clusters/blocks
//! vtg_proto store <bam> <chr> [threads]   # pass A, then stream the chromosome once more:
//!                                         #   record -> process_alignment_into -> ReadRecord (stub
//!                                         #   Noninformative assignment) -> owner block's BlockStore
//! ```
//! Each block's store is built incrementally (one reused `ReadRecord`) and dropped as soon as
//! the stream passes the block's end, so only open blocks are resident. RSS is sampled from
//! `/proc/self/statm` every 100 ms by a background thread.
//!
//! Stubs (U1/U6 not integrated): no gene DB (`NoGenes`: gaps always split clusters, all blocks
//! intergenic), `corrected_exons = exons`, `strand = mapped_strand`, empty assignment (real
//! assignments add ~40 B/read of matches/events).

use iq_bam::{fetch_with, BamHandle, FetchOpts};
use iq_read::polya::{ExternalPolyA, PolyAParams};
use iq_regions::{owner_block, PassA, PassAConfig};
use iq_types::{
    AlignedRead, Assignment, AssignmentType, Block, BlockStore, Iv, Iv0, Params, ReadRecord, Strand, NO_POLYA,
};
use std::path::Path;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Arc;
use std::time::{Duration, Instant};

// ---------------------------------------------------------------------------------------------
// RSS sampler
// ---------------------------------------------------------------------------------------------

fn rss_bytes() -> u64 {
    let s = std::fs::read_to_string("/proc/self/statm").unwrap_or_default();
    let pages: u64 = s.split_whitespace().nth(1).and_then(|x| x.parse().ok()).unwrap_or(0);
    pages * 4096
}

fn vm_hwm() -> u64 {
    let s = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
    s.lines()
        .find(|l| l.starts_with("VmHWM:"))
        .and_then(|l| l.split_whitespace().nth(1))
        .and_then(|x| x.parse::<u64>().ok())
        .map_or(0, |kb| kb * 1024)
}

struct Sampler {
    peak: Arc<AtomicU64>,
    stop: Arc<AtomicBool>,
    h: Option<std::thread::JoinHandle<()>>,
}

impl Sampler {
    fn start() -> Self {
        let peak = Arc::new(AtomicU64::new(rss_bytes()));
        let stop = Arc::new(AtomicBool::new(false));
        let (p, s) = (peak.clone(), stop.clone());
        let h = std::thread::spawn(move || {
            while !s.load(Ordering::Relaxed) {
                p.fetch_max(rss_bytes(), Ordering::Relaxed);
                std::thread::sleep(Duration::from_millis(100));
            }
        });
        Sampler { peak, stop, h: Some(h) }
    }
    fn peak(&self) -> u64 {
        self.peak.fetch_max(rss_bytes(), Ordering::Relaxed).max(rss_bytes())
    }
    /// Resets the peak to the current RSS (start of a new phase).
    fn reset(&self) {
        self.peak.store(rss_bytes(), Ordering::Relaxed);
    }
}

impl Drop for Sampler {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::Relaxed);
        if let Some(h) = self.h.take() {
            let _ = h.join();
        }
    }
}

fn mb(b: u64) -> f64 {
    b as f64 / 1e6
}

// ---------------------------------------------------------------------------------------------
// Pass A
// ---------------------------------------------------------------------------------------------

fn run_pass_a(bam: &BamHandle, chr: &str, smp: &Sampler) -> Vec<Block> {
    let tid = bam
        .ref_id(chr)
        .unwrap_or_else(|| panic!("unknown chr {chr}; have {:?}", &bam.references()[..bam.references().len().min(30)]));
    let len = bam.reference_lengths()[tid];
    let cfg = PassAConfig::for_chr(chr, len, 1_000_000, -1); // 4.1 defaults
    let rss0 = rss_bytes();
    smp.reset();
    let t = Instant::now();
    let mut pa = PassA::new(iq_regions::NoGenes, cfg);
    let opts = FetchOpts { seq: false, batch_records: 4096, qual: false };
    for b in fetch_with(bam, chr, Iv0 { s: 0, e: len as i64 - 1 }, opts) {
        for r in b.iter() {
            pa.push_lite(&r);
        }
    }
    let res = pa.finish();
    let dt = t.elapsed().as_secs_f64();
    let st = res.stats;
    println!("== pass A {chr}");
    println!(
        "records {} (primary {}, secondary {}, supplementary {}, skipped {})",
        st.records, st.primary, st.secondary, st.supplementary, st.skipped
    );
    println!("wall {dt:.2} s, {:.0} rec/s", st.records as f64 / dt);
    println!(
        "RSS before {:.1} MB, peak {:.1} MB (delta {:.1} MB); pass-A internal peak_cluster_bytes {:.1} MB; dup candidates {}",
        mb(rss0),
        mb(smp.peak()),
        mb(smp.peak().saturating_sub(rss0)),
        mb(res.peak_cluster_bytes as u64),
        res.dup_candidates.len()
    );
    let nblocks: usize = res.clusters.iter().map(|c| c.blocks.len()).sum();
    println!("clusters {}, blocks {}", res.clusters.len(), nblocks);
    let mut big: Vec<&iq_types::Cluster> = res.clusters.iter().collect();
    big.sort_by_key(|c| std::cmp::Reverse(c.n_records));
    println!("largest clusters (by records):");
    for c in big.iter().take(5) {
        println!("  cluster {}-{} n_records {} blocks {}", c.region.s, c.region.e, c.n_records, c.blocks.len());
        for b in &c.blocks {
            println!(
                "    block {}-{} len {} max_cov {} owned_primary {} skip {}",
                b.region.s,
                b.region.e,
                b.region.e - b.region.s + 1,
                b.max_cov,
                b.owned_reads,
                b.skip
            );
        }
    }
    let mut owned: Vec<u64> = res.clusters.iter().flat_map(|c| c.blocks.iter().map(|b| b.owned_reads)).collect();
    owned.sort_unstable();
    let tot: u64 = owned.iter().sum();
    if !owned.is_empty() {
        println!(
            "block owned_primary: total {tot}, median {}, p99 {}, max {} ({:.1}% of total)",
            owned[owned.len() / 2],
            owned[owned.len() * 99 / 100],
            owned[owned.len() - 1],
            100.0 * owned[owned.len() - 1] as f64 / tot.max(1) as f64
        );
    }
    res.clusters.into_iter().flat_map(|c| c.blocks).collect()
}

// ---------------------------------------------------------------------------------------------
// Store building
// ---------------------------------------------------------------------------------------------

#[derive(Default, Clone)]
struct BlockStats {
    region: (i64, i64),
    reads: u64,
    introns: u64,    // Σ introns over reads
    chains: u64,     // distinct chains (incl. empty)
    chain_ivs: u64,  // Σ introns over distinct chains
    singletons: u64, // chains used by exactly one read
    mono: u64,       // mono-exon reads
    names: u64,      // Σ qname bytes
    dels: u64,       // boundary deletions stored
    heap: u64,       // BlockStore::heap_bytes() (capacity based, includes Vec slack + names)
    tight_hot: u64,  // len-based hot columns (no names, no chain table)
    chain_tbl: u64,  // chain table estimate (flat + offs + pure + hash index)
    top100: u64,     // reads covered by the 100 most frequent chains
}

struct Open {
    bi: usize,
    st: BlockStore,
    chain_cnt: Vec<u32>,
    stats: BlockStats,
}

/// Hot columns per read actually present in BlockStore v1 (store.rs), len-based:
/// name_hash 8, name_off 4 (cold arena offset), chain 4, start/end 8, raw 4, polya 16, flags 2,
/// mapq 1, group 4, atype+gene_atype 2, match_off 4, dels_off 4.
const HOT_FIXED: u64 = 8 + 4 + 4 + 8 + 4 + 16 + 2 + 1 + 4 + 2 + 4 + 4;

fn close(o: Open) -> BlockStats {
    let mut s = o.stats;
    let st = &o.st;
    s.chains = st.chains.len() as u64;
    s.chain_ivs = st.chains.n_intervals() as u64;
    s.heap = st.heap_bytes() as u64;
    s.tight_hot = s.reads * HOT_FIXED + s.dels * 8;
    // flat Iv 8 B + offs 4 B + pure 1 B per chain + FxHashMap<u64, SmallVec<[u32;1]>> entry
    // (8 + 16 B + 1 ctrl) at ~7/8 load -> ~29 B per chain.
    s.chain_tbl = s.chain_ivs * 8 + s.chains * (4 + 1 + 29);
    let mut c = o.chain_cnt;
    s.singletons = c.iter().filter(|&&x| x == 1).count() as u64;
    c.sort_unstable_by(|a, b| b.cmp(a));
    s.top100 = c.iter().take(100).map(|&x| x as u64).sum();
    s
}

fn run_store(bam: &BamHandle, chr: &str, blocks: &[Block], smp: &Sampler, params: &Params) -> Vec<BlockStats> {
    let tid = bam.ref_id(chr).unwrap();
    let len = bam.reference_lengths()[tid];
    let pp = PolyAParams::from_params(params);
    let mut ar = AlignedRead {
        exons: Vec::new(),
        read_blocks: Vec::new(),
        cigar_blocks: Vec::new(),
        polya: NO_POLYA,
        exons_changed: false,
        deletions: Vec::new(),
    };
    let mut rr = ReadRecord {
        idx: 0,
        name: String::new(),
        exons: Vec::new(),
        corrected_exons: Vec::new(),
        boundary_deletions: Vec::new(),
        polya: NO_POLYA,
        polya_found: false,
        strand: Strand::Unknown,
        mapped_strand: Strand::Unknown,
        mapq: 0,
        multimapper: false,
        group: 0,
        assignment: Assignment::empty(AssignmentType::Noninformative),
        genomic_region: Iv { s: 0, e: 0 },
    };
    let rss0 = rss_bytes();
    smp.reset();
    let t = Instant::now();
    let mut open: Vec<Open> = Vec::new(); // sorted by block index
    let mut done: Vec<BlockStats> = Vec::new();
    let (mut n_rec, mut n_kept) = (0u64, 0u64);
    let mut peak_open_bytes = 0u64;
    let mut last_report = Instant::now();
    for b in fetch_with(bam, chr, Iv0 { s: 0, e: len as i64 - 1 }, FetchOpts::default()) {
        for r in b.iter() {
            n_rec += 1;
            // retire blocks the stream has passed (coordinate-sorted input)
            while let Some(o) = open.first() {
                if blocks[o.bi].region.e < r.pos0 {
                    let o = open.remove(0);
                    done.push(close(o));
                } else {
                    break;
                }
            }
            // 4.1 collect filter: no unmapped / supplementary / secondary (use_secondary off)
            if r.is_unmapped() || r.is_supplementary() || r.is_secondary() || r.cigar.is_empty() {
                continue;
            }
            if !iq_read::process_alignment_into(&r, &pp, ExternalPolyA::default(), &mut ar) {
                continue;
            }
            let span = Iv0 { s: r.pos0, e: r.end0() - 1 };
            let bi = owner_block(span, blocks);
            let pos = match open.binary_search_by_key(&bi, |o| o.bi) {
                Ok(p) => p,
                Err(p) => {
                    let reg = blocks[bi].region;
                    open.insert(
                        p,
                        Open {
                            bi,
                            st: BlockStore::new(Iv::from(reg)),
                            chain_cnt: Vec::new(),
                            stats: BlockStats { region: (reg.s, reg.e), ..Default::default() },
                        },
                    );
                    p
                }
            };
            let o = &mut open[pos];
            rr.idx = o.st.len() as u32;
            rr.name.clear();
            rr.name.push_str(std::str::from_utf8(r.qname).unwrap_or("?"));
            rr.exons.clear();
            rr.exons.extend_from_slice(&ar.exons);
            rr.corrected_exons.clear();
            rr.corrected_exons.extend_from_slice(&ar.exons);
            rr.boundary_deletions.clear();
            rr.boundary_deletions.extend_from_slice(&ar.deletions);
            rr.polya = ar.polya;
            rr.polya_found =
                ar.polya.ext_a != -1 || ar.polya.ext_t != -1 || ar.polya.int_a != -1 || ar.polya.int_t != -1;
            rr.mapped_strand = if r.is_reverse() { Strand::Minus } else { Strand::Plus };
            rr.strand = rr.mapped_strand;
            rr.mapq = r.mapq;
            rr.genomic_region = o.st.region;
            o.st.push(&rr);
            n_kept += 1;
            // chain usage count (chain id of the row just pushed)
            let cid = o.st.view(o.st.len() - 1).chain_id() as usize;
            if cid >= o.chain_cnt.len() {
                o.chain_cnt.resize(cid + 1, 0);
            }
            o.chain_cnt[cid] += 1;
            o.stats.reads += 1;
            o.stats.introns += ar.exons.len() as u64 - 1;
            o.stats.mono += (ar.exons.len() == 1) as u64;
            o.stats.names += r.qname.len() as u64;
            o.stats.dels += ar.deletions.len() as u64;
        }
        if last_report.elapsed() > Duration::from_secs(30) {
            last_report = Instant::now();
            let ob: u64 = open.iter().map(|o| o.st.heap_bytes() as u64).sum();
            peak_open_bytes = peak_open_bytes.max(ob);
            eprintln!(
                "  [{:.0}s] {n_rec} rec, {n_kept} stored, open blocks {}, open store heap {:.1} MB, RSS {:.1} MB",
                t.elapsed().as_secs_f64(),
                open.len(),
                mb(ob),
                mb(rss_bytes())
            );
        }
    }
    let ob: u64 = open.iter().map(|o| o.st.heap_bytes() as u64).sum();
    peak_open_bytes = peak_open_bytes.max(ob);
    let peak_rss = smp.peak();
    for o in open.drain(..) {
        done.push(close(o));
    }
    let dt = t.elapsed().as_secs_f64();
    println!("== store build {chr}");
    println!(
        "records {n_rec}, stored {n_kept}; wall {dt:.2} s ({:.0} rec/s, {:.0} stored reads/s)",
        n_rec as f64 / dt,
        n_kept as f64 / dt
    );
    println!(
        "RSS before {:.1} MB, peak {:.1} MB (delta {:.1} MB); peak open-store heap seen {:.1} MB; VmHWM {:.1} MB",
        mb(rss0),
        mb(peak_rss),
        mb(peak_rss.saturating_sub(rss0)),
        mb(peak_open_bytes),
        mb(vm_hwm())
    );
    done
}

fn print_stats(label: &str, s: &BlockStats) {
    let n = s.reads.max(1) as f64;
    let names_b = s.names; // qname bytes (offsets are in HOT_FIXED)
    println!("-- {label}: region {}-{}", s.region.0, s.region.1);
    println!(
        "   reads {} (mono-exon {:.1}%), introns/read {:.2}, boundary dels/read {:.3}",
        s.reads,
        100.0 * s.mono as f64 / n,
        s.introns as f64 / n,
        s.dels as f64 / n
    );
    println!(
        "   distinct chains {} (reads/chain {:.1}; singletons {} = {:.1}% of reads; top-100 chains cover {:.1}% of reads); introns/distinct chain {:.2}",
        s.chains,
        n / s.chains.max(1) as f64,
        s.singletons,
        100.0 * s.singletons as f64 / n,
        100.0 * s.top100 as f64 / n,
        s.chain_ivs as f64 / s.chains.max(1) as f64
    );
    println!(
        "   B/read: hot columns {:.1} + chain table {:.1} = {:.1} (tight, no names); cold names {:.1}; heap_bytes() {:.1} (capacity, incl. names, excl. hash index)",
        s.tight_hot as f64 / n,
        s.chain_tbl as f64 / n,
        (s.tight_hot + s.chain_tbl) as f64 / n,
        names_b as f64 / n,
        s.heap as f64 / n
    );
    println!(
        "   totals: hot+chains {:.1} MB, names {:.1} MB, heap_bytes {:.1} MB, chain table {:.2} MB",
        mb(s.tight_hot + s.chain_tbl),
        mb(names_b),
        mb(s.heap),
        mb(s.chain_tbl)
    );
    println!(
        "   without interning the chains would cost {:.1} B/read (8 B/intron + 4 B offset)",
        (s.introns * 8 + s.reads * 4) as f64 / n
    );
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        eprintln!("usage: vtg_proto passa|store <bam> <chr> [inflate_threads]");
        std::process::exit(2);
    }
    let (mode, path, chr) = (&args[1], &args[2], &args[3]);
    let threads: usize = args.get(4).and_then(|x| x.parse().ok()).unwrap_or(4);
    let smp = Sampler::start();
    let mut bam = BamHandle::open(Path::new(path)).expect("open bam");
    bam.set_threads(threads);
    println!("bam {path} chr {chr} inflate threads {threads}; RSS at start {:.1} MB", mb(rss_bytes()));
    let blocks = run_pass_a(&bam, chr, &smp);
    if mode == "passa" {
        return;
    }
    let params = Params::pacbio();
    let stats = run_store(&bam, chr, &blocks, &smp, &params);
    let mut agg = BlockStats::default();
    for s in &stats {
        agg.reads += s.reads;
        agg.introns += s.introns;
        agg.chains += s.chains;
        agg.chain_ivs += s.chain_ivs;
        agg.singletons += s.singletons;
        agg.mono += s.mono;
        agg.names += s.names;
        agg.dels += s.dels;
        agg.heap += s.heap;
        agg.tight_hot += s.tight_hot;
        agg.chain_tbl += s.chain_tbl;
        agg.top100 += s.top100;
    }
    let mut sorted: Vec<&BlockStats> = stats.iter().collect();
    sorted.sort_by_key(|s| std::cmp::Reverse(s.reads));
    println!("blocks with stored reads: {}", stats.len());
    for (i, s) in sorted.iter().take(3).enumerate() {
        print_stats(&format!("block #{} by reads", i + 1), s);
    }
    print_stats("ALL blocks (sum; top-100 = Σ per block)", &agg);
}
