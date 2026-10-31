//! U2 throughput benchmark and real-data cross-check.
//!
//! ```text
//! bam_bench full    <bam> <chr> [reps]       full-chromosome fetch, noodles (+ htslib) x threads x seq
//! bam_bench regions <bam> <chr> [n] [seed]   n random regions of 10-200 kb (pass-B pattern)
//! bam_bench check   <bam> <regions.tsv>      per region "chr:s-e<TAB>count<TAB>crc32(qname\tpos\tflag\tend\tseq\tHP\n...)"
//! bam_bench fasta   <fa> <regions.tsv>       per region "chr:s-e<TAB>fnv64(slice)<TAB>len"
//! ```
//! regions.tsv: `chr<TAB>s<TAB>e` (1-based closed). Build with `--features htslib` to add the
//! htslib rows.

use iq_bam::{try_fetch, BamHandle, FetchOpts, FastaIndex, RecordBatch};
use iq_types::{Iv, Iv0};
use std::path::Path;
use std::time::Instant;

struct XorShift(u64);
impl XorShift {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
}

/// zlib-compatible CRC-32 (IEEE), incremental: `crc = crc32(crc, bytes)`.
fn crc32(crc: u32, bytes: &[u8]) -> u32 {
    static TABLE: [u32; 256] = {
        let mut t = [0u32; 256];
        let mut i = 0;
        while i < 256 {
            let mut c = i as u32;
            let mut k = 0;
            while k < 8 {
                c = if c & 1 != 0 { 0xedb88320 ^ (c >> 1) } else { c >> 1 };
                k += 1;
            }
            t[i] = c;
            i += 1;
        }
        t
    };
    let mut c = !crc;
    for &b in bytes {
        c = TABLE[((c ^ b as u32) & 0xff) as usize] ^ (c >> 8);
    }
    !c
}

fn fnv(h: &mut u64, bytes: &[u8]) {
    for &b in bytes {
        *h ^= b as u64;
        *h = h.wrapping_mul(0x100000001b3);
    }
}

/// Touches what pass A/B read (qname, flag, pos, mapq, cigar → end; seq if decoded).
#[derive(Default)]
struct Sink {
    records: u64,
    acc: u64,
}
impl Sink {
    fn eat(&mut self, b: &RecordBatch) {
        for r in b.iter() {
            self.records += 1;
            self.acc = self.acc.wrapping_add(r.qname.len() as u64 + r.flag as u64 + r.mapq as u64 + r.end0() as u64 + r.seq.len() as u64);
        }
    }
}

fn regions_of(len: u64, n: usize, seed: u64) -> Vec<Iv0> {
    let mut rng = XorShift(seed.max(1));
    (0..n)
        .map(|_| {
            let w = 10_000 + rng.next() % 190_001;
            let s = rng.next() % len.saturating_sub(w).max(1);
            Iv0 { s: s as i64, e: (s + w - 1) as i64 }
        })
        .collect()
}

fn row(backend: &str, mode: &str, threads: usize, seq: bool, recs: u64, mb: f64, secs: f64, extra: &str) {
    println!(
        "{backend}\t{mode}\t{threads}\t{}\t{recs}\t{secs:.3}\t{:.0}\t{:.1}\t{extra}",
        if seq { "seq" } else { "noseq" },
        recs as f64 / secs,
        mb / secs
    );
}

#[allow(unused_variables)]
fn htslib_full(bam: &Path, chr: &str, region: Iv0, threads: usize, seq: bool) -> Option<(u64, f64)> {
    #[cfg(feature = "htslib")]
    {
        let mut h = iq_bam::HtsBam::open(bam).unwrap();
        h.threads = threads;
        h.tag = Some(*b"HP");
        let t = Instant::now();
        let mut sink = Sink::default();
        for b in h.try_fetch(chr, region, FetchOpts { seq, ..Default::default() }).unwrap() {
            sink.eat(&b.unwrap());
        }
        std::hint::black_box(sink.acc);
        return Some((sink.records, t.elapsed().as_secs_f64()));
    }
    #[allow(unreachable_code)]
    None
}

#[allow(unused_variables)]
fn htslib_regions(bam: &Path, chr: &str, regs: &[Iv0]) -> Option<(u64, f64)> {
    #[cfg(feature = "htslib")]
    {
        let mut h = iq_bam::HtsBam::open(bam).unwrap();
        h.tag = Some(*b"HP");
        let t = Instant::now();
        let mut sink = Sink::default();
        for &r in regs {
            for b in h.try_fetch(chr, r, FetchOpts::default()).unwrap() {
                sink.eat(&b.unwrap());
            }
        }
        std::hint::black_box(sink.acc);
        return Some((sink.records, t.elapsed().as_secs_f64()));
    }
    #[allow(unreachable_code)]
    None
}

fn read_regions(p: &str) -> Vec<(String, Iv)> {
    std::fs::read_to_string(p)
        .unwrap()
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (c[0].to_owned(), Iv { s: c[1].parse().unwrap(), e: c[2].parse().unwrap() })
        })
        .collect()
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let usage = "usage: bam_bench full|regions|check|fasta ...";
    let cmd = a.get(1).expect(usage).as_str();
    match cmd {
        "full" => {
            let bam = Path::new(&a[2]);
            let chr = &a[3];
            let reps: usize = a.get(4).map_or(2, |x| x.parse().unwrap());
            let mb = std::fs::metadata(bam).unwrap().len() as f64 / 1e6;
            let h0 = BamHandle::open(bam).unwrap();
            let len = h0.reference_lengths()[h0.ref_id(chr).expect("chr")];
            let region = Iv0 { s: 0, e: len as i64 - 1 };
            println!("backend\tmode\tthreads\tseq\trecords\tsecs\trec_per_s\tMB_per_s\tnote");
            for seq in [false, true] {
                for threads in [1usize, 2, 4, 8] {
                    let mut best: Option<(u64, f64)> = None;
                    for _ in 0..reps {
                        let mut h = BamHandle::open(bam).unwrap();
                        h.set_tag(Some(*b"HP"));
                        h.set_threads(threads);
                        let t = Instant::now();
                        let mut sink = Sink::default();
                        for b in try_fetch(&h, chr, region, FetchOpts { seq, ..Default::default() }).unwrap() {
                            sink.eat(&b.unwrap());
                        }
                        std::hint::black_box(sink.acc);
                        let r = (sink.records, t.elapsed().as_secs_f64());
                        best = Some(best.map_or(r, |b: (u64, f64)| if r.1 < b.1 { r } else { b }));
                    }
                    let (n, s) = best.unwrap();
                    row("noodles", "full", threads, seq, n, mb, s, "");
                    let mut hb: Option<(u64, f64)> = None;
                    for _ in 0..reps {
                        if let Some(r) = htslib_full(bam, chr, region, threads, seq) {
                            hb = Some(hb.map_or(r, |b: (u64, f64)| if r.1 < b.1 { r } else { b }));
                        }
                    }
                    if let Some((n, s)) = hb {
                        row("htslib", "full", threads, seq, n, mb, s, "");
                    }
                }
            }
        }
        "regions" => {
            let bam = Path::new(&a[2]);
            let chr = &a[3];
            let n: usize = a.get(4).map_or(2000, |x| x.parse().unwrap());
            let seed: u64 = a.get(5).map_or(42, |x| x.parse().unwrap());
            let mut h = BamHandle::open(bam).unwrap();
            h.set_tag(Some(*b"HP"));
            let len = h.reference_lengths()[h.ref_id(chr).expect("chr")];
            let regs = regions_of(len, n, seed);
            let kb: f64 = regs.iter().map(|r| (r.e - r.s + 1) as f64).sum::<f64>() / 1e3;
            println!("backend\tmode\tthreads\tseq\trecords\tsecs\trec_per_s\tregions_per_s\tnote");
            for _rep in 0..2 {
                let t = Instant::now();
                let mut sink = Sink::default();
                for &r in &regs {
                    for b in try_fetch(&h, chr, r, FetchOpts::default()).unwrap() {
                        sink.eat(&b.unwrap());
                    }
                }
                std::hint::black_box(sink.acc);
                let s = t.elapsed().as_secs_f64();
                println!("noodles\tregions\t1\tseq\t{}\t{s:.3}\t{:.0}\t{:.1}\t{n} regions, {kb:.0} kb", sink.records, sink.records as f64 / s, n as f64 / s);
                if let Some((m, s)) = htslib_regions(bam, chr, &regs) {
                    println!("htslib\tregions\t1\tseq\t{m}\t{s:.3}\t{:.0}\t{:.1}\t{n} regions, {kb:.0} kb", m as f64 / s, n as f64 / s);
                }
            }
        }
        "check" => {
            let bam = Path::new(&a[2]);
            let mut h = BamHandle::open(bam).unwrap();
            h.set_tag(Some(*b"HP"));
            for (chr, iv) in read_regions(&a[3]) {
                let mut hash = 0u32;
                let mut count = 0u64;
                for b in try_fetch(&h, &chr, Iv0 { s: iv.s as i64 - 1, e: iv.e as i64 - 1 }, FetchOpts::default()).unwrap() {
                    let b = b.unwrap();
                    for r in b.iter() {
                        count += 1;
                        hash = crc32(hash, r.qname);
                        hash = crc32(hash, format!("\t{}\t{}\t{}\t", r.pos0, r.flag, r.end0()).as_bytes());
                        hash = crc32(hash, r.seq);
                        hash = crc32(hash, format!("\t{:?}\n", r.tag).as_bytes());
                    }
                }
                println!("{chr}:{}-{}\t{count}\t{hash:08x}", iv.s, iv.e);
            }
        }
        "fasta" => {
            let fa = FastaIndex::open(Path::new(&a[2])).unwrap();
            for (chr, iv) in read_regions(&a[3]) {
                let s = fa.slice(&chr, iv).unwrap();
                let mut hash = 0xcbf29ce484222325u64;
                fnv(&mut hash, &s);
                println!("{chr}:{}-{}\t{hash:016x}\t{}", iv.s, iv.e, s.len());
            }
        }
        _ => panic!("{usage}"),
    }
}
