//! Mock collect/construct stages with configurable cost and memory (U21 stage 1).
//!
//! The mock reproduces every *coupling* the real pipeline has, so scheduling mistakes change
//! its outputs:
//! - `detected_known`: a known model of gene `g` is printed only by the first block (in
//!   genomic order) containing `g`; splitting blocks that share a gene prints it twice;
//! - duplicate read names (chromosome-local, resolved at the barrier): the first copy (block,
//!   row order) stays and is counted as ambiguous, later copies are suspended;
//! - the sample-wide polyA decision: with `requires_polya` only polyA reads support novel
//!   models;
//! - novel IDs: block-local, renumbered per chromosome in Phase W in genomic order.
//!
//! Costs are busy-spins (CPU time, so utilisation is measurable from `/proc/self/stat`);
//! memory is real touched heap (`ballast`), so RSS is measurable. Blocks with more than
//! `par_threshold` reads do their per-read work as nested rayon tasks, like a giant block.

use crate::pipeline::{Pipeline, ReadTotals};
use iq_types::{Block, Iv0};
use rayon::prelude::*;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fmt::Write as _;
use std::path::PathBuf;
use std::sync::Mutex;
use std::time::{Duration, Instant};

#[derive(Clone, Debug)]
pub struct MockBlock {
    pub genes: Vec<u32>,
    pub reads: u32,
}

#[derive(Clone, Debug)]
pub struct MockChr {
    pub name: String,
    pub blocks: Vec<MockBlock>,
}

#[derive(Clone, Debug)]
pub struct MockSpec {
    pub chrs: Vec<MockChr>,
    pub seed: u64,
    pub passa_ns_per_read: u64,
    pub collect_ns_per_read: u64,
    pub construct_ns_per_read: u64,
    /// Touched heap per stored read (the `BlockStore`).
    pub store_bytes_per_read: u64,
    /// Touched transient heap per read while a block is constructed.
    pub construct_bytes_per_read: u64,
    /// Per-read polyA probability.
    pub polya_rate: f64,
    /// Per-read probability that the name duplicates a read of an earlier block (same chr).
    pub dup_rate: f64,
    /// Blocks with more reads than this use nested rayon parallelism.
    pub par_threshold: u32,
    /// Distinct "intron chains" per 64 reads (novel-model diversity).
    pub chains_per_64: u32,
    /// If set, `finish` writes `OUT.gtf` and `OUT.counts.tsv` here.
    pub out_dir: Option<PathBuf>,
    /// If set, stores can be spilled: the ballast is written to a file here and freed.
    pub spill_dir: Option<PathBuf>,
    /// Scheduler estimates use these bytes/read (tests make them match the ballast).
    pub est_collect_per_read: u64,
    pub est_construct_per_read: u64,
}

impl Default for MockSpec {
    fn default() -> Self {
        MockSpec {
            chrs: vec![],
            seed: 1,
            passa_ns_per_read: 0,
            collect_ns_per_read: 0,
            construct_ns_per_read: 0,
            store_bytes_per_read: 0,
            construct_bytes_per_read: 0,
            polya_rate: 0.8,
            dup_rate: 0.002,
            par_threshold: 20_000,
            chains_per_64: 3,
            out_dir: None,
            spill_dir: None,
            est_collect_per_read: 160,
            est_construct_per_read: 48,
        }
    }
}

/// SplitMix64.
fn mix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9e3779b97f4a7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}

fn unit_f(h: u64) -> f64 {
    (h >> 11) as f64 / (1u64 << 53) as f64
}

impl MockSpec {
    /// A genome-like layout: `n_chr` chromosomes, chromosome `i` gets `weights[i]` blocks-ish;
    /// genes span 1-3 consecutive blocks with probability `p_span`; block sizes are
    /// log-uniform up to `max_reads`.
    pub fn synthetic(seed: u64, n_chr: usize, blocks_per_chr: usize, max_reads: u32, p_share: f64) -> Self {
        let mut chrs = Vec::new();
        let mut gene = 0u32;
        for c in 0..n_chr {
            let mut blocks = Vec::new();
            let mut prev_genes: Vec<u32> = vec![];
            for b in 0..blocks_per_chr {
                let h = mix(seed ^ mix((c as u64) << 32 | b as u64));
                let reads = (((max_reads as f64).ln() * unit_f(mix(h))).exp()) as u32;
                let mut genes = Vec::new();
                let r = unit_f(mix(h ^ 7));
                if r < p_share && !prev_genes.is_empty() {
                    genes.push(*prev_genes.last().unwrap()); // continues the previous gene
                }
                if unit_f(mix(h ^ 9)) < 0.7 {
                    gene += 1;
                    genes.push(gene);
                }
                if unit_f(mix(h ^ 11)) < 0.1 {
                    gene += 1;
                    genes.push(gene);
                }
                prev_genes = genes.clone();
                blocks.push(MockBlock { genes, reads });
            }
            chrs.push(MockChr { name: format!("{}", c + 1), blocks });
        }
        MockSpec { chrs, seed, ..Default::default() }
    }

    fn read(&self, chr: usize, block: usize, i: u32) -> MockRead {
        let h = mix(self.seed ^ mix(((chr as u64) << 48) ^ ((block as u64) << 24) ^ i as u64));
        let mut name = ((chr as u64) << 48) | ((block as u64) << 24) | i as u64;
        if block > 0 && unit_f(mix(h ^ 1)) < self.dup_rate {
            let ob = (mix(h ^ 2) % block as u64) as usize;
            let n = self.chrs[chr].blocks[ob].reads;
            if n > 0 {
                let oi = (mix(h ^ 3) % n as u64) as u32;
                name = ((chr as u64) << 48) | ((ob as u64) << 24) | oi as u64;
            }
        }
        let n = self.chrs[chr].blocks[block].reads.max(1);
        let k = (n as u64 * self.chains_per_64.max(1) as u64 / 64).max(1);
        MockRead {
            name,
            polya: unit_f(mix(h ^ 4)) < self.polya_rate,
            // skewed chain distribution: a few frequent chains, a long tail
            chain: ((unit_f(mix(h ^ 5)).powi(3) * k as f64) as u32).min(k as u32 - 1),
            has_gene: unit_f(mix(h ^ 6)) < 0.5,
            status: 0,
        }
    }
}

fn spin(d: Duration) {
    if d.is_zero() {
        return;
    }
    let t = Instant::now();
    while t.elapsed() < d {
        std::hint::spin_loop();
    }
}

/// Spins `ns_per_read × n` in chunks (so nested parallel chunks spread the cost).
fn spin_reads(ns_per_read: u64, n: usize) {
    spin(Duration::from_nanos(ns_per_read * n as u64));
}

#[derive(Clone, Copy, Debug)]
pub struct MockRead {
    pub name: u64,
    pub polya: bool,
    pub chain: u32,
    pub has_gene: bool,
    /// 0 normal, 1 duplicate kept (multimapper → ambiguous), 2 suspended.
    pub status: u8,
}

pub struct MockStore {
    pub reads: Vec<MockRead>,
    pub ballast: Vec<u8>,
    pub spilled: Option<PathBuf>,
}

pub struct MockBlockOut {
    start: i64,
    known: Vec<u32>,
    novel: Vec<(u32, f64)>,
    known_counts: Vec<(u32, f64)>,
    ambiguous: u64,
    no_feature: u64,
}

pub struct MockChrOut {
    gtf: String,
    counts: Vec<(String, f64)>,
    ambiguous: u64,
    no_feature: u64,
}

/// Final texts (`gtf`, `counts`) of the last `finish`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MockResult {
    pub gtf: String,
    pub counts: String,
}

pub struct MockPipeline {
    pub spec: MockSpec,
    pub result: Mutex<Option<MockResult>>,
    starts: Vec<Vec<i64>>,
}

impl MockPipeline {
    pub fn new(spec: MockSpec) -> Self {
        let starts = spec
            .chrs
            .iter()
            .map(|c| {
                let mut s = 0i64;
                c.blocks
                    .iter()
                    .map(|_| {
                        s += 100_000;
                        s
                    })
                    .collect()
            })
            .collect();
        MockPipeline { spec, result: Mutex::new(None), starts }
    }

    pub fn take_result(&self) -> MockResult {
        self.result.lock().unwrap().take().expect("finish ran")
    }
}

pub struct MockGlobal {
    pub requires_polya: bool,
}

impl Pipeline for MockPipeline {
    type Store = MockStore;
    type UnitState = BTreeSet<u32>;
    type Global = MockGlobal;
    type BlockOut = MockBlockOut;
    type ChrOut = MockChrOut;

    fn n_chromosomes(&self) -> usize {
        self.spec.chrs.len()
    }
    fn chr_name(&self, chr: usize) -> &str {
        &self.spec.chrs[chr].name
    }
    fn chr_weight(&self, chr: usize) -> u64 {
        self.spec.chrs[chr].blocks.iter().map(|b| b.reads as u64).sum()
    }

    fn pass_a(&self, chr: usize) -> Result<Vec<Block>, String> {
        spin_reads(self.spec.passa_ns_per_read, self.chr_weight(chr) as usize);
        Ok(self.spec.chrs[chr]
            .blocks
            .iter()
            .enumerate()
            .map(|(i, b)| Block {
                region: Iv0 { s: self.starts[chr][i], e: self.starts[chr][i] + 50_000 },
                max_cov: b.reads,
                skip: 1,
                gene_ids: b.genes.clone(),
                owned_reads: b.reads as u64,
            })
            .collect())
    }

    fn est_pass_a(&self, _chr: usize) -> u64 {
        1 << 20
    }

    fn collect_block(&self, chr: usize, bi: usize, block: &Block) -> Result<MockStore, String> {
        let n = block.owned_reads as u32;
        let s = &self.spec;
        let reads: Vec<MockRead> = if n > s.par_threshold {
            (0..n)
                .into_par_iter()
                .chunks(4096)
                .flat_map_iter(|ch| {
                    spin_reads(s.collect_ns_per_read, ch.len());
                    ch.into_iter().map(|i| s.read(chr, bi, i)).collect::<Vec<_>>()
                })
                .collect()
        } else {
            spin_reads(s.collect_ns_per_read, n as usize);
            (0..n).map(|i| s.read(chr, bi, i)).collect()
        };
        let ballast = vec![1u8; (n as u64 * s.store_bytes_per_read) as usize];
        Ok(MockStore { reads, ballast, spilled: None })
    }

    fn store_bytes(&self, st: &MockStore) -> u64 {
        (st.reads.len() * std::mem::size_of::<MockRead>() + st.ballast.len()) as u64
    }

    fn can_spill(&self) -> bool {
        self.spec.spill_dir.is_some()
    }

    fn spill(&self, chr: usize, bi: usize, st: &mut MockStore) -> Result<u64, String> {
        let dir = self.spec.spill_dir.as_ref().expect("spill_dir");
        let path = dir.join(format!("c{chr}_b{bi}.spill"));
        std::fs::write(&path, &st.ballast).map_err(|e| format!("spill {}: {e}", path.display()))?;
        st.ballast = Vec::new();
        st.spilled = Some(path);
        Ok(self.store_bytes(st))
    }

    fn barrier_chr(&self, _chr: usize, _blocks: &[Block], stores: &mut [MockStore]) -> Result<ReadTotals, String> {
        // first copy wins (block, row order); later copies suspended; the kept copy of a
        // duplicated name becomes a multimapper (ambiguous)
        let mut first: HashMap<u64, (usize, usize)> = HashMap::new();
        let mut dup_first: Vec<(usize, usize)> = Vec::new();
        for (b, st) in stores.iter_mut().enumerate() {
            for (r, rd) in st.reads.iter_mut().enumerate() {
                match first.entry(rd.name) {
                    std::collections::hash_map::Entry::Vacant(e) => {
                        e.insert((b, r));
                    }
                    std::collections::hash_map::Entry::Occupied(e) => {
                        rd.status = 2;
                        dup_first.push(*e.get());
                    }
                }
            }
        }
        for (b, r) in dup_first {
            stores[b].reads[r].status = 1;
        }
        let mut t = ReadTotals::default();
        for st in stores.iter() {
            for rd in &st.reads {
                if rd.status != 2 {
                    t.total += 1;
                    t.polya += rd.polya as u64;
                }
            }
        }
        Ok(t)
    }

    fn decide(&self, t: ReadTotals) -> MockGlobal {
        let frac = if t.total > 0 { t.polya as f64 / t.total as f64 } else { 0.0 };
        MockGlobal { requires_polya: frac >= 0.7 }
    }

    fn construct_block(
        &self,
        chr: usize,
        bi: usize,
        block: &Block,
        store: MockStore,
        known_set: &mut BTreeSet<u32>,
        g: &MockGlobal,
    ) -> Result<MockBlockOut, String> {
        let s = &self.spec;
        let n = store.reads.len();
        let transient = vec![1u8; (n as u64 * s.construct_bytes_per_read) as usize];
        let k = store.reads.iter().map(|r| r.chain as usize + 1).max().unwrap_or(0);
        let support = |r: &MockRead| r.status == 0 && (!g.requires_polya || r.polya);
        let chain_counts: Vec<u32> = if n > s.par_threshold as usize {
            store
                .reads
                .par_chunks(4096)
                .map(|ch| {
                    spin_reads(s.construct_ns_per_read, ch.len());
                    let mut c = vec![0u32; k];
                    for r in ch.iter().filter(|r| support(r)) {
                        c[r.chain as usize] += 1;
                    }
                    c
                })
                .reduce(
                    || vec![0u32; k],
                    |mut a, b| {
                        a.iter_mut().zip(b).for_each(|(x, y)| *x += y);
                        a
                    },
                )
        } else {
            spin_reads(s.construct_ns_per_read, n);
            let mut c = vec![0u32; k];
            for r in store.reads.iter().filter(|r| support(r)) {
                c[r.chain as usize] += 1;
            }
            c
        };
        let mut known = Vec::new();
        for &gene in &block.gene_ids {
            if known_set.insert(gene) {
                known.push(gene);
            }
        }
        let novel_idx: Vec<Option<usize>> = {
            let mut i = 0;
            chain_counts
                .iter()
                .map(|&c| {
                    (c >= 3).then(|| {
                        i += 1;
                        i - 1
                    })
                })
                .collect()
        };
        let mut novel: Vec<(u32, f64)> =
            chain_counts.iter().enumerate().filter(|(_, &c)| c >= 3).map(|(ch, _)| (ch as u32, 0.0)).collect();
        let mut kc: BTreeMap<u32, f64> = BTreeMap::new();
        let (mut ambiguous, mut no_feature) = (0, 0);
        for r in &store.reads {
            match r.status {
                2 => continue,
                1 => ambiguous += 1,
                _ => {
                    if let Some(Some(j)) = novel_idx.get(r.chain as usize) {
                        novel[*j].1 += 1.0;
                    } else if r.has_gene && !block.gene_ids.is_empty() {
                        *kc.entry(block.gene_ids[0]).or_default() += 1.0;
                    } else {
                        no_feature += 1;
                    }
                }
            }
        }
        drop(transient);
        if let Some(p) = &store.spilled {
            let _ = std::fs::remove_file(p);
        }
        drop(store);
        Ok(MockBlockOut {
            start: self.starts[chr][bi],
            known,
            novel,
            known_counts: kc.into_iter().collect(),
            ambiguous,
            no_feature,
        })
    }

    fn merge_chr(&self, chr: usize, _blocks: &[Block], outs: Vec<MockBlockOut>) -> Result<MockChrOut, String> {
        let name = &self.spec.chrs[chr].name;
        let mut gtf = String::new();
        let mut counts: BTreeMap<String, f64> = BTreeMap::new();
        let mut next = 0u64;
        let (mut amb, mut nf) = (0, 0);
        for o in outs {
            for g in &o.known {
                let _ = writeln!(gtf, "{name}\t{}\tknown\tgene_{g}", o.start);
            }
            for &(ch, c) in &o.novel {
                next += 1;
                let id = format!("transcript{next}.{name}.nic");
                let _ = writeln!(gtf, "{name}\t{}\tnovel\t{id}\tchain{ch}", o.start);
                if c > 0.0 {
                    *counts.entry(id).or_default() += c;
                }
            }
            for &(g, c) in &o.known_counts {
                *counts.entry(format!("known_{g}")).or_default() += c;
            }
            amb += o.ambiguous;
            nf += o.no_feature;
        }
        Ok(MockChrOut { gtf, counts: counts.into_iter().collect(), ambiguous: amb, no_feature: nf })
    }

    fn finish(&self, outs: Vec<MockChrOut>, totals: ReadTotals) -> Result<(), String> {
        let mut order: Vec<usize> = (0..outs.len()).collect();
        order.sort_by(|&a, &b| iq_quant::output::natural_key_cmp(&self.spec.chrs[a].name, &self.spec.chrs[b].name));
        let frac = if totals.total > 0 { totals.polya as f64 / totals.total as f64 } else { 0.0 };
        let mut gtf = format!("# mock\n# requires_polya {}\n", frac >= 0.7);
        let mut counts = String::from("feature_id\tcount\n");
        let (mut amb, mut nf) = (0, 0);
        for &c in &order {
            gtf.push_str(&outs[c].gtf);
            for (id, v) in &outs[c].counts {
                let _ = writeln!(counts, "{id}\t{}", iq_quant::fmt_fixed2(*v));
            }
            amb += outs[c].ambiguous;
            nf += outs[c].no_feature;
        }
        let _ = write!(counts, "__ambiguous\t{amb}\n__no_feature\t{nf}\n__totals\t{}\t{}\n", totals.total, totals.polya);
        if let Some(d) = &self.spec.out_dir {
            std::fs::create_dir_all(d).map_err(|e| e.to_string())?;
            std::fs::write(d.join("OUT.gtf"), &gtf).map_err(|e| e.to_string())?;
            std::fs::write(d.join("OUT.counts.tsv"), &counts).map_err(|e| e.to_string())?;
        }
        *self.result.lock().unwrap() = Some(MockResult { gtf, counts });
        Ok(())
    }

    fn est_collect(&self, block: &Block) -> u64 {
        block.owned_reads * self.spec.est_collect_per_read
    }
    fn est_construct(&self, block: &Block) -> u64 {
        block.owned_reads * self.spec.est_construct_per_read
    }
}
