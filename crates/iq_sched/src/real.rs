//! The real pipeline (U21 stage 2): [`Pipeline`] over the done crates.
//!
//! | step | crate |
//! |---|---|
//! | setup | iq_annot (GFF3/GTF), iq_bam (BAM header/BAI, FASTA mmap) |
//! | pass A | iq_bam fetch (no seq) → iq_regions `PassA` with the annotation gene oracle |
//! | collect | iq_annot block index → iq_bam fetch → iq_collect `BlockCollector` (U8) |
//! | barrier | U9 `iq_collect::multimap::resolve_chr_duplicates` → polyA fraction → `requires_polya_for_construction` |
//! | construct | [`crate::construct::construct_block`] (U10b-U19 + U20 counter) |
//! | merge / write | iq_quant `ChrPhaseW`, `write_outputs` |

use crate::pipeline::{Pipeline, ReadTotals};
use crate::Config;
use iq_annot::{AnnotationIndex, GeneTypes};
use iq_bam::{BamHandle, BamOptions, FastaIndex, FetchOpts};
use iq_collect::{BlockCollector, ClusterOwner, CollectStats, GroupIds, NoGroups, SharedGroups};
use iq_model::DetectedKnown;
use iq_quant::{BlockOutput, ChrOutput, ChrPhaseW, OutputSpec, TxCounter};
use iq_regions::{GeneOracle, PassA, PassAConfig};
use iq_types::{
    Block, BlockStore, ChrIdx, GeneIdx, IsoformIndex, Iv, Iv0, Params, PolyAUsage, RefSlice, TranscriptModel,
};
use std::io::Write as _;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;

/// Bases added on each side of a cluster region for the reference window of its blocks.
pub const REF_MARGIN: i64 = 1000;
/// Blocks with at least this many owned reads fetch with multithreaded BGZF inflate.
pub const BIG_BLOCK_READS: u64 = 20_000;

struct ChrInfo {
    name: String,
    ann: Option<ChrIdx>,
    len: u64,
}

/// What pass A leaves for the later phases of one chromosome.
#[derive(Default)]
struct ChrMeta {
    regions: Vec<Iv0>,
    cluster_of_block: Vec<Iv0>,
    dup_candidates: Vec<u64>,
}

pub struct RealPipeline {
    pub cfg: Config,
    pub params: Params,
    pub ann: AnnotationIndex,
    fasta: FastaIndex,
    bam_passa: BamHandle,
    bam_collect: BamHandle,
    /// Same with multithreaded inflate, for blocks of ≥ [`BIG_BLOCK_READS`] reads.
    bam_big: BamHandle,
    chrs: Vec<ChrInfo>,
    meta: Vec<Mutex<ChrMeta>>,
    groups: Option<SharedGroups>,
    pub cmd_line: String,
    /// `--emit-read-assignments` sink (dev); chromosomes append as they merge.
    reads_out: Option<Mutex<std::io::BufWriter<std::fs::File>>>,
    /// Optional phasing step (`--phase-vcf`).
    phase: Option<PhaseState>,
}

/// State of the optional phasing step.
struct PhaseState {
    sites: iq_phase::HetSites,
    cfg: iq_phase::Cfg,
    seed: u64,
    /// Loci per chromosome (index into `chrs`), filled at merge.
    loci: Mutex<Vec<(usize, Vec<iq_phase::Locus>)>>,
    vcf_s: f64,
    /// Het-evidence collection: summed nanoseconds, records scanned, reads with calls, calls.
    ev_ns: AtomicU64,
    ev_records: AtomicU64,
    ev_reads: AtomicU64,
    ev_calls: AtomicU64,
    /// Locus building at merge (summed nanoseconds).
    loci_ns: AtomicU64,
}

pub struct RealStore {
    pub store: BlockStore,
    pub index: IsoformIndex,
    pub stats: CollectStats,
    /// Name hashes of reads this block shares with a neighbouring block (4.1 double
    /// processing of valley-straddling reads); U9 candidates at the barrier.
    pub shared: Vec<u64>,
    /// Het evidence of the block's reads (phasing only).
    pub phase_ev: Option<iq_phase::BlockEvidence>,
}

pub struct RealGlobal {
    pub params: Params,
    pub polya_fraction: f64,
}

pub struct RealBlockOut {
    index: IsoformIndex,
    models: Vec<TranscriptModel>,
    counter: TxCounter,
    reads: Vec<(String, u32, iq_types::FinalAssign)>,
    phase: Option<iq_phase::PhaseReads>,
}

struct AnnOracle<'a> {
    ann: &'a AnnotationIndex,
    chr: Option<ChrIdx>,
}

impl GeneOracle for AnnOracle<'_> {
    fn genes_for_region0(&self, r: Iv0, out: &mut Vec<GeneIdx>) {
        if let Some(c) = self.chr {
            out.extend(self.ann.genes_for_region0(c, r));
        }
    }
}

fn err<E: std::fmt::Display>(ctx: &str) -> impl Fn(E) -> String + '_ {
    move |e| format!("{ctx}: {e}")
}

/// Parameters for a run: preset or hidden `--params-json`, plus the CLI flags.
pub fn run_params(cfg: &Config) -> Result<Params, String> {
    let mut p = match &cfg.params_json {
        Some(path) => {
            let s = std::fs::read_to_string(path).map_err(err(&path.display().to_string()))?;
            let v: iq_types::serde_json::Value = iq_types::serde_json::from_str(&s).map_err(err("params json"))?;
            let m = v.as_object().ok_or("params json: not an object")?;
            Params::from_params_json(m)?
        }
        None => Params::preset(cfg.data_type),
    };
    p.report_novel_unspliced = cfg.report_novel_unspliced;
    Ok(p)
}

impl RealPipeline {
    /// Phase 0: annotation, FASTA, BAM header/index, chromosome list.
    pub fn open(cfg: &Config) -> Result<Self, String> {
        // the phased VCF loads while the annotation is parsed
        let vcf_job = cfg.phase.as_ref().map(|pc| {
            let (path, sample) = (pc.vcf.clone(), pc.sample.clone());
            std::thread::spawn(move || {
                let t = std::time::Instant::now();
                iq_phase::vcf::load(&path, sample.as_deref()).map(|s| (s, t.elapsed().as_secs_f64()))
            })
        });
        let params = run_params(cfg)?;
        let gt = if cfg.gene_types_isoquant { GeneTypes::Isoquant } else { GeneTypes::All };
        let a = cfg.annotation.to_string_lossy();
        let ann = if a.ends_with(".gtf") || a.ends_with(".gtf.gz") {
            iq_annot::parse_gtf(&cfg.annotation, gt)
        } else {
            iq_annot::parse_gff3(&cfg.annotation, gt)
        }
        .map_err(err(&a))?;
        let fasta = FastaIndex::open(&cfg.reference).map_err(err(&cfg.reference.display().to_string()))?;
        let tag = match &cfg.read_group_tag {
            Some(t) if t.len() == 2 => Some([t.as_bytes()[0], t.as_bytes()[1]]),
            Some(t) => return Err(format!("read-group tag must be 2 characters, got {t:?}")),
            None => None,
        };
        let bam_path = cfg.bam.display().to_string();
        let bam_passa = BamHandle::open_with(&cfg.bam, BamOptions { tag: None, threads: cfg.threads.clamp(1, 4) })
            .map_err(err(&bam_path))?;
        let bam_collect = BamHandle::open_with(&cfg.bam, BamOptions { tag, threads: 1 }).map_err(err(&bam_path))?;
        let bam_big = BamHandle::open_with(&cfg.bam, BamOptions { tag, threads: cfg.threads.clamp(1, 4) })
            .map_err(err(&bam_path))?;

        // BAM ∩ FASTA ∩ annotation, longest FASTA record first (dataset_processor.py:342-391)
        let mut chrs: Vec<ChrInfo> = bam_passa
            .references()
            .iter()
            .filter_map(|n| {
                let len = fasta.chrom_len(n)?;
                let a = ann.chrom(n)?;
                Some(ChrInfo { name: n.clone(), ann: Some(a), len })
            })
            .collect();
        if chrs.is_empty() {
            return Err("no chromosome is shared by the BAM, the reference and the annotation".into());
        }
        chrs.sort_by(|a, b| b.len.cmp(&a.len).then_with(|| a.name.cmp(&b.name)));
        let meta = chrs.iter().map(|_| Mutex::new(ChrMeta::default())).collect();
        let reads_out = match &cfg.emit_read_assignments {
            Some(p) => {
                let f = std::fs::File::create(p).map_err(err(&p.display().to_string()))?;
                let mut w = std::io::BufWriter::new(f);
                writeln!(w, "read_id\tchr\tgene_ids\tassignment_type\tcompatible_transcripts\tunique_transcript\tcounted_weight\tgroup")
                    .map_err(err(&p.display().to_string()))?;
                Some(Mutex::new(w))
            }
            None => None,
        };
        let phase = match (vcf_job, &cfg.phase) {
            (Some(j), Some(pc)) => {
                let (sites, vcf_s) = j.join().map_err(|_| "phase VCF reader panicked".to_string())??;
                if cfg.verbose || cfg.profile {
                    eprintln!(
                        "[iqrs] phase VCF {:.2}s: sample {}, {} phased het SNVs on {} contigs",
                        vcf_s,
                        sites.sample,
                        sites.n_sites(),
                        sites.chrs.len()
                    );
                }
                Some(PhaseState {
                    sites,
                    cfg: iq_phase::Cfg { burnin: pc.burnin, sweeps: pc.sweeps, ..Default::default() },
                    seed: pc.seed,
                    loci: Mutex::new(Vec::new()),
                    vcf_s,
                    ev_ns: AtomicU64::new(0),
                    ev_records: AtomicU64::new(0),
                    ev_reads: AtomicU64::new(0),
                    ev_calls: AtomicU64::new(0),
                    loci_ns: AtomicU64::new(0),
                })
            }
            _ => None,
        };
        Ok(RealPipeline {
            cfg: cfg.clone(),
            params,
            ann,
            fasta,
            bam_passa,
            bam_collect,
            bam_big,
            chrs,
            meta,
            groups: tag.map(|_| SharedGroups::new()),
            cmd_line: String::new(),
            reads_out,
            phase,
        })
    }

    fn block_index(&self, chr: usize, block: &Block) -> IsoformIndex {
        let c = &self.chrs[chr];
        let delta = self.params.matching.delta;
        if block.gene_ids.is_empty() {
            iq_annot::index_from_region(c.ann.unwrap_or(0), block.region.s as i32, block.region.e as i32, delta)
        } else {
            iq_annot::index_from_genes(&self.ann, &block.gene_ids, delta)
        }
    }

    /// Reference window (1-based) for block `bi`: its cluster ∪ its index span, ±[`REF_MARGIN`]
    /// (U7 strand detection can read isoform introns outside the read span).
    fn ref_window(&self, chr: usize, bi: usize, index: &IsoformIndex) -> Iv {
        let cl = self.meta[chr].lock().unwrap().cluster_of_block[bi];
        let len = self.chrs[chr].len as i64;
        let (mut lo, mut hi) = (cl.s + 1, cl.e + 1);
        if index.start <= index.end && index.start > 0 {
            lo = lo.min(index.start as i64);
            hi = hi.max(index.end as i64);
        }
        let s = (lo - REF_MARGIN).max(1);
        let e = (hi + REF_MARGIN).min(len.max(1));
        Iv { s: s as i32, e: e as i32 }
    }
}

impl RealPipeline {
    /// Phasing: samples every locus and writes the two phased tables.
    fn finish_phase(&self, p: &PhaseState, spec: &OutputSpec) -> Result<(), String> {
        let mut per_chr = std::mem::take(&mut *p.loci.lock().unwrap_or_else(|e| e.into_inner()));
        per_chr.sort_by_key(|(c, _)| *c);
        let loci: Vec<iq_phase::Locus> = per_chr.into_iter().flat_map(|(_, l)| l).collect();
        let (res, st) = iq_phase::run_loci(&loci, &p.cfg, p.seed);
        let t = std::time::Instant::now();
        iq_phase::write_outputs(&spec.dir, &spec.prefix, &loci, &res).map_err(err("write phased outputs"))?;
        if self.cfg.verbose || self.cfg.profile {
            let sweeps = (p.cfg.burnin + p.cfg.sweeps) as f64;
            eprintln!(
                "[iqrs] phasing: vcf {:.2}s; het evidence {:.2}s cpu over {} records ({} reads with calls, {} calls); \
                 loci build {:.2}s cpu; sampling {:.2}s wall / {:.2}s cpu ({} + {} sweeps, {:.3} ms cpu per sweep \
                 genome-wide); write {:.2}s; {} loci, {} genes, {} transcripts, {} reads, {} classes",
                p.vcf_s,
                p.ev_ns.load(Ordering::Relaxed) as f64 * 1e-9,
                p.ev_records.load(Ordering::Relaxed),
                p.ev_reads.load(Ordering::Relaxed),
                p.ev_calls.load(Ordering::Relaxed),
                p.loci_ns.load(Ordering::Relaxed) as f64 * 1e-9,
                st.wall_s,
                st.cpu_s,
                p.cfg.burnin,
                p.cfg.sweeps,
                st.cpu_s * 1e3 / sweeps,
                t.elapsed().as_secs_f64(),
                st.loci,
                st.genes,
                st.transcripts,
                st.reads,
                st.classes
            );
            for &(sec, i) in &st.slowest {
                let l = &loci[i];
                eprintln!(
                    "[iqrs] phasing: slow locus {} {:.2}s: {} reads, {} transcripts, {} classes, {} class-transcript pairs",
                    l.name(),
                    sec,
                    l.n_reads(),
                    l.tx.len(),
                    l.n_classes(),
                    l.cls_tx.len()
                );
            }
        }
        Ok(())
    }

    /// `--emit-read-assignments` rows of one block: model indices resolved to final ids.
    /// `counted_weight` is the read's `unique_only` weight (what the counter added per
    /// transcript); `group` the read-group value (NA without one).
    fn format_reads(&self, chr: usize, o: &RealBlockOut, names: &iq_quant::BlockNames, out: &mut String) {
        use std::fmt::Write as _;
        let pos: std::collections::HashMap<u32, usize> = o.models.iter().enumerate().map(|(i, m)| (m.key, i)).collect();
        let pool =
            self.groups.as_ref().map(|g| g.0.read().unwrap_or_else(|e| e.into_inner()).clone()).unwrap_or_default();
        for (name, group, fa) in &o.reads {
            let mut tx: Vec<&str> = Vec::new();
            let mut genes: Vec<&str> = Vec::new();
            for m in &fa.tx {
                if let Some(&i) = pos.get(m) {
                    if !tx.contains(&names.tx[i].as_str()) {
                        tx.push(&names.tx[i]);
                    }
                    if !genes.contains(&names.gene[i].as_str()) {
                        genes.push(&names.gene[i]);
                    }
                }
            }
            let fate = iq_quant::counter::read_fate(iq_quant::Strategy::UniqueOnly, fa.atype, !tx.is_empty(), tx.len());
            let w = match fate {
                iq_quant::ReadFate::Unique => 1.0,
                iq_quant::ReadFate::Ambiguous { weight_one: true } => 1.0,
                _ => 0.0,
            };
            let uniq = if w > 0.0 && tx.len() == 1 { tx[0] } else { "." };
            let _ = writeln!(
                out,
                "{name}\t{}\t{}\t{}\t{}\t{uniq}\t{w}\t{}",
                self.chrs[chr].name,
                if genes.is_empty() { ".".to_string() } else { genes.join(",") },
                fa.atype,
                if tx.is_empty() { ".".to_string() } else { tx.join(",") },
                iq_quant::output::group_name(&pool, *group),
            );
        }
    }
}

impl Pipeline for RealPipeline {
    type Store = RealStore;
    type UnitState = DetectedKnown;
    type Global = RealGlobal;
    type BlockOut = RealBlockOut;
    type ChrOut = ChrOutput;

    fn n_chromosomes(&self) -> usize {
        self.chrs.len()
    }
    fn chr_name(&self, chr: usize) -> &str {
        &self.chrs[chr].name
    }
    fn chr_weight(&self, chr: usize) -> u64 {
        self.chrs[chr].len
    }
    fn est_pass_a(&self, _chr: usize) -> u64 {
        // vtg_proto: 85 MB peak for a 17.8 M-record chromosome slice (mostly the dup filter)
        96 << 20
    }

    fn pass_a(&self, chr: usize) -> Result<Vec<Block>, String> {
        let c = &self.chrs[chr];
        let f = &self.params.filter;
        let cfg = PassAConfig::for_chr(&c.name, c.len, f.max_coverage_small_chr, f.max_coverage_normal_chr);
        let mut pa = PassA::new(AnnOracle { ann: &self.ann, chr: c.ann }, cfg);
        let region = Iv0 { s: 0, e: c.len as i64 - 1 };
        let it = iq_bam::try_fetch(
            &self.bam_passa,
            &c.name,
            region,
            FetchOpts { seq: false, batch_records: 4096, qual: false },
        )
        .map_err(err(&c.name))?;
        let mut last_pos0 = i64::MIN;
        for b in it {
            let b = b.map_err(err(&c.name))?;
            for r in b.iter() {
                if r.pos0 < last_pos0 {
                    return Err(format!(
                        "{}: records out of order on {} near position {}; the BAM must be coordinate-sorted \
                         (`samtools sort` then `samtools index`)",
                        self.cfg.bam.display(),
                        c.name,
                        last_pos0 + 1
                    ));
                }
                last_pos0 = r.pos0;
                pa.push_lite(&r);
            }
        }
        let res = pa.finish();
        let mut meta = ChrMeta { dup_candidates: res.dup_candidates, ..Default::default() };
        let mut blocks = Vec::new();
        for cl in res.clusters {
            for b in cl.blocks {
                meta.regions.push(b.region);
                meta.cluster_of_block.push(cl.region);
                blocks.push(b);
            }
        }
        *self.meta[chr].lock().unwrap() = meta;
        Ok(blocks)
    }

    fn collect_block(&self, chr: usize, bi: usize, block: &Block) -> Result<RealStore, String> {
        let c = &self.chrs[chr];
        let index = self.block_index(chr, block);
        let win = self.ref_window(chr, bi, &index);
        let seq = self.fasta.slice(&c.name, win).ok_or_else(|| format!("{}: not in reference", c.name))?;
        let reference = RefSlice { start: win.s, seq: &seq };
        let regions = self.meta[chr].lock().unwrap().regions.clone();
        let owner = ClusterOwner { regions: &regions, me: bi };
        let groups: &dyn GroupIds = match &self.groups {
            Some(g) => g,
            None => &NoGroups,
        };
        let mut col = BlockCollector::new(block, &index, &reference, &self.params, &owner, groups);
        let bam = if block.owned_reads >= BIG_BLOCK_READS { &self.bam_big } else { &self.bam_collect };
        // phasing: het calls of every record, inline in this read pass
        let sites = self.phase.as_ref().and_then(|p| p.sites.get(&c.name));
        let mut ev = sites.map(|_| (iq_phase::EvidenceBuf::new(), 0u64, 0u64));
        let opts = FetchOpts { seq: true, batch_records: 4096, qual: ev.is_some() };
        let it = iq_bam::try_fetch(bam, &c.name, block.region, opts).map_err(err(&c.name))?;
        for b in it {
            let b = b.map_err(err(&c.name))?;
            let recs: Vec<_> = b.iter().collect();
            col.push(&recs);
            if let (Some((buf, ns, n)), Some(sites)) = (ev.as_mut(), sites) {
                let t = std::time::Instant::now();
                for (i, r) in recs.iter().enumerate() {
                    buf.push(sites, r, b.qual(i));
                }
                *n += recs.len() as u64;
                *ns += t.elapsed().as_nanos() as u64;
            }
        }
        let (store, stats, shared) = col.finish_shared();
        let phase_ev = ev.map(|(buf, ns, n)| {
            let t = std::time::Instant::now();
            let keep: rustc_hash::FxHashSet<u64> = (0..store.len()).map(|i| store.view(i).name_hash()).collect();
            let e = buf.finish(|h| keep.contains(&h));
            let p = self.phase.as_ref().unwrap();
            p.ev_ns.fetch_add(ns + t.elapsed().as_nanos() as u64, Ordering::Relaxed);
            p.ev_records.fetch_add(n, Ordering::Relaxed);
            p.ev_reads.fetch_add(e.n_reads() as u64, Ordering::Relaxed);
            p.ev_calls.fetch_add(e.code.len() as u64, Ordering::Relaxed);
            e
        });
        Ok(RealStore { store, index, stats, shared, phase_ev })
    }

    fn store_bytes(&self, s: &RealStore) -> u64 {
        // length-based resident bytes (vtg_proto #9: heap_bytes overstates by 1.35-2x)
        s.store.len_bytes() as u64 + s.phase_ev.as_ref().map_or(0, |e| e.bytes() as u64)
    }

    fn barrier_chr(&self, chr: usize, _blocks: &[Block], stores: &mut [RealStore]) -> Result<ReadTotals, String> {
        // U9 candidates: pass-A duplicate names ∪ reads collected by two blocks (4.1 double
        // processing of block-straddling reads; such a record is not a pass-A duplicate)
        let mut cands = std::mem::take(&mut self.meta[chr].lock().unwrap().dup_candidates);
        for s in stores.iter_mut() {
            cands.append(&mut s.shared);
        }
        cands.sort_unstable();
        cands.dedup();
        // U9: suspended copies are removed from their stores (rows re-indexed densely)
        let mut bs: Vec<BlockStore> = stores.iter_mut().map(|s| std::mem::take(&mut s.store)).collect();
        let o = iq_collect::multimap::resolve_chr_duplicates(&mut bs, Some(&cands), self.params.multimap_strategy);
        for (s, b) in stores.iter_mut().zip(bs) {
            s.store = b;
        }
        Ok(ReadTotals { total: o.total, polya: o.polya })
    }

    fn decide(&self, t: ReadTotals) -> RealGlobal {
        let frac = if t.total > 0 { t.polya as f64 / t.total as f64 } else { 0.0 };
        let mut p = self.params.clone();
        let set = |flag: bool| match p.polya_requirement {
            PolyAUsage::Auto => flag,
            PolyAUsage::Never => false,
            _ => true,
        };
        // dataset_processor.py:283-300
        let req = set(frac >= p.polya_percentage_threshold);
        let mono_int = set(p.model.require_monointronic_polya || req);
        let mono_ex = set(p.model.require_monoexonic_polya || req);
        p.requires_polya_for_construction = req;
        p.model.require_monointronic_polya = mono_int;
        p.model.require_monoexonic_polya = mono_ex;
        if self.cfg.verbose || self.cfg.profile {
            eprintln!(
                "[iqrs] assignments {}, polyA {} ({:.1}%); requires_polya_for_construction = {req}",
                t.total,
                t.polya,
                frac * 100.0
            );
        }
        if frac < p.low_polya_percentage_threshold && p.polya_requirement != PolyAUsage::Never {
            eprintln!("[iqrs] warning: polyA percentage is suspiciously low ({:.1}%)", frac * 100.0);
        }
        RealGlobal { params: p, polya_fraction: frac }
    }

    fn construct_block(
        &self,
        chr: usize,
        bi: usize,
        _block: &Block,
        s: RealStore,
        dk: &mut DetectedKnown,
        g: &RealGlobal,
    ) -> Result<RealBlockOut, String> {
        let c = &self.chrs[chr];
        let win = self.ref_window(chr, bi, &s.index);
        let seq = self.fasta.slice(&c.name, win).ok_or_else(|| format!("{}: not in reference", c.name))?;
        let reference = RefSlice { start: win.s, seq: &seq };
        let emit = self.reads_out.is_some();
        let phase = self.phase.as_ref().map(|_| s.phase_ev.as_ref());
        let bm =
            crate::construct::construct_block(&self.ann, &g.params, &s.index, &s.store, reference, dk, emit, phase);
        Ok(RealBlockOut { index: s.index, models: bm.models, counter: bm.counter, reads: bm.reads, phase: bm.phase })
    }

    fn merge_chr(&self, chr: usize, _blocks: &[Block], outs: Vec<RealBlockOut>) -> Result<ChrOutput, String> {
        let mut w = ChrPhaseW::new(&self.ann, &self.chrs[chr].name);
        let mut text = String::new();
        let mut acc = self.phase.as_ref().map(|_| iq_phase::ChrAcc::new());
        let mut phase_ns = 0u64;
        for o in &outs {
            let names = w.add_block(&BlockOutput { index: &o.index, models: &o.models, counter: &o.counter });
            if self.reads_out.is_some() {
                self.format_reads(chr, o, &names, &mut text);
            }
            if let (Some(acc), Some(pr)) = (acc.as_mut(), o.phase.as_ref()) {
                let t = std::time::Instant::now();
                let pos: rustc_hash::FxHashMap<u32, usize> =
                    o.models.iter().enumerate().map(|(i, m)| (m.key, i)).collect();
                for r in 0..pr.len() {
                    let tx = pr
                        .models_of(r)
                        .iter()
                        .filter_map(|m| pos.get(m).map(|&i| (names.tx[i].as_str(), names.gene[i].as_str())));
                    acc.add_read(Some(pr.hash[r]), tx, pr.d[r]);
                }
                phase_ns += t.elapsed().as_nanos() as u64;
            }
        }
        if let (Some(acc), Some(p)) = (acc, self.phase.as_ref()) {
            let t = std::time::Instant::now();
            let loci = acc.build(&self.chrs[chr].name);
            p.loci.lock().unwrap_or_else(|e| e.into_inner()).push((chr, loci));
            p.loci_ns.fetch_add(phase_ns + t.elapsed().as_nanos() as u64, Ordering::Relaxed);
        }
        if let Some(out) = &self.reads_out {
            let mut f = out.lock().unwrap_or_else(|e| e.into_inner());
            f.write_all(text.as_bytes()).map_err(err("emit read assignments"))?;
        }
        Ok(w.finish())
    }

    fn finish(&self, outs: Vec<ChrOutput>, _t: ReadTotals) -> Result<(), String> {
        if let Some(out) = &self.reads_out {
            out.lock().unwrap_or_else(|e| e.into_inner()).flush().map_err(err("emit read assignments"))?;
        }
        let spec = OutputSpec {
            dir: self.cfg.out_dir.clone(),
            prefix: self.cfg.prefix.clone(),
            cmd_line: self.cmd_line.clone(),
            version: format!("4.1.0-iqrs-{}", env!("CARGO_PKG_VERSION")),
            group_tag: self.cfg.read_group_tag.clone(),
            grouped_out: self.cfg.grouped_out.clone(),
            gzip: false,
        };
        std::fs::create_dir_all(&spec.dir).map_err(err(&spec.dir.display().to_string()))?;
        let pool = match &self.groups {
            Some(g) => g.0.read().unwrap_or_else(|e| e.into_inner()).clone(),
            None => Default::default(),
        };
        iq_quant::write_outputs(&spec, &outs, &pool, self.bam_collect.unmapped()).map_err(err("write outputs"))?;
        if let Some(p) = &self.phase {
            self.finish_phase(p, &spec)?;
        }
        Ok(())
    }

    fn fixed_bytes(&self) -> u64 {
        // annotation + per-thread scratch (not measured precisely)
        self.ann.heap_bytes() as u64 + self.cfg.threads as u64 * (16 << 20)
    }
}
