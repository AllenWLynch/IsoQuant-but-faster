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
}

pub struct RealStore {
    pub store: BlockStore,
    pub index: IsoformIndex,
    pub stats: CollectStats,
    /// Name hashes of reads this block shares with a neighbouring block (4.1 double
    /// processing of valley-straddling reads); U9 candidates at the barrier.
    pub shared: Vec<u64>,
}

pub struct RealGlobal {
    pub params: Params,
    pub polya_fraction: f64,
}

pub struct RealBlockOut {
    index: IsoformIndex,
    models: Vec<TranscriptModel>,
    counter: TxCounter,
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
        let it = iq_bam::try_fetch(&self.bam_passa, &c.name, region, FetchOpts { seq: false, batch_records: 4096 })
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
        let it = iq_bam::try_fetch(bam, &c.name, block.region, FetchOpts { seq: true, batch_records: 4096 })
            .map_err(err(&c.name))?;
        for b in it {
            let b = b.map_err(err(&c.name))?;
            let recs: Vec<_> = b.iter().collect();
            col.push(&recs);
        }
        let (store, stats, shared) = col.finish_shared();
        Ok(RealStore { store, index, stats, shared })
    }

    fn store_bytes(&self, s: &RealStore) -> u64 {
        // length-based resident bytes (vtg_proto #9: heap_bytes overstates by 1.35-2x)
        s.store.len_bytes() as u64
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
        let bm = crate::construct::construct_block(&self.ann, &g.params, &s.index, &s.store, reference, dk);
        Ok(RealBlockOut { index: s.index, models: bm.models, counter: bm.counter })
    }

    fn merge_chr(&self, chr: usize, _blocks: &[Block], outs: Vec<RealBlockOut>) -> Result<ChrOutput, String> {
        let mut w = ChrPhaseW::new(&self.ann, &self.chrs[chr].name);
        for o in &outs {
            w.add_block(&BlockOutput { index: &o.index, models: &o.models, counter: &o.counter });
        }
        Ok(w.finish())
    }

    fn finish(&self, outs: Vec<ChrOutput>, _t: ReadTotals) -> Result<(), String> {
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
        Ok(())
    }

    fn fixed_bytes(&self) -> u64 {
        // annotation + per-thread scratch (not measured precisely)
        self.ann.heap_bytes() as u64 + self.cfg.threads as u64 * (16 << 20)
    }
}
