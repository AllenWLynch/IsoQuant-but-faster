//! Fixture runners: U8 `collect_block` and U9 `resolve`.
//!
//! U8 runs from real BAM records. It fetches the block region from the capture's slice BAM
//! (`<capdir>/slice/*.bam`, HP tag decoded) and builds the block index live from the GFF3 (`iq_annot`,
//! D-R1 gate in `genes_for_region0`; `IQRS_U8_INDEX=dump` uses the captured dump instead).
//! The reference is the whole chromosome from `data/ref/genome.fa` (`$IQRS_REFERENCE`).
//! Ownership (4.1 fetch window, as in both builds since D-R2 was retired) is decided among
//! all U8 block regions of the same chromosome in the capture. `assignment_id` (a 4.1 global counter, also advanced by dropped reads) is
//! not compared.

use crate::multimap::{resolve, BasicRead};
use crate::{BlockCollector, ClusterOwner, CollectStats, FnGroups, StageTimes};
use iq_annot::fixture::{cached_annotation, gff3_path, remap_ref_index};
use iq_annot::{index_from_genes, index_from_region, GeneTypes};
use iq_bam::{BamHandle, BamOptions, FastaIndex};
use iq_types::fixture::schema::{
    BasicAssignmentFixture, ReadRecordFixture, U8CollectBlockIn, U8CollectBlockOut, U9ResolveIn, U9ResolveOut,
};
use iq_types::fixture::{
    diff_values, read_lines, unit_files, Diff, FixError, FixResult, FixtureCtx, ReadRecordExtras, TxSpace, UnitFixture,
    Verdict,
};
use iq_types::serde::Deserialize;
use iq_types::{Block, BlockStore, DevId, IsoformIndex, Iv0, Params, RecordLite, RefSlice};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, OnceLock};

// ------------------------------------------------------------------------------------------
// shared resources (cached per process)
// ------------------------------------------------------------------------------------------

fn cache<K: std::hash::Hash + Eq + Clone, V: Clone>(
    m: &'static OnceLock<Mutex<HashMap<K, V>>>,
    k: K,
    f: impl FnOnce() -> FixResult<V>,
) -> FixResult<V> {
    let m = m.get_or_init(Default::default);
    if let Some(v) = m.lock().unwrap_or_else(|e| e.into_inner()).get(&k) {
        return Ok(v.clone());
    }
    let v = f()?;
    m.lock().unwrap_or_else(|e| e.into_inner()).insert(k, v.clone());
    Ok(v)
}

/// Reference genome: `$IQRS_REFERENCE` or the first `data/ref/genome.fa` above `root`.
pub fn genome_path(root: &Path) -> FixResult<PathBuf> {
    if let Some(p) = std::env::var_os("IQRS_REFERENCE") {
        return Ok(PathBuf::from(p));
    }
    let root = root.canonicalize().unwrap_or_else(|_| root.to_owned());
    for a in root.ancestors() {
        let p = a.join("data/ref/genome.fa");
        if p.is_file() {
            return Ok(p);
        }
    }
    Err(FixError("U8: reference genome not found (set IQRS_REFERENCE)".into()))
}

/// Whole-chromosome sequence (1-based start 1), cached.
pub fn chrom_seq(fa: &Path, chr: &str) -> FixResult<Arc<Vec<u8>>> {
    static C: OnceLock<Mutex<HashMap<(PathBuf, String), Arc<Vec<u8>>>>> = OnceLock::new();
    cache(&C, (fa.to_owned(), chr.to_owned()), || {
        let idx = FastaIndex::open(fa).map_err(|e| FixError(e.to_string()))?;
        let len = idx.chrom_len(chr).ok_or_else(|| FixError(format!("{chr} not in {}", fa.display())))?;
        let s = idx.slice(chr, iq_types::Iv::new(1, len as i32)).unwrap();
        Ok(Arc::new(s.into_owned()))
    })
}

/// The capture's slice BAM, opened with the read-group tag of `params.json`.
pub fn slice_bam(ctx: &FixtureCtx) -> FixResult<Arc<BamHandle>> {
    static C: OnceLock<Mutex<HashMap<PathBuf, Arc<BamHandle>>>> = OnceLock::new();
    let dir = ctx.root.join("slice");
    let pj = ctx.params_json()?;
    cache(&C, dir.clone(), || {
        let mut bams: Vec<PathBuf> = std::fs::read_dir(&dir)
            .map_err(|e| FixError(format!("{}: {e}", dir.display())))?
            .filter_map(|e| e.ok().map(|e| e.path()))
            .filter(|p| p.extension().is_some_and(|x| x == "bam"))
            .collect();
        bams.sort();
        let bam = bams.first().ok_or_else(|| FixError(format!("no BAM in {}", dir.display())))?;
        let tag = read_group_tag(&pj);
        let h = BamHandle::open_with(bam, BamOptions { tag, threads: 0 }).map_err(|e| FixError(e.to_string()))?;
        Ok(Arc::new(h))
    })
}

/// `--read_group tag:XX` → `XX` (other groupers are out of scope).
pub fn read_group_tag(pj: &iq_types::serde_json::Map<String, iq_types::serde_json::Value>) -> Option<[u8; 2]> {
    let v = pj.get("read_group")?;
    let s = match v {
        iq_types::serde_json::Value::Array(a) => a.first()?.as_str()?,
        iq_types::serde_json::Value::String(s) => s.as_str(),
        _ => return None,
    };
    let t = s.strip_prefix("tag:")?.as_bytes();
    (t.len() == 2).then(|| [t[0], t[1]])
}

#[derive(Deserialize)]
#[serde(crate = "iq_types::serde")]
struct RegionOnly {
    #[serde(rename = "in")]
    input: RegionIn,
}
#[derive(Deserialize)]
#[serde(crate = "iq_types::serde")]
struct RegionIn {
    chr: String,
    block_no: u32,
    region: Iv0,
}

/// All U8 block regions of `chr` in the capture, sorted (block order).
fn sibling_regions(ctx: &FixtureCtx, chr: &str) -> FixResult<Arc<Vec<(u32, Iv0)>>> {
    static C: OnceLock<Mutex<HashMap<(PathBuf, String), Arc<Vec<(u32, Iv0)>>>>> = OnceLock::new();
    cache(&C, (ctx.root.clone(), chr.to_owned()), || {
        let mut v = Vec::new();
        for f in unit_files(&ctx.root, "U8")? {
            for c in read_lines::<RegionOnly>(&f)? {
                let c = c?;
                if c.input.chr == chr {
                    v.push((c.input.block_no, c.input.region));
                }
            }
        }
        v.sort_by_key(|x| x.1);
        Ok(Arc::new(v))
    })
}

/// Block index: live from the annotation (default) or the captured dump.
pub fn block_index(ctx: &FixtureCtx, chr: &str, region: Iv0, dump_sha: &str, delta: i32) -> FixResult<Arc<IsoformIndex>> {
    if std::env::var("IQRS_U8_INDEX").is_ok_and(|v| v == "dump") {
        return ctx.index(dump_sha);
    }
    let ann = cached_annotation(&gff3_path(ctx)?, GeneTypes::Isoquant)?;
    let genes = match ann.chrom(chr) {
        Some(c) => ann.genes_for_region0(c, region),
        None => Default::default(),
    };
    let idx = if genes.is_empty() {
        index_from_region(0, region.s as i32, region.e as i32, delta)
    } else {
        index_from_genes(&ann, &genes, delta)
    };
    Ok(Arc::new(remap_ref_index(idx, &ann, chr, &mut ctx.ids())))
}

/// Everything a U8 case needs, loaded outside any timing: fetched batches, index, params.
pub struct U8Prepared {
    pub input: U8CollectBlockIn,
    pub params: Params,
    pub index: Arc<IsoformIndex>,
    pub seq: Arc<Vec<u8>>,
    pub regions: Vec<Iv0>,
    pub me: usize,
    pub batches: Vec<iq_bam::RecordBatch>,
}

impl U8Prepared {
    pub fn load(i: &U8CollectBlockIn, ctx: &FixtureCtx) -> FixResult<Self> {
        let params: Params = ctx.params()?;
        let bam = slice_bam(ctx)?;
        let seq = chrom_seq(&genome_path(&ctx.root)?, &i.chr)?;
        let index = block_index(ctx, &i.chr, i.region, &i.index, params.matching.delta)?;
        let sibs = sibling_regions(ctx, &i.chr)?;
        let regions: Vec<Iv0> = sibs.iter().map(|x| x.1).collect();
        let me = sibs
            .iter()
            .position(|x| x.0 == i.block_no && x.1 == i.region)
            .ok_or_else(|| FixError("block not among U8 regions".into()))?;
        // the 4.1 fetch window
        let batches = iq_bam::fetch(&bam, &i.chr, i.region).collect();
        Ok(U8Prepared { input: i.clone(), params, index, seq, regions, me, batches })
    }

    pub fn n_records(&self) -> usize {
        self.batches.iter().map(|b| b.len()).sum()
    }

    /// Collects the block; returns the store, stats and (if `timed`) stage times.
    pub fn collect(&self, ctx: &FixtureCtx, parallel: bool, timed: bool) -> (BlockStore, CollectStats, Option<StageTimes>) {
        let i = &self.input;
        let reference = RefSlice { start: 1, seq: &self.seq };
        let owner = ClusterOwner { regions: &self.regions, me: self.me };
        let block =
            Block { region: i.region, max_cov: i.max_coverage as u32, skip: i.skip_read_fraction, gene_ids: Vec::new(), owned_reads: 0 };
        let groups = FnGroups(|t: Option<&iq_types::TagValue>| {
            ctx.ids().group_id(&[iq_read::groups::group_string(t).into_owned()])
        });
        let mut c = BlockCollector::new(&block, &self.index, &reference, &self.params, &owner, &groups).parallel(parallel);
        if timed {
            c = c.timed();
        }
        let mut lites: Vec<RecordLite<'_>> = Vec::new();
        for b in &self.batches {
            lites.clear();
            lites.extend(b.iter());
            c.push(&lites);
        }
        let times = c.times();
        let (store, stats) = c.finish();
        (store, stats, times)
    }

    /// Fixture DTOs of a collected store.
    pub fn to_fixture(&self, store: &BlockStore, ctx: &FixtureCtx) -> Vec<ReadRecordFixture> {
        let ids = ctx.ids();
        let ex = ReadRecordExtras { chr: self.input.chr.clone(), ..Default::default() };
        (0..store.len()).map(|k| ReadRecordFixture::from_core(&store.record(k), &ids, &ex)).collect()
    }
}

// ------------------------------------------------------------------------------------------
// U8
// ------------------------------------------------------------------------------------------

pub struct U8Collect;

impl UnitFixture for U8Collect {
    const UNIT: &'static str = "U8";
    const FN: &'static str = "collect_block";
    const DEVIATIONS: &'static [DevId] = &[DevId::R1, DevId::P1, DevId::P3, DevId::A3, DevId::A4];
    type In = U8CollectBlockIn;
    type Out = U8CollectBlockOut;

    fn run(i: &U8CollectBlockIn, ctx: &FixtureCtx) -> FixResult<U8CollectBlockOut> {
        let par = !std::env::var("IQRS_U8_SERIAL").is_ok_and(|v| v == "1");
        let p = U8Prepared::load(i, ctx)?;
        let (store, _, _) = p.collect(ctx, par, false);
        Ok(U8CollectBlockOut { reads: p.to_fixture(&store, ctx) })
    }

    fn compare(got: &U8CollectBlockOut, want: &U8CollectBlockOut, ctx: &FixtureCtx) -> Verdict {
        let mut w = want.clone();
        for (k, r) in w.reads.iter_mut().enumerate() {
            r.assignment_id = got.reads.get(k).map_or(0, |g| g.assignment_id);
        }
        if *got == w {
            return Verdict::Exact;
        }
        match diff_values(got, &w, &ctx.diff_opts) {
            Some(d) => Verdict::Mismatch(d),
            None => Verdict::Mismatch(Diff::error("PartialEq differs")),
        }
    }
}

// ------------------------------------------------------------------------------------------
// U9
// ------------------------------------------------------------------------------------------

pub struct U9Resolve;

impl UnitFixture for U9Resolve {
    const UNIT: &'static str = "U9";
    const FN: &'static str = "resolve";
    const DEVIATIONS: &'static [DevId] = &[DevId::C1];
    type In = U9ResolveIn;
    type Out = U9ResolveOut;

    fn run(i: &U9ResolveIn, ctx: &FixtureCtx) -> FixResult<U9ResolveOut> {
        let mut group: Vec<BasicRead> = {
            let mut ids = ctx.ids();
            i.assignments
                .iter()
                .enumerate()
                .map(|(k, a)| {
                    let mut genes: Vec<_> = a.genes.iter().map(|g| ids.gene_key(g)).collect();
                    genes.sort_unstable();
                    genes.dedup();
                    let mut isoforms: Vec<_> = a.isoforms.iter().map(|t| ids.iso_ref(t, TxSpace::Ref)).collect();
                    isoforms.sort_unstable();
                    isoforms.dedup();
                    BasicRead {
                        block: 0,
                        row: k as u32,
                        chr: ids.chroms.intern(&a.chr),
                        start: a.start,
                        end: a.end,
                        genomic_region: a.region,
                        multimapper: a.multimapper,
                        polya_found: a.polya_found,
                        atype: a.atype,
                        gene_atype: a.gene_atype,
                        penalty: a.penalty,
                        genes,
                        isoforms,
                    }
                })
                .collect()
        };
        resolve(&mut group, i.strategy);
        let out = i
            .assignments
            .iter()
            .zip(&group)
            .map(|(a, r)| BasicAssignmentFixture {
                multimapper: r.multimapper,
                atype: r.atype,
                gene_atype: r.gene_atype,
                ..a.clone()
            })
            .collect();
        Ok(U9ResolveOut { assignments: out })
    }
}

pub(crate) fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    reg.add::<U8Collect>().add::<U9Resolve>();
}
