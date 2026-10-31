//! U3 fixture runners.
//!
//! - `split_coverage_regions`: direct call on the recorded coverage dict.
//! - `cluster`, `process_alignments_in_region`: the captured cases carry no records, so the
//!   runners stream the capture's own slice BAM (`<capdir>/slice/*.bam`, the exact input of
//!   the 4.1 run) through the real [`PassA`] and look up the cluster/block by region. This
//!   checks cluster boundaries, read counts, block numbering, max coverage, the cap and gene
//!   lists end to end. A minimal sequential BAM decoder lives here (test-only; production
//!   BAM IO is `iq_bam`).
//! - `genes_in_region`: the gene-query double below on the recorded region.
//!
//! Gene queries use a double built from `fixtures/annotation/U1.annotation.jsonl.gz`
//! ([`AnnGenes`]); it mirrors `iq_annot::AnnotationIndex::genes_for_region0`, including the
//! D-R1 gate, so default-mode mismatches can be attributed.

use std::collections::HashMap;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, OnceLock};

use iq_types::fixture::schema::*;
use iq_types::fixture::{open_lines, FixError, FixResult, FixtureCtx, Registry, UnitFixture, Verdict};
use iq_types::serde::Deserialize;
use iq_types::{quirk, Cluster, DevId, GeneIdx, Iv0, RecordLite};

use crate::pass_a::{GeneOracle, PassA, PassAConfig, PassAResult};
use crate::split::{coverage_cutoff, is_small_chr, skip_read_fraction, split_coverage_regions, CoverageBins};

pub(crate) fn register_fixtures(reg: &mut Registry) {
    reg.add::<U3Split>();
    reg.add::<U3Cluster>();
    reg.add::<U3Process>();
    reg.add::<U3Genes>();
}

const DEVS: &[DevId] = &[DevId::R1];

// ---------------------------------------------------------------- gene query double

/// Genes of one chromosome from the U1 annotation dump, in dump order (= sorted by start,
/// file order).
pub(crate) struct AnnGenes {
    ids: Vec<String>,
    iv: Vec<(i64, i64)>,
}

impl GeneOracle for AnnGenes {
    fn genes_for_region0(&self, r: Iv0, out: &mut Vec<GeneIdx>) {
        for (i, &(s, e)) in self.iv.iter().enumerate() {
            let hit = if quirk(DevId::R1) {
                // 4.1: 0-based tuple straight into gffutils region(); a falsy start drops the
                // lower bound and leaves `start < end` (same as iq_annot)
                if r.s == 0 {
                    s < r.e
                } else {
                    s <= r.e && e >= r.s
                }
            } else {
                // 1-based query (s+1, e+1)
                let (qs, qe) = (r.s + 1, r.e + 1);
                s <= qe && e >= qs
            };
            if hit {
                out.push(i as GeneIdx);
            }
        }
    }
}

#[derive(Deserialize)]
#[serde(crate = "iq_types::serde")]
struct AnnLine {
    out: AnnOut,
}
#[derive(Deserialize)]
#[serde(crate = "iq_types::serde")]
struct AnnOut {
    chr: String,
    genes: Vec<AnnGene>,
}
#[derive(Deserialize)]
#[serde(crate = "iq_types::serde")]
struct AnnGene {
    id: String,
    start: i64,
    end: i64,
}

fn annotation_path(ctx: &FixtureCtx) -> PathBuf {
    ctx.root.join("../../annotation/U1.annotation.jsonl.gz")
}

fn ann_genes(ctx: &FixtureCtx, chr: &str) -> FixResult<Arc<AnnGenes>> {
    static CACHE: OnceLock<Mutex<HashMap<(PathBuf, String), Arc<AnnGenes>>>> = OnceLock::new();
    let path = annotation_path(ctx);
    let key = (path.clone(), chr.to_owned());
    let mut g = CACHE.get_or_init(Default::default).lock().unwrap();
    if let Some(a) = g.get(&key) {
        return Ok(a.clone());
    }
    let mut found = None;
    use std::io::BufRead;
    for line in open_lines(&path)?.lines() {
        let line = line?;
        // cheap pre-filter before parsing a whole chromosome's features
        if !line.contains(&format!("\"chr\":\"{chr}\"")) {
            continue;
        }
        let l: AnnLine = iq_types::serde_json::from_str(&line).map_err(|e| FixError(format!("annotation: {e}")))?;
        if l.out.chr == chr {
            found = Some(AnnGenes {
                ids: l.out.genes.iter().map(|g| g.id.clone()).collect(),
                iv: l.out.genes.iter().map(|g| (g.start, g.end)).collect(),
            });
            break;
        }
    }
    let a = Arc::new(found.unwrap_or(AnnGenes { ids: vec![], iv: vec![] }));
    g.insert(key, a.clone());
    Ok(a)
}

// ---------------------------------------------------------------- minimal BAM reader

fn read_u32(r: &mut impl Read) -> std::io::Result<u32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}

/// Streams the records of reference `chr` of a BAM file (BGZF read as multi-member gzip)
/// into `f`. Returns the reference length from the header.
fn for_each_bam_record(path: &Path, chr: &str, mut f: impl FnMut(u64, &RecordLite<'_>)) -> FixResult<u64> {
    let file = std::fs::File::open(path).map_err(|e| FixError(format!("{}: {e}", path.display())))?;
    let mut r = std::io::BufReader::with_capacity(1 << 20, flate2::read::MultiGzDecoder::new(file));
    let mut magic = [0u8; 4];
    r.read_exact(&mut magic)?;
    if &magic != b"BAM\x01" {
        return Err(FixError(format!("{}: not a BAM file", path.display())));
    }
    let l_text = read_u32(&mut r)? as usize;
    std::io::copy(&mut (&mut r).take(l_text as u64), &mut std::io::sink())?;
    let n_ref = read_u32(&mut r)?;
    let mut tid = None;
    let mut chr_len = 0u64;
    for i in 0..n_ref {
        let l_name = read_u32(&mut r)? as usize;
        let mut name = vec![0u8; l_name];
        r.read_exact(&mut name)?;
        let l_ref = read_u32(&mut r)?;
        if &name[..l_name.saturating_sub(1)] == chr.as_bytes() {
            tid = Some(i as i32);
            chr_len = l_ref as u64;
        }
    }
    let Some(tid) = tid else { return Err(FixError(format!("{}: no reference {chr}", path.display()))) };
    let mut buf = Vec::new();
    let mut cigar = Vec::new();
    loop {
        let mut bs = [0u8; 4];
        match r.read_exact(&mut bs) {
            Ok(()) => {}
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
            Err(e) => return Err(e.into()),
        }
        let n = u32::from_le_bytes(bs) as usize;
        buf.resize(n, 0);
        r.read_exact(&mut buf)?;
        let i32_at = |o: usize| i32::from_le_bytes(buf[o..o + 4].try_into().unwrap());
        let u16_at = |o: usize| u16::from_le_bytes(buf[o..o + 2].try_into().unwrap());
        if i32_at(0) != tid {
            continue;
        }
        let pos0 = i32_at(4) as i64;
        let l_read_name = buf[8] as usize;
        let mapq = buf[9];
        let n_cigar = u16_at(12) as usize;
        let flag = u16_at(14);
        let name_end = 32 + l_read_name;
        cigar.clear();
        for k in 0..n_cigar {
            let o = name_end + 4 * k;
            cigar.push(u32::from_le_bytes(buf[o..o + 4].try_into().unwrap()));
        }
        // a `<l_seq>S<ref_len>N` CG placeholder has the true reference length: fine for end0
        let rec = RecordLite {
            qname: &buf[32..name_end - 1],
            flag,
            pos0,
            mapq,
            cigar: &cigar,
            seq: &[],
            tag: None,
        };
        f(chr_len, &rec);
    }
    Ok(chr_len)
}

fn slice_bam(ctx: &FixtureCtx) -> FixResult<PathBuf> {
    let dir = ctx.root.join("slice");
    let mut v: Vec<PathBuf> = std::fs::read_dir(&dir)
        .map_err(|e| FixError(format!("{}: {e}", dir.display())))?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|x| x == "bam"))
        .collect();
    v.sort();
    v.into_iter().next().ok_or_else(|| FixError(format!("{}: no slice BAM", dir.display())))
}

struct ChrPassA {
    chr_len: u64,
    cutoff: i64,
    res: PassAResult,
    genes: Arc<AnnGenes>,
}

/// Pass A over the capture's slice BAM, cached per (capture, chr, active quirks).
fn chr_pass_a(ctx: &FixtureCtx, chr: &str) -> FixResult<Arc<ChrPassA>> {
    type Key = (PathBuf, String, bool);
    static CACHE: OnceLock<Mutex<HashMap<Key, Arc<ChrPassA>>>> = OnceLock::new();
    let key = (ctx.root.clone(), chr.to_owned(), quirk(DevId::R1));
    let mut g = CACHE.get_or_init(Default::default).lock().unwrap();
    if let Some(x) = g.get(&key) {
        return Ok(x.clone());
    }
    let genes = ann_genes(ctx, chr)?;
    let params = ctx.params()?;
    let bam = slice_bam(ctx)?;
    let cutoff_of = |len| coverage_cutoff(chr, len, params.filter.max_coverage_small_chr, params.filter.max_coverage_normal_chr);
    let mut p = None;
    let chr_len = for_each_bam_record(&bam, chr, |len, r| {
        p.get_or_insert_with(|| PassA::new(&*genes, PassAConfig { coverage_cutoff: cutoff_of(len), expected_primary: 0 }))
            .push_lite(r)
    })?;
    let cutoff = cutoff_of(chr_len);
    let res = p.map(|p| p.finish()).unwrap_or_default();
    let x = Arc::new(ChrPassA { chr_len, cutoff, res, genes: genes.clone() });
    g.insert(key, x.clone());
    Ok(x)
}

/// Sentinel for "no such cluster/block in iqrs' pass A".
const MISSING: u32 = u32::MAX;

fn find_cluster(clusters: &[Cluster], region: Iv0) -> Option<(usize, u32)> {
    let mut first = 0u32;
    for (i, c) in clusters.iter().enumerate() {
        if c.region == region {
            return Some((i, first));
        }
        first += c.blocks.len() as u32;
    }
    None
}

// ---------------------------------------------------------------- runners

pub(crate) struct U3Split;
impl UnitFixture for U3Split {
    const UNIT: &'static str = "U3";
    const FN: &'static str = "split_coverage_regions";
    const DEVIATIONS: &'static [DevId] = DEVS;
    type In = U3SplitCoverageIn;
    type Out = U3SplitCoverageOut;
    fn run(i: &U3SplitCoverageIn, _ctx: &FixtureCtx) -> FixResult<U3SplitCoverageOut> {
        let (first, counts) = match (i.coverage.first(), i.coverage.last()) {
            (Some(a), Some(b)) => {
                let mut v = vec![0u32; (b.0 - a.0 + 1) as usize];
                for &(k, c) in i.coverage.iter() {
                    v[(k - a.0) as usize] = c as u32;
                }
                (a.0, v)
            }
            _ => (0, vec![]),
        };
        let splits = split_coverage_regions(i.region, &CoverageBins { first_bin: first, counts: &counts }, i.read_count);
        Ok(U3SplitCoverageOut { splits })
    }
}

pub(crate) struct U3Cluster;
impl UnitFixture for U3Cluster {
    const UNIT: &'static str = "U3";
    const FN: &'static str = "cluster";
    const DEVIATIONS: &'static [DevId] = DEVS;
    type In = U3ClusterIn;
    type Out = U3ClusterOut;
    fn run(i: &U3ClusterIn, ctx: &FixtureCtx) -> FixResult<U3ClusterOut> {
        let pa = chr_pass_a(ctx, &i.chr)?;
        Ok(match find_cluster(&pa.res.clusters, i.region) {
            Some((k, first)) if pa.res.clusters[k].n_records == i.read_count => {
                U3ClusterOut { n_blocks: pa.res.clusters[k].blocks.len() as u32, first_block_no: first }
            }
            _ => U3ClusterOut { n_blocks: 0, first_block_no: MISSING },
        })
    }
}

pub(crate) struct U3Process;
impl UnitFixture for U3Process {
    const UNIT: &'static str = "U3";
    const FN: &'static str = "process_alignments_in_region";
    const DEVIATIONS: &'static [DevId] = DEVS;
    type In = U3ProcessRegionIn;
    type Out = U3ProcessRegionOut;
    fn run(i: &U3ProcessRegionIn, ctx: &FixtureCtx) -> FixResult<U3ProcessRegionOut> {
        let pa = chr_pass_a(ctx, &i.chr)?;
        if pa.chr_len != i.chr_len || is_small_chr(&i.chr, pa.chr_len) != i.small_chr || pa.cutoff != i.coverage_cutoff {
            return Err(FixError(format!(
                "chromosome setup differs: len {} small {} cutoff {}",
                pa.chr_len,
                is_small_chr(&i.chr, pa.chr_len),
                pa.cutoff
            )));
        }
        let found = pa.res.clusters.iter().flat_map(|c| c.blocks.iter()).enumerate().find(|(_, b)| b.region == i.region);
        let (block_no, skip, genes) = match found {
            Some((n, b)) if b.max_cov as u64 == i.max_coverage => (n as u32, b.skip, b.gene_ids.clone()),
            _ => {
                let mut g = Vec::new();
                pa.genes.genes_for_region0(i.region, &mut g);
                (MISSING, skip_read_fraction(i.max_coverage, i.coverage_cutoff), g)
            }
        };
        let gene_ids: Vec<String> = genes.iter().map(|&g| pa.genes.ids[g as usize].clone()).collect();
        Ok(U3ProcessRegionOut {
            block_no,
            skip_read_fraction: skip,
            genic: !gene_ids.is_empty(),
            gene_ids,
            index: String::new(),
        })
    }
    /// `index` (the IsoformIndex sha1) is U1's output; not judged here.
    fn compare(got: &U3ProcessRegionOut, want: &U3ProcessRegionOut, ctx: &FixtureCtx) -> Verdict {
        let mut g = got.clone();
        g.index = want.index.clone();
        if g == *want {
            return Verdict::Exact;
        }
        match iq_types::fixture::diff_values(&g, want, &ctx.diff_opts) {
            Some(d) => Verdict::Mismatch(d),
            None => Verdict::Mismatch(iq_types::fixture::Diff::error("PartialEq differs")),
        }
    }
}

pub(crate) struct U3Genes;
impl UnitFixture for U3Genes {
    const UNIT: &'static str = "U3";
    const FN: &'static str = "genes_in_region";
    const DEVIATIONS: &'static [DevId] = DEVS;
    type In = U3GenesInRegionIn;
    type Out = GeneIdsOut;
    fn run(i: &U3GenesInRegionIn, ctx: &FixtureCtx) -> FixResult<GeneIdsOut> {
        let genes = ann_genes(ctx, &i.chr)?;
        let mut g = Vec::new();
        genes.genes_for_region0(i.region, &mut g);
        let mut gene_ids: Vec<String> = g.iter().map(|&k| genes.ids[k as usize].clone()).collect();
        gene_ids.sort();
        Ok(GeneIdsOut { gene_ids })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pass_a::{pass_a, PassARecord};
    use iq_types::{with_quirks, QuirkSet};

    /// D-R1 on both sides: a gene whose 1-based start is the base right after the 0-based
    /// region end bridges a gap only with the 1-based query (iqrs).
    #[test]
    fn d_r1_gap_bridge_both_sides() {
        let genes = AnnGenes { ids: vec!["g".into()], iv: vec![(301, 401)] };
        let rec = |s: i64, e: i64| PassARecord { pos0: s, end0: e + 1, flag: 0, name_hash: s as u64 };
        let cfg = PassAConfig { coverage_cutoff: -1, expected_primary: 0 };
        let run = || pass_a([rec(100, 300), rec(400, 500)], &genes, cfg).clusters.len();
        assert_eq!(with_quirks(QuirkSet::NONE, run), 1);
        if iq_types::QUIRKS41 {
            // 4.1: query (100, 300) misses a gene starting at 301
            assert_eq!(with_quirks(QuirkSet::only(DevId::R1), run), 2);
        }
    }
}
