//! Readers for the IsoQuant / iqrs output files compared by `iqrs-compare`.
//!
//! Everything is streamed line by line through a 1 MB buffer; `.gz` inputs are decoded with
//! `MultiGzDecoder` (also handles bgzip). Only the fields the metrics need are kept.

use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

/// Opens a possibly gzip-compressed text file.
pub fn open(p: &Path) -> Result<Box<dyn BufRead + Send>> {
    let f = File::open(p).with_context(|| format!("cannot open {}", p.display()))?;
    let gz = p.extension().is_some_and(|e| e == "gz");
    Ok(if gz {
        Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(f)))
    } else {
        Box::new(BufReader::with_capacity(1 << 20, f))
    })
}

/// Calls `f` for every line (without the line terminator).
pub fn for_each_line(p: &Path, mut f: impl FnMut(&str) -> Result<()>) -> Result<()> {
    let mut r = open(p)?;
    let mut buf = String::new();
    let mut n = 0usize;
    loop {
        buf.clear();
        if r.read_line(&mut buf).with_context(|| format!("{}: read error after line {n}", p.display()))? == 0 {
            return Ok(());
        }
        n += 1;
        f(buf.trim_end_matches(['\n', '\r'])).with_context(|| format!("{}:{n}", p.display()))?;
    }
}

/// Files of `dir` whose name satisfies `pred`, sorted by name.
pub fn find_files(dir: &Path, pred: impl Fn(&str) -> bool) -> Result<Vec<PathBuf>> {
    let mut v = Vec::new();
    for e in std::fs::read_dir(dir).with_context(|| format!("cannot list {}", dir.display()))? {
        let e = e?;
        let name = e.file_name();
        let Some(name) = name.to_str() else { continue };
        if pred(name) && e.path().is_file() {
            v.push(e.path());
        }
    }
    v.sort();
    Ok(v)
}

/// Exactly one file matching `pred` (error if several), or `None`.
pub fn find_one(dir: &Path, what: &str, pred: impl Fn(&str) -> bool) -> Result<Option<PathBuf>> {
    let v = find_files(dir, pred)?;
    if v.len() > 1 {
        bail!("{}: several {what} files: {:?}", dir.display(), v);
    }
    Ok(v.into_iter().next())
}

/// Value of a GTF attribute `key "value";`.
pub fn gtf_attr<'a>(s: &'a str, key: &str) -> Option<&'a str> {
    let b = s.as_bytes();
    let mut i = 0;
    while let Some(p) = s[i..].find(key) {
        let a = i + p;
        let before_ok = a == 0 || b[a - 1] == b' ' || b[a - 1] == b';' || b[a - 1] == b'\t';
        let rest = &s[a + key.len()..];
        if before_ok {
            if let Some(r) = rest.strip_prefix(" \"") {
                return r.find('"').map(|e| &r[..e]);
            }
        }
        i = a + key.len();
    }
    None
}

/// One transcript model from `transcript_models.gtf`.
#[derive(Clone, Debug)]
pub struct Model {
    pub id: String,
    pub gene: String,
    /// Chromosome index (side-local while parsing, global after [`Side::remap_chroms`]).
    pub chr: u32,
    /// `b'+'`, `b'-'` or `b'.'`.
    pub strand: u8,
    /// Sorted exons, 1-based closed.
    pub exons: Vec<(i32, i32)>,
    /// Transcript `source` column (`IsoQuant` for novel models).
    pub source: String,
    /// Count from `discovered_transcript_counts.tsv` (0 when absent).
    pub count: f64,
}

impl Model {
    pub fn span(&self) -> (i32, i32) {
        (self.exons.first().map_or(0, |e| e.0), self.exons.last().map_or(0, |e| e.1))
    }
    pub fn is_mono(&self) -> bool {
        self.exons.len() == 1
    }
    /// Introns as `(last base of exon k, first base of exon k+1)`.
    pub fn introns(&self) -> impl Iterator<Item = (i32, i32)> + '_ {
        self.exons.windows(2).map(|w| (w[0].1, w[1].0))
    }
}

/// Stats rows of `discovered_transcript_counts.tsv`.
#[derive(Clone, Debug, Default)]
pub struct CountStats {
    pub assigned: f64,
    pub ambiguous: f64,
    pub no_feature: f64,
    pub not_aligned: f64,
    pub n_features: usize,
    /// Count rows whose feature is not a model of the GTF.
    pub unknown_features: usize,
    pub unknown_weight: f64,
}

/// One grouped-counts file.
#[derive(Clone, Debug)]
pub struct Grouped {
    /// File name with the leading `<prefix>.` removed, used to pair files across sides.
    pub key: String,
    /// `(model index or u32::MAX, group index, count)`.
    pub rows: Vec<(u32, u32, f64)>,
    pub groups: Vec<String>,
}

/// Fine read fates (4.1 `ReadAssignmentType` names as printed).
pub const FINE_FATES: [&str; 10] = [
    "unique",
    "unique_minor_difference",
    "ambiguous",
    "inconsistent",
    "inconsistent_ambiguous",
    "inconsistent_non_intronic",
    "noninformative",
    "intergenic",
    "discarded",
    "other",
];
/// Coarse fates gated by the acceptance metric (plus `absent` in the confusion matrix).
pub const COARSE_FATES: [&str; 5] = ["unique", "ambiguous", "inconsistent", "no_feature", "other"];

pub fn fine_fate(s: &str) -> u8 {
    FINE_FATES.iter().position(|&f| f == s).unwrap_or(9) as u8
}
pub fn coarse_of(fine: u8) -> usize {
    match fine {
        0 | 1 => 0,
        2 => 1,
        3..=5 => 2,
        6..=8 => 3,
        _ => 4,
    }
}

/// No model (`ReadRec::model`).
pub const NO_MODEL: u32 = 0x0FFF_FFFF;

/// One read of a per-read fate file, packed as `(read-name hash, fine << 28 | first model)`.
#[derive(Copy, Clone, Debug)]
pub struct ReadRec(pub u64, pub u32);
impl ReadRec {
    pub fn hash(&self) -> u64 {
        self.0
    }
    pub fn fine(&self) -> u8 {
        (self.1 >> 28) as u8
    }
    /// First model (side-local index) or [`NO_MODEL`].
    pub fn model(&self) -> u32 {
        self.1 & NO_MODEL
    }
}

#[derive(Debug, Default)]
pub struct Reads {
    pub path: PathBuf,
    /// Sorted by hash, one entry per read (16 B/read; a hash map needed ~3x that at 10^8 reads).
    pub recs: Vec<ReadRec>,
    pub n_rows: usize,
    pub unknown_models: usize,
    /// Reads with records at more than one locus (non-consecutive rows).
    pub multi_locus: usize,
}

/// Everything read from one output directory.
#[derive(Debug, Default)]
pub struct Side {
    pub label: String,
    pub dir: PathBuf,
    pub gtf: PathBuf,
    pub chroms: Vec<String>,
    pub models: Vec<Model>,
    pub id_map: FxHashMap<String, u32>,
    pub stats: CountStats,
    pub grouped: Vec<Grouped>,
    pub reads: Option<Reads>,
    pub warnings: Vec<String>,
}

pub fn is_gtf(n: &str) -> bool {
    n.ends_with("transcript_models.gtf") || n.ends_with("transcript_models.gtf.gz")
}
pub fn is_counts(n: &str) -> bool {
    n.ends_with("discovered_transcript_counts.tsv") || n.ends_with("discovered_transcript_counts.tsv.gz")
}
pub fn is_grouped(n: &str) -> bool {
    n.contains("discovered_transcript_grouped_") && (n.ends_with("counts.linear.tsv") || n.ends_with("counts.linear.tsv.gz"))
}
pub fn is_reads(n: &str) -> bool {
    ["transcript_model_reads.tsv", "read_fates.tsv"].iter().any(|s| n.ends_with(s) || n.ends_with(&format!("{s}.gz")))
}

fn intern(chroms: &mut Vec<String>, map: &mut FxHashMap<String, u32>, name: &str) -> u32 {
    if let Some(&i) = map.get(name) {
        return i;
    }
    let i = chroms.len() as u32;
    chroms.push(name.to_string());
    map.insert(name.to_string(), i);
    i
}

impl Side {
    /// Reads GTF, counts and grouped counts of `dir`, and the per-read file `reads` if given.
    pub fn load(label: &str, dir: &Path, reads: Option<&Path>) -> Result<Side> {
        let mut s = Side { label: label.to_string(), dir: dir.to_path_buf(), ..Default::default() };
        let gtf = find_one(dir, "transcript_models.gtf", is_gtf)?
            .with_context(|| format!("{}: no *transcript_models.gtf", dir.display()))?;
        s.gtf = gtf.clone();
        s.read_gtf(&gtf)?;
        match find_one(dir, "discovered_transcript_counts.tsv", is_counts)? {
            Some(p) => s.read_counts(&p)?,
            None => s.warnings.push("no *discovered_transcript_counts.tsv: all counts are 0".into()),
        }
        for p in find_files(dir, is_grouped)? {
            let g = s.read_grouped(&p)?;
            s.grouped.push(g);
        }
        if let Some(r) = reads {
            s.reads = Some(s.read_reads(r)?);
        }
        Ok(s)
    }

    fn read_gtf(&mut self, p: &Path) -> Result<()> {
        let mut chr_map = FxHashMap::default();
        let (chroms, models, id_map) = (&mut self.chroms, &mut self.models, &mut self.id_map);
        let mut last: Option<(String, u32)> = None; // fast path: exons follow their transcript
        for_each_line(p, |l| {
            if l.starts_with('#') || l.is_empty() {
                return Ok(());
            }
            let mut f = l.splitn(9, '\t');
            let (Some(chr), Some(src), Some(ty), Some(s), Some(e), _, Some(strand), _, Some(attrs)) =
                (f.next(), f.next(), f.next(), f.next(), f.next(), f.next(), f.next(), f.next(), f.next())
            else {
                bail!("malformed GTF line");
            };
            if ty != "exon" && ty != "transcript" {
                return Ok(());
            }
            let Some(tid) = gtf_attr(attrs, "transcript_id") else { return Ok(()) };
            let idx = match &last {
                Some((id, i)) if id == tid => *i,
                _ => {
                    let i = match id_map.get(tid) {
                        Some(&i) => i,
                        None => {
                            let i = models.len() as u32;
                            let c = intern(chroms, &mut chr_map, chr);
                            models.push(Model {
                                id: tid.to_string(),
                                gene: gtf_attr(attrs, "gene_id").unwrap_or("").to_string(),
                                chr: c,
                                strand: strand.as_bytes().first().copied().unwrap_or(b'.'),
                                exons: Vec::new(),
                                source: src.to_string(),
                                count: 0.0,
                            });
                            id_map.insert(tid.to_string(), i);
                            i
                        }
                    };
                    last = Some((tid.to_string(), i));
                    i
                }
            };
            let m = &mut models[idx as usize];
            if ty == "transcript" {
                m.source = src.to_string();
            } else {
                let s: i32 = s.parse().context("bad exon start")?;
                let e: i32 = e.parse().context("bad exon end")?;
                m.exons.push((s, e));
            }
            Ok(())
        })?;
        let mut empty = 0;
        for m in models.iter_mut() {
            m.exons.sort_unstable();
            if m.exons.is_empty() {
                empty += 1;
            }
        }
        if empty > 0 {
            self.warnings.push(format!("{empty} transcripts without exons in the GTF (ignored by matching)"));
        }
        Ok(())
    }

    fn read_counts(&mut self, p: &Path) -> Result<()> {
        let (models, id_map, st) = (&mut self.models, &self.id_map, &mut self.stats);
        for_each_line(p, |l| {
            if l.is_empty() || l.starts_with('#') || l.starts_with("feature_id") {
                return Ok(());
            }
            let Some((id, c)) = l.split_once('\t') else { bail!("malformed counts line") };
            let c: f64 = c.split('\t').next().unwrap_or("").trim().parse().context("bad count")?;
            match id {
                "__ambiguous" => st.ambiguous += c,
                "__no_feature" => st.no_feature += c,
                "__not_aligned" => st.not_aligned += c,
                _ if id.starts_with("__") => {}
                _ => {
                    st.assigned += c;
                    st.n_features += 1;
                    match id_map.get(id) {
                        Some(&i) => models[i as usize].count += c,
                        None => {
                            st.unknown_features += 1;
                            st.unknown_weight += c;
                        }
                    }
                }
            }
            Ok(())
        })?;
        if st.unknown_features > 0 {
            self.warnings.push(format!(
                "{} count rows ({} reads) name features absent from the GTF",
                st.unknown_features, st.unknown_weight
            ));
        }
        Ok(())
    }

    fn read_grouped(&self, p: &Path) -> Result<Grouped> {
        let name = p.file_name().and_then(|n| n.to_str()).unwrap_or("").to_string();
        let key = match name.find("discovered_transcript_grouped_") {
            Some(i) => name[i..].trim_end_matches(".gz").to_string(),
            None => name.clone(),
        };
        let mut g = Grouped { key, rows: Vec::new(), groups: Vec::new() };
        let mut gmap: FxHashMap<String, u32> = FxHashMap::default();
        for_each_line(p, |l| {
            if l.is_empty() || l.starts_with('#') || l.starts_with("feature_id") {
                return Ok(());
            }
            let mut f = l.split('\t');
            let (Some(id), Some(grp), Some(c)) = (f.next(), f.next(), f.next()) else {
                bail!("malformed grouped line")
            };
            if id.starts_with("__") {
                return Ok(());
            }
            let c: f64 = c.trim().parse().context("bad count")?;
            let gi = match gmap.get(grp) {
                Some(&i) => i,
                None => {
                    let i = g.groups.len() as u32;
                    g.groups.push(grp.to_string());
                    gmap.insert(grp.to_string(), i);
                    i
                }
            };
            g.rows.push((self.id_map.get(id).copied().unwrap_or(u32::MAX), gi, c));
            Ok(())
        })?;
        Ok(g)
    }

    /// Per-read fate file: either 4.1's `transcript_model_reads.tsv(.gz)` (one row per
    /// read x model, columns `read_id`, `isoform_id`, `isoform_assignment_type`) or the iqrs
    /// `--emit-read-fates` file (`read_id`, `model_ids` comma-joined or `*`, `assignment_type`).
    /// Columns are located by header name; read IDs are stored as 64-bit hashes.
    fn read_reads(&self, p: &Path) -> Result<Reads> {
        let mut r = Reads { path: p.to_path_buf(), ..Default::default() };
        let mut cols: Option<(usize, usize, usize)> = None;
        let find = |hdr: &[&str], names: &[&str]| hdr.iter().position(|h| names.contains(h));
        for_each_line(p, |l| {
            if l.is_empty() {
                return Ok(());
            }
            if cols.is_none() {
                let h = l.trim_start_matches('#');
                let hdr: Vec<&str> = h.split('\t').collect();
                if let Some(ri) = find(&hdr, &["read_id"]) {
                    let mi = find(&hdr, &["isoform_id", "model_ids", "model_id", "transcript_id"])
                        .context("per-read file: no isoform_id/model_ids column")?;
                    let ti = find(&hdr, &["isoform_assignment_type", "assignment_type"])
                        .context("per-read file: no assignment type column")?;
                    cols = Some((ri, mi, ti));
                    return Ok(());
                }
                if l.starts_with('#') {
                    return Ok(());
                }
                bail!("per-read file: no header with a read_id column");
            }
            if l.starts_with('#') {
                return Ok(());
            }
            let (ri, mi, ti) = cols.unwrap();
            let mut rid = "";
            let mut mid = "";
            let mut ty = "";
            for (k, f) in l.split('\t').enumerate() {
                if k == ri {
                    rid = f;
                } else if k == mi {
                    mid = f;
                } else if k == ti {
                    ty = f;
                }
            }
            r.n_rows += 1;
            let fine = fine_fate(ty) as u32;
            let h = hash_str(rid);
            let mut model = NO_MODEL;
            for m in mid.split(',') {
                if m.is_empty() || m == "." || m == "*" {
                    continue;
                }
                match self.id_map.get(m) {
                    Some(&i) => {
                        model = i.min(NO_MODEL - 1);
                        break;
                    }
                    None => r.unknown_models += 1,
                }
            }
            // rows of one read are consecutive in 4.1's file: merge in place
            match r.recs.last_mut() {
                Some(last) if last.0 == h => {
                    if last.model() == NO_MODEL {
                        last.1 = (last.1 & !NO_MODEL) | model;
                    }
                }
                _ => r.recs.push(ReadRec(h, (fine << 28) | model)),
            }
            Ok(())
        })?;
        // A read can occur at several loci (36 reads in SRR32588716 chr8): keep its most
        // assigned record (smallest fine-fate index, then smallest model index), independent
        // of file order.
        r.recs.sort_unstable_by_key(|x| (x.0, x.1));
        let n0 = r.recs.len();
        r.recs.dedup_by(|b, a| a.0 == b.0);
        r.multi_locus = n0 - r.recs.len();
        r.recs.shrink_to_fit();
        Ok(r)
    }

    /// Re-indexes chromosomes into the shared table `global`.
    pub fn remap_chroms(&mut self, global: &mut Vec<String>, gmap: &mut FxHashMap<String, u32>) {
        let remap: Vec<u32> = self.chroms.iter().map(|c| intern(global, gmap, c)).collect();
        for m in &mut self.models {
            m.chr = remap[m.chr as usize];
        }
    }
}

/// FNV-1a 64 over the read name; collisions are negligible at 10^8 reads for this purpose.
pub fn hash_str(s: &str) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in s.as_bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn attrs() {
        let a = r#"gene_id "gene:G1"; transcript_id "transcript7.8.nnic"; exons "5";"#;
        assert_eq!(gtf_attr(a, "transcript_id"), Some("transcript7.8.nnic"));
        assert_eq!(gtf_attr(a, "gene_id"), Some("gene:G1"));
        assert_eq!(gtf_attr(a, "exon_id"), None);
        assert_eq!(gtf_attr(r#"xgene_id "no"; gene_id "yes";"#, "gene_id"), Some("yes"));
    }

    #[test]
    fn fates() {
        assert_eq!(coarse_of(fine_fate("unique_minor_difference")), 0);
        assert_eq!(coarse_of(fine_fate("inconsistent_non_intronic")), 2);
        assert_eq!(coarse_of(fine_fate("discarded")), 3);
        assert_eq!(coarse_of(fine_fate("weird")), 4);
    }
}
