//! U20 fixture runners: `gff_dump`, `counter_add_read_info`, `counter_dump`.
//!
//! State the cases do not carry is reconstructed from the capture itself:
//! - `gff_dump` records the exon-id counter but not the `FeatureIdStorage` dict. The dict
//!   before a block = annotation `exon_id`s (parsed GFF3, `iq_annot`) + every novel exon id
//!   *created* (printed bare, as 4.1 does) by earlier `gff_dump` cases of the chromosome.
//!   [`Prescan`] reads those from the recorded lines of the unit file.
//! - `counter_add_read_info` records the counter state before/after but not the earlier
//!   reads (sampled at 5%). The runner predicts the read's *delta* on a fresh counter and
//!   compares deltas: per feature, the groups whose value changed and by how much; groups we
//!   touch must be present afterwards; `__ambiguous` / `__no_feature` increments. The
//!   `confirmed_added` set (newly confirmed features) must be a subset of the features the
//!   read confirms; that check needs the read's `corrected_exons` and the model's introns,
//!   found through the case key ([`Prescan::counter_keys`]) and the block's `gff_dump`.
//! - `counter_dump` cannot reproduce the `.usable` file (`reads_for_tpm` is not recorded and
//!   the file is not an iqrs output); it is left out of the comparison.

use crate::counter::{FeatureCounter, ReadFate, Strategy};
use crate::gtf::{BlockAnnot, GtfChrWriter, PrintModel};
use crate::output::{group_name, ChrCounts};
use iq_types::fixture::schema::{
    BlockRef, Case, IndexDump, U20CounterAddIn, U20CounterAddOut, U20CounterDumpIn, U20CounterDumpOut, U20GffDumpIn,
    U20GffDumpOut,
};
use iq_types::fixture::{
    read_lines, unit_files, Diff, DiffEntry, FixError, FixResult, FixtureCtx, Registry, UnitFixture, Verdict,
};
use iq_types::serde::{Deserialize, Serialize};
use iq_types::{AssignmentType, GroupPool, Iv, Pos, Strand};
use rustc_hash::{FxHashMap, FxHashSet};
use std::borrow::Cow;
use std::cell::RefCell;
use std::hash::{Hash, Hasher};
use std::path::PathBuf;
use std::sync::{Arc, Mutex, OnceLock};

pub(crate) fn register_fixtures(reg: &mut Registry) {
    reg.add::<U20GffDump>();
    reg.add::<U20CounterAdd>();
    reg.add::<U20CounterDump>();
}

// ------------------------------------------------------------------------------------------
// Prescan of the unit file
// ------------------------------------------------------------------------------------------

#[derive(Default)]
struct Prescan {
    /// `(chr, block_no)` → position of the block's `gff_dump` among the chromosome's.
    gff_order: FxHashMap<(String, u32), usize>,
    /// Per chromosome: `(s, e, strand)` → `(gff case order, N)` of the case that created
    /// the novel exon id `N`.
    novel_exons: FxHashMap<String, FxHashMap<(Pos, Pos, char), (usize, u64)>>,
    /// `(chr, block_no, transcript_id)` → number of exons of the printed model.
    model_exons: FxHashMap<(String, u32, String), usize>,
    /// Hash of a counter case's `in` → case keys.
    counter_keys: FxHashMap<u64, Vec<String>>,
}

fn in_hash(v: &iq_types::serde_json::Value) -> u64 {
    let s = iq_types::serde_json::to_string(v).unwrap_or_default();
    let mut h = std::collections::hash_map::DefaultHasher::new();
    s.hash(&mut h);
    h.finish()
}

fn attr_value<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    let pat = format!("{key} \"");
    let i = attrs.find(&pat)? + pat.len();
    let j = attrs[i..].find('"')?;
    Some(&attrs[i..i + j])
}

fn prescan(ctx: &FixtureCtx) -> FixResult<Arc<Prescan>> {
    static CACHE: OnceLock<Mutex<FxHashMap<PathBuf, Arc<Prescan>>>> = OnceLock::new();
    let cache = CACHE.get_or_init(Default::default);
    if let Some(p) = cache.lock().unwrap().get(&ctx.root) {
        return Ok(p.clone());
    }
    let mut p = Prescan::default();
    let mut n_gff: FxHashMap<String, usize> = FxHashMap::default();
    for f in unit_files(&ctx.root, "U20")? {
        for case in read_lines::<Case>(&f)? {
            let case = case?;
            let fn_name = case.input.get("fn").and_then(|v| v.as_str()).unwrap_or("");
            match fn_name {
                "gff_dump" => {
                    let i: U20GffDumpIn = U20GffDumpIn::deserialize(&case.input)?;
                    if i.printer != "transcript_models" {
                        continue;
                    }
                    let chr = i.block.chr.clone();
                    let k = n_gff.entry(chr.clone()).or_default();
                    let order = *k;
                    *k += 1;
                    p.gff_order.insert((chr.clone(), i.block.block_no), order);
                    for m in &i.models {
                        p.model_exons.insert((chr.clone(), i.block.block_no, m.id.clone()), m.exons.len());
                    }
                    let Ok(o) = U20GffDumpOut::deserialize(&case.out) else {
                        continue;
                    };
                    let novel = p.novel_exons.entry(chr.clone()).or_default();
                    for line in &o.lines {
                        let cols: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
                        if cols.len() < 9 || cols[2] == "gene" || cols[2] == "transcript" {
                            continue;
                        }
                        let Some(id) = attr_value(cols[8], "exon_id") else {
                            continue;
                        };
                        if !id.is_empty() && id.bytes().all(|b| b.is_ascii_digit()) {
                            let (Ok(s), Ok(e)) = (cols[3].parse::<Pos>(), cols[4].parse::<Pos>()) else {
                                continue;
                            };
                            let strand = cols[6].chars().next().unwrap_or('.');
                            novel.entry((s, e, strand)).or_insert((order, id.parse().unwrap_or(0)));
                        }
                    }
                }
                "counter_add_read_info" => {
                    p.counter_keys.entry(in_hash(&case.input)).or_default().push(case.case.clone());
                }
                _ => {}
            }
        }
    }
    let p = Arc::new(p);
    cache.lock().unwrap().insert(ctx.root.clone(), p.clone());
    Ok(p)
}

// ------------------------------------------------------------------------------------------
// gff_dump
// ------------------------------------------------------------------------------------------

/// [`BlockAnnot`] over a captured `GeneInfo` dump + the parsed annotation (exon ids).
struct DumpAnnot<'a> {
    dump: &'a IndexDump,
    models: &'a [PrintModel<'a>],
    regions: FxHashMap<&'a str, Iv>,
    sources: FxHashMap<&'a str, &'a str>,
    attrs: FxHashMap<&'a str, &'a str>,
    ann: Option<(&'a iq_annot::AnnotationIndex, iq_types::ChrIdx)>,
}

impl<'a> DumpAnnot<'a> {
    fn new(
        dump: &'a IndexDump,
        models: &'a [PrintModel<'a>],
        ann: Option<(&'a iq_annot::AnnotationIndex, iq_types::ChrIdx)>,
    ) -> Self {
        let regions = if dump.exon_features.is_empty() {
            FxHashMap::default() // gene_info.empty()
        } else {
            dump.gene_regions.iter().map(|(g, iv)| (g.as_str(), *iv)).collect()
        };
        DumpAnnot {
            dump,
            models,
            regions,
            sources: dump.sources.iter().map(|(k, v)| (k.as_str(), v.as_str())).collect(),
            attrs: dump.feature_attributes.iter().map(|(k, v)| (k.as_str(), v.as_str())).collect(),
            ann,
        }
    }
}

impl BlockAnnot for DumpAnnot<'_> {
    fn gene_region(&self, i: usize) -> Option<Iv> {
        self.regions.get(self.models[i].gid.as_ref()).copied()
    }
    fn gene_source(&self, i: usize) -> Option<Cow<'_, str>> {
        self.sources.get(self.models[i].gid.as_ref()).map(|s| Cow::Borrowed(*s))
    }
    fn gene_attrs(&self, i: usize) -> Option<Cow<'_, str>> {
        self.attrs.get(self.models[i].gid.as_ref()).map(|s| Cow::Borrowed(*s))
    }
    fn tx_attrs(&self, i: usize) -> Option<Cow<'_, str>> {
        self.attrs.get(self.models[i].tid.as_ref()).map(|s| Cow::Borrowed(*s))
    }
    fn exon_key(&self, i: usize, s: Pos, e: Pos, strand: Strand) -> bool {
        let _ = self.dump;
        let k = format!("{}_{s}_{e}_{strand}", self.models[i].tid);
        self.attrs.contains_key(k.as_str())
    }
    fn annotated_exon_id(&self, s: Pos, e: Pos, strand: Strand) -> Option<Cow<'_, str>> {
        let (ann, chr) = self.ann?;
        ann.exon_id(chr, Iv::new(s, e), strand).map(Cow::Borrowed)
    }
}

pub struct U20GffDump;

impl UnitFixture for U20GffDump {
    const UNIT: &'static str = "U20";
    const FN: &'static str = "gff_dump";
    type In = U20GffDumpIn;
    type Out = U20GffDumpOut;

    fn run(i: &U20GffDumpIn, ctx: &FixtureCtx) -> FixResult<U20GffDumpOut> {
        if i.printer != "transcript_models" {
            return Err(FixError("only the transcript_models printer is ported".into()));
        }
        let dump = ctx.index_dump(&i.index)?;
        let ann =
            iq_annot::fixture::cached_annotation(&iq_annot::fixture::gff3_path(ctx)?, iq_annot::GeneTypes::Isoquant)?;
        let chr = i.block.chr.as_str();
        let pre = prescan(ctx)?;
        let order = *pre.gff_order.get(&(chr.to_owned(), i.block.block_no)).ok_or("block not in prescan")?;

        let mut w = GtfChrWriter::new(chr);
        w.printed_gene_ids = i.printed_gene_ids.iter().cloned().collect();
        if let Some(nov) = pre.novel_exons.get(chr) {
            for (&(s, e, st), &(k, n)) in nov {
                if k < order {
                    w.exon_ids.seed(s, e, Strand::parse(&st.to_string()).unwrap_or(Strand::Unknown), n);
                }
            }
        }
        w.exon_ids.counter = i.exon_id_counter as u64;

        let models: Vec<PrintModel<'_>> = i
            .models
            .iter()
            .map(|m| PrintModel {
                tid: Cow::Borrowed(m.id.as_str()),
                gid: Cow::Borrowed(m.gene.as_str()),
                strand: m.strand,
                exons: Cow::Borrowed(m.exons.as_slice()),
                source: Cow::Borrowed(m.source.as_str()),
                other: m.other_features.iter().map(|(s, e, t)| (*s, *e, Cow::Borrowed(t.as_str()))).collect(),
                additional: m.additional_info.iter().map(|(k, v)| (Cow::Borrowed(k.as_str()), v.clone())).collect(),
            })
            .collect();
        let annot = DumpAnnot::new(&dump, &models, ann.chrom(chr).map(|c| (&*ann, c)));
        let mut out = String::new();
        w.dump_block(&annot, &models, &mut out);
        Ok(U20GffDumpOut {
            lines: out.split_inclusive('\n').map(str::to_owned).collect(),
            exon_id_counter: w.exon_ids.counter as i64,
        })
    }
}

// ------------------------------------------------------------------------------------------
// counter_add_read_info
// ------------------------------------------------------------------------------------------

/// Delta form of a `counter_add_read_info` output (see module docs).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(crate = "iq_types::serde", from = "U20CounterAddOut")]
pub struct AddDelta {
    /// Sorted by feature: `(feature, changed groups (group, delta), groups present after)`.
    pub features: Vec<(String, Vec<(String, f64)>, Vec<String>)>,
    pub ambiguous: i64,
    pub not_assigned: i64,
    /// Recorded: `confirmed_added`. Predicted: features the read confirms.
    pub confirmed: Vec<String>,
    /// Recorded only: false if a group other than the changed ones changed value.
    pub consistent: bool,
}

impl From<U20CounterAddOut> for AddDelta {
    fn from(o: U20CounterAddOut) -> Self {
        let mut features: Vec<(String, Vec<(String, f64)>, Vec<String>)> = o
            .features
            .into_iter()
            .map(|(f, gl)| {
                let changed = gl.iter().filter(|g| g.1 != g.2).map(|g| (g.0.clone(), g.2 - g.1)).collect();
                let present = gl.iter().map(|g| g.0.clone()).collect();
                (f, changed, present)
            })
            .collect();
        features.sort_by(|a, b| a.0.cmp(&b.0));
        let mut confirmed = o.confirmed_added;
        confirmed.sort();
        AddDelta {
            features,
            ambiguous: o.ambiguous.1 - o.ambiguous.0,
            not_assigned: o.not_assigned.1 - o.not_assigned.0,
            confirmed,
            consistent: true,
        }
    }
}

thread_local! {
    /// The input of the case being judged (`compare` needs it for the confirm check).
    static CURRENT: RefCell<Option<U20CounterAddIn>> = const { RefCell::new(None) };
}

fn is_gene_counter(name: &str) -> bool {
    name.contains("discovered_gene")
}

/// The read's distinct features (match order) and the assignment type the counter uses.
fn read_features(i: &U20CounterAddIn) -> (AssignmentType, Vec<String>) {
    let gene = is_gene_counter(&i.counter.name);
    let mut feats: Vec<String> = Vec::new();
    for m in &i.assignment.matches {
        let f = if gene { m.gene.as_ref() } else { m.tx.as_ref() };
        if let Some(f) = f {
            if !f.is_empty() && !feats.contains(f) {
                feats.push(f.clone());
            }
        }
    }
    (if gene { i.assignment.gene_atype } else { i.assignment.atype }, feats)
}

pub struct U20CounterAdd;

impl UnitFixture for U20CounterAdd {
    const UNIT: &'static str = "U20";
    const FN: &'static str = "counter_add_read_info";
    type In = U20CounterAddIn;
    type Out = AddDelta;

    fn run(i: &U20CounterAddIn, _ctx: &FixtureCtx) -> FixResult<AddDelta> {
        CURRENT.with(|c| *c.borrow_mut() = Some(i.clone()));
        if i.read.is_none() {
            return Err("read is None (not_aligned branch) not expected for model counters".into());
        }
        let (atype, feats) = read_features(i);
        let group = if i.counter.grouped {
            let g = i.read_group.get(i.counter.group_index as usize).ok_or("grouped counter without read group")?;
            g.clone()
        } else {
            "NA".to_owned()
        };
        let ids: Vec<u32> = (0..feats.len() as u32).collect();
        let mut c = FeatureCounter::new();
        let strategy =
            if is_gene_counter(&i.counter.name) { Strategy::UniqueSplicingConsistent } else { Strategy::UniqueOnly };
        let fate = c.add_features(strategy, atype, !i.assignment.matches.is_empty(), &ids, 0);
        let mut features: Vec<(String, Vec<(String, f64)>, Vec<String>)> = feats
            .iter()
            .enumerate()
            .map(|(k, f)| {
                let s = c.feats.get(k);
                let touched = s.is_some_and(|s| s.touched);
                let changed = match s {
                    Some(s) if s.touched && s.total != 0.0 => vec![(group.clone(), s.total)],
                    _ => vec![],
                };
                (f.clone(), changed, if touched { vec![group.clone()] } else { vec![] })
            })
            .collect();
        features.sort_by(|a, b| a.0.cmp(&b.0));
        let confirmed = if fate == ReadFate::Unique { vec![feats[0].clone()] } else { vec![] };
        Ok(AddDelta {
            features,
            ambiguous: c.ambiguous as i64,
            not_assigned: c.no_feature as i64,
            confirmed,
            consistent: true,
        })
    }

    fn compare(got: &AddDelta, want: &AddDelta, ctx: &FixtureCtx) -> Verdict {
        let mut entries = Vec::new();
        let mut e =
            |path: &str, g: String, w: String| entries.push(DiffEntry { path: path.to_owned(), got: g, want: w });
        let gf: Vec<&String> = got.features.iter().map(|f| &f.0).collect();
        let wf: Vec<&String> = want.features.iter().map(|f| &f.0).collect();
        if gf != wf {
            e("features", format!("{gf:?}"), format!("{wf:?}"));
        } else {
            for (g, w) in got.features.iter().zip(&want.features) {
                if g.1 != w.1 {
                    e(&format!("features.{}.delta", g.0), format!("{:?}", g.1), format!("{:?}", w.1));
                }
                if !g.2.iter().all(|x| w.2.contains(x)) {
                    e(&format!("features.{}.present", g.0), format!("{:?}", g.2), format!("{:?}", w.2));
                }
                if w.1.len() > 1 {
                    e(&format!("features.{}.multi_change", g.0), String::new(), format!("{:?}", w.1));
                }
            }
        }
        if got.ambiguous != want.ambiguous {
            e("ambiguous", got.ambiguous.to_string(), want.ambiguous.to_string());
        }
        if got.not_assigned != want.not_assigned {
            e("not_assigned", got.not_assigned.to_string(), want.not_assigned.to_string());
        }
        // confirmed_added ⊆ confirmed-by-this-read (exact condition needs the read and model)
        if !want.confirmed.is_empty() {
            match confirm_check(got, ctx) {
                Ok(conf) => {
                    if !want.confirmed.iter().all(|f| conf.contains(f)) {
                        e("confirmed_added", format!("{conf:?}"), format!("{:?}", want.confirmed));
                    }
                }
                Err(err) => return Verdict::Mismatch(Diff::error(format!("confirm check: {err}"))),
            }
        }
        if entries.is_empty() {
            Verdict::Exact
        } else {
            Verdict::Mismatch(Diff { entries, ..Default::default() })
        }
    }
}

/// `confirms_feature` for the read being judged (`long_read_counter.py:116-121`): genes:
/// unique gene type; transcripts: unique and (model without introns or read with > 1
/// corrected exon).
fn confirm_check(got: &AddDelta, ctx: &FixtureCtx) -> FixResult<Vec<String>> {
    let Some(i) = CURRENT.with(|c| c.borrow().clone()) else {
        return Err("no current input".into());
    };
    if got.confirmed.is_empty() || is_gene_counter(&i.counter.name) {
        return Ok(got.confirmed.clone());
    }
    let pre = prescan(ctx)?;
    let v = iq_types::serde_json::to_value(&i)?;
    let keys = pre.counter_keys.get(&in_hash(&v)).ok_or("case not found in prescan")?;
    let parsed: FxHashSet<(String, String, String, u32)> = keys
        .iter()
        .filter_map(|k| iq_types::fixture::CaseKey::parse(k))
        .filter_map(|k| Some((k.sample.to_owned(), k.slice.to_owned(), k.chr.to_owned(), k.block_no?)))
        .collect();
    if parsed.len() != 1 {
        return Err(format!("ambiguous case identity ({} blocks)", parsed.len()).into());
    }
    let (sample, slice, chr, block_no) = parsed.into_iter().next().unwrap();
    let read_idx = i.read.as_ref().ok_or("no read")?.read_idx as usize;
    let tid = &got.confirmed[0];
    let n_exons = *pre.model_exons.get(&(chr.clone(), block_no, tid.clone())).ok_or("model not printed in block")?;
    if n_exons <= 1 {
        return Ok(got.confirmed.clone());
    }
    let block = ctx.block(&BlockRef { sample, slice, chr, block_no })?;
    let rec = block.records.get(read_idx).ok_or("read_idx out of range")?;
    Ok(if rec.corrected_exons.len() > 1 { got.confirmed.clone() } else { vec![] })
}

// ------------------------------------------------------------------------------------------
// counter_dump
// ------------------------------------------------------------------------------------------

pub struct U20CounterDump;

/// Chromosome counter state of a `counter_dump` input (group names interned into `pool`).
fn counts_from_dump(i: &U20CounterDumpIn, pool: &mut GroupPool) -> ChrCounts {
    let mut c = ChrCounts::new();
    let all: FxHashSet<&str> = i.all_features.iter().map(String::as_str).collect();
    let conf: FxHashSet<&str> = i.confirmed_features.iter().map(String::as_str).collect();
    for (f, groups) in &i.feature_counter {
        let e = c.entry(f);
        for (g, v) in groups {
            if i.counter.grouped {
                e.groups.push((pool.intern(g), *v));
            } else {
                e.total = *v; // the ungrouped counter only has group 0 ("NA")
            }
        }
    }
    for f in &i.all_features {
        c.entry(f);
    }
    for e in &mut c.feats {
        e.in_all = all.contains(e.id.as_str());
        e.confirmed = conf.contains(e.id.as_str());
    }
    c.ambiguous = i.ambiguous as u64;
    c.no_feature = i.not_assigned as u64;
    c
}

impl UnitFixture for U20CounterDump {
    const UNIT: &'static str = "U20";
    const FN: &'static str = "counter_dump";
    type In = U20CounterDumpIn;
    type Out = U20CounterDumpOut;

    fn run(i: &U20CounterDumpIn, _ctx: &FixtureCtx) -> FixResult<U20CounterDumpOut> {
        let mut pool = GroupPool::new();
        let c = counts_from_dump(i, &mut pool);
        let name = i.counter.name.clone();
        let mut files = Vec::new();
        if i.counter.grouped {
            let mut s = String::from("feature_id\tgroup_id\tcount\n");
            c.write_grouped_rows(|g| group_name(&pool, g), &mut s);
            files.push((name, s));
        } else {
            let mut s = String::from("feature_id\tcount\n");
            c.write_ungrouped_rows(&mut s);
            files.push((name.clone(), s));
            files.push((format!("{name}.stats"), c.stats_rows(i.not_aligned as u64)));
        }
        Ok(U20CounterDumpOut { files })
    }

    fn compare(got: &U20CounterDumpOut, want: &U20CounterDumpOut, _ctx: &FixtureCtx) -> Verdict {
        let want: Vec<&(String, String)> = want.files.iter().filter(|f| !f.0.ends_with(".usable")).collect();
        let got: Vec<&(String, String)> = got.files.iter().collect();
        if got == want {
            return Verdict::Exact;
        }
        let mut entries = Vec::new();
        for k in 0..got.len().max(want.len()) {
            let g = got.get(k).map(|f| format!("{}: {}", f.0, f.1)).unwrap_or_default();
            let w = want.get(k).map(|f| format!("{}: {}", f.0, f.1)).unwrap_or_default();
            if g != w {
                entries.push(DiffEntry { path: format!("files[{k}]"), got: g, want: w });
            }
        }
        Verdict::Mismatch(Diff { entries, ..Default::default() })
    }
}

// ------------------------------------------------------------------------------------------
// End-to-end: final files of the fullchr8 captures
// ------------------------------------------------------------------------------------------

#[cfg(test)]
mod e2e {
    use super::*;
    use crate::ids::NovelNamer;
    use crate::output::{assemble, ChrFeat, ChrOutput, OutputSpec};
    use iq_types::{GeneKey, ModelKind, TranscriptModel};
    use std::path::Path;

    struct Capture {
        gff: Vec<U20GffDumpIn>,
        dumps: Vec<U20CounterDumpIn>,
    }

    fn load(root: &Path) -> Option<Capture> {
        let files = unit_files(root, "U20").ok()?;
        if files.is_empty() {
            return None;
        }
        let mut cap = Capture { gff: vec![], dumps: vec![] };
        for f in files {
            for case in read_lines::<Case>(&f).unwrap() {
                let case = case.unwrap();
                match case.input.get("fn").and_then(|v| v.as_str()) {
                    Some("gff_dump") => cap.gff.push(U20GffDumpIn::deserialize(&case.input).unwrap()),
                    Some("counter_dump") => cap.dumps.push(U20CounterDumpIn::deserialize(&case.input).unwrap()),
                    _ => {}
                }
            }
        }
        Some(cap)
    }

    fn rename_attr(line: &str, key: &str, map: &FxHashMap<String, String>) -> String {
        let pat = format!("{key} \"");
        let Some(i) = line.find(&pat) else {
            return line.to_owned();
        };
        let j = i + pat.len();
        let k = j + line[j..].find('"').unwrap();
        match map.get(&line[j..k]) {
            Some(n) => format!("{}{}{}", &line[..j], n, &line[k..]),
            None => line.to_owned(),
        }
    }

    fn run_capture(sample: &str) {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../fixtures").join(sample).join("fullchr8");
        let Some(cap) = load(&root) else {
            eprintln!("skipping: no capture at {}", root.display());
            return;
        };
        let ctx = FixtureCtx::new(&root);
        let ann = iq_annot::fixture::cached_annotation(
            &iq_annot::fixture::gff3_path(&ctx).unwrap(),
            iq_annot::GeneTypes::Isoquant,
        )
        .unwrap();
        let chr = "8";
        let chr_idx = ann.chrom(chr);
        let out_dir = root.join("outputs");
        let want_gtf = std::fs::read_to_string(out_dir.join("OUT.transcript_models.gtf")).unwrap();
        let want_counts = std::fs::read_to_string(out_dir.join("OUT.discovered_transcript_counts.tsv")).unwrap();
        let want_grouped =
            std::fs::read_to_string(out_dir.join("OUT.discovered_transcript_grouped_tag_HP_counts.linear.tsv"))
                .unwrap();
        let cmd = want_gtf.lines().nth(1).unwrap().strip_prefix("# Command line: ").unwrap().to_owned();
        let spec = OutputSpec {
            dir: PathBuf::new(),
            prefix: "OUT".into(),
            cmd_line: cmd,
            version: "4.1.0".into(),
            group_tag: Some("HP".into()),
            grouped_out: None,
            gzip: false,
        };

        // counts state of the chromosome (ungrouped + grouped dumps)
        let mut pool = GroupPool::new();
        let tx_dumps: Vec<&U20CounterDumpIn> =
            cap.dumps.iter().filter(|d| d.counter.name.contains("discovered_transcript")).collect();
        let ung = tx_dumps.iter().find(|d| !d.counter.grouped).unwrap();
        let grp = tx_dumps.iter().find(|d| d.counter.grouped).unwrap();
        let mut counts = counts_from_dump(ung, &mut pool);
        let gcounts = counts_from_dump(grp, &mut pool);
        for f in &gcounts.feats {
            counts.entry(&f.id).groups = f.groups.clone();
        }

        // (a) 4.1 ids: sequential printing over all blocks + merge = the 4.1 files
        let mut gtf = GtfChrWriter::new(chr);
        let mut gtf_ours = GtfChrWriter::new(chr);
        let mut namer = NovelNamer::new(chr, chr_idx.map(|c| ann.forbidden_ids(c)).unwrap_or(&[]));
        let mut map: FxHashMap<String, String> = FxHashMap::default();
        let mut back: FxHashMap<String, String> = FxHashMap::default();
        let mut text41 = String::new();
        let mut text_ours = String::new();
        let mut gene_strs: Vec<String> = Vec::new();
        for g in &cap.gff {
            assert_eq!(g.block.chr, chr);
            let dump = ctx.index_dump(&g.index).unwrap();
            let pms: Vec<PrintModel<'_>> = g
                .models
                .iter()
                .map(|m| PrintModel {
                    tid: Cow::Borrowed(m.id.as_str()),
                    gid: Cow::Borrowed(m.gene.as_str()),
                    strand: m.strand,
                    exons: Cow::Borrowed(m.exons.as_slice()),
                    source: Cow::Borrowed(m.source.as_str()),
                    other: m.other_features.iter().map(|(s, e, t)| (*s, *e, Cow::Borrowed(t.as_str()))).collect(),
                    additional: vec![],
                })
                .collect();
            let annot = DumpAnnot::new(&dump, &pms, chr_idx.map(|c| (&*ann, c)));
            gtf.dump_block(&annot, &pms, &mut text41);

            // (b) our names (D-M1 renumbering) for the same models
            let mut novel_genes: Vec<&str> = Vec::new();
            let models: Vec<TranscriptModel> = g
                .models
                .iter()
                .enumerate()
                .map(|(k, m)| TranscriptModel {
                    key: k as u32,
                    strand: m.strand,
                    gene: if m.gene.starts_with("novel_gene_") {
                        let p = novel_genes.iter().position(|x| *x == m.gene).unwrap_or_else(|| {
                            novel_genes.push(&m.gene);
                            novel_genes.len() - 1
                        });
                        GeneKey::Novel(p as u32)
                    } else {
                        let p = gene_strs.iter().position(|x| *x == m.gene).unwrap_or_else(|| {
                            gene_strs.push(m.gene.clone());
                            gene_strs.len() - 1
                        });
                        GeneKey::Ref(p as u32)
                    },
                    kind: match m.ty {
                        iq_types::TranscriptModelType::Known => ModelKind::Known(k as u32),
                        iq_types::TranscriptModelType::NovelInCatalog => ModelKind::Nic,
                        _ => ModelKind::Nnic,
                    },
                    exons: m.exons.clone(),
                    intron_path: Box::new([]),
                    polya_confirmed: false,
                    tss_confirmed: false,
                    source_known: None,
                })
                .collect();
            let names =
                namer.name_block(&models, |t| g.models[t as usize].id.as_str(), |x| gene_strs[x as usize].as_str());
            for (k, m) in g.models.iter().enumerate() {
                for (a, b) in [(&m.id, &names.tx[k]), (&m.gene, &names.gene[k])] {
                    let prev = map.insert(a.clone(), b.clone());
                    assert!(prev.is_none() || prev.as_ref() == Some(b), "{a} renamed inconsistently");
                    let prevb = back.insert(b.clone(), a.clone());
                    assert!(prevb.is_none() || prevb.as_ref() == Some(a), "{b} assigned to two ids");
                }
            }
            let pms2: Vec<PrintModel<'_>> = pms
                .iter()
                .enumerate()
                .map(|(k, p)| PrintModel {
                    tid: names.tx[k].clone().into(),
                    gid: names.gene[k].clone().into(),
                    ..p.clone()
                })
                .collect();
            let annot2 = DumpAnnot::new(&dump, &pms2, chr_idx.map(|c| (&*ann, c)));
            gtf_ours.dump_block(&annot2, &pms2, &mut text_ours);
        }
        let t = assemble(&spec, &[ChrOutput { chr: chr.into(), gtf: text41, counts: counts.clone() }], &pool, 0);
        let want_gtf_n = want_gtf.clone();
        assert!(String::from_utf8(t.gtf).unwrap() == want_gtf_n, "{sample}: GTF differs");
        assert!(t.counts == want_counts, "{sample}: counts differ");
        assert!(t.grouped.as_deref() == Some(want_grouped.as_str()), "{sample}: grouped differs");

        // (b) renumbered output == 4.1 output with ids mapped
        let mut rcounts = ChrCounts::new();
        for f in &counts.feats {
            let id = map.get(&f.id).cloned().unwrap_or_else(|| f.id.clone());
            *rcounts.entry(&id) = ChrFeat { id: id.clone(), ..f.clone() };
        }
        rcounts.ambiguous = counts.ambiguous;
        rcounts.no_feature = counts.no_feature;
        let r = assemble(&spec, &[ChrOutput { chr: chr.into(), gtf: text_ours, counts: rcounts }], &pool, 0);
        let renamed_gtf: String = want_gtf_n
            .split_inclusive('\n')
            .map(|l| rename_attr(&rename_attr(l, "transcript_id", &map), "gene_id", &map))
            .collect();
        assert!(String::from_utf8(r.gtf).unwrap() == renamed_gtf, "{sample}: renumbered GTF differs");
        let resort = |text: &str, ncol: usize| -> String {
            let mut lines: Vec<String> = text.lines().map(|l| l.to_owned()).collect();
            let head = lines.remove(0);
            let mut stats: Vec<String> = lines.iter().filter(|l| l.starts_with("__")).cloned().collect();
            let mut rows: Vec<(String, String)> = lines
                .iter()
                .filter(|l| !l.starts_with("__"))
                .map(|l| {
                    let (id, rest) = l.split_once('\t').unwrap();
                    let id = map.get(id).cloned().unwrap_or_else(|| id.to_owned());
                    (id, rest.to_owned())
                })
                .collect();
            rows.sort_by(|a, b| a.0.as_bytes().cmp(b.0.as_bytes())); // stable: group order kept
            let _ = ncol;
            let mut out = head + "\n";
            for (a, b) in rows {
                out += &format!("{a}\t{b}\n");
            }
            for s in stats.drain(..) {
                out += &format!("{s}\n");
            }
            out
        };
        assert!(r.counts == resort(&want_counts, 2), "{sample}: renumbered counts differ");
        assert!(r.grouped.unwrap() == resort(&want_grouped, 3), "{sample}: renumbered grouped differs");
        assert!(map.keys().any(|k| k.contains(".nic") || k.contains(".nnic")), "no novel models");
    }

    #[test]
    fn fullchr8_ont() {
        run_capture("SRR33046254");
    }

    #[test]
    fn fullchr8_isoseq() {
        run_capture("SRR32588720");
    }
}
