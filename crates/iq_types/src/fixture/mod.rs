//! Fixture machinery (schema `harness/iqcap/SCHEMA.md`).
//!
//! - [`schema`]: DTOs mirroring the JSON exactly; [`convert`]: DTO ↔ core types;
//!   [`ids`]: id-string maps.
//! - [`read_cases`] / [`CaseReader`]: gzip JSON Lines readers; [`write_case_line`] writes
//!   Python-identical bytes ([`pyjson`]).
//! - [`FixtureCtx`]: one capture directory (`fixtures/<sample>/<slice>/`) with cached side
//!   tables (`index/<sha1>.json`, `block/...`, `params.json`).
//! - [`UnitFixture`]: what every unit implements, one impl per `(unit, fn)`;
//!   [`Verdict`]; [`Registry`] + [`run_path`] power `iqrs-fixture run`.
//! - [`diff`]: generic field-level diff.
//!
//! ## Layout
//! A capture dir holds `units/<UNIT>[.<pid>].jsonl.gz`, `index/`, `block/`, `params.json`.
//! [`FixtureCtx::open`] accepts the capture dir or any file inside `units/`.

pub mod convert;
pub mod diff;
pub mod ids;
pub mod pyjson;
pub mod roundtrip;
pub mod schema;

pub use convert::{store_diff, ModelExtras, ReadRecordExtras, StoreExtras};
pub use diff::{diff_json, diff_values, Diff, DiffEntry, DiffOptions};
pub use ids::{IdMaps, TxSpace};
pub use pyjson::{to_py_string, PyF64};
pub use schema::{BlockHeader, Case, IndexDump, ReadRecordFixture};

use crate::index::IsoformIndex;
use crate::params::Params;
use crate::read::ReadRecord;
use crate::store::BlockStore;
use flate2::read::MultiGzDecoder;
use rustc_hash::FxHashMap;
use serde::de::DeserializeOwned;
use serde::Serialize;
use std::collections::BTreeMap;
use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, MutexGuard};

/// Schema version this crate implements.
pub const SCHEMA_VERSION: u32 = 1;

/// Error type of the fixture layer.
#[derive(Debug)]
pub struct FixError(pub String);

impl fmt::Display for FixError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for FixError {}
impl From<io::Error> for FixError {
    fn from(e: io::Error) -> Self {
        FixError(e.to_string())
    }
}
impl From<serde_json::Error> for FixError {
    fn from(e: serde_json::Error) -> Self {
        FixError(e.to_string())
    }
}
impl From<String> for FixError {
    fn from(e: String) -> Self {
        FixError(e)
    }
}
impl From<&str> for FixError {
    fn from(e: &str) -> Self {
        FixError(e.to_owned())
    }
}
pub type FixResult<T> = Result<T, FixError>;

// ------------------------------------------------------------------------------------------
// Reading / writing JSONL.gz
// ------------------------------------------------------------------------------------------

/// Opens a `.jsonl.gz` (or plain `.jsonl`) file as a line reader.
pub fn open_lines(path: &Path) -> FixResult<Box<dyn BufRead>> {
    let f = File::open(path).map_err(|e| FixError(format!("{}: {e}", path.display())))?;
    let gz = path.extension().is_some_and(|e| e == "gz");
    Ok(if gz {
        Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(f)))
    } else {
        Box::new(BufReader::with_capacity(1 << 20, f))
    })
}

/// Iterator over the non-empty lines of a fixture file, deserialized as `T`.
pub struct CaseReader<T> {
    lines: io::Lines<Box<dyn BufRead>>,
    path: PathBuf,
    line_no: usize,
    _t: std::marker::PhantomData<T>,
}

impl<T: DeserializeOwned> Iterator for CaseReader<T> {
    type Item = FixResult<T>;
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let line = match self.lines.next()? {
                Ok(l) => l,
                Err(e) => return Some(Err(FixError(format!("{}: {e}", self.path.display())))),
            };
            self.line_no += 1;
            if line.trim().is_empty() {
                continue;
            }
            return Some(
                serde_json::from_str(&line)
                    .map_err(|e| FixError(format!("{}:{}: {e}", self.path.display(), self.line_no))),
            );
        }
    }
}

/// Reads every line of `path` as `T` (e.g. `Case<In, Out>` or raw `Case`).
pub fn read_lines<T: DeserializeOwned>(path: &Path) -> FixResult<CaseReader<T>> {
    Ok(CaseReader { lines: open_lines(path)?.lines(), path: path.to_owned(), line_no: 0, _t: std::marker::PhantomData })
}

/// Reads the cases of a unit file with typed payloads.
pub fn read_cases<I: DeserializeOwned, O: DeserializeOwned>(path: &Path) -> FixResult<CaseReader<Case<I, O>>> {
    read_lines(path)
}

/// Writes one value as a Python-identical JSON line (`\n` terminated).
pub fn write_case_line<W: Write, T: Serialize>(mut w: W, v: &T) -> FixResult<()> {
    pyjson::to_py_writer(&mut w, v)?;
    w.write_all(b"\n")?;
    Ok(())
}

/// Unit files for `unit` under `<capdir>/units/`: `<UNIT>.jsonl.gz` or
/// `<UNIT>.<pid>.jsonl.gz` (schema §0.1), sorted.
pub fn unit_files(capdir: &Path, unit: &str) -> FixResult<Vec<PathBuf>> {
    let dir = capdir.join("units");
    let mut out = Vec::new();
    if !dir.is_dir() {
        return Ok(out);
    }
    for e in std::fs::read_dir(&dir)? {
        let p = e?.path();
        let name = p.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if let Some(rest) = name.strip_prefix(unit) {
            if rest == ".jsonl.gz"
                || (rest.starts_with('.')
                    && rest.ends_with(".jsonl.gz")
                    && rest[1..].starts_with(|c: char| c.is_ascii_digit()))
            {
                out.push(p);
            }
        }
    }
    out.sort();
    Ok(out)
}

/// Every `*.jsonl.gz` under `<capdir>/units/`, sorted.
pub fn all_unit_files(capdir: &Path) -> FixResult<Vec<PathBuf>> {
    let dir = capdir.join("units");
    let mut out = Vec::new();
    if dir.is_dir() {
        for e in std::fs::read_dir(&dir)? {
            let p = e?.path();
            if p.to_str().is_some_and(|s| s.ends_with(".jsonl.gz")) {
                out.push(p);
            }
        }
    }
    out.sort();
    Ok(out)
}

/// Parsed case key `<sample>:<slice>:<chr>:<block_no>:<local>` (schema §0.3). `None` for
/// keys of other shapes (U1a `annotation:<chr>`, U2 `<sample>:<region>:<n>`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CaseKey<'a> {
    pub sample: &'a str,
    pub slice: &'a str,
    pub chr: &'a str,
    /// `None` for `-` (outside a block).
    pub block_no: Option<u32>,
    pub local: &'a str,
}

impl<'a> CaseKey<'a> {
    pub fn parse(key: &'a str) -> Option<CaseKey<'a>> {
        let mut it = key.splitn(5, ':');
        let (sample, slice, chr, block, local) = (it.next()?, it.next()?, it.next()?, it.next()?, it.next()?);
        let block_no = if block == "-" { None } else { Some(block.parse().ok()?) };
        Some(CaseKey { sample, slice, chr, block_no, local })
    }
}

// ------------------------------------------------------------------------------------------
// FixtureCtx
// ------------------------------------------------------------------------------------------

/// A loaded block side table (schema §2.2).
#[derive(Debug)]
pub struct LoadedBlock {
    pub header: BlockHeader,
    /// Records as written (DTOs).
    pub fixtures: Vec<ReadRecordFixture>,
    /// Core records (reference transcript space).
    pub records: Vec<ReadRecord>,
    pub extras: Vec<ReadRecordExtras>,
    pub store: BlockStore,
}

impl LoadedBlock {
    pub fn read_names(&self) -> Vec<&str> {
        self.records.iter().map(|r| r.name.as_str()).collect()
    }
}

/// Context of one capture directory. Side tables are loaded lazily and cached; id maps are
/// shared by every conversion in the run (behind a mutex: `ctx.ids()`).
pub struct FixtureCtx {
    pub root: PathBuf,
    pub diff_opts: DiffOptions,
    ids: Mutex<IdMaps>,
    params_json: Mutex<Option<Arc<serde_json::Map<String, serde_json::Value>>>>,
    index_dumps: Mutex<FxHashMap<String, Arc<IndexDump>>>,
    indices: Mutex<FxHashMap<String, Arc<IsoformIndex>>>,
    blocks: Mutex<FxHashMap<(String, u32), Arc<LoadedBlock>>>,
}

impl FixtureCtx {
    /// Context rooted at `root` (a capture dir).
    pub fn new(root: impl Into<PathBuf>) -> Self {
        FixtureCtx {
            root: root.into(),
            diff_opts: DiffOptions::default(),
            ids: Mutex::new(IdMaps::new()),
            params_json: Mutex::new(None),
            index_dumps: Mutex::new(FxHashMap::default()),
            indices: Mutex::new(FxHashMap::default()),
            blocks: Mutex::new(FxHashMap::default()),
        }
    }

    /// Context for a path: a capture dir, or a file under `<capdir>/units/`.
    pub fn open(path: &Path) -> Self {
        Self::new(capture_dir_of(path))
    }

    /// Shared id maps.
    pub fn ids(&self) -> MutexGuard<'_, IdMaps> {
        self.ids.lock().expect("id maps poisoned")
    }

    /// Raw `params.json` (schema §2.3).
    pub fn params_json(&self) -> FixResult<Arc<serde_json::Map<String, serde_json::Value>>> {
        let mut g = self.params_json.lock().unwrap();
        if let Some(p) = &*g {
            return Ok(p.clone());
        }
        let path = self.root.join("params.json");
        let txt = std::fs::read_to_string(&path).map_err(|e| FixError(format!("{}: {e}", path.display())))?;
        let m: serde_json::Map<String, serde_json::Value> = serde_json::from_str(&txt)?;
        let a = Arc::new(m);
        *g = Some(a.clone());
        Ok(a)
    }

    /// Effective [`Params`] from `params.json`.
    pub fn params(&self) -> FixResult<Params> {
        Ok(Params::from_params_json(&*self.params_json()?)?)
    }

    /// The raw `index/<sha1>.json` dump.
    pub fn index_dump(&self, sha1: &str) -> FixResult<Arc<IndexDump>> {
        if let Some(d) = self.index_dumps.lock().unwrap().get(sha1) {
            return Ok(d.clone());
        }
        let path = self.root.join("index").join(format!("{sha1}.json"));
        let txt = std::fs::read_to_string(&path).map_err(|e| FixError(format!("{}: {e}", path.display())))?;
        let d: Arc<IndexDump> =
            Arc::new(serde_json::from_str(&txt).map_err(|e| FixError(format!("{}: {e}", path.display())))?);
        self.index_dumps.lock().unwrap().insert(sha1.to_owned(), d.clone());
        Ok(d)
    }

    /// `IsoformIndex` for `index/<sha1>.json` (converted with the shared id maps).
    pub fn index(&self, sha1: &str) -> FixResult<Arc<IsoformIndex>> {
        if let Some(i) = self.indices.lock().unwrap().get(sha1) {
            return Ok(i.clone());
        }
        let d = self.index_dump(sha1)?;
        let idx = Arc::new(IsoformIndex::from_dump(&d, &mut self.ids())?);
        self.indices.lock().unwrap().insert(sha1.to_owned(), idx.clone());
        Ok(idx)
    }

    /// Path of a block side table.
    pub fn block_path(&self, sample: &str, chr: &str, block_no: u32) -> PathBuf {
        self.root.join("block").join(sample).join(chr).join(format!("{block_no}.jsonl.gz"))
    }

    /// Loads (and caches) `block/<sample>/<chr>/<block_no>.jsonl.gz`.
    pub fn block(&self, b: &schema::BlockRef) -> FixResult<Arc<LoadedBlock>> {
        let key = (format!("{}\0{}", b.sample, b.chr), b.block_no);
        if let Some(x) = self.blocks.lock().unwrap().get(&key) {
            return Ok(x.clone());
        }
        let path = self.block_path(&b.sample, &b.chr, b.block_no);
        let lb = Arc::new(load_block(&path, &mut self.ids())?);
        self.blocks.lock().unwrap().insert(key, lb.clone());
        Ok(lb)
    }
}

/// Loads a block file: header line + `ReadRecord`s.
pub fn load_block(path: &Path, ids: &mut IdMaps) -> FixResult<LoadedBlock> {
    let mut lines = open_lines(path)?.lines();
    let first = lines.next().ok_or_else(|| FixError(format!("{}: empty block file", path.display())))??;
    let header: BlockHeader =
        serde_json::from_str(&first).map_err(|e| FixError(format!("{}:1: {e}", path.display())))?;
    let mut fixtures = Vec::with_capacity(header.n_reads as usize);
    for (k, l) in lines.enumerate() {
        let l = l?;
        if l.trim().is_empty() {
            continue;
        }
        fixtures.push(
            serde_json::from_str::<ReadRecordFixture>(&l)
                .map_err(|e| FixError(format!("{}:{}: {e}", path.display(), k + 2)))?,
        );
    }
    let (records, extras): (Vec<_>, Vec<_>) = fixtures.iter().map(|f| f.to_core(ids)).unzip();
    let mut store = BlockStore::from_records(&records);
    if records.is_empty() {
        store.region = header.region;
    }
    Ok(LoadedBlock { header, fixtures, records, extras, store })
}

/// Capture dir for a path: the dir itself if it has `units/`, `index/`, or `params.json`;
/// for a file under `units/`, its grandparent; otherwise the file's parent.
pub fn capture_dir_of(path: &Path) -> PathBuf {
    if path.is_dir() {
        return path.to_owned();
    }
    let parent = path.parent().unwrap_or(Path::new("."));
    if parent.file_name().is_some_and(|n| n == "units") {
        parent.parent().unwrap_or(Path::new(".")).to_owned()
    } else {
        parent.to_owned()
    }
}

// ------------------------------------------------------------------------------------------
// UnitFixture / Verdict
// ------------------------------------------------------------------------------------------

/// Result of judging one case.
#[derive(Clone, Debug, PartialEq)]
pub enum Verdict {
    Exact,
    /// Differs from the capture only for a known, documented reason (e.g. `"D-T1"`: peak tie
    /// order on the SIMD CPU that recorded the capture).
    Explained(&'static str),
    Mismatch(Diff),
}

impl Verdict {
    pub fn is_ok(&self) -> bool {
        !matches!(self, Verdict::Mismatch(_))
    }
}

/// A unit's fixture runner for one hooked function (`(UNIT, FN)`).
///
/// Implement `run` (call the unit's pure entry function and convert its result to the
/// fixture output DTO). Override `compare` for tolerances or set semantics; the default is
/// exact structural equality with a field-level diff.
pub trait UnitFixture {
    /// Unit id as in the envelope (`"U6"`, `"U10a"`, ...).
    const UNIT: &'static str;
    /// `in.fn` value handled by this runner.
    const FN: &'static str;
    type In: DeserializeOwned;
    type Out: DeserializeOwned + Serialize + PartialEq + fmt::Debug;

    fn run(input: &Self::In, ctx: &FixtureCtx) -> FixResult<Self::Out>;

    fn compare(got: &Self::Out, want: &Self::Out, ctx: &FixtureCtx) -> Verdict {
        if got == want {
            return Verdict::Exact;
        }
        match diff_values(got, want, &ctx.diff_opts) {
            Some(d) => Verdict::Mismatch(d),
            // equal as JSON but not by PartialEq (e.g. -0.0 vs 0.0): report as mismatch
            None => Verdict::Mismatch(Diff::error(format!("PartialEq differs: got {got:?}, want {want:?}"))),
        }
    }

    /// Runs and judges one case (panics and run errors become mismatches).
    fn judge(input: &Self::In, want: &Self::Out, ctx: &FixtureCtx) -> Verdict {
        match catch_unwind(AssertUnwindSafe(|| Self::run(input, ctx))) {
            Ok(Ok(got)) => Self::compare(&got, want, ctx),
            Ok(Err(e)) => Verdict::Mismatch(Diff::error(format!("run failed: {e}"))),
            Err(p) => Verdict::Mismatch(Diff::error(format!("panic: {}", panic_msg(&p)))),
        }
    }
}

fn panic_msg(p: &Box<dyn std::any::Any + Send>) -> String {
    if let Some(s) = p.downcast_ref::<&str>() {
        (*s).to_owned()
    } else if let Some(s) = p.downcast_ref::<String>() {
        s.clone()
    } else {
        "<non-string panic>".into()
    }
}

// ------------------------------------------------------------------------------------------
// Registry and runner
// ------------------------------------------------------------------------------------------

/// Outcome of one raw case.
#[derive(Clone, Debug, PartialEq)]
pub enum Outcome {
    Verdict(Verdict),
    /// `out` is a recorded Python exception; not run.
    PyException(String),
    /// Payload did not deserialize into the runner's DTOs.
    Decode(String),
}

/// Type-erased runner.
pub trait ErasedRunner: Send + Sync {
    fn unit(&self) -> &'static str;
    fn fn_name(&self) -> &'static str;
    fn judge_raw(&self, case: &Case, ctx: &FixtureCtx) -> Outcome;
}

struct Typed<U>(std::marker::PhantomData<fn() -> U>);

impl<U: UnitFixture> ErasedRunner for Typed<U> {
    fn unit(&self) -> &'static str {
        U::UNIT
    }
    fn fn_name(&self) -> &'static str {
        U::FN
    }
    fn judge_raw(&self, case: &Case, ctx: &FixtureCtx) -> Outcome {
        if let Some(e) = case.out.get("exception").and_then(|v| v.as_str()) {
            return Outcome::PyException(e.to_owned());
        }
        let input: U::In = match U::In::deserialize(&case.input) {
            Ok(x) => x,
            Err(e) => return Outcome::Decode(format!("in: {e}")),
        };
        let want: U::Out = match U::Out::deserialize(&case.out) {
            Ok(x) => x,
            Err(e) => return Outcome::Decode(format!("out: {e}")),
        };
        Outcome::Verdict(U::judge(&input, &want, ctx))
    }
}

/// Registry of fixture runners. Each crate exposes `pub fn register_fixtures(reg: &mut
/// Registry)` and calls `reg.add::<MyRunner>()` for each `(unit, fn)` it implements;
/// `iqrs-fixture` calls every crate's function.
#[derive(Default)]
pub struct Registry {
    runners: BTreeMap<(String, String), Box<dyn ErasedRunner>>,
}

impl Registry {
    pub fn new() -> Self {
        Self::default()
    }
    /// Registers `U`. Panics on a duplicate `(UNIT, FN)`.
    pub fn add<U: UnitFixture + 'static>(&mut self) -> &mut Self {
        let key = (U::UNIT.to_owned(), U::FN.to_owned());
        assert!(!self.runners.contains_key(&key), "duplicate fixture runner {key:?}");
        self.runners.insert(key, Box::new(Typed::<U>(std::marker::PhantomData)));
        self
    }
    pub fn get(&self, unit: &str, fn_name: &str) -> Option<&dyn ErasedRunner> {
        self.runners.get(&(unit.to_owned(), fn_name.to_owned())).map(|b| &**b)
    }
    /// Registered `(unit, fn)` pairs.
    pub fn keys(&self) -> impl Iterator<Item = (&str, &str)> {
        self.runners.keys().map(|(u, f)| (u.as_str(), f.as_str()))
    }
    pub fn len(&self) -> usize {
        self.runners.len()
    }
    pub fn is_empty(&self) -> bool {
        self.runners.is_empty()
    }
}

/// Per-`(unit, fn)` counts.
#[derive(Clone, Debug, Default)]
pub struct Tally {
    pub exact: u64,
    pub explained: BTreeMap<&'static str, u64>,
    pub mismatch: u64,
    /// Run errors / panics (`todo!()` included), counted inside `mismatch` too.
    pub errors: u64,
    pub py_exception: u64,
    pub decode_error: u64,
    pub unregistered: u64,
    /// First few mismatches: (case key, diff).
    pub examples: Vec<(String, Diff)>,
}

impl Tally {
    pub fn explained_total(&self) -> u64 {
        self.explained.values().sum()
    }
    pub fn total(&self) -> u64 {
        self.exact + self.explained_total() + self.mismatch + self.py_exception + self.decode_error + self.unregistered
    }
}

/// Summary of a run: `(unit, fn)` → tally.
#[derive(Clone, Debug, Default)]
pub struct Summary {
    pub by_key: BTreeMap<(String, String), Tally>,
    /// Max examples kept per key.
    pub max_examples: usize,
}

impl Summary {
    pub fn new(max_examples: usize) -> Self {
        Summary { by_key: BTreeMap::new(), max_examples }
    }
    pub fn add(&mut self, unit: &str, fn_name: &str, case: &str, o: Outcome) {
        let max = self.max_examples;
        let t = self.by_key.entry((unit.to_owned(), fn_name.to_owned())).or_default();
        let example = |d: Diff, t: &mut Tally| {
            if t.examples.len() < max {
                t.examples.push((case.to_owned(), d));
            }
        };
        match o {
            Outcome::Verdict(Verdict::Exact) => t.exact += 1,
            Outcome::Verdict(Verdict::Explained(d)) => *t.explained.entry(d).or_default() += 1,
            Outcome::Verdict(Verdict::Mismatch(d)) => {
                t.mismatch += 1;
                if d.error.is_some() {
                    t.errors += 1;
                }
                example(d, t);
            }
            Outcome::PyException(_) => t.py_exception += 1,
            Outcome::Decode(e) => {
                t.decode_error += 1;
                example(Diff::error(format!("decode: {e}")), t);
            }
        }
    }
    pub fn add_unregistered(&mut self, unit: &str, fn_name: &str) {
        self.by_key.entry((unit.to_owned(), fn_name.to_owned())).or_default().unregistered += 1;
    }
    /// True if no mismatches or decode errors (unregistered cases do not count as failures).
    pub fn all_ok(&self) -> bool {
        self.by_key.values().all(|t| t.mismatch == 0 && t.decode_error == 0)
    }
}

impl fmt::Display for Summary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "{:<6} {:<38} {:>8} {:>9} {:>8} {:>6} {:>6} {:>6} {:>7}",
            "unit", "fn", "exact", "explained", "mismatch", "error", "pyexc", "decode", "unreg"
        )?;
        for ((u, func), t) in &self.by_key {
            writeln!(
                f,
                "{:<6} {:<38} {:>8} {:>9} {:>8} {:>6} {:>6} {:>6} {:>7}",
                u,
                func,
                t.exact,
                t.explained_total(),
                t.mismatch,
                t.errors,
                t.py_exception,
                t.decode_error,
                t.unregistered
            )?;
            if !t.explained.is_empty() {
                let s: Vec<String> = t.explained.iter().map(|(d, n)| format!("{d}={n}")).collect();
                writeln!(f, "{:<6} {:<38}   explained by: {}", "", "", s.join(", "))?;
            }
        }
        for ((u, func), t) in &self.by_key {
            for (case, d) in &t.examples {
                writeln!(f, "--- {u} {func} {case}")?;
                write!(f, "{d}")?;
            }
        }
        Ok(())
    }
}

/// Runs every case of one unit file, adding to `summary`. Cases whose `unit` is not
/// `only_unit` (when given) are skipped.
pub fn run_file(
    reg: &Registry,
    path: &Path,
    ctx: &FixtureCtx,
    only_unit: Option<&str>,
    summary: &mut Summary,
) -> FixResult<()> {
    for case in read_lines::<Case>(path)? {
        let case = case?;
        if only_unit.is_some_and(|u| u != case.unit) {
            continue;
        }
        let fn_name = case.input.get("fn").and_then(|v| v.as_str()).unwrap_or("?").to_owned();
        match reg.get(&case.unit, &fn_name) {
            Some(r) => {
                let o = r.judge_raw(&case, ctx);
                summary.add(&case.unit, &fn_name, &case.case, o);
            }
            None => summary.add_unregistered(&case.unit, &fn_name),
        }
    }
    Ok(())
}

/// `iqrs-fixture run <unit> <file|dir>`: a file is run directly; a directory is treated as
/// a capture dir (its `units/<unit>*.jsonl.gz`) or, if it has no `units/`, searched
/// recursively for capture dirs. `unit = "all"` runs every unit file.
pub fn run_path(reg: &Registry, unit: &str, path: &Path, summary: &mut Summary) -> FixResult<()> {
    let only = if unit == "all" { None } else { Some(unit) };
    if path.is_file() {
        let ctx = FixtureCtx::open(path);
        return run_file(reg, path, &ctx, only, summary);
    }
    let mut caps = Vec::new();
    find_capture_dirs(path, &mut caps)?;
    for cap in caps {
        let ctx = FixtureCtx::new(&cap);
        let files = match only {
            Some(u) => unit_files(&cap, u)?,
            None => all_unit_files(&cap)?,
        };
        for f in files {
            run_file(reg, &f, &ctx, only, summary)?;
        }
    }
    Ok(())
}

/// Directories at or below `root` that contain a `units/` subdirectory.
pub fn find_capture_dirs(root: &Path, out: &mut Vec<PathBuf>) -> FixResult<()> {
    if root.join("units").is_dir() {
        out.push(root.to_owned());
        return Ok(());
    }
    if root.is_dir() {
        let mut subs: Vec<PathBuf> =
            std::fs::read_dir(root)?.filter_map(|e| e.ok().map(|e| e.path())).filter(|p| p.is_dir()).collect();
        subs.sort();
        for s in subs {
            find_capture_dirs(&s, out)?;
        }
    }
    Ok(())
}

/// Fixtures root from `$IQRS_FIXTURES`, if set.
pub fn fixtures_root() -> Option<PathBuf> {
    std::env::var_os("IQRS_FIXTURES").map(PathBuf::from)
}

use serde::Deserialize as _;
