//! U15a: `ModelFilter.pre_filter_transcripts/filter_transcripts/detect_similar_isoforms`
//! (`model_construction/model_filter.py:49-141, 194-232`), plus the graph coverage queries
//! (`intron_graph.py:204-209, 752-781`) and the `ModelStore` bookkeeping they use
//! (`context.py:95-108`).
//!
//! Notes:
//! - **defaultdict accesses.** 4.1 reads `transcript_read_ids[id]` (and `internal_counter[id]`)
//!   through `defaultdict`s, which inserts empty entries. Accesses that can leave a visible entry
//!   (`transcript_read_ids`, compared in order) are mirrored with `entry().or_default()`.
//!   `internal_counter` zero entries are invisible (the store diff drops them, like U17+).
//! - **Graph defaultdict accesses.** `get_max_component_coverage` and `is_monointron` index
//!   `outgoing_edges`/`incoming_edges` by the model's path introns. A path intron that is not
//!   a key in either map becomes an edgeless vertex of the live graph, which every later
//!   `get_overlapping_component_max_coverage` scan then visits (its component is itself, with
//!   its `clustered_introns` count). [`GraphView`] tracks these vertices and adds them to the
//!   [`iq_graph::query`] result; the captured graph itself is not modified.
//! - **detect_similar_isoforms.** 4.1 builds a `GeneInfo` + assigner per outer model and
//!   re-profiles every candidate inner model against it, and `filter_transcripts` calls it
//!   twice over (mostly) the same models. Here the outer index and profile builder are built
//!   at most once per outer model and only when an uncached pair needs it; each pair's
//!   assignment class is cached in a [`SimilarCache`] keyed by both models, valid while the
//!   inputs that determine it (both models' exons/strand/gene and the inner model's polyA
//!   info) are unchanged. The second call therefore re-assigns only pairs touching a model whose
//!   ends `correct_novel_transcript_ends` changed. Pairs whose spans do not overlap are not
//!   assigned at all: no read block overlaps a split exon, so 4.1's assignment is
//!   `noninformative`, which neither substitutes nor counts as unconfirmed (exact).
//!   A single profile per inner model shared across outer models is *not* exact: the
//!   profile depends on the outer index's gene region and split exons.
//! - **mapping_quality** over a model without reads raises `ZeroDivisionError` in 4.1; here it
//!   returns 0.0 (the model is then dropped as low quality).

use crate::ctx::Ctx;
use crate::fl::is_matching_assignment;
use iq_assign::profile::{ProfileBuilder, Scratch};
use iq_types::indexmap::IndexMap;
use iq_types::{
    vertex_code, AssignmentType, CombinedProfile, GeneKey, IntronGraph, Iv, ModelIdx, ModelKind, ModelStore, PolyAInfo, Pos,
    ReadIdx, Strand, TranscriptModel, Vertex, NO_POLYA,
};
use std::collections::{HashMap, HashSet};

/// The end-processor calls `filter_transcripts` makes (U15b `TranscriptEndProcessor`).
/// Replayed from captured U15b cases in the U15a fixture runner.
pub trait EndProcessor {
    /// `correct_novel_transcript_ends(model, store.transcript_read_ids[model.id])`.
    fn correct_novel_transcript_ends(&self, model: &mut TranscriptModel, reads: &[ReadIdx], ctx: &Ctx<'_>);
    /// `_mark_terminal_confirmation(model)`; `reads` = `store.transcript_read_ids[model.id]`.
    fn mark_terminal_confirmation(&self, model: &mut TranscriptModel, reads: &[ReadIdx], ctx: &Ctx<'_>);
    /// `_add_known_alternative_end_models()`.
    fn add_known_alternative_end_models(&self, ms: &mut ModelStore, ctx: &Ctx<'_>);
    /// `_drop_duplicate_alt_end_models()`.
    fn drop_duplicate_alt_end_models(&self, ms: &mut ModelStore, ctx: &Ctx<'_>);
}

// ------------------------------------------------------------------------------------------
// Store bookkeeping
// ------------------------------------------------------------------------------------------

/// `internal_counter[id]` (defaultdict read; the zero entry is not materialised).
#[inline]
fn counter(ms: &ModelStore, key: ModelIdx) -> i64 {
    ms.internal_counter.get(&key).copied().unwrap_or(0) as i64
}

/// `ModelStore.mapping_quality` (`context.py:104-108`): mean MAPQ of the bound reads, in f64.
/// Materialises the `transcript_read_ids` entry like the 4.1 defaultdict read. Empty → 0.0
/// (4.1: ZeroDivisionError).
pub fn mapping_quality(ms: &mut ModelStore, ctx: &Ctx<'_>, key: ModelIdx) -> f64 {
    let reads = ms.transcript_read_ids.entry(key).or_default();
    if reads.is_empty() {
        return 0.0;
    }
    let sum: i64 = reads.iter().map(|&r| ctx.store.view(r as usize).mapq() as i64).sum();
    sum as f64 / reads.len() as f64
}

// ------------------------------------------------------------------------------------------
// pre_filter_transcripts
// ------------------------------------------------------------------------------------------

/// `pre_filter_transcripts` (`model_filter.py:49-78`), incl. the chrM hack (KEEP, design §7.6).
pub fn pre_filter(ms: &mut ModelStore, ctx: &Ctx<'_>) {
    let mut vals: Vec<i64> = ms.models.iter().filter(|m| m.exons.len() <= 2).map(|m| counter(ms, m.key)).collect();
    vals.sort_unstable_by(|a, b| b.cmp(a));
    // dirty hack to avoid slow assignment in chrM (4.1 comment)
    let cutoff = if vals.len() > 50 { vals[50] } else { ctx.params.model.min_novel_count as i64 };
    let mapq_cutoff = ctx.params.simple_models_mapq_cutoff as f64;

    let models = std::mem::take(&mut ms.models);
    let mut kept = Vec::with_capacity(models.len());
    for m in models {
        if m.exons.len() > 2 || matches!(m.kind, ModelKind::Known(_)) {
            kept.push(m);
            continue;
        }
        if counter(ms, m.key) < cutoff || mapping_quality(ms, ctx, m.key) < mapq_cutoff {
            ms.delete_from_storage(m.key);
            continue;
        }
        kept.push(m);
    }
    ms.models = kept;
}

// ------------------------------------------------------------------------------------------
// Graph queries on the live graph
// ------------------------------------------------------------------------------------------

#[inline]
fn is_terminal_vertex(v: &Vertex) -> bool {
    matches!(v.0, vertex_code::POLYA | vertex_code::READ_END | vertex_code::TSS_RIGHT)
}
#[inline]
fn is_starting_vertex(v: &Vertex) -> bool {
    matches!(v.0, vertex_code::POLYT | vertex_code::READ_START | vertex_code::TSS_LEFT)
}

/// The captured graph plus the edgeless vertices 4.1's live graph gains from `defaultdict`
/// reads during filtering (see module docs).
pub struct GraphView<'g> {
    pub g: &'g IntronGraph,
    /// Vertices touched that are keys of neither adjacency map, in first-touch order.
    extra: Vec<Vertex>,
}

impl<'g> GraphView<'g> {
    pub fn new(g: &'g IntronGraph) -> Self {
        GraphView { g, extra: Vec::new() }
    }

    fn count(&self, v: &Vertex) -> u32 {
        v.as_intron().and_then(|iv| self.g.clustered_introns.get(&iv).copied()).unwrap_or(0)
    }

    /// Records a `outgoing_edges[v]` / `incoming_edges[v]` read.
    fn touch(&mut self, v: Vertex) {
        if !self.g.outgoing.contains_key(&v) && !self.g.incoming.contains_key(&v) && !self.extra.contains(&v) {
            self.extra.push(v);
        }
    }

    /// `get_max_component_coverage(path)` (`intron_graph.py:752-756`).
    pub fn max_component_coverage(&mut self, path: &[Vertex]) -> u32 {
        let comp = iq_graph::query::connected_component(self.g, path);
        let mut best = 0;
        for v in &comp {
            best = best.max(self.count(v));
        }
        // only the path's own introns can be absent from both maps; the rest were reached by edges
        for &v in path {
            if v.0 >= 0 {
                self.touch(v);
            }
        }
        best
    }

    /// `is_monointron(v)` (`intron_graph.py:204-209`).
    pub fn is_monointron(&mut self, v: Vertex) -> bool {
        self.touch(v);
        let no_out = self.g.outgoing.get(&v).is_none_or(|o| o.iter().all(is_terminal_vertex));
        let no_in = self.g.incoming.get(&v).is_none_or(|o| o.iter().all(is_starting_vertex));
        no_in && no_out
    }

    /// `get_overlapping_component_max_coverage((start, end))` (`intron_graph.py:758-781`).
    pub fn overlapping_component_max_coverage(&self, iv: Iv) -> u32 {
        let mut best = iq_graph::query::overlapping_component_max_coverage(self.g, iv);
        for v in &self.extra {
            // edgeless: span (v0, v1), component {v} (skipped for terminal vertices)
            if v.0 >= 0 && !((iv.e as i64) < v.0 || (iv.s as i64) > v.1) {
                best = best.max(self.count(v));
            }
        }
        best
    }
}

// ------------------------------------------------------------------------------------------
// detect_similar_isoforms
// ------------------------------------------------------------------------------------------

/// The assignment-type facts `detect_similar_isoforms` uses.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
struct PairVerdict {
    matching: bool,
    /// `unique_minor_difference` or `inconsistent_non_intronic`.
    truncation_like: bool,
}

const NONINFORMATIVE: PairVerdict = PairVerdict { matching: false, truncation_like: false };

/// Inputs that determine a pair's assignment.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Snap {
    exons: Vec<Iv>,
    strand: Strand,
    gene: GeneKey,
    polya: PolyAInfo,
}

/// Per-pair assignment cache shared by the `detect_similar_isoforms` calls of one
/// `filter_transcripts` (see module docs).
#[derive(Default, Debug)]
pub struct SimilarCache {
    pairs: HashMap<(ModelIdx, ModelIdx), PairVerdict>,
    snaps: HashMap<ModelIdx, Snap>,
    /// Number of assignments computed (for tests / stats).
    pub n_assigned: usize,
}

impl SimilarCache {
    /// Drops cached pairs of models whose inputs changed since they were cached.
    fn refresh(&mut self, models: &[TranscriptModel]) {
        let mut dirty: HashSet<ModelIdx> = HashSet::new();
        for m in models {
            let s = Snap { exons: m.exons.clone(), strand: m.strand, gene: m.gene, polya: inner_polya(m) };
            match self.snaps.get(&m.key) {
                Some(old) if *old == s => {}
                Some(_) => {
                    dirty.insert(m.key);
                    self.snaps.insert(m.key, s);
                }
                None => {
                    self.snaps.insert(m.key, s);
                }
            }
        }
        if !dirty.is_empty() {
            self.pairs.retain(|(a, b), _| !dirty.contains(a) && !dirty.contains(b));
        }
    }
}

/// The inner model's `PolyAInfo` (`model_filter.py:209-214`). Dead in practice (D-M4:
/// `intron_path` holds introns only), ported because it costs nothing.
fn inner_polya(m: &TranscriptModel) -> PolyAInfo {
    match (m.intron_path.first(), m.intron_path.last()) {
        (Some(f), _) if f.0 == vertex_code::POLYT => PolyAInfo { ext_t: f.1 as Pos, ..NO_POLYA },
        (_, Some(l)) if l.0 == vertex_code::POLYA => PolyAInfo { ext_a: l.1 as Pos, ..NO_POLYA },
        _ => NO_POLYA,
    }
}

#[inline]
fn spans_overlap(a: &TranscriptModel, b: &TranscriptModel) -> bool {
    !(a.end() < b.start() || a.start() > b.end())
}

/// `detect_similar_isoforms(models, remove_unconfirmed)` → `to_substitute` (inner → outer) in
/// 4.1 dict order.
pub fn detect_similar_isoforms(models: &[TranscriptModel], ctx: &Ctx<'_>, remove_unconfirmed: bool) -> IndexMap<ModelIdx, ModelIdx> {
    detect_similar_cached(models, ctx, remove_unconfirmed, &mut SimilarCache::default())
}

/// [`detect_similar_isoforms`] with a cache shared across calls.
pub fn detect_similar_cached(
    models: &[TranscriptModel],
    ctx: &Ctx<'_>,
    remove_unconfirmed: bool,
    cache: &mut SimilarCache,
) -> IndexMap<ModelIdx, ModelIdx> {
    cache.refresh(models);
    let use_tss = ctx.params.use_tss_model;
    let mut ts: IndexMap<ModelIdx, ModelIdx> = IndexMap::new();
    let mut cands: Vec<usize> = Vec::new();
    let mut scratch = Scratch::default();
    let mut prof = CombinedProfile::default();
    for outer in models {
        if outer.exons.len() <= 2 || ts.contains_key(&outer.key) {
            continue;
        }
        // `to_substitute` only grows by inner models already visited in this loop, so testing
        // membership up front is the same as testing it on visit.
        cands.clear();
        let mut need = false;
        for (i, m) in models.iter().enumerate() {
            if matches!(m.kind, ModelKind::Known(_))
                || m.key == outer.key
                || ts.contains_key(&m.key)
                || m.exons.len() == 1
                || m.intron_path.is_empty()
                || m.exons.len() > outer.exons.len()
            {
                continue;
            }
            cands.push(i);
            need |= spans_overlap(m, outer) && !cache.pairs.contains_key(&(outer.key, m.key));
        }
        if need {
            let index = iq_annot::index_from_models(ctx.index.chr, std::slice::from_ref(outer), ctx.params.matching.delta);
            let builder = ProfileBuilder::new(&index, &ctx.params.matching);
            for &i in &cands {
                let m = &models[i];
                if !spans_overlap(m, outer) || cache.pairs.contains_key(&(outer.key, m.key)) {
                    continue;
                }
                builder.build_into(&m.exons, inner_polya(m), &mut scratch, &mut prof);
                let a = ctx.assigner.assign(&index, &prof, false, m.key as u64);
                cache.n_assigned += 1;
                let v = match a {
                    Some(a) => PairVerdict {
                        matching: is_matching_assignment(&a),
                        truncation_like: matches!(
                            a.atype,
                            AssignmentType::UniqueMinorDifference | AssignmentType::InconsistentNonIntronic
                        ),
                    },
                    None => NONINFORMATIVE,
                };
                cache.pairs.insert((outer.key, m.key), v);
            }
        }
        for &i in &cands {
            let m = &models[i];
            let v = if spans_overlap(m, outer) { cache.pairs[&(outer.key, m.key)] } else { NONINFORMATIVE };
            if v.matching
                || (remove_unconfirmed && v.truncation_like && !(m.polya_confirmed || (use_tss && m.tss_confirmed)))
            {
                ts.insert(m.key, outer.key);
            }
        }
    }
    ts
}

// ------------------------------------------------------------------------------------------
// filter_transcripts
// ------------------------------------------------------------------------------------------

/// `filter_transcripts` (`model_filter.py:80-141`).
pub fn filter(ms: &mut ModelStore, ctx: &Ctx<'_>, ends: &dyn EndProcessor) {
    let p = ctx.params;
    let mut cache = SimilarCache::default();
    let models = std::mem::take(&mut ms.models);
    let ts = detect_similar_cached(&models, ctx, false, &mut cache);
    let mut gv = GraphView::new(ctx.graph);
    let min_novel = p.model.min_novel_count as f64;
    let mapq_cutoff = p.simple_models_mapq_cutoff as f64;

    let mut pre = Vec::with_capacity(models.len());
    for mut model in models {
        if matches!(model.kind, ModelKind::Known(_)) {
            pre.push(model);
            continue;
        }
        let key = model.key;
        let path = &model.intron_path;
        let cov = gv.max_component_coverage(path);
        let cutoff = if cov == 0 || path.is_empty() || (path.len() == 1 && gv.is_monointron(path[0])) {
            let cov = gv.overlapping_component_max_coverage(model.span());
            min_novel.max(p.model.min_mono_count_rel * cov as f64)
        } else {
            min_novel.max(p.model.min_novel_count_rel * cov as f64)
        };
        if ts.contains_key(&key) || (counter(ms, key) as f64) < cutoff {
            ms.delete_from_storage(key);
            continue;
        }
        if model.exons.len() <= 2 && mapping_quality(ms, ctx, key) < mapq_cutoff {
            ms.delete_from_storage(key);
            continue;
        }
        let reads = ms.transcript_read_ids.entry(key).or_default().clone();
        ends.correct_novel_transcript_ends(&mut model, &reads, ctx);
        let reads = ms.transcript_read_ids.entry(key).or_default().clone();
        ends.mark_terminal_confirmation(&mut model, &reads, ctx);
        pre.push(model);
    }

    let ts = detect_similar_cached(&pre, ctx, true, &mut cache);
    let mut kept = Vec::with_capacity(pre.len());
    for m in pre {
        if !matches!(m.kind, ModelKind::Known(_)) && ts.contains_key(&m.key) {
            ms.delete_from_storage(m.key);
            continue;
        }
        kept.push(m);
    }
    ms.models = kept;
    if p.create_alt_ts_nics || p.novel_apa {
        ends.add_known_alternative_end_models(ms, ctx);
        ends.drop_duplicate_alt_end_models(ms, ctx);
    }
}

// ------------------------------------------------------------------------------------------
// Fixture runners
// ------------------------------------------------------------------------------------------

pub(crate) fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    reg.add::<fixture::U15aPreFilter>();
    reg.add::<fixture::U15aFilter>();
    reg.add::<fixture::U15aDetectSimilar>();
}

pub(crate) mod fixture {
    use super::*;
    use crate::fixture_util::{compare_stores, StageInput};
    use crate::fl::LiveAssigner;
    use iq_assign::AssignParams;
    use iq_types::fixture::ids::IdMaps;
    use iq_types::fixture::schema::{
        BlockRef, Case, ModelOut, StageModelIn, StageStoreIn, StoreOut, U15aDetectSimilarIn, U15aDetectSimilarOut,
    };
    use iq_types::fixture::{read_lines, unit_files, FixError, FixResult, FixtureCtx, UnitFixture, Verdict};
    use iq_types::{serde_json, Params, Predictions};
    use std::cell::RefCell;
    use std::path::PathBuf;
    use std::sync::{Arc, Mutex, OnceLock};

    fn live_assigner(p: &Params) -> LiveAssigner {
        LiveAssigner { params: AssignParams::new(p.matching.clone(), p.splice.clone(), false) }
    }

    fn with_ctx<R>(si: &StageInput, f: impl FnOnce(&Ctx<'_>) -> R) -> R {
        let preds = Predictions::default();
        let assigner = live_assigner(&si.params);
        let c = Ctx {
            params: &si.params,
            index: &si.index,
            store: &si.block.store,
            graph: &si.graph,
            paths: &si.paths,
            predictions: &preds,
            assigner: &assigner,
        };
        f(&c)
    }

    pub struct U15aPreFilter;
    impl UnitFixture for U15aPreFilter {
        const UNIT: &'static str = "U15a";
        const FN: &'static str = "pre_filter_transcripts";
        type In = StageStoreIn;
        type Out = StoreOut;
        fn run(input: &StageStoreIn, ctx: &FixtureCtx) -> FixResult<StoreOut> {
            let si = StageInput::load(input, ctx)?;
            let mut ms = si.store.clone();
            with_ctx(&si, |c| pre_filter(&mut ms, c));
            Ok(StoreOut { store: si.store_out(&ms, &si.detected_known) })
        }
        fn compare(got: &StoreOut, want: &StoreOut, ctx: &FixtureCtx) -> Verdict {
            compare_stores(&got.store, &want.store, ctx)
        }
    }

    // ---- U15b replay ------------------------------------------------------------------

    type CaseMap = Arc<HashMap<String, Case>>;

    /// All U15b cases of a capture dir by case key (cached).
    fn u15b_cases(ctx: &FixtureCtx) -> FixResult<CaseMap> {
        static CACHE: OnceLock<Mutex<HashMap<PathBuf, CaseMap>>> = OnceLock::new();
        let mut g = CACHE.get_or_init(Default::default).lock().unwrap();
        if let Some(m) = g.get(&ctx.root) {
            return Ok(m.clone());
        }
        let mut m = HashMap::new();
        for f in unit_files(&ctx.root, "U15b")? {
            for case in read_lines::<Case>(&f)? {
                let case = case?;
                m.insert(case.case.clone(), case);
            }
        }
        let m = Arc::new(m);
        g.insert(ctx.root.clone(), m.clone());
        Ok(m)
    }

    /// Replays the captured U15b calls nested in a `filter_transcripts` case (by case key
    /// `<block prefix>:<fn>:<model id>.<n>` / `:n.<n>`), checking that each call's input equals
    /// the captured one. Input mismatches and missing cases are collected in `errors`.
    pub(crate) struct ReplayEnds {
        prefix: String,
        cases: CaseMap,
        pub(crate) ids: RefCell<IdMaps>,
        names: Vec<String>,
        seen: RefCell<HashMap<String, u32>>,
        pub(crate) errors: RefCell<Vec<String>>,
    }

    impl ReplayEnds {
        pub(crate) fn new(b: &BlockRef, cases: CaseMap, ids: IdMaps, names: Vec<String>) -> Self {
            ReplayEnds {
                prefix: format!("{}:{}:{}:{}", b.sample, b.slice, b.chr, b.block_no),
                cases,
                ids: RefCell::new(ids),
                names,
                seen: RefCell::default(),
                errors: RefCell::default(),
            }
        }

        fn case(&self, fn_: &str, local: &str) -> Option<&Case> {
            let base = format!("{}:{fn_}:{local}", self.prefix);
            let n = {
                let mut s = self.seen.borrow_mut();
                let e = s.entry(base.clone()).or_insert(0);
                *e += 1;
                *e - 1
            };
            let key = format!("{base}.{n}");
            let c = self.cases.get(&key);
            if c.is_none() {
                self.errors.borrow_mut().push(format!("U15b case {key} not found"));
            }
            c
        }

        fn model_id(&self, m: &TranscriptModel) -> String {
            let ids = self.ids.borrow();
            match m.kind {
                ModelKind::Known(t) => ids.txs.resolve(t).to_owned(),
                _ => ids.models.resolve(m.key).to_owned(),
            }
        }

        fn replay_model(&self, fn_: &str, model: &mut TranscriptModel, reads: &[ReadIdx]) {
            let id = self.model_id(model);
            let Some(c) = self.case(fn_, &id) else { return };
            let parsed = serde_json::from_value::<StageModelIn>(c.input.clone())
                .and_then(|i| serde_json::from_value::<ModelOut>(c.out.clone()).map(|o| (i, o)));
            let (i, o) = match parsed {
                Ok(x) => x,
                Err(e) => return self.errors.borrow_mut().push(format!("{}: {e}", c.case)),
            };
            let mut ids = self.ids.borrow_mut();
            let (want_in, _) = i.model.to_core(&mut ids);
            if want_in != *model {
                self.errors.borrow_mut().push(format!("{}: input model differs: got {model:?}, want {want_in:?}", c.case));
            }
            if i.reads != reads {
                self.errors.borrow_mut().push(format!("{}: input reads differ: got {reads:?}, want {:?}", c.case, i.reads));
            }
            *model = o.model.to_core(&mut ids).0;
        }

        fn replay_store(&self, fn_: &str, ms: &mut ModelStore) {
            let Some(c) = self.case(fn_, "n") else { return };
            let parsed = serde_json::from_value::<StageStoreIn>(c.input.clone())
                .and_then(|i| serde_json::from_value::<StoreOut>(c.out.clone()).map(|o| (i, o)));
            let (i, o) = match parsed {
                Ok(x) => x,
                Err(e) => return self.errors.borrow_mut().push(format!("{}: {e}", c.case)),
            };
            let names: Vec<&str> = self.names.iter().map(String::as_str).collect();
            let mut ids = self.ids.borrow_mut();
            match i.store.to_core(&mut ids, &names) {
                Ok((want_in, _)) => {
                    if let Some(d) = store_core_diff(ms, &want_in) {
                        self.errors.borrow_mut().push(format!("{}: input store differs: {d}", c.case));
                    }
                }
                Err(e) => self.errors.borrow_mut().push(format!("{}: {e}", c.case)),
            }
            match o.store.to_core(&mut ids, &names) {
                Ok((out, _)) => *ms = out,
                Err(e) => self.errors.borrow_mut().push(format!("{}: {e}", c.case)),
            }
        }
    }

    /// Store comparison modulo `internal_counter` zero entries and `next_key`.
    fn store_core_diff(got: &ModelStore, want: &ModelStore) -> Option<String> {
        let nz = |s: &ModelStore| s.internal_counter.iter().filter(|e| *e.1 != 0).map(|(k, v)| (*k, *v)).collect::<Vec<_>>();
        if got.models != want.models {
            return Some(format!("models: got {:?}, want {:?}", got.models, want.models));
        }
        if got.transcript_read_ids != want.transcript_read_ids {
            return Some(format!("transcript_read_ids: got {:?}, want {:?}", got.transcript_read_ids, want.transcript_read_ids));
        }
        if nz(got) != nz(want) {
            return Some(format!("internal_counter: got {:?}, want {:?}", nz(got), nz(want)));
        }
        if got.read_assignment_counts != want.read_assignment_counts {
            return Some("read_assignment_counts".into());
        }
        if got.construction_assignment != want.construction_assignment {
            return Some("construction_assignment".into());
        }
        if got.nonfl_known_reads != want.nonfl_known_reads {
            return Some("nonfl_known_reads".into());
        }
        None
    }

    impl EndProcessor for ReplayEnds {
        fn correct_novel_transcript_ends(&self, model: &mut TranscriptModel, reads: &[ReadIdx], _: &Ctx<'_>) {
            self.replay_model("correct_novel_transcript_ends", model, reads)
        }
        fn mark_terminal_confirmation(&self, model: &mut TranscriptModel, reads: &[ReadIdx], _: &Ctx<'_>) {
            self.replay_model("_mark_terminal_confirmation", model, reads)
        }
        fn add_known_alternative_end_models(&self, ms: &mut ModelStore, _: &Ctx<'_>) {
            self.replay_store("_add_known_alternative_end_models", ms)
        }
        fn drop_duplicate_alt_end_models(&self, ms: &mut ModelStore, _: &Ctx<'_>) {
            self.replay_store("_drop_duplicate_alt_end_models", ms)
        }
    }

    pub struct U15aFilter;
    impl UnitFixture for U15aFilter {
        const UNIT: &'static str = "U15a";
        const FN: &'static str = "filter_transcripts";
        type In = StageStoreIn;
        type Out = StoreOut;
        fn run(input: &StageStoreIn, ctx: &FixtureCtx) -> FixResult<StoreOut> {
            let mut si = StageInput::load(input, ctx)?;
            let mut ms = si.store.clone();
            let names = si.block.read_names().into_iter().map(str::to_owned).collect();
            let replay = ReplayEnds::new(&input.block, u15b_cases(ctx)?, si.ids.clone(), names);
            with_ctx(&si, |c| filter(&mut ms, c, &replay));
            let errors = replay.errors.take();
            if !errors.is_empty() {
                return Err(FixError(errors.join("; ")));
            }
            // models created by the replayed U15b calls carry their 4.1 ids
            si.ids = replay.ids.into_inner();
            Ok(StoreOut { store: si.store_out(&ms, &si.detected_known) })
        }
        fn compare(got: &StoreOut, want: &StoreOut, ctx: &FixtureCtx) -> Verdict {
            compare_stores(&got.store, &want.store, ctx)
        }
    }

    pub struct U15aDetectSimilar;
    impl UnitFixture for U15aDetectSimilar {
        const UNIT: &'static str = "U15a";
        const FN: &'static str = "detect_similar_isoforms";
        type In = U15aDetectSimilarIn;
        type Out = U15aDetectSimilarOut;
        fn run(input: &U15aDetectSimilarIn, ctx: &FixtureCtx) -> FixResult<U15aDetectSimilarOut> {
            let params = ctx.params()?;
            let block = ctx.block(&input.block)?;
            let index = ctx.index(&input.index)?;
            let mut ids = ctx.ids().clone();
            ids.reset_block();
            let models: Vec<TranscriptModel> = input.models.iter().map(|m| m.to_core(&mut ids).0).collect();
            let name: HashMap<ModelIdx, &str> = models.iter().zip(&input.models).map(|(m, f)| (m.key, f.id.as_str())).collect();
            let graph = IntronGraph::default();
            let paths = iq_types::PathStorage::default();
            let preds = Predictions::default();
            let assigner = live_assigner(&params);
            let c = Ctx {
                params: &params,
                index: &index,
                store: &block.store,
                graph: &graph,
                paths: &paths,
                predictions: &preds,
                assigner: &assigner,
            };
            let ts = detect_similar_isoforms(&models, &c, input.remove_unconfirmed);
            Ok(U15aDetectSimilarOut { to_substitute: ts.iter().map(|(a, b)| (name[a].to_owned(), name[b].to_owned())).collect() })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fl::LiveAssigner;
    use iq_assign::AssignParams;
    use iq_types::{Assignment, BlockStore, IsoformIndex, Params, PathStorage, Predictions, ReadRecord};

    fn read(idx: u32, mapq: u8) -> ReadRecord {
        ReadRecord {
            idx,
            name: format!("r{idx}"),
            exons: vec![Iv::new(10, 20)],
            corrected_exons: vec![Iv::new(10, 20)],
            boundary_deletions: vec![],
            polya: NO_POLYA,
            polya_found: false,
            strand: Strand::Plus,
            mapped_strand: Strand::Plus,
            mapq,
            multimapper: false,
            group: 0,
            assignment: Assignment::new(AssignmentType::Noninformative, Default::default()),
            genomic_region: Iv::new(1, 100000),
        }
    }

    fn model(key: ModelIdx, kind: ModelKind, exons: &[(Pos, Pos)]) -> TranscriptModel {
        let exons: Vec<Iv> = exons.iter().map(|&(s, e)| Iv::new(s, e)).collect();
        let intron_path = exons.windows(2).map(|w| Vertex(w[0].e as i64 + 1, w[1].s as i64 - 1)).collect();
        TranscriptModel {
            key,
            strand: Strand::Plus,
            gene: GeneKey::Novel(0),
            kind,
            exons,
            intron_path,
            polya_confirmed: false,
            tss_confirmed: false,
            source_known: None,
        }
    }

    struct Env {
        params: Params,
        index: IsoformIndex,
        store: BlockStore,
        graph: IntronGraph,
        paths: PathStorage,
        preds: Predictions,
        assigner: LiveAssigner,
    }
    impl Env {
        fn new(recs: &[ReadRecord], graph: IntronGraph) -> Self {
            let params = Params::ont();
            let assigner = LiveAssigner { params: AssignParams::new(params.matching.clone(), params.splice.clone(), false) };
            Env {
                params,
                index: IsoformIndex::default(),
                store: BlockStore::from_records(recs),
                graph,
                paths: PathStorage::default(),
                preds: Predictions::default(),
                assigner,
            }
        }
        fn ctx(&self) -> Ctx<'_> {
            Ctx {
                params: &self.params,
                index: &self.index,
                store: &self.store,
                graph: &self.graph,
                paths: &self.paths,
                predictions: &self.preds,
                assigner: &self.assigner,
            }
        }
    }

    #[test]
    fn chrm_hack_cutoff() {
        // 52 novel mono-exon models with counts 1..=52, one read each in transcript_read_ids
        let n = 52u32;
        let recs: Vec<ReadRecord> = (0..n).map(|i| read(i, 60)).collect();
        let env = Env::new(&recs, IntronGraph::default());
        let mut ms = ModelStore::new(n as usize);
        for i in 0..n {
            let k = ms.new_key();
            ms.models.push(model(k, ModelKind::Nnic, &[(100, 200)]));
            ms.save_assigned_read(i, k);
            ms.internal_counter[&k] = i as i32 + 1;
        }
        assert_eq!(env.params.model.min_novel_count, 3);
        pre_filter(&mut ms, &env.ctx());
        // cutoff = 51st largest count = 2 (not min_novel_count = 3): only count 1 is dropped
        assert_eq!(ms.models.len(), 51);
        assert_eq!(ms.models[0].key, 1);
        assert!(!ms.internal_counter.contains_key(&0));
        assert_eq!(ms.read_assignment_counts[0], 0);

        // 50 models: no hack, min_novel_count applies
        let mut ms2 = ModelStore::new(n as usize);
        for i in 0..50 {
            let k = ms2.new_key();
            ms2.models.push(model(k, ModelKind::Nnic, &[(100, 200)]));
            ms2.save_assigned_read(i, k);
            ms2.internal_counter[&k] = i as i32 + 1;
        }
        pre_filter(&mut ms2, &env.ctx());
        assert_eq!(ms2.models.len(), 48);
    }

    #[test]
    fn low_mapq_simple_model_dropped() {
        let recs = [read(0, 10), read(1, 20), read(2, 60)];
        let env = Env::new(&recs, IntronGraph::default());
        let mut ms = ModelStore::new(3);
        let a = ms.new_key();
        ms.models.push(model(a, ModelKind::Nnic, &[(100, 200)]));
        for r in 0..3 {
            ms.save_assigned_read(r, a);
        }
        // mean 30.0 is not < 30
        assert_eq!(mapping_quality(&mut ms, &env.ctx(), a), 30.0);
        pre_filter(&mut ms, &env.ctx());
        assert_eq!(ms.models.len(), 1);
        // no reads: 0.0 (4.1 would divide by zero); the access materialises the entry
        assert_eq!(mapping_quality(&mut ms, &env.ctx(), 7), 0.0);
        assert_eq!(ms.transcript_read_ids.get(&7), Some(&vec![]));
    }

    #[test]
    fn edgeless_vertices_from_defaultdict_reads() {
        let (a, b, c) = (Vertex(100, 200), Vertex(300, 400), Vertex(5000, 5100));
        let mut g = IntronGraph::default();
        g.clustered_introns.insert(Iv::new(100, 200), 5);
        g.clustered_introns.insert(Iv::new(300, 400), 7);
        g.clustered_introns.insert(Iv::new(5000, 5100), 9);
        g.outgoing.insert(a, vec![b]);
        g.incoming.insert(b, vec![a]);
        let mut gv = GraphView::new(&g);
        assert_eq!(gv.overlapping_component_max_coverage(Iv::new(4900, 5200)), 0);
        assert_eq!(gv.max_component_coverage(&[a]), 7);
        assert!(gv.extra.is_empty());
        // c is a key of neither map: the read adds it as an edgeless vertex
        assert_eq!(gv.max_component_coverage(&[c]), 9);
        assert!(gv.is_monointron(c));
        assert_eq!(gv.extra, vec![c]);
        assert_eq!(gv.overlapping_component_max_coverage(Iv::new(4900, 5200)), 9);
        assert_eq!(gv.overlapping_component_max_coverage(Iv::new(150, 160)), 7);
        assert!(!gv.is_monointron(a));
        // terminal neighbours only -> monointron
        let mut g2 = IntronGraph::default();
        g2.outgoing.insert(b, vec![Vertex(vertex_code::POLYA, 500)]);
        g2.incoming.insert(b, vec![Vertex(vertex_code::READ_START, 250)]);
        assert!(GraphView::new(&g2).is_monointron(b));
    }

    #[test]
    fn similar_cache_matches_uncached() {
        let env = Env::new(&[], IntronGraph::default());
        let c = env.ctx();
        let ex = [(100, 200), (300, 400), (500, 600), (700, 800)];
        let mut models = vec![
            model(0, ModelKind::Nnic, &ex),
            model(1, ModelKind::Nnic, &ex),                                          // identical -> matching
            model(2, ModelKind::Nnic, &ex[1..]),                                    // truncated
            model(3, ModelKind::Nnic, &[(100, 200), (300, 400), (700, 800)]),       // exon skip
            model(4, ModelKind::Nnic, &[(10000, 10100), (10200, 10300), (10400, 10500)]), // far away
        ];
        let mut cache = SimilarCache::default();
        for ru in [false, true] {
            assert_eq!(detect_similar_cached(&models, &c, ru, &mut cache), detect_similar_isoforms(&models, &c, ru));
        }
        let ts = detect_similar_isoforms(&models, &c, false);
        assert_eq!(ts.get(&1), Some(&0));
        assert!(!ts.contains_key(&3));
        // the second call reuses every pair
        let before = cache.n_assigned;
        detect_similar_cached(&models, &c, true, &mut cache);
        assert_eq!(cache.n_assigned, before);
        // model 2's ends change: only its pairs are recomputed
        models[2].exons[0].s = 290;
        let got = detect_similar_cached(&models, &c, true, &mut cache);
        assert_eq!(got, detect_similar_isoforms(&models, &c, true));
        assert!(cache.n_assigned > before);
        // confirmed polyA protects a truncation from remove_unconfirmed, never a match
        models[2].polya_confirmed = true;
        models[1].polya_confirmed = true;
        let got = detect_similar_cached(&models, &c, true, &mut cache);
        assert_eq!(got, detect_similar_isoforms(&models, &c, true));
        assert!(got.contains_key(&1));
    }
}
