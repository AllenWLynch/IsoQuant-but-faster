//! Shared helpers for the U13-U19 store-snapshot stages (owned by U14; imported by the
//! other model units). Everything here except [`next_novel_gene`] is fixture-only.
//!
//! 1. [`StageInput::load`]: loads the case's block, index, U11 graph and U12 paths, and
//!    converts the input `Store` with a **private copy** of the id maps (block-local model
//!    ids start at 0, so `ModelIdx` = creation order within the case).
//! 2. [`next_novel_gene`]: allocation of a fresh block-local novel gene key.
//! 3. Id remapping: models and novel genes created in Rust have no 4.1 id. [`StageInput::store_out`]
//!    names them with placeholders (`PLACEHOLDER` prefix), and [`compare_stores`] renames
//!    the placeholders to the 4.1 ids of the models/genes the 4.1 run created, by creation
//!    order, before calling [`store_diff`]. (4.1 ids depend on the id distributor, which
//!    iqrs does not reproduce: D-M1.)

use crate::ctx::DetectedKnown;
use iq_graph::paths::{known_spliced_isoforms, PathProcessor};
use iq_types::fixture::convert::{store_diff, ModelExtras, StoreExtras};
use iq_types::fixture::ids::IdMaps;
use iq_types::fixture::schema::{IndexDump, StageStoreIn, StoreFixture, U12FillOut};
use iq_types::fixture::{read_lines, unit_files, Case, FixError, FixResult, FixtureCtx, LoadedBlock, Verdict};
use iq_types::{serde_json, GeneKey, IntronGraph, IsoformIndex, ModelIdx, ModelKind, ModelStore, Params, PathStorage};
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::sync::{Arc, Mutex, OnceLock};

/// `GeneKey::Novel(1 + max)` over the store's models (0 if none).
///
/// Block-local numbering (renumbered at output). Caveat: if the model holding
/// the largest novel gene number was deleted, its number is reused. That is harmless in
/// production (genes exist only through models), but the fixture remap then treats the
/// gene as pre-existing.
pub fn next_novel_gene(ms: &ModelStore) -> GeneKey {
    let next = ms
        .models
        .iter()
        .filter_map(|m| match m.gene {
            GeneKey::Novel(g) => Some(g + 1),
            GeneKey::Ref(_) => None,
        })
        .max()
        .unwrap_or(0);
    GeneKey::Novel(next)
}

/// Prefix of placeholder ids given to Rust-created models (`<P>model:<key>`) and novel genes
/// (`<P>gene:<n>`) in [`StageInput::store_out`]. [`compare_stores`] renames them.
pub const PLACEHOLDER: &str = "\u{1}iqrs:";

/// Loads the U11 graph stored under case key `key` (cached per capture dir; see
/// `iq_graph::paths::fixture::graph`).
pub fn graph(ctx: &FixtureCtx, key: &str) -> FixResult<Arc<IntronGraph>> {
    iq_graph::paths::fixture::graph(ctx, key)
}

/// Loads the U12 `IntronPathStorage.fill` output stored under case key `key` (cached per
/// capture dir). `known_isoforms_in_graph` is left empty; see [`StageInput::load`].
pub fn fill_paths(ctx: &FixtureCtx, key: &str) -> FixResult<Arc<PathStorage>> {
    type Cache = Mutex<HashMap<PathBuf, Arc<HashMap<String, Arc<PathStorage>>>>>;
    static CACHE: OnceLock<Cache> = OnceLock::new();
    let cache = CACHE.get_or_init(Default::default);
    let all = {
        let mut g = cache.lock().unwrap();
        match g.get(&ctx.root) {
            Some(m) => m.clone(),
            None => {
                let mut m = HashMap::new();
                for f in unit_files(&ctx.root, "U12")? {
                    for case in read_lines::<Case>(&f)? {
                        let case = case?;
                        if case.input.get("fn").and_then(|v| v.as_str()) != Some("IntronPathStorage.fill") {
                            continue;
                        }
                        if let Ok(out) = serde_json::from_value::<U12FillOut>(case.out) {
                            m.insert(case.case, Arc::new(out.to_core()));
                        }
                    }
                }
                let m = Arc::new(m);
                g.insert(ctx.root.clone(), m.clone());
                m
            }
        }
    };
    all.get(key).cloned().ok_or_else(|| FixError(format!("U12 fill case {key} not found")))
}

/// Everything a store-snapshot stage runner needs, converted.
pub struct StageInput {
    /// Private copy of the shared id maps, with block-local model ids.
    pub ids: IdMaps,
    pub params: Params,
    pub block: Arc<LoadedBlock>,
    pub index: Arc<IsoformIndex>,
    pub index_dump: Arc<IndexDump>,
    /// Empty when the case has no graph.
    pub graph: Arc<IntronGraph>,
    /// U12 fill output plus `known_isoforms_in_graph` recomputed live from the graph (U12
    /// `get_known_spliced_isoforms` is exact). Empty when the case has no paths.
    pub paths: PathStorage,
    /// The input store.
    pub store: ModelStore,
    pub extras: StoreExtras,
    /// Input `detected_known_isoforms` (sorted ids → insertion order of that list).
    pub detected_known: DetectedKnown,
    /// Input models' extras, by key.
    model_extras: HashMap<ModelIdx, ModelExtras>,
    /// Novel gene keys present in the input store.
    input_genes: HashSet<GeneKey>,
}

impl StageInput {
    pub fn load(input: &StageStoreIn, ctx: &FixtureCtx) -> FixResult<Self> {
        let params = ctx.params()?;
        let block = ctx.block(&input.block)?;
        let index = ctx.index(&input.index)?;
        let index_dump = ctx.index_dump(&input.index)?;
        let graph = match &input.graph {
            Some(k) => graph(ctx, k)?,
            None => Arc::new(IntronGraph::default()),
        };
        let mut paths = match &input.paths {
            Some(k) => (*fill_paths(ctx, k)?).clone(),
            None => PathStorage::default(),
        };
        if input.graph.is_some() {
            let proc = PathProcessor::new(&graph, &params);
            paths.known_isoforms_in_graph = known_spliced_isoforms(&index, &proc);
        }
        let mut ids = ctx.ids().clone();
        ids.reset_block();
        let names = block.read_names();
        let (store, extras) = input.store.to_core(&mut ids, &names)?;
        let detected_known = extras.detected_known_isoforms.iter().map(|s| ids.txs.intern(s)).collect();
        let model_extras = store.models.iter().zip(&extras.models).map(|(m, e)| (m.key, e.clone())).collect();
        let input_genes = store.models.iter().map(|m| m.gene).collect();
        Ok(StageInput {
            ids,
            params,
            block,
            index,
            index_dump,
            graph,
            paths,
            store,
            extras,
            detected_known,
            model_extras,
            input_genes,
        })
    }

    /// Fixture-only fields of a model created in Rust: a known model takes `source` and
    /// `other_features` from the index (`from_reference_transcript`); a novel one gets the
    /// `TranscriptModel` defaults (`source = "IsoQuant"`).
    pub fn new_model_extras(&self, kind: ModelKind) -> ModelExtras {
        let chr = self.block.header.chr.clone();
        match kind {
            ModelKind::Known(t) => {
                let id = self.ids.txs.resolve(t);
                let d = &self.index_dump;
                let source = d.sources.iter().find(|(k, _)| k == id).map(|(_, v)| v.clone()).unwrap_or_default();
                let other_features =
                    d.other_features.iter().find(|(k, _)| k == id).map(|(_, v)| v.clone()).unwrap_or_default();
                ModelExtras { chr, source, other_features, additional_info: Vec::new() }
            }
            _ => ModelExtras { chr, source: "IsoQuant".to_owned(), ..Default::default() },
        }
    }

    /// Output snapshot of `ms`. Input models keep their ids and extras; models created in
    /// Rust and their new novel genes get placeholder ids (see [`compare_stores`]).
    /// `detected` is written as sorted ids. `id_counter` is copied from the input (ignored by
    /// [`store_diff`]).
    pub fn store_out(&self, ms: &ModelStore, detected: &DetectedKnown) -> StoreFixture {
        let mut ids = self.ids.clone();
        for k in ids.models.len() as ModelIdx..ms.next_key {
            let got = ids.models.intern(&format!("{PLACEHOLDER}model:{k}"));
            debug_assert_eq!(got, k);
        }
        let mut ms2;
        let ms = if ms.models.iter().any(|m| !self.input_genes.contains(&m.gene) && matches!(m.gene, GeneKey::Novel(_)))
        {
            ms2 = ms.clone();
            for m in &mut ms2.models {
                if let GeneKey::Novel(g) = m.gene {
                    if !self.input_genes.contains(&m.gene) {
                        m.gene = GeneKey::Novel(ids.novel_genes.intern(&format!("{PLACEHOLDER}gene:{g}")));
                    }
                }
            }
            &ms2
        } else {
            ms
        };
        let mut dk: Vec<String> = detected.iter().map(|&t| ids.txs.resolve(t).to_owned()).collect();
        dk.sort();
        let ex = StoreExtras {
            models: ms
                .models
                .iter()
                .map(|m| self.model_extras.get(&m.key).cloned().unwrap_or_else(|| self.new_model_extras(m.kind)))
                .collect(),
            detected_known_isoforms: dk,
            id_counter: self.extras.id_counter,
        };
        StoreFixture::from_core(ms, &ids, &ex, &self.block.read_names())
    }
}

/// Renames placeholder ids in `got` to the ids 4.1 gave the corresponding new objects:
/// placeholder models (in `got` storage order) ↔ `want` model ids absent from `got` (in
/// `want` storage order); likewise for novel genes (first appearance in model order).
/// Unmatched placeholders stay as they are (and show up in the diff).
pub fn align_placeholders(got: &mut StoreFixture, want: &StoreFixture) {
    fn pairing<'a>(got: impl Iterator<Item = &'a str>, want: impl Iterator<Item = &'a str>) -> HashMap<String, String> {
        let mut ph = Vec::new();
        let mut known = HashSet::new();
        for g in got {
            if g.starts_with(PLACEHOLDER) {
                if !ph.contains(&g) {
                    ph.push(g);
                }
            } else {
                known.insert(g);
            }
        }
        let mut fresh: Vec<&str> = Vec::new();
        for w in want {
            if !known.contains(w) && !fresh.contains(&w) {
                fresh.push(w);
            }
        }
        ph.into_iter().zip(fresh).map(|(a, b)| (a.to_owned(), b.to_owned())).collect()
    }
    let mm = pairing(got.models.iter().map(|m| m.id.as_str()), want.models.iter().map(|m| m.id.as_str()));
    let gm = pairing(got.models.iter().map(|m| m.gene.as_str()), want.models.iter().map(|m| m.gene.as_str()));
    if mm.is_empty() && gm.is_empty() {
        return;
    }
    let ren = |s: &mut String, map: &HashMap<String, String>| {
        if let Some(n) = map.get(s.as_str()) {
            *s = n.clone();
        }
    };
    for m in &mut got.models {
        ren(&mut m.id, &mm);
        ren(&mut m.gene, &gm);
    }
    for (k, _) in &mut got.transcript_read_ids {
        ren(k, &mm);
    }
    for (k, _) in &mut got.internal_counter {
        ren(k, &mm);
    }
    for (_, v) in &mut got.construction_assignment {
        ren(v, &mm);
    }
    if let Some(v) = &mut got.construction_assignment_idx {
        for (_, m) in v {
            ren(m, &mm);
        }
    }
}

/// Store-snapshot verdict: [`align_placeholders`], then [`store_diff`].
pub fn compare_stores(got: &StoreFixture, want: &StoreFixture, ctx: &FixtureCtx) -> Verdict {
    let mut got = got.clone();
    align_placeholders(&mut got, want);
    match store_diff(&got, want, &ctx.diff_opts) {
        None => Verdict::Exact,
        Some(d) if d.is_empty() => Verdict::Exact,
        Some(d) => Verdict::Mismatch(d),
    }
}
