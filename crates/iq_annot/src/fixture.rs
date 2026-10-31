//! U1 fixture runners: `dump_annotation`, `get_gene_info_for_region`, `from_models`.
//!
//! Indices are compared in the fixture id space: the builder output (keys = annotation
//! `TxIdx`/`GeneIdx`, or `ModelIdx`) is remapped through the context's `IdMaps` by id string
//! and must then equal `IsoformIndex::from_dump` of the captured dump field for field.
//!
//! The GFF3 for (b) is `$IQRS_GFF3`, or the first `data/ref/annotation.gff3` found in an
//! ancestor of the capture directory. Parsed annotations are cached per process.

use crate::builder::{index_from_genes, index_from_models_named, index_from_region};
use crate::index::{AnnotationIndex, GeneTypes};
use crate::parse::parse_gff3;
use iq_types::fixture::schema::{
    GeneFixture, TxFixture, U1DumpAnnotationIn, U1DumpAnnotationOut, U1FromModelsIn, U1GeneInfoForRegionIn,
};
use iq_types::fixture::{Diff, DiffEntry, FixError, FixResult, FixtureCtx, IdMaps, TxSpace, UnitFixture, Verdict};
use iq_types::serde::{Deserialize, Serialize};
use iq_types::{GeneKey, IsoRef, IsoformIndex};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, OnceLock};

type Cache = Mutex<HashMap<(PathBuf, GeneTypes), Arc<AnnotationIndex>>>;

/// Parses (once per process) the annotation at `path`.
pub fn cached_annotation(path: &Path, gt: GeneTypes) -> FixResult<Arc<AnnotationIndex>> {
    static CACHE: OnceLock<Cache> = OnceLock::new();
    let cache = CACHE.get_or_init(Default::default);
    let mut g = cache.lock().unwrap_or_else(|e| e.into_inner());
    let key = (path.to_owned(), gt);
    if let Some(a) = g.get(&key) {
        return Ok(a.clone());
    }
    let a = Arc::new(parse_gff3(path, gt).map_err(|e| FixError(e.0))?);
    g.insert(key, a.clone());
    Ok(a)
}

/// GFF3 used for region cases: `$IQRS_GFF3` or `<ancestor>/data/ref/annotation.gff3`.
pub fn gff3_path(ctx: &FixtureCtx) -> FixResult<PathBuf> {
    if let Some(p) = std::env::var_os("IQRS_GFF3") {
        return Ok(PathBuf::from(p));
    }
    let root = std::fs::canonicalize(&ctx.root).unwrap_or_else(|_| ctx.root.clone());
    for a in root.ancestors() {
        let p = a.join("data/ref/annotation.gff3");
        if p.is_file() {
            return Ok(p);
        }
    }
    Err(FixError(format!("no data/ref/annotation.gff3 above {} (set IQRS_GFF3)", ctx.root.display())))
}

/// Remaps an index built against `ann` (reference keys) into the fixture id space.
pub fn remap_ref_index(mut idx: IsoformIndex, ann: &AnnotationIndex, chr_name: &str, ids: &mut IdMaps) -> IsoformIndex {
    idx.chr = ids.chroms.intern(chr_name) as u16;
    let gene = |k: GeneKey, ids: &mut IdMaps| match k {
        GeneKey::Ref(g) => ids.gene_key(ann.gene_id(g)),
        n => n,
    };
    for r in &mut idx.iso {
        if let IsoRef::Ref(t) = r.key {
            r.key = ids.iso_ref(ann.tx_id(t), TxSpace::Ref);
        }
        r.gene = gene(r.gene, ids);
    }
    for g in &mut idx.genes {
        g.key = gene(g.key, ids);
    }
    idx
}

/// Field-by-field diff of two indices (empty if equal).
pub fn index_diff(got: &IsoformIndex, want: &IsoformIndex) -> Vec<DiffEntry> {
    let mut out = Vec::new();
    let mut f = |name: &str, a: String, b: String| {
        if a != b {
            let cut = |s: String| if s.len() > 600 { format!("{}…", &s[..600]) } else { s };
            out.push(DiffEntry { path: format!("index.{name}"), got: cut(a), want: cut(b) });
        }
    };
    macro_rules! fld {
        ($($n:ident),*) => { $( f(stringify!($n), format!("{:?}", got.$n), format!("{:?}", want.$n)); )* };
    }
    fld!(chr, start, end, delta, iso, id_rank, iso_exons, iso_introns, introns, exons, split_exons);
    fld!(intron_prof, exon_prof, split_prof, prof_range, span_feat, genes, span_tree);
    out
}

fn judge_index(got: Option<&IsoformIndex>, sha: &str, ctx: &FixtureCtx, mut extra: Vec<DiffEntry>) -> Verdict {
    let Some(got) = got else { return Verdict::Mismatch(Diff::error("runner produced no index")) };
    let want = match ctx.index_dump(sha).map_err(|e| e.0).and_then(|d| IsoformIndex::from_dump(&d, &mut ctx.ids())) {
        Ok(w) => w,
        Err(e) => return Verdict::Mismatch(Diff::error(format!("index {sha}: {e}"))),
    };
    extra.extend(index_diff(got, &want));
    if extra.is_empty() {
        Verdict::Exact
    } else {
        Verdict::Mismatch(Diff { entries: extra, ..Default::default() })
    }
}

// ---- (a) dump_annotation -----------------------------------------------------------------

/// The per-chromosome annotation dump (`iqcap.dump_annotation`) from an [`AnnotationIndex`].
pub fn dump_chromosome(ann: &AnnotationIndex, chr: &str) -> U1DumpAnnotationOut {
    let mut genes = Vec::new();
    if let Some(c) = ann.chrom(chr) {
        for g in ann.chrom_genes(c) {
            let (_, iv, strand) = ann.gene_span(g);
            let transcripts = ann
                .gene_txs(g)
                .map(|t| {
                    let tiv = ann.tx_row_span(t);
                    TxFixture {
                        id: ann.tx_id(t).to_owned(),
                        ty: ann.tx_type(t).to_owned(),
                        start: tiv.s,
                        end: tiv.e,
                        strand: ann.tx_strand(t),
                        source: ann.tx_source(t).to_owned(),
                        attributes: ann.tx_attributes(t),
                        attrs: ann.tx_attrs(t),
                        exons: ann.tx_exons(t).to_vec(),
                        other_features: ann.tx_other_features(t).map(|(iv, ty)| (iv.s, iv.e, ty.to_owned())).collect(),
                        children_order_index: ann.tx_file_rank(t),
                    }
                })
                .collect();
            genes.push(GeneFixture {
                id: ann.gene_id(g).to_owned(),
                start: iv.s,
                end: iv.e,
                strand,
                source: ann.gene_source(g).to_owned(),
                attributes: ann.gene_attributes(g),
                attrs: ann.gene_attrs(g),
                transcripts,
            });
        }
    }
    U1DumpAnnotationOut { chr: chr.to_owned(), genes }
}

pub struct U1Dump;
impl UnitFixture for U1Dump {
    const UNIT: &'static str = "U1";
    const FN: &'static str = "dump_annotation";
    type In = U1DumpAnnotationIn;
    type Out = U1DumpAnnotationOut;
    fn run(input: &Self::In, _ctx: &FixtureCtx) -> FixResult<Self::Out> {
        let gt = match input.gene_types.as_str() {
            "isoquant" => GeneTypes::Isoquant,
            "all" => GeneTypes::All,
            o => return Err(FixError(format!("unknown gene_types {o:?}"))),
        };
        let ann = cached_annotation(Path::new(&input.gff3), gt)?;
        Ok(dump_chromosome(&ann, &input.chr))
    }
}

// ---- (b) get_gene_info_for_region --------------------------------------------------------

/// `out` of (b) plus the built index (not serialized; compared against the dump).
#[derive(Debug, Serialize, Deserialize)]
#[serde(crate = "iq_types::serde")]
pub struct RegionOut {
    pub gene_ids: Vec<String>,
    pub index: String,
    pub empty: bool,
    #[serde(skip)]
    pub built: Option<Box<IsoformIndex>>,
}
impl PartialEq for RegionOut {
    fn eq(&self, o: &Self) -> bool {
        self.gene_ids == o.gene_ids && self.empty == o.empty
    }
}

pub struct U1Region;
impl UnitFixture for U1Region {
    const UNIT: &'static str = "U1";
    const FN: &'static str = "get_gene_info_for_region";
    type In = U1GeneInfoForRegionIn;
    type Out = RegionOut;
    fn run(input: &Self::In, ctx: &FixtureCtx) -> FixResult<Self::Out> {
        let ann = cached_annotation(&gff3_path(ctx)?, GeneTypes::Isoquant)?;
        let genes = match ann.chrom(&input.chr) {
            Some(c) => ann.genes_for_region0(c, input.region),
            None => Default::default(),
        };
        let idx = if genes.is_empty() {
            index_from_region(0, input.region.s as i32, input.region.e as i32, input.delta)
        } else {
            index_from_genes(&ann, &genes, input.delta)
        };
        let idx = remap_ref_index(idx, &ann, &input.chr, &mut ctx.ids());
        Ok(RegionOut {
            gene_ids: genes.iter().map(|&g| ann.gene_id(g).to_owned()).collect(),
            index: String::new(),
            empty: idx.exons.is_empty(),
            built: Some(Box::new(idx)),
        })
    }
    fn compare(got: &Self::Out, want: &Self::Out, ctx: &FixtureCtx) -> Verdict {
        let mut extra = Vec::new();
        if got.gene_ids != want.gene_ids {
            extra.push(DiffEntry {
                path: "gene_ids".into(),
                got: format!("{:?}", got.gene_ids),
                want: format!("{:?}", want.gene_ids),
            });
        }
        if got.empty != want.empty {
            extra.push(DiffEntry { path: "empty".into(), got: got.empty.to_string(), want: want.empty.to_string() });
        }
        judge_index(got.built.as_deref(), &want.index, ctx, extra)
    }
}

// ---- (c) from_models ---------------------------------------------------------------------

#[derive(Debug, Serialize, Deserialize)]
#[serde(crate = "iq_types::serde")]
pub struct ModelsOut {
    pub index: String,
    #[serde(skip)]
    pub built: Option<Box<IsoformIndex>>,
}
impl PartialEq for ModelsOut {
    fn eq(&self, _: &Self) -> bool {
        false // always judged by `compare`
    }
}

pub struct U1FromModels;
impl UnitFixture for U1FromModels {
    const UNIT: &'static str = "U1";
    const FN: &'static str = "from_models";
    type In = U1FromModelsIn;
    type Out = ModelsOut;
    fn run(input: &Self::In, ctx: &FixtureCtx) -> FixResult<Self::Out> {
        let mut ids = ctx.ids();
        let models: Vec<_> = input.models.iter().map(|m| m.to_core(&mut ids).0).collect();
        let chr = input.models.first().map_or(0, |m| ids.chroms.intern(&m.chr) as u16);
        let names: Vec<&str> = input.models.iter().map(|m| m.id.as_str()).collect();
        let idx = index_from_models_named(chr, &models, input.delta, &names);
        Ok(ModelsOut { index: String::new(), built: Some(Box::new(idx)) })
    }
    fn compare(got: &Self::Out, want: &Self::Out, ctx: &FixtureCtx) -> Verdict {
        judge_index(got.built.as_deref(), &want.index, ctx, Vec::new())
    }
}

pub(crate) fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    reg.add::<U1Dump>().add::<U1Region>().add::<U1FromModels>();
}
