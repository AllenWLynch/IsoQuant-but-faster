//! U13: FL-path models (`model_construction/fl_graph_constructor.py:42-206`).
//!
//! Every full-length intron-graph path becomes a known reference model, an alternative-end
//! NIC, or a novel (NIC/NNIC) model; the path's reads are bound to the new model
//! (`ModelStore.save_assigned_read`, `context.py:76-82`).
//!
//! Structure (design §7.4):
//! 1. FL paths are ordered `(len desc, Python tuple order)` (`:94-95`).
//! 2. For every path that could produce a model, the path is profiled and assigned (U5/U6),
//!    and `_reference_isoform_for_path` (`:42-90`) is evaluated. This part is pure, so it runs
//!    in parallel ([`Decision`] per path).
//! 3. The 4.1 decision loop runs serially in path order, because it mutates
//!    `detected_known` and the store.
//!
//! Differences from 4.1 that cannot change outputs:
//! - **Cheap checks first.** 4.1 profiles and assigns every FL path before testing its read
//!   count (`:113` vs `:125,158`). A path with `count < min(min_known_count, min_novel_count)`
//!   can produce no model on either branch (known: `count < min_known_count`; novel:
//!   `count < min_novel_count`; the `create_alt_ts_nics = false` `continue` produces nothing
//!   either), and the skipped work has no other effect (strand-cache fills, the provisional
//!   transcript id), so such paths are dropped before step 2.
//! - **D-M3** (dead branch): the FL-path `PolyAInfo` tests `intron_path[0][0] == polyt`
//!   *after* the terminals were stripped (`:107-112`), so it is always empty: `NO_POLYA`.
//! - **D-M1**: no transcript/gene id is drawn for rejected paths. Models get `ModelIdx`
//!   keys in creation order, novel genes `GeneKey::Novel(k)` in creation order; final ids
//!   are assigned at output (design §5.6).
//! - Technical replicas (`file_name` grouping, `:166-176`) are off in scope; not ported.

use crate::ctx::{Ctx, DetectedKnown, ReadAssigner};
use iq_assign::assigner::{assign_with, AssignScratch};
use iq_assign::profile::{ProfileBuilder, Scratch};
use iq_assign::AssignParams;
use iq_types::{
    vertex_code, Assignment, AssignmentType, CombinedProfile, EventType, GeneKey, IsoRef, IsoformIndex, Iv,
    ModelKind, ModelStore, Pos, RefSlice, Strand, StrandnessReporting, TranscriptModel, TxIdx, Vertex, NO_POLYA,
};
use rayon::prelude::*;
use std::cell::RefCell;

/// Inputs of U13 that [`Ctx`] does not carry.
#[derive(Clone, Copy, Debug)]
pub struct FlAux<'a> {
    /// Chromosome sequence covering every intron of the block's FL paths (4.1 `chr_record`,
    /// read by `get_intron_strand` for introns without a single annotated strand).
    pub reference: RefSlice<'a>,
    /// Rank of each `ctx.index.genes[k]` id in Python string order (the
    /// `select_reference_gene` tie break, `context.py:196`). See [`gene_id_ranks`].
    pub gene_rank: &'a [u32],
}

/// Ranks for [`FlAux::gene_rank`] from the gene id strings of `index.genes` (same order).
pub fn gene_id_ranks(gene_ids: &[&str]) -> Vec<u32> {
    iq_types::index::id_ranks(gene_ids)
}

/// Live U6 assigner over `iq_assign` (one scratch per thread). `quick_mode` per call
/// overrides `params.quick_mode`.
#[derive(Clone, Debug)]
pub struct LiveAssigner {
    pub params: AssignParams,
}

impl ReadAssigner for LiveAssigner {
    fn assign(&self, index: &IsoformIndex, profile: &CombinedProfile, quick_mode: bool, _read_key: u64) -> Option<Assignment> {
        thread_local! {
            static SCRATCH: RefCell<AssignScratch> = RefCell::new(AssignScratch::default());
        }
        SCRATCH.with(|s| {
            let s = &mut *s.borrow_mut();
            if quick_mode == self.params.quick_mode {
                assign_with(index, &self.params, profile, profile.polya, s)
            } else {
                let p = AssignParams { quick_mode, ..self.params.clone() };
                assign_with(index, &p, profile, profile.polya, s)
            }
        })
    }
}

// ------------------------------------------------------------------------------------------
// Strand detection (gene_info.py:924-967, context.py:171-184)
// ------------------------------------------------------------------------------------------

/// `get_intron_strand(intron, chr_record)` (`common.py:854-864`): donor `ref[s..=s+1]`,
/// acceptor `ref[e-1..=e]`, upper-cased; missing bases make a site non-canonical.
pub fn intron_strand_ref(intron: Iv, reference: &RefSlice<'_>) -> Strand {
    let b = |p: Pos| reference.base(p).map(|c| c.to_ascii_uppercase());
    let (Some(l0), Some(l1), Some(r0), Some(r1)) = (b(intron.s), b(intron.s + 1), b(intron.e - 1), b(intron.e)) else {
        return Strand::Unknown;
    };
    let (l, r) = ([l0, l1], [r0, r1]);
    let fwd = matches!((&l, &r), (b"GT", b"AG") | (b"GC", b"AG") | (b"AT", b"AC"));
    let rev = matches!((&l, &r), (b"CT", b"AC") | (b"CT", b"GC") | (b"GT", b"AT"));
    match (fwd, rev) {
        (true, false) => Strand::Plus,
        (false, true) => Strand::Minus,
        _ => Strand::Unknown,
    }
}

/// The construction-context `StrandDetector`: annotated introns with a single isoform
/// strand use it (`set_gene_properties`); every other intron is read from the reference.
#[derive(Clone, Debug)]
pub struct StrandDetector<'a> {
    /// Sorted, unique annotated introns with their strand (`None` = several strands →
    /// reference).
    annotated: Vec<(Iv, Option<Strand>)>,
    reference: RefSlice<'a>,
}

impl<'a> StrandDetector<'a> {
    pub fn new(index: &IsoformIndex, reference: RefSlice<'a>) -> Self {
        let mut all: Vec<(Iv, Strand)> = Vec::with_capacity(index.iso_introns.len());
        for (i, r) in index.iso.iter().enumerate() {
            all.extend(index.iso_introns(i).iter().map(|&iv| (iv, r.strand)));
        }
        all.sort_unstable_by_key(|&(iv, s)| (iv, s as u8));
        let mut annotated: Vec<(Iv, Option<Strand>)> = Vec::with_capacity(all.len());
        for (iv, s) in all {
            match annotated.last_mut() {
                Some((last, st)) if *last == iv => {
                    if *st != Some(s) {
                        *st = None;
                    }
                }
                _ => annotated.push((iv, Some(s))),
            }
        }
        StrandDetector { annotated, reference }
    }

    /// Strand of one intron (`strand_dict` lookup or `get_intron_strand`).
    pub fn intron_strand(&self, iv: Iv) -> Strand {
        match self.annotated.binary_search_by_key(&iv, |x| x.0) {
            Ok(k) => match self.annotated[k].1 {
                Some(s) => s,
                None => intron_strand_ref(iv, &self.reference),
            },
            Err(_) => intron_strand_ref(iv, &self.reference),
        }
    }

    /// `count_canonical_sites` → `(fwd, rev)`.
    pub fn count_canonical_sites(&self, introns: impl IntoIterator<Item = Iv>) -> (u32, u32) {
        let (mut f, mut r) = (0, 0);
        for iv in introns {
            match self.intron_strand(iv) {
                Strand::Plus => f += 1,
                Strand::Minus => r += 1,
                Strand::Unknown => {}
            }
        }
        (f, r)
    }
}

/// `get_clean_strand` from `(fwd, rev)` counts.
pub fn clean_strand(fwd: u32, rev: u32) -> Strand {
    match (fwd, rev) {
        (0, r) if r > 0 => Strand::Minus,
        (f, 0) if f > 0 => Strand::Plus,
        _ => Strand::Unknown,
    }
}

/// `get_strand(introns, has_polya, has_polyt)` from `(fwd, rev)` counts.
pub fn majority_strand(fwd: u32, rev: u32, has_polya: bool, has_polyt: bool) -> Strand {
    if fwd == rev {
        return match (has_polya, has_polyt) {
            (true, false) => Strand::Plus,
            (false, true) => Strand::Minus,
            _ => Strand::Unknown,
        };
    }
    if rev < fwd {
        Strand::Plus
    } else {
        Strand::Minus
    }
}

// ------------------------------------------------------------------------------------------
// Helpers
// ------------------------------------------------------------------------------------------

/// `is_matching_assignment` (`isoform_assignment.py:1100-1118`).
pub fn is_matching_assignment(a: &Assignment) -> bool {
    match a.atype {
        AssignmentType::Unique => true,
        AssignmentType::UniqueMinorDifference => a.matches.first().is_some_and(|m| {
            m.events.iter().all(|e| {
                matches!(
                    e.ty,
                    EventType::None
                        | EventType::Fsm
                        | EventType::ExonMisalignment
                        | EventType::IntronShift
                        | EventType::TerminalSiteMatchLeft
                        | EventType::TerminalSiteMatchRight
                        | EventType::TerminalSiteMatchLeftPrecise
                        | EventType::TerminalSiteMatchRightPrecise
                        | EventType::CorrectPolyaSiteRight
                        | EventType::CorrectPolyaSiteLeft
                        | EventType::ExonElongationLeft
                        | EventType::ExonElongationRight
                )
            })
        }),
        _ => false,
    }
}

#[inline]
fn intron_of(v: &Vertex) -> Iv {
    Iv::new(v.0 as Pos, v.1 as Pos)
}

/// `get_exons(transcript_range, introns)` (`common.py:506-507`) =
/// `junctions_from_blocks([(-inf, s-1)] + introns + [(e+1, inf)])`: gaps that are empty
/// are dropped, as in 4.1.
pub fn get_exons_into(range: (Pos, Pos), introns: &[Vertex], out: &mut Vec<Iv>) {
    out.clear();
    let mut prev_end = range.0 as i64 - 1;
    for v in introns {
        if prev_end + 1 < v.0 {
            out.push(Iv::new((prev_end + 1) as Pos, (v.0 - 1) as Pos));
        }
        prev_end = v.1;
    }
    let last_start = range.1 as i64 + 1;
    if prev_end + 1 < last_start {
        out.push(Iv::new((prev_end + 1) as Pos, (last_start - 1) as Pos));
    }
}

/// Outcome of the pure part of the per-path loop.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Decision {
    /// `reference_isoform` set: a known reference model candidate.
    Ref(TxIdx),
    /// Novel branch.
    Novel,
    /// `continue` (`create_alt_ts_nics = false` and the chain is a known isoform).
    Skip,
}

/// Block-level tables of one `construct_fl` call.
struct Tables<'a> {
    /// Sorted unique `(intron, gene)` pairs (`ctx.intron_genes`).
    intron_genes: Vec<(Iv, GeneKey)>,
    strands: StrandDetector<'a>,
}

/// `_reference_isoform_for_path` + the branch selection of `:115-125`.
fn decide(ctx: &Ctx<'_>, path: &[Vertex], assignment: Option<&Assignment>) -> Decision {
    let introns = &path[1..path.len() - 1];
    let index = ctx.index;
    let known = &ctx.paths.known_isoforms_in_graph;
    let matching = assignment.filter(|a| is_matching_assignment(a));
    let assigned_tx = || matching.and_then(|a| a.matches.first()).and_then(|m| m.tx);
    let as_ref = |tx: Option<IsoRef>| match tx {
        Some(IsoRef::Ref(t)) => index.find_iso(IsoRef::Ref(t)).map(|p| (t, p)),
        _ => None,
    };
    if !ctx.params.create_alt_ts_nics {
        return if matching.is_some() {
            as_ref(assigned_tx()).map_or(Decision::Novel, |(t, _)| Decision::Ref(t))
        } else if known.contains_key(introns) {
            Decision::Skip
        } else {
            Decision::Novel
        };
    }
    let tx = if matching.is_some() {
        assigned_tx()
    } else if let Some(&t) = known.get(introns) {
        Some(IsoRef::Ref(t))
    } else {
        return Decision::Novel;
    };
    let Some((t, p)) = as_ref(tx) else { return Decision::Novel };
    let ex = index.iso_exons(p);
    let apa = ctx.params.matching.apa_delta as i64;
    let (first, last) = (path[0], path[path.len() - 1]);
    let left_diff = (first.1 - ex[0].s as i64).abs() > apa;
    let right_diff = (last.1 - ex[ex.len() - 1].e as i64).abs() > apa;
    if (first.0 == vertex_code::POLYT && left_diff) || (last.0 == vertex_code::POLYA && right_diff) {
        return Decision::Novel;
    }
    if ctx.params.use_tss_model {
        let strand = index.iso[p].strand;
        if (strand == Strand::Plus && first.0 == vertex_code::TSS_LEFT && left_diff)
            || (strand == Strand::Minus && last.0 == vertex_code::TSS_RIGHT && right_diff)
        {
            return Decision::Novel;
        }
    }
    Decision::Ref(t)
}

/// `select_reference_gene` (`context.py:186-200`).
fn select_reference_gene(ctx: &Ctx<'_>, aux: &FlAux<'_>, tables: &Tables<'_>, introns: &[Vertex], strand: Strand) -> Option<GeneKey> {
    let index = ctx.index;
    if index.exons.is_empty() {
        return None;
    }
    // gene -> count, then sort by (count, gene id) descending
    let mut counts: Vec<(GeneKey, u32)> = Vec::new();
    for v in introns {
        let iv = intron_of(v);
        let lo = tables.intron_genes.partition_point(|x| x.0 < iv);
        for &(_, g) in tables.intron_genes[lo..].iter().take_while(|x| x.0 == iv) {
            match counts.iter_mut().find(|c| c.0 == g) {
                Some(c) => c.1 += 1,
                None => counts.push((g, 1)),
            }
        }
    }
    let gene_pos = |g: GeneKey| index.genes.iter().position(|r| r.key == g);
    let mut ranked: Vec<(u32, u32, Strand, GeneKey)> = counts
        .iter()
        .map(|&(g, c)| {
            let p = gene_pos(g);
            let rank = p.and_then(|p| aux.gene_rank.get(p).copied()).unwrap_or(0);
            let st = p.map_or(Strand::Unknown, |p| index.genes[p].strand);
            (c, rank, st, g)
        })
        .collect();
    ranked.sort_unstable_by_key(|r| std::cmp::Reverse((r.0, r.1)));
    ranked.into_iter().find(|r| strand == Strand::Unknown || r.2 == strand).map(|r| r.3)
}

/// Novel branch (`:142-203`): the new model, or `None` when a gate rejects the path.
#[allow(clippy::too_many_arguments)]
fn novel_model(
    ms: &mut ModelStore,
    ctx: &Ctx<'_>,
    aux: &FlAux<'_>,
    tables: &Tables<'_>,
    path: &[Vertex],
    count: u32,
    exons: &[Iv],
    next_novel_gene: &mut u32,
) -> Option<TranscriptModel> {
    let p = ctx.params;
    if count < p.model.min_novel_count {
        return None;
    }
    let introns = &path[1..path.len() - 1];
    let has_polyt = path[0].0 == vertex_code::POLYT;
    let has_polya = path[path.len() - 1].0 == vertex_code::POLYA;
    let (fwd, rev) = tables.strands.count_canonical_sites(introns.iter().map(intron_of));
    let mut strand = majority_strand(fwd, rev, has_polya, has_polyt);
    let clean = clean_strand(fwd, rev);
    if exons.len() == 2 && ((p.model.require_monointronic_polya && !(has_polya || has_polyt)) || clean == Strand::Unknown) {
        return None;
    }
    if (p.report_canonical == StrandnessReporting::OnlyCanonical && clean == Strand::Unknown)
        || (p.report_canonical == StrandnessReporting::OnlyStranded && strand == Strand::Unknown)
    {
        return None;
    }
    let gene = match select_reference_gene(ctx, aux, tables, introns, strand) {
        None => {
            let g = GeneKey::Novel(*next_novel_gene);
            *next_novel_gene += 1;
            g
        }
        Some(g) => {
            if strand == Strand::Unknown {
                strand = ctx.index.genes.iter().find(|r| r.key == g).map_or(Strand::Unknown, |r| r.strand);
            }
            g
        }
    };
    let known_introns = &ctx.index.introns;
    let kind = if introns.iter().all(|v| known_introns.binary_search(&intron_of(v)).is_ok()) {
        ModelKind::Nic
    } else {
        ModelKind::Nnic
    };
    Some(TranscriptModel {
        key: ms.new_key(),
        strand,
        gene,
        kind,
        exons: exons.to_vec(),
        intron_path: introns.into(),
        polya_confirmed: false,
        tss_confirmed: false,
        source_known: None,
    })
}

/// `TranscriptModel.from_reference_transcript` (`gene_info.py:53-68`).
fn known_model(ms: &mut ModelStore, index: &IsoformIndex, t: TxIdx) -> Option<TranscriptModel> {
    let p = index.find_iso(IsoRef::Ref(t))?;
    let r = &index.iso[p];
    Some(TranscriptModel {
        key: ms.new_key(),
        strand: r.strand,
        gene: r.gene,
        kind: ModelKind::Known(t),
        exons: index.iso_exons(p).to_vec(),
        intron_path: Box::new([]),
        polya_confirmed: false,
        tss_confirmed: false,
        source_known: None,
    })
}

/// Paths below this count are assigned serially.
const PAR_MIN_PATHS: usize = 32;

/// `construct_fl_isoforms` (`fl_graph_constructor.py:92-206`). Threads the unit's
/// `detected_known` set (D-M2). Assignments go through `ctx.assigner` (live:
/// [`LiveAssigner`] with the run's matching/splice params, `quick_mode = false`).
pub fn construct_fl(ms: &mut ModelStore, detected_known: &mut DetectedKnown, ctx: &Ctx<'_>, aux: &FlAux<'_>) {
    let p = ctx.params;
    let paths = ctx.paths;
    let (mkc, mnc) = (p.model.min_known_count, p.model.min_novel_count);
    let min_any = mkc.min(mnc);
    // (path, count) of paths that can produce a model; see the module doc ("cheap checks
    // first"). A FL path missing from `paths` (4.1 KeyError; never happens) counts as 0.
    let mut order: Vec<(&[Vertex], u32)> = paths
        .fl_paths
        .iter()
        .filter(|pa| pa.len() > 2)
        .map(|pa| (&pa[..], paths.paths.get(&pa[..]).copied().unwrap_or(0)))
        .filter(|&(_, c)| c >= min_any)
        .collect();
    order.sort_by(|a, b| b.0.len().cmp(&a.0.len()).then_with(|| a.0.cmp(b.0)));
    if order.is_empty() {
        return;
    }

    // pure part: profile + assign + reference isoform, in parallel
    let builder = ProfileBuilder::new(ctx.index, &p.matching);
    let eval = |k: usize, path: &[Vertex], scratch: &mut Scratch, exons: &mut Vec<Iv>| -> Decision {
        let range = (path[0].1 as Pos, path[path.len() - 1].1 as Pos);
        get_exons_into(range, &path[1..path.len() - 1], exons);
        let prof = builder.build(exons, NO_POLYA, scratch);
        let a = ctx.assigner.assign(ctx.index, &prof, false, k as u64);
        decide(ctx, path, a.as_ref())
    };
    let decisions: Vec<Decision> = if order.len() >= PAR_MIN_PATHS {
        order
            .par_iter()
            .enumerate()
            .map_init(|| (Scratch::default(), Vec::new()), |(s, ex), (k, &(pa, _))| eval(k, pa, s, ex))
            .collect()
    } else {
        let (mut s, mut ex) = (Scratch::default(), Vec::new());
        order.iter().enumerate().map(|(k, &(pa, _))| eval(k, pa, &mut s, &mut ex)).collect()
    };

    // serial decision loop
    let mut intron_genes: Vec<(Iv, GeneKey)> = Vec::new();
    for (i, r) in ctx.index.iso.iter().enumerate() {
        intron_genes.extend(ctx.index.iso_introns(i).iter().map(|&iv| (iv, r.gene)));
    }
    intron_genes.sort_unstable();
    intron_genes.dedup();
    let tables = Tables { intron_genes, strands: StrandDetector::new(ctx.index, aux.reference) };
    let mut next_novel_gene = match crate::fixture_util::next_novel_gene(ms) {
        GeneKey::Novel(k) => k,
        GeneKey::Ref(_) => unreachable!("next_novel_gene returns a novel key"),
    };
    let mut exons = Vec::new();
    for (&(path, count), dec) in order.iter().zip(decisions) {
        let model = match dec {
            Decision::Skip => continue,
            Decision::Ref(t) => {
                if detected_known.contains(&t) || count < mkc {
                    None
                } else {
                    let m = known_model(ms, ctx.index, t);
                    detected_known.insert(t);
                    m
                }
            }
            Decision::Novel => {
                let range = (path[0].1 as Pos, path[path.len() - 1].1 as Pos);
                get_exons_into(range, &path[1..path.len() - 1], &mut exons);
                novel_model(ms, ctx, aux, &tables, path, count, &exons, &mut next_novel_gene)
            }
        };
        if let Some(m) = model {
            let key = m.key;
            ms.models.push(m);
            if let Some(reads) = paths.paths_to_reads.get(path) {
                for &r in reads {
                    ms.save_assigned_read(r, key);
                }
            }
        }
    }
}

/// Registers the U13 fixture runner.
pub(crate) fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    reg.add::<fixture::U13Fl>();
}

pub(crate) mod fixture {
    //! U13 runner (plumbing: [`crate::fixture_util`]). The fixture `in` has no reference
    //! windows, so the runner reads the capture's genome (`data/ref/genome.fa`, the file
    //! `capture_slice.sh:21` staged; `$IQRS_REFERENCE` overrides).
    use super::*;
    use crate::fixture_util::{compare_stores, StageInput};
    use iq_types::fixture::schema::{StageStoreIn, StoreOut};
    use iq_types::fixture::{FixError, FixResult, FixtureCtx, UnitFixture, Verdict};
    use iq_types::{DevId, Predictions};
    use std::io::{Read, Seek, SeekFrom};
    use std::path::PathBuf;

    fn genome_path(ctx: &FixtureCtx) -> FixResult<PathBuf> {
        if let Some(p) = std::env::var_os("IQRS_REFERENCE") {
            return Ok(PathBuf::from(p));
        }
        let root = ctx.root.canonicalize().unwrap_or_else(|_| ctx.root.clone());
        for a in root.ancestors() {
            let p = a.join("data/ref/genome.fa");
            if p.is_file() {
                return Ok(p);
            }
        }
        Err(FixError("U13: reference genome not found (set IQRS_REFERENCE)".into()))
    }

    /// Reads `chr:[lo, hi]` (1-based closed, clipped) from an faidx-indexed FASTA.
    fn read_region(fa: &std::path::Path, chr: &str, lo: Pos, hi: Pos) -> FixResult<Vec<u8>> {
        let mut fai_path = fa.as_os_str().to_owned();
        fai_path.push(".fai");
        let fai = std::fs::read_to_string(&fai_path)?;
        let line = fai
            .lines()
            .find(|l| l.split('\t').next() == Some(chr))
            .ok_or_else(|| FixError(format!("{chr} not in fai")))?;
        let f: Vec<u64> = line.split('\t').skip(1).map(|x| x.parse().unwrap_or(0)).collect();
        if f.len() < 4 || f[2] == 0 {
            return Err(FixError(format!("bad fai line for {chr}")));
        }
        let (len, off, lb, lw) = (f[0], f[1], f[2], f[3]);
        let lo0 = (lo.max(1) - 1) as u64;
        let hi0 = (hi.max(0) as u64).min(len);
        if lo0 >= hi0 {
            return Ok(Vec::new());
        }
        let n = hi0 - lo0;
        let mut file = std::fs::File::open(fa)?;
        file.seek(SeekFrom::Start(off + (lo0 / lb) * lw + lo0 % lb))?;
        let mut raw = Vec::with_capacity(((n / lb + 2) * lw) as usize);
        file.take((n / lb + 2) * lw).read_to_end(&mut raw)?;
        let mut seq: Vec<u8> = raw.into_iter().filter(|c| !c.is_ascii_whitespace()).collect();
        seq.truncate(n as usize);
        Ok(seq)
    }

    pub struct U13Fl;

    impl UnitFixture for U13Fl {
        const UNIT: &'static str = "U13";
        const FN: &'static str = "construct_fl_isoforms";
        const DEVIATIONS: &'static [DevId] = &[DevId::M1];
        type In = StageStoreIn;
        type Out = StoreOut;

        fn run(input: &StageStoreIn, ctx: &FixtureCtx) -> FixResult<StoreOut> {
            let si = StageInput::load(input, ctx)?;
            // reference window over every FL-path intron
            let (mut lo, mut hi) = (Pos::MAX, Pos::MIN);
            for pa in &si.paths.fl_paths {
                for v in pa.iter().filter(|v| !v.is_terminal()) {
                    lo = lo.min(v.0 as Pos);
                    hi = hi.max(v.1 as Pos);
                }
            }
            let seq = if lo <= hi { read_region(&genome_path(ctx)?, &input.block.chr, lo, hi)? } else { Vec::new() };
            let reference = RefSlice { start: lo.max(1), seq: &seq };
            let names: Vec<&str> = si.index.genes.iter().map(|g| si.ids.gene_name(g.key)).collect();
            let gene_rank = gene_id_ranks(&names);
            let p = &si.params;
            let assigner = LiveAssigner { params: AssignParams::new(p.matching.clone(), p.splice.clone(), false) };
            let preds = Predictions::default();
            let c = Ctx {
                params: p,
                index: &si.index,
                store: &si.block.store,
                graph: &si.graph,
                paths: &si.paths,
                predictions: &preds,
                assigner: &assigner,
            };
            let mut ms = si.store.clone();
            let mut dk = si.detected_known.clone();
            construct_fl(&mut ms, &mut dk, &c, &FlAux { reference, gene_rank: &gene_rank });
            Ok(StoreOut { store: si.store_out(&ms, &dk) })
        }

        fn compare(got: &StoreOut, want: &StoreOut, ctx: &FixtureCtx) -> Verdict {
            compare_stores(&got.store, &want.store, ctx)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exons_from_path() {
        let introns = [Vertex(101, 199), Vertex(301, 399)];
        let mut out = Vec::new();
        get_exons_into((50, 500), &introns, &mut out);
        assert_eq!(out, vec![Iv::new(50, 100), Iv::new(200, 300), Iv::new(400, 500)]);
        // 4.1 drops empty gaps (junctions_from_blocks): start at the intron start
        get_exons_into((101, 500), &introns, &mut out);
        assert_eq!(out, vec![Iv::new(200, 300), Iv::new(400, 500)]);
    }

    #[test]
    fn strands() {
        assert_eq!(clean_strand(2, 0), Strand::Plus);
        assert_eq!(clean_strand(2, 1), Strand::Unknown);
        assert_eq!(clean_strand(0, 0), Strand::Unknown);
        assert_eq!(majority_strand(2, 1, false, false), Strand::Plus);
        assert_eq!(majority_strand(1, 1, true, false), Strand::Plus);
        assert_eq!(majority_strand(1, 1, true, true), Strand::Unknown);
        assert_eq!(majority_strand(0, 0, false, true), Strand::Minus);
        //            pos 10: G T ... A G at 17,18
        let r = RefSlice { start: 10, seq: b"gtaaaaag" };
        assert_eq!(intron_strand_ref(Iv::new(10, 17), &r), Strand::Plus);
        let r = RefSlice { start: 10, seq: b"CTaaaaAC" };
        assert_eq!(intron_strand_ref(Iv::new(10, 17), &r), Strand::Minus);
        assert_eq!(intron_strand_ref(Iv::new(10, 30), &r), Strand::Unknown);
    }

    #[test]
    fn matching_assignment() {
        use iq_types::{IsoformMatch, MatchClass, MatchEvent};
        let mk = |atype, ev: EventType| {
            let mut a = Assignment::empty(atype);
            a.matches.push(IsoformMatch::new_single(
                MatchClass::FullSpliceMatch,
                None,
                Some(IsoRef::Ref(1)),
                MatchEvent::new(ev),
                Strand::Plus,
                0.0,
            ));
            a
        };
        assert!(is_matching_assignment(&mk(AssignmentType::Unique, EventType::IntronRetention)));
        assert!(is_matching_assignment(&mk(AssignmentType::UniqueMinorDifference, EventType::IntronShift)));
        assert!(!is_matching_assignment(&mk(AssignmentType::UniqueMinorDifference, EventType::IntronRetention)));
        assert!(!is_matching_assignment(&mk(AssignmentType::Ambiguous, EventType::Fsm)));
    }

    #[test]
    fn strand_detector_mixed_annotation_uses_reference() {
        use iq_types::{GeneRec, IsoRec};
        // two isoforms share intron (10,17) on opposite strands -> reference decides
        let mut idx = IsoformIndex::default();
        for (k, st) in [(0u32, Strand::Plus), (1, Strand::Minus)] {
            idx.iso.push(IsoRec {
                key: IsoRef::Ref(k),
                gene: GeneKey::Ref(k),
                strand: st,
                span: Iv::new(1, 30),
                exon_off: 0,
                n_exon: 0,
                intron_off: idx.iso_introns.len() as u32,
                n_intron: 1,
            });
            idx.iso_introns.push(Iv::new(10, 17));
            idx.genes.push(GeneRec { key: GeneKey::Ref(k), strand: st, region: Iv::new(1, 30) });
        }
        idx.iso.push(IsoRec {
            key: IsoRef::Ref(2),
            gene: GeneKey::Ref(0),
            strand: Strand::Minus,
            span: Iv::new(1, 30),
            exon_off: 0,
            n_exon: 0,
            intron_off: idx.iso_introns.len() as u32,
            n_intron: 1,
        });
        idx.iso_introns.push(Iv::new(20, 25));
        let r = RefSlice { start: 10, seq: b"GTaaaaAG" };
        let sd = StrandDetector::new(&idx, r);
        assert_eq!(sd.intron_strand(Iv::new(10, 17)), Strand::Plus); // from reference
        assert_eq!(sd.intron_strand(Iv::new(20, 25)), Strand::Minus); // annotated, single strand
        assert_eq!(sd.count_canonical_sites([Iv::new(10, 17), Iv::new(20, 25)]), (1, 1));
    }
}
