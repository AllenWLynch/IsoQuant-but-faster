//! U15b: `TranscriptEndProcessor` (`model_construction/end_processor.py:30-394`): terminal
//! (polyA/TSS) end refinement of novel models, terminal confirmation flags, and
//! alternative-end NIC models derived from known (or, with `--novel_apa`, novel) models.
//!
//! [`Ends`] implements U15a's [`crate::filter::EndProcessor`]; `derive_alternative_end_models`
//! is an inherent method (only `_add_known_alternative_end_models` calls it).
//!
//! Notes:
//! - **Histograms** (`_terminal_histograms`): per side, `(position, read)` entries in read
//!   order; the peak histogram is `position → count` in first-occurrence (4.1 dict) order,
//!   which fixes the float summation order inside [`iq_terminal::peaks`].
//! - **D-T1.** Peak ranking ties go through [`Ends::ties`] (iqrs: position ascending). The
//!   fixture runners accept a 4.1 output reproduced by the reverse tie order as `D-T1`.
//! - **Batching.** 4.1 builds a one-row DataFrame and calls `XGBClassifier.predict` per
//!   model per side. Here a peak call is a direct tree walk per candidate peak with no
//!   per-call setup, `_mark_terminal_confirmation` needs no peak call at all (presence test
//!   only), and `_add_known_alternative_end_models` computes the alternative ends of all
//!   eligible models in parallel (rayon) before the sequential dedup/read-moving pass. That
//!   is exact: a source's read list is only changed by moves onto its own NICs, after its
//!   own derivation (storage with a repeated model key falls back to the sequential order).
//! - **D-E1 (mono-exon end overwrite, dropped):** 4.1 `correct_novel_transcript_ends`
//!   captures `last_exon_left` before editing; for a one-exon model whose start *and* end are
//!   both corrected, the end edit writes `(old_start, new_end)` and undoes the start
//!   correction (regression from the 4.1 refactor fe84ccdc; the 2021 original read it live).
//!   iqrs reads the last exon start live (clamp and write); `quirk(D-E1)` restores 4.1.
//! - `source_known` of a derived NIC stays `None` (4.1 models have no such field and the
//!   fixture conversion never sets it).

use crate::ctx::Ctx;
use crate::filter::EndProcessor;
use iq_terminal::peaks::{detect_peaks_ties, Histogram, Ties};
use iq_terminal::{bundled_models, XgbModel};
use iq_types::{junctions_from_blocks, quirk, DevId, Iv, ModelIdx, ModelKind, ModelStore, Peak, Pos, ReadIdx, Strand, TranscriptModel};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

/// The live end processor (stateless apart from the D-T1 tie order).
#[derive(Debug, Clone, Copy)]
pub struct Ends {
    /// Peak tie order (D-T1). [`Ties::Numpy`] (numpy argsort emulation) in production.
    pub ties: Ties,
}

impl Default for Ends {
    fn default() -> Self {
        Ends { ties: Ties::Numpy }
    }
}

/// One alternative end of a source model: the replaced terminal and the reads routed to it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AltEnd {
    /// `Some` for a start-side (genomic left) alternative end.
    pub start: Option<Pos>,
    /// `Some` for an end-side (genomic right) alternative end.
    pub end: Option<Pos>,
    pub reads: Vec<ReadIdx>,
}

/// Read-terminus entries of one model, per side, in read order (`_terminal_histograms`).
#[derive(Default)]
struct Terminals {
    start: Vec<(Pos, ReadIdx)>,
    end: Vec<(Pos, ReadIdx)>,
}

/// `_confirmed_polya_pos`: polyA-confirmed reads only, at `external_polya_pos` (`+`) /
/// `external_polyt_pos` (`-`).
#[inline]
fn confirmed_polya_pos(r: &iq_types::ReadView<'_>, strand: Strand) -> Option<Pos> {
    if !r.polya_found() {
        return None;
    }
    let p = r.polya();
    match strand {
        Strand::Plus if p.ext_a != -1 => Some(p.ext_a),
        Strand::Minus if p.ext_t != -1 => Some(p.ext_t),
        _ => None,
    }
}

/// `_terminal_histograms(source_model, reads)` as `(position, read)` entries.
fn terminals(model: &TranscriptModel, reads: &[ReadIdx], ctx: &Ctx<'_>) -> Terminals {
    let fer = model.exons[0].e;
    let lel = model.exons[model.exons.len() - 1].s;
    let mut t = Terminals { start: Vec::with_capacity(reads.len()), end: Vec::with_capacity(reads.len()) };
    for &ri in reads {
        let r = ctx.store.view(ri as usize);
        let (s, e) = match model.strand {
            Strand::Plus => (Some(r.start()), confirmed_polya_pos(&r, Strand::Plus)),
            Strand::Minus => (confirmed_polya_pos(&r, Strand::Minus), Some(r.end())),
            Strand::Unknown => (Some(r.start()), Some(r.end())),
        };
        if let Some(s) = s.filter(|&s| s < fer) {
            t.start.push((s, ri));
        }
        if let Some(e) = e.filter(|&e| e > lel) {
            t.end.push((e, ri));
        }
    }
    t
}

/// `{p: len(rs) for p, rs in reads_map.items()}`.
fn histogram(entries: &[(Pos, ReadIdx)]) -> Histogram {
    let mut h = Histogram::with_capacity(entries.len().min(64));
    for &(p, _) in entries {
        *h.entry(p).or_insert(0) += 1;
    }
    h
}

#[inline]
fn within(a: Pos, b: Pos, d: i32) -> bool {
    (a as i64 - b as i64).abs() <= d as i64
}

/// Intron chain key (`junctions_from_blocks`).
fn chain_key(exons: &[Iv]) -> Vec<Iv> {
    junctions_from_blocks(exons)
}

/// `(chain, [(start, end)])` pairs of the block's reference isoforms
/// (`gene_info.all_isoforms_exons`).
fn reference_pairs(ctx: &Ctx<'_>) -> HashMap<Vec<Iv>, Vec<(Pos, Pos)>> {
    let idx = ctx.index;
    let mut m: HashMap<Vec<Iv>, Vec<(Pos, Pos)>> = HashMap::with_capacity(idx.n_iso());
    for i in 0..idx.n_iso() {
        let ex = idx.iso_exons(i);
        if ex.is_empty() {
            continue;
        }
        m.entry(chain_key(ex)).or_default().push((ex[0].s, ex[ex.len() - 1].e));
    }
    m
}

/// `for r in reads: if len(transcript_read_ids[src]) <= 1: break; move_read_to_model(r, src, dst)`
/// in one pass over the source list. `ModelStore::move_read_to_model` does a linear
/// `list.remove` per read, which is quadratic when an alternative end claims a large share of a
/// deep known model (fbxw7 on Revio: ~10^5 reads moved, 60 s). Same final state: each move
/// shrinks the source by one, so the loop moves `min(len(reads), len0 - 1)` reads; removing
/// the first occurrence of each (distinct, `claimed` dedups) keeps the remaining order;
/// `internal_counter[src]` is touched (inserted) before the target, as in the first move.
fn move_reads_batched(ms: &mut ModelStore, reads: &[ReadIdx], src: ModelIdx, dst: ModelIdx) {
    let len0 = ms.transcript_read_ids.get(&src).map_or(0, Vec::len);
    let m = reads.len().min(len0.saturating_sub(1));
    if m == 0 {
        return;
    }
    let moved = &reads[..m];
    let mut pending: HashSet<ReadIdx> = moved.iter().copied().collect();
    debug_assert_eq!(pending.len(), m, "alternative-end reads are distinct");
    let v = ms.transcript_read_ids.get_mut(&src).expect("source list");
    v.retain(|r| !pending.remove(r));
    assert!(pending.is_empty(), "list.remove(x): x not in list");
    *ms.internal_counter.entry(src).or_default() -= m as i32;
    for &r in moved {
        ms.read_assignment_counts[r as usize] -= 1;
        ms.save_assigned_read(r, dst);
    }
}

fn any_close(pairs: Option<&Vec<(Pos, Pos)>>, s: Pos, e: Pos, d: i32) -> bool {
    pairs.is_some_and(|v| v.iter().any(|&(ps, pe)| within(s, ps, d) && within(e, pe, d)))
}

impl Ends {
    /// `_terminal_model(strand, left)`: polyA side always; TSS side only with
    /// `use_tss_model`; unstranded → `None` (positional fallback).
    fn terminal_model(&self, strand: Strand, left: bool, ctx: &Ctx<'_>) -> Option<&'static XgbModel> {
        let polya_side = match strand {
            Strand::Plus => !left,
            Strand::Minus => left,
            Strand::Unknown => return None,
        };
        let m = bundled_models();
        if polya_side {
            Some(&m.polya)
        } else if ctx.params.use_tss_model {
            Some(&m.tss)
        } else {
            None
        }
    }

    fn peaks(&self, hist: &Histogram, model: &XgbModel) -> Vec<Peak> {
        detect_peaks_ties(hist, model, self.ties).map(|d| d.peaks).unwrap_or_default()
    }

    /// `_peak_boundary`: position of the first highest-count accepted peak passing `valid`.
    fn peak_boundary(&self, hist: &Histogram, model: &XgbModel, valid: impl Fn(Pos) -> bool) -> Option<Pos> {
        if hist.is_empty() {
            return None;
        }
        let mut best: Option<Peak> = None;
        for p in self.peaks(hist, model).into_iter().filter(|p| valid(p.position)) {
            if best.is_none_or(|b| p.count > b.count) {
                best = Some(p);
            }
        }
        best.map(|p| p.position)
    }

    /// `_alternative_end_positions`.
    fn alternative_end_positions(
        &self,
        hist: &Histogram,
        model: &XgbModel,
        annotated: Pos,
        clamp: impl Fn(Pos) -> bool,
        ctx: &Ctx<'_>,
    ) -> Vec<Pos> {
        if hist.is_empty() {
            return Vec::new();
        }
        let peaks: Vec<Peak> = self.peaks(hist, model).into_iter().filter(|p| clamp(p.position)).collect();
        let Some(dominant) = peaks.iter().map(|p| p.count).max() else { return Vec::new() };
        let p = ctx.params;
        let cutoff = (p.model.min_novel_count as f64).max(p.model.terminal_position_rel * dominant as f64);
        peaks
            .iter()
            .filter(|pk| !within(pk.position, annotated, p.matching.apa_delta) && pk.count as f64 >= cutoff)
            .map(|pk| pk.position)
            .collect()
    }

    /// `correct_novel_transcript_ends(model, reads)`.
    pub fn correct_ends(&self, model: &mut TranscriptModel, reads: &[ReadIdx], ctx: &Ctx<'_>) {
        let d = ctx.params.matching.apa_delta;
        let ts = model.start();
        let te = model.end();
        let fer = model.exons[0].e;
        let last = model.exons.len() - 1;
        let lel = model.exons[last].s;
        let t = terminals(model, reads, ctx);
        let sh = histogram(&t.start);
        let eh = histogram(&t.end);

        if !sh.keys().any(|&p| within(p, ts, d)) {
            let new = match self.terminal_model(model.strand, true, ctx) {
                Some(m) => self.peak_boundary(&sh, m, |p| p < fer),
                // `_closest_inward(sorted(keys), start, greater=True)`
                None => sh.keys().copied().filter(|&p| p > ts).min(),
            };
            if let Some(ns) = new {
                if ns != ts && ns < fer {
                    model.exons[0] = Iv::new(ns, fer);
                }
            }
        }
        if !eh.keys().any(|&p| within(p, te, d)) {
            let new = match self.terminal_model(model.strand, false, ctx) {
                Some(m) => self.peak_boundary(&eh, m, |p| p > lel),
                None => eh.keys().copied().filter(|&p| p < te).max(),
            };
            if let Some(ne) = new {
                // D-E1: 4.1 uses `last_exon_left` captured before the start edit, so a
                // mono-exon model loses its start correction; iqrs reads it live (as the
                // 2021 original, 0bd9b1e2). Identical for multi-exon models.
                let lel_now = if quirk(DevId::E1) { lel } else { model.exons[last].s };
                if ne != te && ne > lel_now {
                    model.exons[last] = Iv::new(lel_now, ne);
                }
            }
        }
    }

    /// `_mark_terminal_confirmation(model)`; `reads` = `store.transcript_read_ids[model.id]`.
    pub fn mark_confirmation(&self, model: &mut TranscriptModel, reads: &[ReadIdx], ctx: &Ctx<'_>) {
        if reads.is_empty() {
            model.polya_confirmed = false;
            model.tss_confirmed = false;
            return;
        }
        let d = ctx.params.matching.apa_delta;
        let t = terminals(model, reads, ctx);
        let (ms, me) = (model.start(), model.end());
        let ss = t.start.iter().any(|&(p, _)| within(p, ms, d));
        let es = t.end.iter().any(|&(p, _)| within(p, me, d));
        match model.strand {
            Strand::Plus => (model.tss_confirmed, model.polya_confirmed) = (ss, es),
            Strand::Minus => (model.polya_confirmed, model.tss_confirmed) = (ss, es),
            Strand::Unknown => (model.polya_confirmed, model.tss_confirmed) = (false, false),
        }
    }

    /// The alternative ends `derive_alternative_end_models` would build (pure; no ids).
    pub fn alternative_ends(&self, source: &TranscriptModel, reads: &[ReadIdx], ctx: &Ctx<'_>) -> Vec<AltEnd> {
        if reads.is_empty() {
            return Vec::new();
        }
        let d = ctx.params.matching.apa_delta;
        let ann_start = source.start();
        let ann_end = source.end();
        let fer = source.exons[0].e;
        let lel = source.exons[source.exons.len() - 1].s;
        let sm = self.terminal_model(source.strand, true, ctx);
        let em = self.terminal_model(source.strand, false, ctx);
        if sm.is_none() && em.is_none() {
            return Vec::new();
        }
        let t = terminals(source, reads, ctx);
        let mut claimed: HashSet<ReadIdx> = HashSet::new();
        let mut out = Vec::new();

        // assign_to_peaks: per position (dict order), the nearest of alt peaks + annotated
        // end (ties → smaller coordinate); reads go to an alt peak within apa_delta, once.
        // assign_to_peaks: per position (dict order), the nearest of alt peaks + annotated
        // end (ties → smaller coordinate); reads go to an alt peak within apa_delta, once.
        // Peak positions are distinct (find_peaks locations), so a slot index stands for
        // 4.1's position-keyed `per_peak` dict.
        let mut side = |entries: &[(Pos, ReadIdx)], peaks: &[Pos], annotated: Pos, start_side: bool, out: &mut Vec<AltEnd>| {
            if peaks.is_empty() {
                return;
            }
            let mut per_peak: Vec<Vec<ReadIdx>> = vec![Vec::new(); peaks.len()];
            let hist = histogram(entries);
            let mut order: Vec<(usize, ReadIdx)> = entries.iter().map(|&(p, r)| (hist.get_index_of(&p).unwrap(), r)).collect();
            order.sort_by_key(|x| x.0); // stable: read order within a position
            let mut k = 0;
            for (pi, &p) in hist.keys().enumerate() {
                let k0 = k;
                while k < order.len() && order[k].0 == pi {
                    k += 1;
                }
                let dist = |c: Pos| (c as i64 - p as i64).abs();
                let mut best = (dist(annotated), annotated);
                let mut slot = None;
                for (i, &c) in peaks.iter().enumerate() {
                    if (dist(c), c) < best {
                        best = (dist(c), c);
                        slot = Some(i);
                    }
                }
                if let Some(i) = slot.filter(|_| best.0 <= d as i64) {
                    for &(_, r) in &order[k0..k] {
                        if claimed.insert(r) {
                            per_peak[i].push(r);
                        }
                    }
                }
            }
            for (&pos, reads) in peaks.iter().zip(per_peak) {
                out.push(if start_side {
                    AltEnd { start: Some(pos), end: None, reads }
                } else {
                    AltEnd { start: None, end: Some(pos), reads }
                });
            }
        };
        if let Some(m) = sm {
            let peaks = self.alternative_end_positions(&histogram(&t.start), m, ann_start, |p| p < fer, ctx);
            side(&t.start, &peaks, ann_start, true, &mut out);
        }
        if let Some(m) = em {
            let peaks = self.alternative_end_positions(&histogram(&t.end), m, ann_end, |p| p > lel, ctx);
            side(&t.end, &peaks, ann_end, false, &mut out);
        }
        out
    }

    /// `_nic_model_with_boundary`: copy of `source` with one terminal replaced.
    pub fn nic_model(source: &TranscriptModel, alt: &AltEnd, key: ModelIdx) -> TranscriptModel {
        let kind = match source.kind {
            ModelKind::Known(_) => ModelKind::Nic,
            k => k,
        };
        let mut exons = source.exons.clone();
        if let Some(s) = alt.start {
            exons[0] = Iv::new(s, exons[0].e);
        }
        if let Some(e) = alt.end {
            let l = exons.len() - 1;
            exons[l] = Iv::new(exons[l].s, e);
        }
        TranscriptModel {
            key,
            strand: source.strand,
            gene: source.gene,
            kind,
            exons,
            intron_path: source.intron_path.clone(),
            polya_confirmed: false,
            tss_confirmed: false,
            source_known: None,
        }
    }

    /// `derive_alternative_end_models(source, reads)`: new models (keys from `next_key`, in
    /// 4.1 id order) with their supporting reads.
    pub fn derive_alternative_end_models(
        &self,
        source: &TranscriptModel,
        reads: &[ReadIdx],
        ctx: &Ctx<'_>,
        next_key: &mut ModelIdx,
    ) -> Vec<(TranscriptModel, Vec<ReadIdx>)> {
        self.alternative_ends(source, reads, ctx)
            .into_iter()
            .map(|a| {
                let k = *next_key;
                *next_key += 1;
                (Self::nic_model(source, &a, k), a.reads)
            })
            .collect()
    }

    /// `_add_known_alternative_end_models()`.
    pub fn add_known_alt_end_models(&self, ms: &mut ModelStore, ctx: &Ctx<'_>) {
        let d = ctx.params.matching.apa_delta;
        let novel_apa = ctx.params.novel_apa;
        let mut existing = reference_pairs(ctx);
        for m in &ms.models {
            existing.entry(chain_key(&m.exons)).or_default().push((m.start(), m.end()));
        }
        let eligible: Vec<usize> =
            (0..ms.models.len()).filter(|&i| novel_apa || matches!(ms.models[i].kind, ModelKind::Known(_))).collect();
        let derive = |i: usize, ms: &ModelStore| {
            let m = &ms.models[i];
            let reads = ms.transcript_read_ids.get(&m.key).map_or(&[][..], |v| v.as_slice());
            self.alternative_ends(m, reads, ctx)
        };
        let unique_keys = {
            let mut seen = HashSet::with_capacity(eligible.len());
            eligible.iter().all(|&i| seen.insert(ms.models[i].key))
        };
        let mut pre: Vec<Option<Vec<AltEnd>>> = if unique_keys {
            let ms_ref = &*ms;
            eligible.par_iter().map(|&i| Some(derive(i, ms_ref))).collect()
        } else {
            vec![None; eligible.len()]
        };

        let mut new_models = Vec::new();
        for (j, &i) in eligible.iter().enumerate() {
            let src = ms.models[i].key;
            // `self.store.transcript_read_ids[source_id]` (defaultdict entry)
            ms.transcript_read_ids.entry(src).or_default();
            let alts = pre[j].take().unwrap_or_else(|| derive(i, ms));
            for alt in alts {
                let key = ms.new_key();
                let nic = Self::nic_model(&ms.models[i], &alt, key);
                let ck = chain_key(&nic.exons);
                let (ns, ne) = (nic.start(), nic.end());
                if any_close(existing.get(&ck), ns, ne, d) {
                    continue;
                }
                existing.entry(ck).or_default().push((ns, ne));
                new_models.push(nic);
                move_reads_batched(ms, &alt.reads, src, key);
            }
        }
        ms.models.extend(new_models);
    }

    /// `_drop_duplicate_alt_end_models()`.
    pub fn drop_duplicate_alt_end_models(&self, ms: &mut ModelStore, ctx: &Ctx<'_>) {
        let d = ctx.params.matching.apa_delta;
        let refs = reference_pairs(ctx);
        let mut kept_pairs: HashMap<Vec<Iv>, Vec<(Pos, Pos)>> = HashMap::new();
        let models = std::mem::take(&mut ms.models);
        let mut kept = Vec::with_capacity(models.len());
        for m in models {
            if matches!(m.kind, ModelKind::Known(_)) {
                kept.push(m);
                continue;
            }
            let ck = chain_key(&m.exons);
            let (s, e) = (m.start(), m.end());
            if any_close(refs.get(&ck), s, e, d) || any_close(kept_pairs.get(&ck), s, e, d) {
                if ms.internal_counter.contains_key(&m.key) {
                    ms.delete_from_storage(m.key);
                }
                continue;
            }
            kept_pairs.entry(ck).or_default().push((s, e));
            kept.push(m);
        }
        ms.models = kept;
    }
}

impl EndProcessor for Ends {
    fn correct_novel_transcript_ends(&self, model: &mut TranscriptModel, reads: &[ReadIdx], ctx: &Ctx<'_>) {
        self.correct_ends(model, reads, ctx)
    }
    fn mark_terminal_confirmation(&self, model: &mut TranscriptModel, reads: &[ReadIdx], ctx: &Ctx<'_>) {
        self.mark_confirmation(model, reads, ctx)
    }
    fn add_known_alternative_end_models(&self, ms: &mut ModelStore, ctx: &Ctx<'_>) {
        self.add_known_alt_end_models(ms, ctx)
    }
    fn drop_duplicate_alt_end_models(&self, ms: &mut ModelStore, ctx: &Ctx<'_>) {
        Ends::drop_duplicate_alt_end_models(self, ms, ctx)
    }
}

// ------------------------------------------------------------------------------------------
// Fixture runners
// ------------------------------------------------------------------------------------------

pub(crate) fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    reg.add::<fixture::U15bCorrectEnds>();
    reg.add::<fixture::U15bMarkConfirmation>();
    reg.add::<fixture::U15bDerive>();
    reg.add::<fixture::U15bAddKnown>();
    reg.add::<fixture::U15bDropDuplicates>();
}

pub(crate) mod fixture {
    use super::*;
    use crate::ctx::ReadAssigner;
    use crate::fixture_util::{compare_stores, StageInput};
    use iq_types::fixture::convert::ModelExtras;
    use iq_types::fixture::schema::{ModelFixture, ModelOut, StageModelIn, StageStoreIn, StoreOut, U15bDeriveAltEndOut};
    use iq_types::fixture::{FixResult, FixtureCtx, UnitFixture, Verdict};
    use iq_types::{
        Assignment, CombinedProfile, DevId, IntronGraph, IsoformIndex, PathStorage, Predictions, TranscriptModelType,
    };

    /// End processing never assigns reads.
    pub(crate) struct NoAssign;
    impl ReadAssigner for NoAssign {
        fn assign(&self, _: &IsoformIndex, _: &CombinedProfile, _: bool, _: u64) -> Option<Assignment> {
            None
        }
    }

    /// `UnitFixture::judge` (parity = all quirks; default = none, then one DevId of `devs` at
    /// a time) with a D-T1 fallback: a mismatch that the reverse peak tie order reproduces
    /// exactly is `Explained(D-T1)`.
    pub(crate) fn judge_ties<I, O>(
        run: impl Fn(&I, &FixtureCtx, Ends) -> FixResult<O>,
        compare: impl Fn(&O, &O, &FixtureCtx) -> Verdict,
        devs: &[DevId],
        input: &I,
        want: &O,
        ctx: &FixtureCtx,
    ) -> Verdict {
        use iq_types::fixture::diff::Diff;
        use iq_types::fixture::RunMode;
        use iq_types::{with_quirks, QuirkSet};
        use std::panic::{catch_unwind, AssertUnwindSafe};
        let attempt = |set: QuirkSet, e: Ends| match catch_unwind(AssertUnwindSafe(|| with_quirks(set, || run(input, ctx, e)))) {
            Ok(Ok(got)) => compare(&got, want, ctx),
            Ok(Err(err)) => Verdict::Mismatch(Diff::error(format!("run failed: {err}"))),
            Err(_) => Verdict::Mismatch(Diff::error("panic")),
        };
        let reverse = Ends { ties: Ties::PositionDesc };
        let base = if ctx.mode == RunMode::Parity { QuirkSet::all() } else { QuirkSet::NONE };
        let v = attempt(base, Ends::default());
        if v.is_ok() {
            return v;
        }
        // fixtures come from a SIMD CPU, where large tied peak arrays may sort differently
        if attempt(base, reverse) == Verdict::Exact || attempt(base, Ends { ties: Ties::PositionAsc }) == Verdict::Exact {
            return Verdict::Explained(DevId::T1);
        }
        if ctx.mode == RunMode::Default && iq_types::QUIRKS41 {
            for &d in devs.iter().filter(|&&d| d != DevId::T1) {
                if attempt(QuirkSet::only(d), Ends::default()) == Verdict::Exact {
                    return Verdict::Explained(d);
                }
            }
        }
        v
    }

    /// Per-model case: converted model + its extras, block, index, params.
    fn model_case<R>(input: &StageModelIn, ctx: &FixtureCtx, f: impl FnOnce(TranscriptModel, ModelExtras, &mut iq_types::fixture::ids::IdMaps, &Ctx<'_>) -> R) -> FixResult<R> {
        let params = ctx.params()?;
        let block = ctx.block(&input.block)?;
        let index = ctx.index(&input.index)?;
        let mut ids = ctx.ids().clone();
        ids.reset_block();
        let (model, ex) = input.model.to_core(&mut ids);
        let graph = IntronGraph::default();
        let paths = PathStorage::default();
        let preds = Predictions::default();
        let c = Ctx {
            params: &params,
            index: &index,
            store: &block.store,
            graph: &graph,
            paths: &paths,
            predictions: &preds,
            assigner: &NoAssign,
        };
        Ok(f(model, ex, &mut ids, &c))
    }

    fn run_model(input: &StageModelIn, ctx: &FixtureCtx, e: Ends, mark: bool) -> FixResult<ModelOut> {
        model_case(input, ctx, |mut m, ex, ids, c| {
            if mark {
                e.mark_confirmation(&mut m, &input.reads, c);
            } else {
                e.correct_ends(&mut m, &input.reads, c);
            }
            ModelOut { model: ModelFixture::from_core(&m, ids, &ex) }
        })
    }

    fn default_compare<O: PartialEq + iq_types::serde::Serialize>(got: &O, want: &O, ctx: &FixtureCtx) -> Verdict {
        if got == want {
            return Verdict::Exact;
        }
        match iq_types::fixture::diff::diff_values(got, want, &ctx.diff_opts) {
            Some(d) => Verdict::Mismatch(d),
            None => Verdict::Mismatch(iq_types::fixture::diff::Diff::error("PartialEq differs")),
        }
    }

    pub struct U15bCorrectEnds;
    impl UnitFixture for U15bCorrectEnds {
        const UNIT: &'static str = "U15b";
        const FN: &'static str = "correct_novel_transcript_ends";
        const DEVIATIONS: &'static [DevId] = &[DevId::T1, DevId::E1];
        type In = StageModelIn;
        type Out = ModelOut;
        fn run(input: &StageModelIn, ctx: &FixtureCtx) -> FixResult<ModelOut> {
            run_model(input, ctx, Ends::default(), false)
        }
        fn judge(input: &StageModelIn, want: &ModelOut, ctx: &FixtureCtx) -> Verdict {
            judge_ties(|i, c, e| run_model(i, c, e, false), default_compare, Self::DEVIATIONS, input, want, ctx)
        }
    }

    pub struct U15bMarkConfirmation;
    impl UnitFixture for U15bMarkConfirmation {
        const UNIT: &'static str = "U15b";
        const FN: &'static str = "_mark_terminal_confirmation";
        type In = StageModelIn;
        type Out = ModelOut;
        fn run(input: &StageModelIn, ctx: &FixtureCtx) -> FixResult<ModelOut> {
            run_model(input, ctx, Ends::default(), true)
        }
    }

    fn run_derive(input: &StageModelIn, ctx: &FixtureCtx, e: Ends) -> FixResult<U15bDeriveAltEndOut> {
        let counter = input.id_counter.unwrap_or(0);
        model_case(input, ctx, |m, ex, ids, c| {
            let mut next = ids.models.len() as ModelIdx;
            let res = e.derive_alternative_end_models(&m, &input.reads, c, &mut next);
            let chr = &input.model.chr;
            let result = res
                .into_iter()
                .enumerate()
                .map(|(i, (nic, reads))| {
                    // 4.1 id: `transcript<get_transcript_id()>.<chr>.<nic|nnic>` (D-M1: iqrs
                    // keys stand in for these ids; the name is rebuilt here for comparison)
                    let suffix = if nic.kind.model_type() == TranscriptModelType::NovelInCatalog { "nic" } else { "nnic" };
                    let name = format!("transcript{}.{chr}.{suffix}", counter + 1 + i as i64);
                    let k = ids.models.intern(&name);
                    debug_assert_eq!(k, nic.key);
                    let nex = ModelExtras { chr: ex.chr.clone(), source: "IsoQuant".to_owned(), ..Default::default() };
                    (ModelFixture::from_core(&nic, ids, &nex), reads)
                })
                .collect::<Vec<_>>();
            U15bDeriveAltEndOut { id_counter: counter + result.len() as i64, result }
        })
    }

    pub struct U15bDerive;
    impl UnitFixture for U15bDerive {
        const UNIT: &'static str = "U15b";
        const FN: &'static str = "derive_alternative_end_models";
        const DEVIATIONS: &'static [DevId] = &[DevId::T1];
        type In = StageModelIn;
        type Out = U15bDeriveAltEndOut;
        fn run(input: &StageModelIn, ctx: &FixtureCtx) -> FixResult<U15bDeriveAltEndOut> {
            run_derive(input, ctx, Ends::default())
        }
        fn judge(input: &StageModelIn, want: &U15bDeriveAltEndOut, ctx: &FixtureCtx) -> Verdict {
            judge_ties(run_derive, default_compare, Self::DEVIATIONS, input, want, ctx)
        }
    }

    fn run_store(input: &StageStoreIn, ctx: &FixtureCtx, e: Ends) -> FixResult<StoreOut> {
        let si = StageInput::load(input, ctx)?;
        let mut ms = si.store.clone();
        let preds = Predictions::default();
        let c = Ctx {
            params: &si.params,
            index: &si.index,
            store: &si.block.store,
            graph: &si.graph,
            paths: &si.paths,
            predictions: &preds,
            assigner: &NoAssign,
        };
        if input.fn_ == "_add_known_alternative_end_models" {
            e.add_known_alt_end_models(&mut ms, &c);
        } else {
            Ends::drop_duplicate_alt_end_models(&e, &mut ms, &c);
        }
        Ok(StoreOut { store: si.store_out(&ms, &si.detected_known) })
    }

    fn store_compare(got: &StoreOut, want: &StoreOut, ctx: &FixtureCtx) -> Verdict {
        compare_stores(&got.store, &want.store, ctx)
    }

    pub struct U15bAddKnown;
    impl UnitFixture for U15bAddKnown {
        const UNIT: &'static str = "U15b";
        const FN: &'static str = "_add_known_alternative_end_models";
        const DEVIATIONS: &'static [DevId] = &[DevId::T1];
        type In = StageStoreIn;
        type Out = StoreOut;
        fn run(input: &StageStoreIn, ctx: &FixtureCtx) -> FixResult<StoreOut> {
            run_store(input, ctx, Ends::default())
        }
        fn compare(got: &StoreOut, want: &StoreOut, ctx: &FixtureCtx) -> Verdict {
            store_compare(got, want, ctx)
        }
        fn judge(input: &StageStoreIn, want: &StoreOut, ctx: &FixtureCtx) -> Verdict {
            judge_ties(run_store, store_compare, Self::DEVIATIONS, input, want, ctx)
        }
    }

    pub struct U15bDropDuplicates;
    impl UnitFixture for U15bDropDuplicates {
        const UNIT: &'static str = "U15b";
        const FN: &'static str = "_drop_duplicate_alt_end_models";
        type In = StageStoreIn;
        type Out = StoreOut;
        fn run(input: &StageStoreIn, ctx: &FixtureCtx) -> FixResult<StoreOut> {
            run_store(input, ctx, Ends::default())
        }
        fn compare(got: &StoreOut, want: &StoreOut, ctx: &FixtureCtx) -> Verdict {
            store_compare(got, want, ctx)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::fixture::NoAssign;
    use super::*;
    use iq_types::{Assignment, AssignmentType, BlockStore, GeneKey, IntronGraph, IsoformIndex, Params, PathStorage, PolyAInfo, Predictions, ReadRecord, NO_POLYA};

    #[test]
    fn batched_moves_match_sequential_list_remove() {
        // Source list with duplicates; cases cover the `len <= 1` break and a missing target key.
        let src_reads: Vec<ReadIdx> = vec![3, 1, 4, 1, 5, 9, 2, 6, 5, 3];
        for alt in [vec![1, 5, 3], vec![9], vec![6, 2, 9, 4, 1, 5, 3], (0..7).collect::<Vec<_>>()[1..].to_vec()] {
            let alt: Vec<ReadIdx> = alt.into_iter().filter(|r| src_reads.contains(r)).collect();
            for keep in [src_reads.len(), 3, 2, 1] {
                let mk = || {
                    let mut ms = ModelStore::new(10);
                    ms.transcript_read_ids.insert(7, Vec::new());
                    for &r in &src_reads[..keep] {
                        ms.save_assigned_read(r, 0);
                    }
                    ms.internal_counter.shift_remove(&0); // exercise the counter insertion order
                    ms
                };
                let alt: Vec<ReadIdx> = alt.iter().copied().filter(|r| src_reads[..keep].contains(r)).collect();
                let mut a = mk();
                for &r in &alt {
                    if a.transcript_read_ids.get(&0).map_or(0, Vec::len) <= 1 {
                        break;
                    }
                    a.move_read_to_model(r, 0, 8);
                }
                let mut b = mk();
                move_reads_batched(&mut b, &alt, 0, 8);
                assert_eq!(a.transcript_read_ids.iter().collect::<Vec<_>>(), b.transcript_read_ids.iter().collect::<Vec<_>>(), "alt {alt:?} keep {keep}");
                assert_eq!(a.internal_counter.iter().collect::<Vec<_>>(), b.internal_counter.iter().collect::<Vec<_>>());
                assert_eq!(a.read_assignment_counts, b.read_assignment_counts);
                assert_eq!(a.construction_assignment, b.construction_assignment);
            }
        }
    }

    fn read(idx: u32, s: Pos, e: Pos, polya: PolyAInfo) -> ReadRecord {
        ReadRecord {
            idx,
            name: format!("r{idx}"),
            exons: vec![Iv::new(s, e)],
            corrected_exons: vec![Iv::new(s, e)],
            boundary_deletions: vec![],
            polya,
            polya_found: polya != NO_POLYA,
            strand: Strand::Plus,
            mapped_strand: Strand::Plus,
            mapq: 60,
            multimapper: false,
            group: 0,
            assignment: Assignment::new(AssignmentType::Noninformative, Default::default()),
            genomic_region: Iv::new(1, 1_000_000),
        }
    }

    fn model(key: ModelIdx, kind: ModelKind, strand: Strand, exons: &[(Pos, Pos)]) -> TranscriptModel {
        TranscriptModel {
            key,
            strand,
            gene: GeneKey::Novel(0),
            kind,
            exons: exons.iter().map(|&(s, e)| Iv::new(s, e)).collect(),
            intron_path: Box::new([]),
            polya_confirmed: false,
            tss_confirmed: false,
            source_known: None,
        }
    }

    fn with_ctx<R>(reads: &[ReadRecord], tss: bool, f: impl FnOnce(&Ctx<'_>) -> R) -> R {
        let mut params = Params::pacbio();
        params.use_tss_model = tss;
        let store = BlockStore::from_records(reads);
        let index = IsoformIndex::default();
        let graph = IntronGraph::default();
        let paths = PathStorage::default();
        let preds = Predictions::default();
        let c = Ctx { params: &params, index: &index, store: &store, graph: &graph, paths: &paths, predictions: &preds, assigner: &NoAssign };
        f(&c)
    }

    #[test]
    fn mark_confirmation_by_strand() {
        let pa = PolyAInfo { ext_a: 1990, ext_t: -1, int_a: -1, int_t: -1 };
        let reads = [read(0, 1010, 1990, pa), read(1, 1300, 1700, NO_POLYA)];
        with_ctx(&reads, true, |c| {
            let e = Ends::default();
            let mut m = model(0, ModelKind::Nnic, Strand::Plus, &[(1000, 1200), (1500, 2000)]);
            e.mark_confirmation(&mut m, &[0, 1], c);
            assert!(m.tss_confirmed && m.polya_confirmed);
            // minus strand: polyA side is the left, needs ext_t → not confirmed; TSS (right end) is
            m.strand = Strand::Minus;
            e.mark_confirmation(&mut m, &[0, 1], c);
            assert!(!m.polya_confirmed && m.tss_confirmed);
            m.strand = Strand::Unknown;
            e.mark_confirmation(&mut m, &[0, 1], c);
            assert!(!m.polya_confirmed && !m.tss_confirmed);
            m.strand = Strand::Plus;
            m.tss_confirmed = true;
            e.mark_confirmation(&mut m, &[], c);
            assert!(!m.polya_confirmed && !m.tss_confirmed);
        });
    }

    #[test]
    fn closest_inward_fallback_without_tss_model() {
        // '+' strand, no TSS model: start moves to the nearest read start inside (> start)
        let reads = [read(0, 1100, 1900, NO_POLYA), read(1, 1150, 1900, NO_POLYA), read(2, 900, 1900, NO_POLYA)];
        with_ctx(&reads, false, |c| {
            let mut m = model(0, ModelKind::Nnic, Strand::Plus, &[(1000, 1200), (1500, 2000)]);
            Ends::default().correct_ends(&mut m, &[0, 1, 2], c);
            // 900 is outward; 1100 is the first start > 1000 (and |1100-1000| > apa_delta)
            assert_eq!(m.exons[0], Iv::new(1100, 1200));
            // end side: polyA model, no confirmed reads → empty histogram → unchanged
            assert_eq!(m.exons[1], Iv::new(1500, 2000));
        });
    }

    #[test]
    fn mono_exon_end_edit_d_e1() {
        // unstranded → positional fallback on both sides
        let reads = [read(0, 1100, 1800, NO_POLYA)];
        with_ctx(&reads, false, |c| {
            let run = |set: iq_types::QuirkSet| {
                let mut m = model(0, ModelKind::Nnic, Strand::Unknown, &[(1000, 2000)]);
                iq_types::with_quirks(set, || Ends::default().correct_ends(&mut m, &[0], c));
                m.exons
            };
            // iqrs: both corrections kept
            assert_eq!(run(iq_types::QuirkSet::NONE), vec![Iv::new(1100, 1800)]);
            if iq_types::QUIRKS41 {
                // 4.1: the end edit restores the old start
                assert_eq!(run(iq_types::QuirkSet::only(DevId::E1)), vec![Iv::new(1000, 1800)]);
            }
        });
    }

    #[test]
    fn drop_duplicates_keeps_first_and_frees_reads() {
        let reads = [read(0, 1000, 2000, NO_POLYA), read(1, 1000, 2000, NO_POLYA)];
        with_ctx(&reads, false, |c| {
            let mut ms = ModelStore::new(2);
            let a = ms.new_key();
            let b = ms.new_key();
            let k = ms.new_key();
            ms.models.push(model(k, ModelKind::Known(0), Strand::Plus, &[(1000, 1200), (1500, 2000)]));
            ms.models.push(model(a, ModelKind::Nnic, Strand::Plus, &[(1010, 1200), (1500, 2040)]));
            ms.models.push(model(b, ModelKind::Nic, Strand::Plus, &[(1030, 1200), (1500, 1990)]));
            ms.save_assigned_read(0, a);
            ms.save_assigned_read(1, b);
            Ends::default().drop_duplicate_alt_end_models(&mut ms, c);
            assert_eq!(ms.models.iter().map(|m| m.key).collect::<Vec<_>>(), vec![k, a]);
            assert!(!ms.transcript_read_ids.contains_key(&b));
            assert_eq!(ms.read_assignment_counts, vec![1, 0]);
        });
    }

    /// Integration: U15a `filter_transcripts` fixtures with the live [`Ends`] instead of the
    /// replay stub. `IQRS_FIXTURES` (default `../fixtures` next to the workspace).
    /// `cargo test -p iq_model --release --features quirks41 -- --ignored filter_with_live_ends --nocapture`
    #[test]
    #[ignore]
    fn filter_with_live_ends() {
        use crate::filter::filter;
        use crate::fixture_util::{compare_stores, StageInput};
        use crate::fl::LiveAssigner;
        use iq_assign::AssignParams;
        use iq_types::fixture::schema::{StageStoreIn, StoreOut};
        use iq_types::fixture::{fixtures_root, run_path, FixResult, FixtureCtx, Registry, RunMode, Summary, UnitFixture, Verdict};

        struct LiveFilter;
        impl UnitFixture for LiveFilter {
            const UNIT: &'static str = "U15a";
            const FN: &'static str = "filter_transcripts";
            type In = StageStoreIn;
            type Out = StoreOut;
            fn run(input: &StageStoreIn, ctx: &FixtureCtx) -> FixResult<StoreOut> {
                let si = StageInput::load(input, ctx)?;
                let mut ms = si.store.clone();
                let preds = Predictions::default();
                let assigner = LiveAssigner { params: AssignParams::new(si.params.matching.clone(), si.params.splice.clone(), false) };
                let c = Ctx {
                    params: &si.params,
                    index: &si.index,
                    store: &si.block.store,
                    graph: &si.graph,
                    paths: &si.paths,
                    predictions: &preds,
                    assigner: &assigner,
                };
                filter(&mut ms, &c, &Ends::default());
                Ok(StoreOut { store: si.store_out(&ms, &si.detected_known) })
            }
            fn compare(got: &StoreOut, want: &StoreOut, ctx: &FixtureCtx) -> Verdict {
                compare_stores(&got.store, &want.store, ctx)
            }
        }

        let root = fixtures_root().unwrap_or_else(|| std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../fixtures"));
        let mut reg = Registry::new();
        reg.add::<LiveFilter>();
        let mode = if iq_types::QUIRKS41 { RunMode::Parity } else { RunMode::Default };
        let mut summary = Summary::new(3);
        run_path(&reg, "U15a", &root, mode, &mut summary).unwrap();
        print!("{summary}");
        let t = &summary.by_key[&("U15a".to_owned(), "filter_transcripts".to_owned())];
        assert!(t.exact > 0);
        assert_eq!(t.mismatch, 0);
    }

    /// U15b fixtures through this crate's registry (independent of the `iqrs-fixture`
    /// binary), parity and default mode.
    /// `cargo test -p iq_model --release --features quirks41 -- --ignored u15b_fixtures --nocapture`
    #[test]
    #[ignore]
    fn u15b_fixtures() {
        use iq_types::fixture::{fixtures_root, run_path, Registry, RunMode, Summary};
        let root = fixtures_root().unwrap_or_else(|| std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../fixtures"));
        let mut reg = Registry::new();
        crate::register_fixtures(&mut reg);
        let modes: &[RunMode] = if iq_types::QUIRKS41 { &[RunMode::Parity, RunMode::Default] } else { &[RunMode::Default] };
        for &mode in modes {
            let mut summary = Summary::new(3);
            run_path(&reg, "U15b", &root, mode, &mut summary).unwrap();
            println!("mode {mode:?}");
            print!("{summary}");
            assert!(summary.by_key.values().all(|t| t.mismatch == 0 && t.decode_error == 0 && t.exact > 0));
        }
    }
}
