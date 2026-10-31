//! U6: `LongReadAssigner` (`assignment/long_read_assigner.py:58-830`).
//!
//! Layout:
//! - [`assign`] is the entry (`assign_to_isoform`, `:438-507`). It builds an [`Assigner`]
//!   (index, params, polyA info, scratch) and dispatches.
//! - This file holds the **consistent** path (`match_consistent*`, `:509-618`), the
//!   classification helpers shared with the inconsistent path (nucleotide scores, ambiguity
//!   resolution, exon elongation, splice-match categorisation, `classify_assignment`,
//!   `assign_to_overlapping_genes`) and the U6 fixture runner.
//! - `inconsistent.rs` (U6b) implements `Assigner::match_inconsistent` (`:620-830`).
//!
//! Isoforms are referred to by their position in `IsoformIndex::iso` (`u32`). Python sets
//! and `sorted(isoform_ids)` become `iso`-ordered candidate lists and sorts by `id_rank`
//! (design §6.4.1).

use crate::polya_verify::PolyAVerifier;
use crate::profile::Scratch;
use iq_types::fixture::schema::{AssignParamsFixture, AssignmentFixture, U6AssignIn, U6AssignOut};
use iq_types::fixture::{FixResult, FixtureCtx, UnitFixture};
use iq_types::index::prof_kind;
use iq_types::params::{MatchingParams, SpliceCorrection, MINIMAL_SCORE, TOP_SCORED_FACTOR};
use iq_types::{
    AmbiguityResolving, Assignment, AssignmentType, CombinedProfile, EventType, Events, IsoRef, IsoformIndex,
    IsoformMatch, Iv, MappedReadProfile, MatchClass, MatchEvent, Matches, PolyAInfo, Strand,
};

/// Parameters the assigner is built with (the 4.1 `params` object; schema §2.3
/// `AssignParams`) plus `quick_mode`.
#[derive(Clone, Debug, PartialEq)]
pub struct AssignParams {
    pub matching: MatchingParams,
    pub splice: SpliceCorrection,
    pub quick_mode: bool,
}

impl AssignParams {
    pub fn new(matching: MatchingParams, splice: SpliceCorrection, quick_mode: bool) -> Self {
        AssignParams { matching, splice, quick_mode }
    }
    /// From a fixture's inline `params` (fields not in the fixture keep preset values).
    pub fn from_fixture(p: &AssignParamsFixture, base: &MatchingParams, quick_mode: bool) -> Self {
        let mut m = base.clone();
        m.delta = p.delta;
        m.minor_exon_extension = p.minor_exon_extension;
        m.major_exon_extension = p.major_exon_extension;
        m.min_abs_exon_overlap = p.min_abs_exon_overlap;
        m.min_rel_exon_overlap = p.min_rel_exon_overlap;
        m.minimal_exon_overlap = p.minimal_exon_overlap;
        m.minimal_intron_absence_overlap = p.minimal_intron_absence_overlap;
        m.micro_intron_length = p.micro_intron_length;
        m.max_fake_terminal_exon_len = p.max_fake_terminal_exon_len;
        m.max_intron_shift = p.max_intron_shift;
        m.max_missed_exon_len = p.max_missed_exon_len;
        m.max_intron_abs_diff = p.max_intron_abs_diff;
        m.max_intron_rel_diff = p.max_intron_rel_diff;
        m.max_suspicious_intron_abs_len = p.max_suspicious_intron_abs_len;
        m.max_suspicious_intron_rel_len = p.max_suspicious_intron_rel_len;
        m.resolve_ambiguous = p.resolve_ambiguous;
        m.apa_delta = p.apa_delta;
        m.count_exons = p.count_exons;
        let splice = SpliceCorrection {
            fuzzy_junctions: p.correct_fuzzy_junctions,
            intron_shifts: p.correct_intron_shifts,
            skipped_exons: p.correct_skipped_exons,
            terminal_exons: p.correct_terminal_exons,
            fake_terminal_exons: p.correct_fake_terminal_exons,
            microintron_retention: p.correct_microintron_retention,
        };
        AssignParams { matching: m, splice, quick_mode }
    }
}

/// Reusable buffers of the assigner (candidate isoform lists, scores). Owned by the caller
/// (one per thread) so steady-state assignment does not allocate for candidate sets.
#[derive(Debug, Default)]
pub struct AssignScratch {
    /// Candidate lists for the consistent path (containing / overlapping / consistent /
    /// matched).
    pub(crate) c0: Vec<u32>,
    pub(crate) c1: Vec<u32>,
    pub(crate) c2: Vec<u32>,
    pub(crate) c3: Vec<u32>,
    pub(crate) c4: Vec<u32>,
    /// Free buffers for the inconsistent path (U6b); never touched by the consistent path.
    pub(crate) inc0: Vec<u32>,
    pub(crate) inc1: Vec<u32>,
    pub(crate) inc2: Vec<u32>,
    /// `(iso, score)` scratch of [`Assigner::resolve_by_nucleotide_score`].
    pub(crate) scores: Vec<(u32, f64)>,
}

/// Nucleotide similarity used by [`Assigner::resolve_by_nucleotide_score`].
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub(crate) enum NucScore {
    /// `coverage_based_nucleotide_score` (`:292-298`).
    Coverage,
    /// `jaccard_based_nucleotide_score` (`:302-308`).
    Jaccard,
}

/// How `assign_to_isoform` routed a read (used by the fixture breakdown).
#[derive(Debug)]
pub(crate) enum Route {
    /// Empty index, or no read block matched (`:462-488`): final.
    Early(Assignment),
    /// `match_consistent` produced the final assignment.
    Consistent(Assignment),
    /// Falls through to `match_inconsistent`.
    Inconsistent,
}

/// One `LongReadAssigner` call context: everything both paths need.
pub(crate) struct Assigner<'a> {
    /// `self.gene_info`.
    pub index: &'a IsoformIndex,
    /// `self.params` (+ `quick_mode`).
    pub params: &'a AssignParams,
    /// `combined_read_profile.polya_info`.
    pub polya: PolyAInfo,
    pub scratch: &'a mut AssignScratch,
}

/// `assign_to_isoform` (`:438-507`). `None` where 4.1 returns `None` (schema U6 `out`
/// `assignment: null`; 4.1 always returns an assignment here).
pub fn assign(
    index: &IsoformIndex,
    params: &AssignParams,
    profile: &CombinedProfile,
    polya: PolyAInfo,
    scratch: &mut Scratch,
) -> Option<Assignment> {
    assign_with(index, params, profile, polya, &mut scratch.assign)
}

/// [`assign`] with an explicit [`AssignScratch`].
pub fn assign_with(
    index: &IsoformIndex,
    params: &AssignParams,
    profile: &CombinedProfile,
    polya: PolyAInfo,
    scratch: &mut AssignScratch,
) -> Option<Assignment> {
    let mut a = Assigner { index, params, polya, scratch };
    Some(match a.route(profile) {
        Route::Early(x) | Route::Consistent(x) => x,
        Route::Inconsistent => a.match_inconsistent(profile),
    })
}

// ------------------------------------------------------------------------------------------
// Free helpers (common.py ports)
// ------------------------------------------------------------------------------------------

#[inline]
fn ilen(iv: &Iv) -> i64 {
    iv.e as i64 - iv.s as i64 + 1
}

/// `jaccard_similarity` (`common.py:311-361`). The 4.1 `included` arrays only ever mark
/// the current position of each list, so one flag per list suffices.
pub(crate) fn jaccard_similarity(a: &[Iv], b: &[Iv]) -> f64 {
    let (mut union, mut inter) = (0i64, 0i64);
    let (mut p1, mut p2) = (0usize, 0usize);
    let (mut inc1, mut inc2) = (false, false);
    while p1 < a.len() && p2 < b.len() {
        let (b1, b2) = (a[p1], b[p2]);
        if b1.overlaps(&b2) {
            inter += (b1.e.min(b2.e) as i64) - (b1.s.max(b2.s) as i64) + 1;
            if !inc2 && !inc1 {
                union += (b1.e.max(b2.e) as i64) - (b1.s.min(b2.s) as i64) + 1;
            } else if inc2 {
                union += 0.max(b1.e as i64 - b2.e as i64);
            } else {
                union += 0.max(b2.e as i64 - b1.e as i64);
            }
            inc1 = true;
            inc2 = true;
            if b2.e < b1.e {
                p2 += 1;
                inc2 = false;
            } else {
                p1 += 1;
                inc1 = false;
            }
        } else if b2.e < b1.s {
            if !inc2 {
                union += ilen(&b2);
            }
            p2 += 1;
            inc2 = false;
        } else {
            if !inc1 {
                union += ilen(&b1);
            }
            p1 += 1;
            inc1 = false;
        }
    }
    if p1 < a.len() {
        if !inc1 {
            union += ilen(&a[p1]);
        }
        union += a[p1 + 1..].iter().map(ilen).sum::<i64>();
    }
    if p2 < b.len() {
        if !inc2 {
            union += ilen(&b[p2]);
        }
        union += b[p2 + 1..].iter().map(ilen).sum::<i64>();
    }
    inter as f64 / union as f64
}

/// `read_coverage_fraction` (`common.py:425-445`).
pub(crate) fn read_coverage_fraction(read: &[Iv], iso: &[Iv]) -> f64 {
    let mut inter = 0i64;
    let (mut p1, mut p2) = (0usize, 0usize);
    while p1 < read.len() && p2 < iso.len() {
        let (b1, b2) = (read[p1], iso[p2]);
        if b1.overlaps(&b2) {
            inter += (b1.e.min(b2.e) as i64) - (b1.s.max(b2.s) as i64) + 1;
            if b2.e < b1.e {
                p2 += 1;
            } else {
                p1 += 1;
            }
        } else if b2.e < b1.s {
            p2 += 1;
        } else {
            p1 += 1;
        }
    }
    let read_len: i64 = read.iter().map(ilen).sum();
    inter as f64 / read_len as f64
}

/// `extra_exon_percentage` (`common.py:448-459`), `region` given as `(s, e)` in i64.
pub(crate) fn extra_exon_percentage(region: (i64, i64), read: &[Iv]) -> f64 {
    let (mut total, mut outside) = (0i64, 0i64);
    for e in read {
        let (s, t) = (e.s as i64, e.e as i64);
        total += t - s + 1;
        if s < region.0 {
            outside += t.min(region.0 - 1) - s + 1;
        }
        if t > region.1 {
            outside += t - s.max(region.1 + 1) + 1;
        }
    }
    outside as f64 / total as f64
}

/// `MatchClassification.get_mono_exon_classification` (IA:106-125).
pub(crate) fn mono_exon_classification(events: &[MatchEvent]) -> MatchClass {
    use EventType as E;
    if events.iter().any(|e| {
        matches!(
            e.ty,
            E::AlternativePolyaSiteLeft | E::AlternativePolyaSiteRight | E::InternalPolyaLeft | E::InternalPolyaRight
        )
    }) {
        return MatchClass::NovelNotInCatalog;
    }
    if events.iter().any(|e| e.ty == E::UnsplicedIntronRetention) {
        MatchClass::NovelInCatalog
    } else if events
        .iter()
        .any(|e| matches!(e.ty, E::IncompleteIntronRetentionLeft | E::IncompleteIntronRetentionRight))
    {
        MatchClass::Genic
    } else {
        match events.first().map(|e| e.ty) {
            Some(E::FakeMicroIntronRetention) => MatchClass::IncompleteSpliceMatch,
            Some(E::MonoExonMatch) => MatchClass::MonoExonMatch,
            Some(E::MonoExonic) => MatchClass::IncompleteSpliceMatch,
            _ => MatchClass::Undefined,
        }
    }
}

/// `MatchClassification.get_inconsistency_classification` (IA:98-103).
pub(crate) fn inconsistency_classification(events: &[MatchEvent]) -> MatchClass {
    if events.iter().any(|e| e.ty.is_nnic()) {
        MatchClass::NovelNotInCatalog
    } else if events.iter().any(|e| e.ty.is_nic()) {
        MatchClass::NovelInCatalog
    } else {
        MatchClass::Undefined
    }
}

/// `classify_assignment` (`:657-681`). `is_ambiguous` = more than one selected isoform;
/// `events` = the events of all selected isoforms (a set in 4.1; only membership matters).
pub(crate) fn classify_assignment<'e>(
    is_ambiguous: bool,
    events: impl Iterator<Item = &'e MatchEvent> + Clone,
) -> AssignmentType {
    use AssignmentType as A;
    if events.clone().all(|e| e.ty.is_consistent()) {
        if is_ambiguous {
            A::Ambiguous
        } else {
            A::Unique
        }
    } else if events.clone().any(|e| e.ty.is_major_inconsistency()) {
        if is_ambiguous {
            A::InconsistentAmbiguous
        } else if events.clone().any(|e| e.ty.is_intronic_inconsistency()) {
            A::Inconsistent
        } else {
            A::InconsistentNonIntronic
        }
    } else if events.clone().any(|e| e.ty.is_minor_error()) {
        if is_ambiguous {
            A::Ambiguous
        } else {
            A::UniqueMinorDifference
        }
    } else {
        A::Noninformative
    }
}

/// `equal_profiles_in_range` (`common.py:713-721`).
#[inline]
fn equal_profiles_in_range(iso: &[i8], read: &[i8], range: (i32, i32)) -> bool {
    let (lo, hi) = (range.0.max(0) as usize, (range.1.max(0) as usize).min(read.len()));
    (lo..hi).all(|i| read[i] == 0 || iso[i] == read[i])
}

/// `has_overlapping_features` (`common.py:750-758`) over `[lo, hi)`.
#[inline]
fn has_overlapping_features(iso: &[i8], read: &[i8], lo: i64, hi: i64) -> bool {
    let lo = lo.max(0) as usize;
    let hi = (hi.max(0) as usize).min(read.len());
    (lo..hi).any(|i| iso[i] == 1 && read[i] == 1)
}

// ------------------------------------------------------------------------------------------
// Assigner: shared helpers
// ------------------------------------------------------------------------------------------

impl<'a> Assigner<'a> {
    #[inline]
    pub(crate) fn m(&self) -> &'a MatchingParams {
        &self.params.matching
    }

    /// `(read_exons[0][0], read_exons[-1][1])` of the split-exon profile.
    #[inline]
    pub(crate) fn read_region(&self, p: &CombinedProfile) -> Iv {
        let f = &p.split_exon.read_features;
        Iv::new(f[0].s, f[f.len() - 1].e)
    }

    /// `gene_info.transcript_region(iso)`.
    #[inline]
    pub(crate) fn iso_region(&self, iso: u32) -> Iv {
        self.index.iso[iso as usize].span
    }

    #[inline]
    pub(crate) fn iso_exons(&self, iso: u32) -> &'a [Iv] {
        self.index.iso_exons(iso as usize)
    }

    #[inline]
    pub(crate) fn iso_introns(&self, iso: u32) -> &'a [Iv] {
        self.index.iso_introns(iso as usize)
    }

    /// Sorts isoform positions like 4.1 `sorted(isoform_ids)` (Python string order).
    #[inline]
    pub(crate) fn sort_by_id(&self, v: &mut [u32]) {
        let r = &self.index.id_rank;
        v.sort_unstable_by_key(|&i| r[i as usize]);
    }

    /// The verifier for this index / params.
    #[inline]
    pub(crate) fn polya_verifier(&self) -> PolyAVerifier<'a> {
        PolyAVerifier::new(self.index, &self.params.matching)
    }

    /// `PolyAVerifier.verify_read_ends(combined_read_profile, iso, events)` (in place).
    pub(crate) fn verify_read_ends(&self, p: &CombinedProfile, iso: u32, events: &mut Events) {
        self.polya_verifier().verify_read_ends(&p.split_exon.read_features, self.polya, iso as usize, events);
    }

    /// `IsoformMatch(class, gene=gene_id_map[iso], transcript=iso, events, strand, penalty)`
    /// with a list of events (`none` events filtered out).
    pub(crate) fn iso_match(
        &self,
        class: MatchClass,
        iso: u32,
        events: impl IntoIterator<Item = MatchEvent>,
        penalty: f64,
    ) -> IsoformMatch {
        let r = &self.index.iso[iso as usize];
        IsoformMatch::new(class, Some(r.gene), Some(r.key), events, r.strand, penalty)
    }

    /// Same with a single (unfiltered) event object.
    pub(crate) fn iso_match_single(&self, class: MatchClass, iso: u32, event: MatchEvent) -> IsoformMatch {
        let r = &self.index.iso[iso as usize];
        IsoformMatch::new_single(class, Some(r.gene), Some(r.key), event, r.strand, 0.0)
    }

    /// Isoforms whose span overlaps `q`, appended to `out` in `iso` order (span tree query
    /// without allocation).
    pub(crate) fn isoforms_overlapping_span(&self, q: Iv, out: &mut Vec<u32>) {
        let t = &self.index.span_tree;
        let hi = t.starts.partition_point(|&s| s <= q.e);
        let lo = t.prefix_max_end[..hi].partition_point(|&e| e < q.s);
        let start = out.len();
        for &i in &t.order[lo..hi] {
            if self.index.iso[i as usize].span.e >= q.s {
                out.push(i);
            }
        }
        out[start..].sort_unstable();
    }

    /// `find_overlapping_isoforms(read_split_exon_profile, split_exon_profiles, hint)`
    /// (`:83-93`): isoforms sharing a split exon with the read. `hint = None` means all
    /// isoforms (only span-overlapping ones can qualify). Appends to `out` in hint / `iso`
    /// order.
    pub(crate) fn find_overlapping_isoforms(&self, p: &CombinedProfile, hint: Option<&[u32]>, out: &mut Vec<u32>) {
        let rp = &p.split_exon;
        let check = |iso: u32| {
            let (s, e) = self.index.prof_range[iso as usize][prof_kind::SPLIT_EXON];
            let lo = (rp.gene_profile_range.0 as i64).max(s as i64);
            let hi = (rp.gene_profile_range.1 as i64).min(e as i64);
            has_overlapping_features(self.index.split_prof.row(iso as usize), &rp.gene_profile, lo, hi)
        };
        match hint {
            Some(h) => out.extend(h.iter().copied().filter(|&i| check(i))),
            None => {
                let start = out.len();
                self.isoforms_overlapping_span(self.read_region(p), out);
                let mut k = start;
                for j in start..out.len() {
                    if check(out[j]) {
                        out[k] = out[j];
                        k += 1;
                    }
                }
                out.truncate(k);
            }
        }
    }

    /// `find_containing_isoforms` (`:95-105`): `contains_approx(span, read_region,
    /// min_abs_exon_overlap)` over all isoforms. Appends in `iso` order.
    pub(crate) fn find_containing_isoforms(&self, p: &CombinedProfile, out: &mut Vec<u32>) {
        let rr = self.read_region(p);
        let d = self.m().min_abs_exon_overlap;
        let start = out.len();
        self.isoforms_overlapping_span(Iv::new(rr.s.saturating_sub(d), rr.e.saturating_add(d)), out);
        let mut k = start;
        for j in start..out.len() {
            let sp = self.index.iso[out[j] as usize].span;
            if sp.e as i64 + d as i64 >= rr.e as i64 && sp.s as i64 - d as i64 <= rr.s as i64 {
                out[k] = out[j];
                k += 1;
            }
        }
        out.truncate(k);
    }

    /// `find_matching_isoforms(read_profile, profiles(kind), hint)` (`:140-158`):
    /// `equal_profiles_in_range` over the read's `gene_profile_range`. Appends in hint order.
    pub(crate) fn find_matching_isoforms(&self, rp: &MappedReadProfile, kind: usize, hint: &[u32], out: &mut Vec<u32>) {
        let prof = self.index.prof(kind);
        for &i in hint {
            if equal_profiles_in_range(prof.row(i as usize), &rp.gene_profile, rp.gene_profile_range) {
                out.push(i);
            }
        }
    }

    /// `coverage_based_nucleotide_score` / `jaccard_based_nucleotide_score` (`:292-308`).
    pub(crate) fn nucleotide_score(&self, kind: NucScore, read_exons: &[Iv], iso: u32) -> f64 {
        let ex = self.iso_exons(iso);
        let base = match kind {
            NucScore::Coverage => read_coverage_fraction(read_exons, ex),
            NucScore::Jaccard => jaccard_similarity(read_exons, ex),
        };
        let mee = self.m().minor_exon_extension as i64;
        let region = (ex[0].s as i64 - mee, ex[ex.len() - 1].e as i64 + mee);
        base - extra_exon_percentage(region, read_exons)
    }

    /// `resolve_by_nucleotide_score` (`:267-290`) with `min_similarity = -0.5`. Writes the
    /// kept isoforms to `out` (cleared first) sorted by id. `top_scored_factor = 0` keeps all
    /// isoforms with score ≥ -0.5; otherwise `score * factor >= best` (also when `best < 0`,
    /// as in 4.1).
    pub(crate) fn resolve_by_nucleotide_score(
        &mut self,
        p: &CombinedProfile,
        isos: &[u32],
        kind: NucScore,
        top_scored_factor: f64,
        out: &mut Vec<u32>,
    ) {
        out.clear();
        if isos.is_empty() {
            return;
        }
        let read_exons = &p.split_exon.read_features;
        let mut scores = std::mem::take(&mut self.scratch.scores);
        scores.clear();
        scores.extend(isos.iter().map(|&i| (i, self.nucleotide_score(kind, read_exons, i))));
        let best = scores.iter().map(|x| x.1).fold(f64::NEG_INFINITY, f64::max);
        let min_sim = MINIMAL_SCORE;
        if top_scored_factor == 0.0 {
            out.extend(scores.iter().filter(|x| x.1 >= min_sim).map(|x| x.0));
        } else {
            // Matches IsoQuant 4.1: with a negative best score `score * factor >= best` can drop
            // every candidate.
            out.extend(scores.iter().filter(|x| x.1 * top_scored_factor >= best && x.1 >= min_sim).map(|x| x.0));
        }
        self.scratch.scores = scores;
        self.sort_by_id(out);
    }

    /// `categorize_exon_elongation_subtype` (`:199-263`): appends the elongation /
    /// terminal-site events of the read vs isoform `iso` to `out`.
    pub(crate) fn categorize_exon_elongation_subtype(&self, p: &CombinedProfile, iso: u32, out: &mut Events) {
        let rp = &p.split_exon;
        let feats = &self.index.split_exons;
        let k = feats.len() as i64;
        let row = self.index.split_prof.row(iso as usize);
        let (pf, pl) = self.index.prof_range[iso as usize][prof_kind::SPLIT_EXON];
        let (r0, r1) = (rp.gene_profile_range.0 as i64, rp.gene_profile_range.1 as i64);
        let gp = &rp.gene_profile;

        let iso_first = pf as i64;
        let mut common_first = -1i64;
        let mut i = iso_first.max(r0);
        while i < k {
            if row[i as usize] == 1 && gp[i as usize] == 1 {
                common_first = i;
                break;
            }
            i += 1;
        }
        let iso_last = pl as i64 - 1;
        let mut common_last = -1i64;
        let mut i = iso_last.min(r1 - 1);
        while i >= 0 {
            if row[i as usize] == 1 && gp[i as usize] == 1 {
                common_last = i;
                break;
            }
            i -= 1;
        }
        if common_first == -1 || common_last == -1 {
            // Matches IsoQuant 4.1: logs "Odd case" and indexes `split_exons[-1]`.
            if k == 0 {
                return;
            }
        }
        // Python negative indexing
        let feat = |c: i64| feats[if c < 0 { (k + c) as usize } else { c as usize }];
        let m = self.m();
        let (minor, delta) = (m.minor_exon_extension as i64, m.delta as i64);

        let first_read = rp.read_features[0];
        let f = feat(common_first);
        if first_read.overlaps(&f) {
            let extra = f.s as i64 - first_read.s as i64;
            if common_first == iso_first {
                if extra.abs() <= minor {
                    let t = if extra.abs() <= delta {
                        EventType::TerminalSiteMatchLeftPrecise
                    } else {
                        EventType::TerminalSiteMatchLeft
                    };
                    out.push(MatchEvent::with_info(t, extra));
                }
                if extra > minor {
                    out.push(MatchEvent::with_info(EventType::MajorExonElongationLeft, extra));
                } else if extra > delta {
                    out.push(MatchEvent::with_info(EventType::ExonElongationLeft, extra));
                }
            } else if minor >= extra && extra > delta {
                out.push(MatchEvent::with_info(EventType::ExonElongationLeft, extra));
            }
        }

        let last_read = rp.read_features[rp.read_features.len() - 1];
        let f = feat(common_last);
        if last_read.overlaps(&f) {
            let extra = last_read.e as i64 - f.e as i64;
            if common_last == iso_last {
                if extra.abs() <= minor {
                    let t = if extra.abs() <= delta {
                        EventType::TerminalSiteMatchRightPrecise
                    } else {
                        EventType::TerminalSiteMatchRight
                    };
                    out.push(MatchEvent::with_info(t, extra));
                }
                if extra > minor {
                    out.push(MatchEvent::with_info(EventType::MajorExonElongationRight, extra));
                } else if extra > delta {
                    out.push(MatchEvent::with_info(EventType::ExonElongationRight, extra));
                }
            } else if minor >= extra && extra > delta {
                out.push(MatchEvent::with_info(EventType::ExonElongationRight, extra));
            }
        }
    }

    /// `is_fsm` (`:357-361`): all isoform introns lie within the read region.
    pub(crate) fn is_fsm(&self, read_region: Iv, iso: u32) -> bool {
        let introns = self.iso_introns(iso);
        let r = Iv::new(introns[0].s, introns[introns.len() - 1].e);
        read_region.contains(&r)
    }

    /// `detect_ism_subtype` (`:335-354`).
    pub(crate) fn detect_ism_subtype(&self, read_region: Iv, iso: u32) -> MatchEvent {
        let introns = self.iso_introns(iso);
        let (is, ie) = (introns[0].s, introns[introns.len() - 1].e);
        let left = is < read_region.s;
        let right = ie > read_region.e;
        MatchEvent::new(match (left, right) {
            (true, true) => EventType::IsmInternal,
            (true, false) => EventType::IsmLeft,
            (false, true) => EventType::IsmRight,
            (false, false) => EventType::None,
        })
    }

    /// `categorize_correct_splice_match` (`:364-385`).
    pub(crate) fn categorize_correct_splice_match(&self, p: &CombinedProfile, iso: u32) -> IsoformMatch {
        let rr = self.read_region(p);
        if p.intron.read_profile.is_empty() || self.iso_introns(iso).is_empty() {
            self.iso_match_single(MatchClass::MonoExonMatch, iso, MatchEvent::new(EventType::MonoExonMatch))
        } else if self.is_fsm(rr, iso) {
            self.iso_match_single(MatchClass::FullSpliceMatch, iso, MatchEvent::new(EventType::Fsm))
        } else {
            self.iso_match_single(MatchClass::IncompleteSpliceMatch, iso, self.detect_ism_subtype(rr, iso))
        }
    }

    /// `categorize_correct_unspliced_match` (`:395-410`).
    pub(crate) fn categorize_correct_unspliced_match(&self, iso: u32) -> IsoformMatch {
        let ev = if self.iso_exons(iso).len() == 1 {
            MatchEvent::new(EventType::MonoExonMatch)
        } else {
            MatchEvent::new(EventType::MonoExonic)
        };
        let class = mono_exon_classification(std::slice::from_ref(&ev));
        self.iso_match(class, iso, [ev], 0.0)
    }

    /// True if gene regions exist for this index: 4.1 `get_gene_regions()` iterates
    /// `gene_db_list`, which is empty for `from_models` / `from_region` indices.
    fn has_gene_regions(&self) -> bool {
        self.index.iso.first().is_none_or(|r| matches!(r.key, IsoRef::Ref(_)))
    }

    /// `assign_to_overlapping_genes` (`:420-436`).
    pub(crate) fn assign_to_overlapping_genes(&self, read_region: Iv, class: MatchClass) -> Assignment {
        let mut matches = Matches::new();
        if self.has_gene_regions() {
            for g in &self.index.genes {
                if read_region.overlaps(&g.region) {
                    matches.push(IsoformMatch::new(class, Some(g.key), None, [], g.strand, 0.0));
                }
            }
        }
        if matches.is_empty() {
            matches.push(IsoformMatch::new(class, None, None, [], Strand::Unknown, 0.0));
        }
        Assignment::new(AssignmentType::Noninformative, matches)
    }

    // --------------------------------------------------------------------------------------
    // Dispatch and consistent path
    // --------------------------------------------------------------------------------------

    /// `assign_to_isoform` (`:438-507`) up to the `match_inconsistent` fall-through.
    pub(crate) fn route(&mut self, p: &CombinedProfile) -> Route {
        if self.index.exons.is_empty() {
            let m = IsoformMatch::new(MatchClass::Intergenic, None, None, [], Strand::Unknown, 0.0);
            return Route::Early(Assignment::new(AssignmentType::Intergenic, Matches::from_elem(m, 1)));
        }
        let se = &p.split_exon;
        let ip = &p.intron;
        if se.read_profile.iter().all(|&v| v != 1) || se.gene_profile.iter().all(|&v| v == 0 || v == -2) {
            let rr = self.read_region(p);
            let f = &self.index.split_exons;
            let gene_region = Iv::new(f[0].s, f[f.len() - 1].e);
            let a = if !rr.overlaps(&gene_region) {
                let m = IsoformMatch::new(MatchClass::Intergenic, None, None, [], Strand::Unknown, 0.0);
                Assignment::new(AssignmentType::Noninformative, Matches::from_elem(m, 1))
            } else if se.gene_profile.iter().all(|&v| v != 1) {
                self.assign_to_overlapping_genes(rr, MatchClass::GenicIntron)
            } else {
                self.assign_to_overlapping_genes(rr, MatchClass::Genic)
            };
            return Route::Early(a);
        }
        if ip.read_profile.iter().any(|&v| v == -1 || v == 0) || se.read_profile.iter().any(|&v| v == -1 || v == 0) {
            return Route::Inconsistent;
        }
        match self.match_consistent(p) {
            Some(a) => Route::Consistent(a),
            None => Route::Inconsistent,
        }
    }

    /// `match_consistent` (`:509-550`).
    pub(crate) fn match_consistent(&mut self, p: &CombinedProfile) -> Option<Assignment> {
        let mut containing = std::mem::take(&mut self.scratch.c0);
        let mut overlapping = std::mem::take(&mut self.scratch.c1);
        let mut consistent = std::mem::take(&mut self.scratch.c2);
        containing.clear();
        overlapping.clear();
        consistent.clear();
        let mut result = None;
        self.find_containing_isoforms(p, &mut containing);
        if !containing.is_empty() {
            self.find_overlapping_isoforms(p, Some(&containing), &mut overlapping);
            if !overlapping.is_empty() {
                self.find_matching_isoforms(&p.intron, prof_kind::INTRON, &overlapping, &mut consistent);
                result = if p.intron.read_profile.is_empty() {
                    self.match_consistent_unspliced(p, &consistent)
                } else {
                    self.match_consistent_spliced(p, &consistent)
                };
            }
        }
        self.scratch.c0 = containing;
        self.scratch.c1 = overlapping;
        self.scratch.c2 = consistent;

        let mut a = result?;
        self.check_read_ends(p, &mut a);
        self.verify_read_ends_for_assignment(p, &mut a);
        if a.atype.is_inconsistent() {
            return None;
        }
        Some(a)
    }

    /// `match_consistent_spliced` (`:552-590`).
    fn match_consistent_spliced(&mut self, p: &CombinedProfile, consistent: &[u32]) -> Option<Assignment> {
        let mut exon_matched = std::mem::take(&mut self.scratch.c3);
        let mut resolved = std::mem::take(&mut self.scratch.c4);
        exon_matched.clear();
        resolved.clear();
        let mut use_exon = false;
        if consistent.len() > 1 {
            self.find_matching_isoforms(&p.split_exon, prof_kind::SPLIT_EXON, consistent, &mut exon_matched);
            use_exon = !exon_matched.is_empty();
        }
        let mut use_resolved = false;
        {
            let matched: &[u32] = if use_exon { &exon_matched } else { consistent };
            if matched.len() > 1 {
                let rr = self.read_region(p);
                let ra = self.m().resolve_ambiguous;
                if ra == AmbiguityResolving::All
                    || (ra == AmbiguityResolving::MonoexonAndFsm && matched.iter().any(|&i| self.is_fsm(rr, i)))
                {
                    self.resolve_by_nucleotide_score(p, matched, NucScore::Jaccard, TOP_SCORED_FACTOR, &mut resolved);
                    use_resolved = true;
                }
            }
        }
        let matched: &[u32] = if use_resolved {
            &resolved
        } else if use_exon {
            &exon_matched
        } else {
            consistent
        };
        let res = match matched.len() {
            0 => None,
            1 => {
                let m = self.categorize_correct_splice_match(p, matched[0]);
                Some(Assignment::new(AssignmentType::Unique, Matches::from_elem(m, 1)))
            }
            _ => {
                let mut ids: iq_types::smallvec::SmallVec<[u32; 8]> = matched.iter().copied().collect();
                self.sort_by_id(&mut ids);
                let ms: Matches = ids.iter().map(|&i| self.categorize_correct_splice_match(p, i)).collect();
                Some(Assignment::new(AssignmentType::Ambiguous, ms))
            }
        };
        self.scratch.c3 = exon_matched;
        self.scratch.c4 = resolved;
        res
    }

    /// `match_consistent_unspliced` (`:592-617`).
    fn match_consistent_unspliced(&mut self, p: &CombinedProfile, consistent: &[u32]) -> Option<Assignment> {
        let mut buf = std::mem::take(&mut self.scratch.c3);
        buf.clear();
        let matched: &[u32] = if consistent.len() > 1 && self.m().resolve_ambiguous != AmbiguityResolving::None {
            self.resolve_by_nucleotide_score(p, consistent, NucScore::Jaccard, TOP_SCORED_FACTOR, &mut buf);
            &buf
        } else {
            consistent
        };
        let res = match matched.len() {
            0 => None,
            1 => {
                let m = self.categorize_correct_unspliced_match(matched[0]);
                Some(Assignment::new(AssignmentType::Unique, Matches::from_elem(m, 1)))
            }
            _ => {
                let mut ids: iq_types::smallvec::SmallVec<[u32; 8]> = matched.iter().copied().collect();
                self.sort_by_id(&mut ids);
                let ms: Matches = ids.iter().map(|&i| self.categorize_correct_unspliced_match(i)).collect();
                Some(Assignment::new(AssignmentType::Ambiguous, ms))
            }
        };
        self.scratch.c3 = buf;
        res
    }

    /// `check_read_ends` (`:315-332`): only the events (the type writes are dead, D-A5).
    fn check_read_ends(&self, p: &CombinedProfile, a: &mut Assignment) {
        let mut evs = Events::new();
        for m in a.matches.iter_mut() {
            let Some(tx) = m.tx else { continue };
            let Some(iso) = self.iso_pos(tx) else { continue };
            evs.clear();
            self.categorize_exon_elongation_subtype(p, iso, &mut evs);
            for &e in &evs {
                m.add_event(e);
            }
        }
    }

    /// `verify_read_ends_for_assignment` (`:819-830`). `assignment_type` is overwritten;
    /// `gene_assignment_type` keeps its constructor value, as in 4.1.
    fn verify_read_ends_for_assignment(&self, p: &CombinedProfile, a: &mut Assignment) {
        let mut n = 0usize;
        for m in a.matches.iter_mut() {
            let Some(tx) = m.tx else { continue };
            let Some(iso) = self.iso_pos(tx) else { continue };
            self.verify_read_ends(p, iso, &mut m.events);
            n += 1;
        }
        let evs = a.matches.iter().filter(|m| m.tx.is_some()).flat_map(|m| m.events.iter());
        a.atype = classify_assignment(n > 1, evs);
    }

    /// Position of an isoform key in the index. Matches built here come from the index, so
    /// this is a search over at most the assignment's isoforms' neighbourhood; ambiguous
    /// assignments are small.
    #[inline]
    pub(crate) fn iso_pos(&self, key: IsoRef) -> Option<u32> {
        self.index.find_iso(key).map(|i| i as u32)
    }
}

// ------------------------------------------------------------------------------------------
// Fixture runner
// ------------------------------------------------------------------------------------------

/// U6 `assign_to_isoform` runner.
pub struct U6Assign;

impl U6Assign {
    fn params(i: &U6AssignIn) -> AssignParams {
        AssignParams::from_fixture(&i.params, &MatchingParams::default_ont(), i.quick_mode)
    }
}

impl UnitFixture for U6Assign {
    const UNIT: &'static str = "U6";
    const FN: &'static str = "assign_to_isoform";
    type In = U6AssignIn;
    type Out = U6AssignOut;

    fn run(i: &U6AssignIn, ctx: &FixtureCtx) -> FixResult<U6AssignOut> {
        let index = ctx.index(&i.index)?;
        let params = Self::params(i);
        let mut s = AssignScratch::default();
        let a = assign_with(&index, &params, &i.profile, i.profile.polya, &mut s);
        Ok(U6AssignOut { assignment: a.map(|a| AssignmentFixture::from_core(&a, &ctx.ids())) })
    }
}

/// Registers the U6 fixture runners.
pub fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    reg.add::<U6Assign>();
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(v: &[(i32, i32)]) -> Vec<Iv> {
        v.iter().map(|&(s, e)| Iv::new(s, e)).collect()
    }

    #[test]
    fn jaccard_matches_python() {
        // python: jaccard_similarity([(1,10),(20,30)], [(5,25)]) = 12/30
        let a = iv(&[(1, 10), (20, 30)]);
        let b = iv(&[(5, 25)]);
        assert_eq!(jaccard_similarity(&a, &b), 12.0 / 30.0);
        assert_eq!(jaccard_similarity(&b, &a), 12.0 / 30.0);
        let c = iv(&[(1, 10)]);
        assert_eq!(jaccard_similarity(&c, &c), 1.0);
        let d = iv(&[(100, 109)]);
        assert_eq!(jaccard_similarity(&c, &d), 0.0);
    }

    #[test]
    fn coverage_and_flank() {
        let r = iv(&[(1, 10), (20, 29)]);
        let i = iv(&[(5, 25)]);
        assert_eq!(read_coverage_fraction(&r, &i), 12.0 / 20.0);
        assert_eq!(extra_exon_percentage((5, 25), &r), 8.0 / 20.0);
    }

    #[test]
    fn classify() {
        let fsm = [MatchEvent::new(EventType::Fsm)];
        assert_eq!(classify_assignment(false, fsm.iter()), AssignmentType::Unique);
        let el = [MatchEvent::new(EventType::Fsm), MatchEvent::new(EventType::ExonElongationLeft)];
        assert_eq!(classify_assignment(false, el.iter()), AssignmentType::UniqueMinorDifference);
        let apa = [MatchEvent::new(EventType::AlternativePolyaSiteRight)];
        assert_eq!(classify_assignment(false, apa.iter()), AssignmentType::InconsistentNonIntronic);
        assert_eq!(classify_assignment(true, apa.iter()), AssignmentType::InconsistentAmbiguous);
    }

    /// U6 breakdown by call_site × route (early / consistent / inconsistent).
    /// `cargo test --release -p iq_assign u6_breakdown -- --ignored --nocapture`
    /// Env: `IQRS_FIXTURES` (default `../fixtures`), `IQRS_U6_CALL_SITE`, `IQRS_U6_ROUTE`
    /// (`early|consistent|inconsistent`, route as computed by *this* code),
    /// `IQRS_U6_ATYPE` (expected atype code), `IQRS_U6_EXAMPLES` (n diffs to print).
    #[test]
    #[ignore]
    fn u6_breakdown() {
        use iq_types::fixture::{find_capture_dirs, read_cases, unit_files, Verdict};
        use std::collections::BTreeMap;
        let root = std::env::var("IQRS_FIXTURES")
            .unwrap_or_else(|_| concat!(env!("CARGO_MANIFEST_DIR"), "/../../../fixtures").into());
        let f_cs = std::env::var("IQRS_U6_CALL_SITE").ok();
        let f_route = std::env::var("IQRS_U6_ROUTE").ok();
        let f_atype: Option<u8> = std::env::var("IQRS_U6_ATYPE").ok().and_then(|v| v.parse().ok());
        let n_ex: usize = std::env::var("IQRS_U6_EXAMPLES").ok().and_then(|v| v.parse().ok()).unwrap_or(5);
        std::panic::set_hook(Box::new(|_| {}));
        let mut caps = Vec::new();
        find_capture_dirs(std::path::Path::new(&root), &mut caps).unwrap();
        // (call_site, route) -> [exact, explained, mismatch]
        let mut tally: BTreeMap<(String, &str), [u64; 3]> = BTreeMap::new();
        let mut shown = 0;
        for cap in caps {
            let ctx = FixtureCtx::new(&cap);
            for f in unit_files(&cap, "U6").unwrap() {
                for c in read_cases::<U6AssignIn, U6AssignOut>(&f).unwrap() {
                    let c = c.unwrap();
                    let (i, want) = (&c.input, &c.out);
                    if f_cs.as_ref().is_some_and(|x| *x != i.call_site) {
                        continue;
                    }
                    if f_atype.is_some_and(|t| want.assignment.as_ref().map(|a| a.atype.code()) != Some(t)) {
                        continue;
                    }
                    let index = ctx.index(&i.index).unwrap();
                    let params = U6Assign::params(i);
                    let mut s = AssignScratch::default();
                    let route = {
                        let mut a =
                            Assigner { index: &index, params: &params, polya: i.profile.polya, scratch: &mut s };
                        match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| a.route(&i.profile))) {
                            Ok(Route::Early(_)) => "early",
                            Ok(Route::Consistent(_)) => "consistent",
                            Ok(Route::Inconsistent) => "inconsistent",
                            Err(_) => "panic",
                        }
                    };
                    if f_route.as_ref().is_some_and(|x| x != route) {
                        continue;
                    }
                    let v = U6Assign::judge(i, want, &ctx);
                    let t = tally.entry((i.call_site.clone(), route)).or_default();
                    match v {
                        Verdict::Exact => t[0] += 1,
                        Verdict::Explained(_) => t[1] += 1,
                        Verdict::Mismatch(d) => {
                            t[2] += 1;
                            if shown < n_ex && route != "inconsistent" {
                                shown += 1;
                                println!("--- {} [{}]\n{d}", c.case, route);
                            }
                        }
                    }
                }
            }
        }
        println!("{:<14} {:<13} {:>8} {:>9} {:>8}", "call_site", "route", "exact", "explained", "mismatch");
        for ((cs, r), t) in &tally {
            println!("{:<14} {:<13} {:>8} {:>9} {:>8}", cs, r, t[0], t[1], t[2]);
        }
    }
}
