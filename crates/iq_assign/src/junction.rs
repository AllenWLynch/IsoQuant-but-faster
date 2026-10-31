//! U6b: `JunctionComparator` (`assignment/junction_comparator.py:34-472`).
//!
//! Compares a read's introns with one isoform's introns and classifies every contradiction
//! into IsoQuant `MatchEventSubtype`s. A port of 4.1, statement by statement, with one
//! algorithmic change:
//!
//! * `are_known_introns` (`:461-472`) builds a full K-length intron profile (K = all introns
//!   of the index) per contradictory region, per candidate, per read. Here the same sweep
//!   (`OverlappingFeaturesProfileConstructor.construct_profile_for_features`,
//!   `long_read_profiles.py:110-142`) runs only over the selected read introns, starting at
//!   the first known intron the 4.1 sweep could stop at. That start is found by binary search
//!   on a prefix maximum of known-intron ends (built once per comparator, lazily). The skipped
//!   prefix only ever executes the `gene_pos += 1` branch in 4.1, so the answer is identical,
//!   including the corner case of tiny read exons, where one known intron can be consumed by
//!   the previous read intron. Cost: O(log K + region + overlapping known introns).
//!
//! Regions are `(i64, i64)` index pairs as in 4.1, with `ABSENT_POSITION` for an absent side.

use iq_types::enums::match_consts::{ABSENT_POSITION, EXTRA_LEFT_REGION, EXTRA_RIGHT_REGION};
use iq_types::params::MatchingParams;
use iq_types::smallvec::SmallVec;
use iq_types::{EventType, Events, Iv, MatchEvent, Pos};
use std::cell::OnceCell;

/// A contradictory region pair `(read_cregion, isoform_cregion)`.
type RegionPair = ((i64, i64), (i64, i64));

const ABSENT: i64 = ABSENT_POSITION;
const ABSENT_REGION: RegionPair = ((ABSENT, ABSENT), (ABSENT, ABSENT));

// ---------------------------------------------------------------------------------------
// `common.py` interval helpers (closed 1-based intervals; lengths in i64)
// ---------------------------------------------------------------------------------------

/// `equal_ranges(a, b, delta)`.
#[inline]
pub(crate) fn equal_ranges(a: Iv, b: Iv, delta: i32) -> bool {
    (a.s as i64 - b.s as i64).abs() <= delta as i64 && (a.e as i64 - b.e as i64).abs() <= delta as i64
}

/// `left_of(a, b)`: `a.e < b.s`.
#[inline]
fn left_of(a: Iv, b: Iv) -> bool {
    a.e < b.s
}

/// `contains_well_inside(bigger, smaller, delta)`.
#[inline]
fn contains_well_inside(bigger: Iv, smaller: Iv, delta: i32) -> bool {
    bigger.e as i64 >= smaller.e as i64 + delta as i64 && bigger.s as i64 <= smaller.s as i64 - delta as i64
}

/// `contains_approx(bigger, smaller, delta)`.
#[inline]
fn contains_approx(bigger: Iv, smaller: Iv, delta: i32) -> bool {
    bigger.e as i64 + delta as i64 >= smaller.e as i64 && bigger.s as i64 - delta as i64 <= smaller.s as i64
}

/// `overlaps_at_least(a, b, delta)` (`common.py:199-208`).
#[inline]
pub(crate) fn overlaps_at_least(a: Iv, b: Iv, delta: i64) -> bool {
    let ovlp1 = a.e as i64 - b.s as i64;
    let ovlp2 = b.e as i64 - a.s as i64;
    if ovlp1 < 0 || ovlp2 < 0 {
        return false;
    }
    let d = delta - 1;
    if a.e < b.e {
        ovlp1 >= d || a.s >= b.s
    } else {
        ovlp2 >= d || a.s <= b.s
    }
}

/// `get_exon(region, junctions, pos)` for a non-negative `pos <= len` (`common.py:493-503`);
/// 4.1's `pos = -1` is `junctions.len()`.
#[inline]
fn get_exon(region: Iv, j: &[Iv], pos: usize) -> Iv {
    debug_assert!(pos <= j.len());
    if pos == 0 {
        Iv::new(region.s, j[0].s - 1)
    } else if pos == j.len() {
        Iv::new(j[pos - 1].e + 1, region.e)
    } else {
        Iv::new(j[pos - 1].e + 1, j[pos].s - 1)
    }
}

/// `get_following_exon_from_junctions(region, introns, pos)` (`common.py:471-477`);
/// 4.1's `pos = -1` is `introns.len() - 1`.
#[inline]
fn following_exon(region: Iv, j: &[Iv], pos: usize) -> Iv {
    let end = if pos == j.len() - 1 { region.e } else { j[pos + 1].s - 1 };
    Iv::new(j[pos].e + 1, end)
}

/// `get_preceding_exon_from_junctions(region, introns, pos)` for `pos <= len`
/// (`common.py:480-490`).
#[inline]
fn preceding_exon(region: Iv, j: &[Iv], pos: usize) -> Iv {
    debug_assert!(pos <= j.len());
    let start = if pos == 0 { region.s } else { j[pos - 1].e + 1 };
    if pos == j.len() {
        Iv::new(start, region.e)
    } else {
        Iv::new(start, j[pos].s - 1)
    }
}

/// `alternative_sites[(side, known)]` (`isoform_assignment.py:406-409`).
#[inline]
fn alternative_site(left: bool, known: bool) -> EventType {
    match (left, known) {
        (true, true) => EventType::AltLeftSiteKnown,
        (true, false) => EventType::AltLeftSiteNovel,
        (false, true) => EventType::AltRightSiteKnown,
        (false, false) => EventType::AltRightSiteNovel,
    }
}

#[inline]
fn us(x: i64) -> usize {
    debug_assert!(x >= 0 && x != ABSENT);
    x as usize
}

// ---------------------------------------------------------------------------------------
// JunctionComparator
// ---------------------------------------------------------------------------------------

/// `JunctionComparator(params, intron_profile_constructor)`: matching parameters plus the
/// index's sorted known introns (`gene_info.intron_profiles.features`, compared with
/// `equal_ranges(delta = params.delta)`).
///
/// Build one per (index, params) and reuse it across reads; the lookup table behind
/// [`Self::are_known_introns`] is built on first use.
#[derive(Debug)]
pub struct JunctionComparator<'a> {
    p: &'a MatchingParams,
    known: &'a [Iv],
    /// `prefix_max_end[k] = max(known[0..=k].e)`.
    prefix_max_end: OnceCell<Vec<Pos>>,
}

impl<'a> JunctionComparator<'a> {
    /// `known_introns` must be the index's intron feature list (`IsoformIndex::introns`).
    pub fn new(params: &'a MatchingParams, known_introns: &'a [Iv]) -> Self {
        JunctionComparator { p: params, known: known_introns, prefix_max_end: OnceCell::new() }
    }

    pub fn params(&self) -> &MatchingParams {
        self.p
    }

    /// `compare_junctions` (`:39-159`). Appends the detected events to `events` (callers
    /// normally pass an empty list). Appends at least one event: `[none]` when nothing
    /// contradicts (and the read is spliced).
    pub fn compare_junctions(&self, read_j: &[Iv], read_region: Iv, iso_j: &[Iv], iso_region: Iv, events: &mut Events) {
        if read_j.is_empty() {
            self.mono_exon_subtype(read_region, iso_j, events);
            return;
        }
        let delta = self.p.delta;
        let mut read_pos = 0usize;
        let mut iso_pos = 0usize;
        let mut read_present: SmallVec<[i8; 48]> = SmallVec::from_elem(0, read_j.len());
        let mut iso_present: SmallVec<[i8; 48]> = SmallVec::from_elem(0, iso_j.len());
        let mut pairs: SmallVec<[RegionPair; 8]> = SmallVec::new();
        let mut cur = ABSENT_REGION;

        while read_pos < read_j.len() && iso_pos < iso_j.len() {
            let (ij, rj) = (iso_j[iso_pos], read_j[read_pos]);
            if equal_ranges(ij, rj, delta) {
                read_present[read_pos] = 1;
                iso_present[iso_pos] = 1;
                if cur != ABSENT_REGION {
                    pairs.push(cur);
                    cur = ABSENT_REGION;
                }
                read_pos += 1;
                iso_pos += 1;
            } else if ij.overlaps(&rj) {
                read_present[read_pos] = -1;
                iso_present[iso_pos] = -1;
                if cur == ABSENT_REGION {
                    cur = ((read_pos as i64, read_pos as i64), (iso_pos as i64, iso_pos as i64));
                } else {
                    cur = ((cur.0 .0, read_pos as i64), (cur.1 .0, iso_pos as i64));
                }
                if rj.e < ij.e {
                    read_pos += 1;
                } else {
                    iso_pos += 1;
                }
            } else if left_of(ij, rj) {
                if cur != ABSENT_REGION {
                    pairs.push(cur);
                    cur = ABSENT_REGION;
                }
                if read_pos > 0 || read_region.overlaps(&ij) {
                    if iso_present[iso_pos] != -1 {
                        pairs.push(((ABSENT, read_pos as i64), (iso_pos as i64, iso_pos as i64)));
                    }
                    iso_present[iso_pos] = -1;
                }
                iso_pos += 1;
            } else {
                if cur != ABSENT_REGION {
                    pairs.push(cur);
                    cur = ABSENT_REGION;
                }
                if iso_pos > 0 || iso_region.overlaps(&rj) {
                    if read_present[read_pos] != -1 {
                        pairs.push(((read_pos as i64, read_pos as i64), (ABSENT, iso_pos as i64)));
                    }
                    read_present[read_pos] = -1;
                }
                read_pos += 1;
            }
        }
        if cur != ABSENT_REGION {
            pairs.push(cur);
        }

        // terminating regions
        while read_pos < read_j.len() {
            if iso_region.overlaps(&read_j[read_pos]) {
                if read_present[read_pos] != -1 {
                    pairs.push(((read_pos as i64, read_pos as i64), (ABSENT, iso_pos as i64)));
                    read_present[read_pos] = -1;
                }
            } else {
                break;
            }
            read_pos += 1;
        }
        while iso_pos < iso_j.len() {
            if read_region.overlaps(&iso_j[iso_pos]) {
                if iso_present[iso_pos] != -1 {
                    pairs.push(((ABSENT, read_pos as i64), (iso_pos as i64, iso_pos as i64)));
                    iso_present[iso_pos] = -1;
                }
            } else {
                break;
            }
            iso_pos += 1;
        }

        let start = events.len();
        if read_present.contains(&-1) || iso_present.contains(&-1) {
            // detect_contradiction_type (:161-181)
            for &(rc, ic) in &pairs {
                if let Some(ev) = self.compare_overlapping_contradictional_regions(read_region, read_j, iso_region, iso_j, rc, ic) {
                    events.push(ev);
                }
            }
        }
        if read_present[0] == 0 || read_present[read_present.len() - 1] == 0 {
            self.add_extra_out_exon_events(events, &read_present, read_region, read_j, iso_region.s);
        }
        if events.len() == start {
            events.push(MatchEvent::new(EventType::None));
        }
    }

    /// `compare_overlapping_contradictional_regions` (`:183-322`).
    fn compare_overlapping_contradictional_regions(
        &self,
        read_region: Iv,
        read_j: &[Iv],
        iso_region: Iv,
        iso_j: &[Iv],
        rc: (i64, i64),
        ic: (i64, i64),
    ) -> Option<MatchEvent> {
        let p = self.p;
        let ev = |t: EventType| Some(MatchEvent::with_regions(t, ic, rc));
        if rc.0 == ABSENT {
            // (4.1 logs a warning when ic.0 != ic.1; only ic.0 is used)
            let intron = iso_j[us(ic.0)];
            if read_region.contains(&intron) {
                if intron.len() <= p.micro_intron_length as i64
                    && contains_well_inside(preceding_exon(read_region, read_j, us(rc.1)), intron, p.minimal_exon_overlap)
                {
                    return ev(EventType::FakeMicroIntronRetention);
                }
                return ev(EventType::IntronRetention);
            } else if overlaps_at_least(read_region, intron, p.minor_exon_extension as i64) {
                if intron.s <= read_region.s {
                    return ev(EventType::IncompleteIntronRetentionLeft);
                }
                return ev(EventType::IncompleteIntronRetentionRight);
            }
            return None;
        } else if ic.0 == ABSENT {
            let n = read_j.len();
            if self.are_known_introns(read_j, us(rc.0), us(rc.1)) {
                return ev(EventType::ExtraIntronKnown);
            } else if self.are_suspicious_introns(read_region, read_j, us(rc.0), us(rc.1)) {
                return ev(EventType::None);
            } else if rc.0 == 0 && get_exon(read_region, read_j, 0).len() <= p.max_fake_terminal_exon_len as i64 {
                return Some(MatchEvent::with_regions(EventType::FakeTerminalExonLeft, EXTRA_LEFT_REGION, rc));
            } else if rc.1 == n as i64 - 1 && get_exon(read_region, read_j, n).len() <= p.max_fake_terminal_exon_len as i64 {
                return Some(MatchEvent::with_regions(EventType::FakeTerminalExonRight, EXTRA_RIGHT_REGION, rc));
            }
            return ev(EventType::ExtraIntronNovel);
        }

        let (r0, r1, i0, i1) = (us(rc.0), us(rc.1), us(ic.0), us(ic.1));
        let read_intron_total_len: i64 = read_j[r0..=r1].iter().map(Iv::len).sum();
        let iso_intron_total_len: i64 = iso_j[i0..=i1].iter().map(Iv::len).sum();
        let total_intron_len_diff = (read_intron_total_len - iso_intron_total_len).abs();
        // `diff <= min(abs_diff, rel_diff * max(...))` (int vs float compare in Python)
        let intron_length_is_similar = (total_intron_len_diff as f64)
            <= (p.max_intron_abs_diff as f64)
                .min(p.max_intron_rel_diff * (read_intron_total_len.max(iso_intron_total_len) as f64));

        let read_introns_known = self.are_known_introns(read_j, r0, r1);

        let left_exons_overlap = get_exon(read_region, read_j, r0).overlaps(&get_exon(iso_region, iso_j, i0));
        let right_exons_overlap = get_exon(read_region, read_j, r1 + 1).overlaps(&get_exon(iso_region, iso_j, i1 + 1));
        let surrounded_by_exons = left_exons_overlap && right_exons_overlap;
        let read_left_site = read_j[r0].s as i64;
        let read_right_site = read_j[r1].e as i64;
        let iso_left_site = iso_j[i0].s as i64;
        let iso_right_site = iso_j[i1].e as i64;
        let delta = p.delta as i64;
        let similar_left_site = (read_left_site - iso_left_site).abs() <= 2 * delta;
        let similar_right_site = (read_right_site - iso_right_site).abs() <= 2 * delta;
        let similar_bounds = similar_left_site && similar_right_site;
        let read_iv = Iv::new(read_left_site as Pos, read_right_site as Pos);
        let iso_iv = Iv::new(iso_left_site as Pos, iso_right_site as Pos);
        let read_introns_inside = contains_approx(iso_iv, read_iv, p.delta);
        let iso_introns_inside = contains_approx(read_iv, iso_iv, p.delta);

        let single_read = r1 == r0;
        let single_iso = i1 == i0;
        let mut event: Option<EventType> = None;
        if surrounded_by_exons && single_read && single_iso {
            // single intron is altered
            event = Some(self.classify_single_intron_alternation(
                read_region,
                read_j,
                iso_region,
                iso_j,
                r0,
                i0,
                intron_length_is_similar,
                read_introns_known,
            ));
        } else if read_j.len() > 1
            && single_read
            && single_iso
            && ((r0 == 0 && i0 == 0 && similar_right_site)
                || (r0 == read_j.len() - 1 && i0 == iso_j.len() - 1 && similar_left_site))
        {
            // terminal exon alternation
            let (read_exon, iso_exon) = if r0 == 0 && i0 == 0 {
                (preceding_exon(read_region, read_j, 0), preceding_exon(iso_region, iso_j, 0))
            } else {
                (following_exon(read_region, read_j, read_j.len() - 1), following_exon(iso_region, iso_j, iso_j.len() - 1))
            };
            event = Some(if (read_exon.len() - iso_exon.len()).abs() < 2 * delta {
                if r0 == 0 {
                    EventType::TerminalExonMisalignmentLeft
                } else {
                    EventType::TerminalExonMisalignmentRight
                }
            } else if read_introns_known {
                EventType::TerminalExonShiftKnown
            } else {
                EventType::TerminalExonShiftNovel
            });
        } else if surrounded_by_exons
            && similar_bounds
            && r1 - r0 == i1 - i0
            && r1 - r0 >= 1
            && total_intron_len_diff <= 2 * delta
        {
            // several introns of the same total length as the reference ones
            event = Some(if read_introns_known {
                EventType::MutuallyExclusiveExonsKnown
            } else {
                EventType::MutuallyExclusiveExonsNovel
            });
        } else if surrounded_by_exons && read_introns_inside && single_read && i1 > i0 {
            // skipped exon(s)
            event = self.classify_skipped_exons(iso_j, i0, i1, intron_length_is_similar, read_introns_known, similar_bounds);
        } else if surrounded_by_exons && similar_bounds && r1 > r0 && single_iso {
            // gain exon
            event = Some(if read_introns_known {
                EventType::ExonGainKnown
            } else if self.are_suspicious_introns(read_region, read_j, r0, r1) {
                EventType::IntronRetention
            } else {
                EventType::ExonGainNovel
            });
        } else if surrounded_by_exons && intron_length_is_similar && iso_introns_inside && r1 > r0 && single_iso {
            // exon detach
            event = Some(if read_introns_known {
                EventType::ExonDetachKnown
            } else {
                EventType::ExonDetachNovel
            });
        }

        let event = match event {
            Some(e) => e,
            // none of the above: complex alternative structure
            None => {
                if read_introns_known {
                    EventType::AlternativeStructureKnown
                } else if surrounded_by_exons && self.are_suspicious_introns(read_region, read_j, r0, r1) {
                    EventType::IntronRetention
                } else {
                    EventType::AlternativeStructureNovel
                }
            }
        };
        ev(event)
    }

    /// `classify_skipped_exons` (`:324-341`).
    fn classify_skipped_exons(
        &self,
        iso_j: &[Iv],
        i0: usize,
        i1: usize,
        intron_length_is_similar: bool,
        read_introns_known: bool,
        similar_bounds: bool,
    ) -> Option<EventType> {
        let total_exon_len: i64 = (i0..i1).map(|i| iso_j[i + 1].s as i64 - iso_j[i].e as i64 + 1).sum();
        if intron_length_is_similar {
            Some(if total_exon_len <= self.p.max_missed_exon_len as i64 {
                EventType::ExonMisalignment
            } else if read_introns_known {
                EventType::ExonMergeKnown
            } else {
                EventType::ExonMergeNovel
            })
        } else if similar_bounds {
            Some(if read_introns_known { EventType::ExonSkippingKnown } else { EventType::ExonSkippingNovel })
        } else {
            None
        }
    }

    /// `classify_single_intron_alternation` (`:343-379`).
    fn classify_single_intron_alternation(
        &self,
        read_region: Iv,
        read_j: &[Iv],
        iso_region: Iv,
        iso_j: &[Iv],
        rp: usize,
        ip: usize,
        intron_length_is_similar: bool,
        read_introns_known: bool,
    ) -> EventType {
        let p = self.p;
        let left_shift = (iso_j[ip].s as i64 - read_j[rp].s as i64).abs();
        if intron_length_is_similar {
            if left_shift <= p.max_intron_shift as i64 {
                EventType::IntronShift
            } else if read_introns_known {
                EventType::IntronMigration
            } else {
                EventType::IntronAlternationNovel
            }
        } else {
            let mut event =
                if read_introns_known { EventType::IntronAlternationKnown } else { EventType::IntronAlternationNovel };
            if left_shift <= p.delta as i64 {
                // left splice sites are similar, check for exon overlap
                let read_exon = following_exon(read_region, read_j, rp);
                let iso_exon = following_exon(iso_region, iso_j, ip);
                if overlaps_at_least(read_exon, iso_exon, self.min_exon_overlap(iso_exon)) {
                    event = alternative_site(false, read_introns_known);
                }
            } else if (iso_j[ip].e as i64 - read_j[rp].e as i64).abs() <= p.delta as i64 {
                // right splice sites are similar, check for exon overlap
                let read_exon = preceding_exon(read_region, read_j, rp);
                let iso_exon = preceding_exon(iso_region, iso_j, ip);
                if overlaps_at_least(read_exon, iso_exon, self.min_exon_overlap(iso_exon)) {
                    event = alternative_site(true, read_introns_known);
                }
            }
            if matches!(
                event,
                EventType::IntronAlternationNovel | EventType::AltLeftSiteNovel | EventType::AltRightSiteNovel
            ) && self.are_suspicious_introns(read_region, read_j, rp, rp)
            {
                event = EventType::IntronRetention;
            }
            event
        }
    }

    /// `max(1, min(min_abs_exon_overlap, round(min_rel_exon_overlap * len(exon))))`, with
    /// Python's round-half-to-even.
    #[inline]
    fn min_exon_overlap(&self, iso_exon: Iv) -> i64 {
        let rel = (self.p.min_rel_exon_overlap * iso_exon.len() as f64).round_ties_even() as i64;
        1.max((self.p.min_abs_exon_overlap as i64).min(rel))
    }

    /// `are_suspicious_introns(read_region, read_junctions, (lo, hi))` (`:381-393`).
    pub fn are_suspicious_introns(&self, read_region: Iv, read_j: &[Iv], lo: usize, hi: usize) -> bool {
        let p = self.p;
        let mut total_intron_len = 0i64;
        for iv in &read_j[lo..=hi] {
            if iv.len() > p.max_suspicious_intron_abs_len as i64 {
                return false;
            }
            total_intron_len += iv.len();
        }
        let mut total_exon_len = 0i64;
        for c in lo..=hi {
            total_exon_len += preceding_exon(read_region, read_j, c).len();
        }
        total_exon_len += following_exon(read_region, read_j, hi).len();
        (total_intron_len as f64) <= (total_exon_len as f64) * p.max_suspicious_intron_rel_len
    }

    /// `get_mono_exon_subtype` (`:395-425`).
    fn mono_exon_subtype(&self, read_region: Iv, iso_j: &[Iv], events: &mut Events)
    {
        let p = self.p;
        let start = events.len();
        let read_reg = (ABSENT, 0);
        if iso_j.is_empty() {
            events.push(MatchEvent::new(EventType::MonoExonMatch));
            return;
        }
        for (i, &intron) in iso_j.iter().enumerate() {
            let ir = (i as i64, i as i64);
            if read_region.contains(&intron) {
                let t = if intron.len() <= p.micro_intron_length as i64
                    && contains_well_inside(read_region, intron, p.minimal_exon_overlap)
                {
                    EventType::FakeMicroIntronRetention
                } else {
                    EventType::UnsplicedIntronRetention
                };
                events.push(MatchEvent::with_regions(t, ir, read_reg));
            } else if overlaps_at_least(read_region, intron, p.minor_exon_extension as i64) && !intron.contains(&read_region) {
                let t = if intron.s <= read_region.s {
                    EventType::IncompleteIntronRetentionLeft
                } else {
                    EventType::IncompleteIntronRetentionRight
                };
                events.push(MatchEvent::with_regions(t, ir, read_reg));
            }
        }
        if events.len() == start {
            events.push(MatchEvent::new(EventType::MonoExonic));
        }
    }

    /// `add_extra_out_exon_events` (`:427-459`).
    fn add_extra_out_exon_events(
        &self,
        events: &mut Events,
        prof: &[i8],
        read_region: Iv,
        read_j: &[Iv],
        iso_start: Pos,
    )
    {
        let max_fake = self.p.max_fake_terminal_exon_len as i64;
        let mut extra_left = prof[0] == 0;
        let mut extra_right = prof[prof.len() - 1] == 0;
        if prof.iter().all(|&x| x == 0) {
            if read_j[0].s < iso_start {
                extra_right = false;
            } else {
                extra_left = false;
            }
        }
        if extra_left {
            let mut rp = 0usize;
            if get_exon(read_region, read_j, rp).len() <= max_fake {
                events.push(MatchEvent::with_regions(EventType::FakeTerminalExonLeft, EXTRA_LEFT_REGION, (0, 0)));
                rp += 1;
            }
            while rp < prof.len() && prof[rp] == 0 {
                events.push(MatchEvent::with_regions(
                    EventType::ExtraIntronFlankingLeft,
                    EXTRA_LEFT_REGION,
                    (rp as i64, rp as i64),
                ));
                rp += 1;
            }
        }
        if extra_right {
            let mut rp = prof.len() as i64 - 1;
            if get_exon(read_region, read_j, prof.len()).len() <= max_fake {
                events.push(MatchEvent::with_regions(EventType::FakeTerminalExonRight, EXTRA_RIGHT_REGION, (rp, rp)));
                rp -= 1;
            }
            while rp >= 0 && prof[rp as usize] == 0 {
                events.push(MatchEvent::with_regions(EventType::ExtraIntronFlankingRight, EXTRA_RIGHT_REGION, (rp, rp)));
                rp -= 1;
            }
        }
    }

    /// `are_known_introns(junctions, (lo, hi))` (`:470-472`): true iff every read intron
    /// `read_j[lo..=hi]` gets `read_profile == 1` in 4.1's
    /// `construct_profile_for_features(read_j[lo..=hi])` against the known introns.
    pub fn are_known_introns(&self, read_j: &[Iv], lo: usize, hi: usize) -> bool {
        let sel = &read_j[lo..=hi];
        let known = self.known;
        let delta = self.p.delta;
        let pm = self.prefix_max_end.get_or_init(|| {
            let mut m = Pos::MIN;
            known
                .iter()
                .map(|f| {
                    m = m.max(f.e);
                    m
                })
                .collect()
        });
        // 4.1 starts at gene_pos = 0; every known intron before the first one with
        // e >= sel[0].s only takes the `gene_pos += 1` branch.
        let mut g = pm.partition_point(|&e| e < sel[0].s);
        let mut r = 0usize;
        let mut matched = false; // read_profile[r] == 1
        while g < known.len() && r < sel.len() {
            let (rf, kf) = (sel[r], known[g]);
            if rf.e < kf.s {
                if !matched {
                    return false;
                }
                r += 1;
                matched = false;
            } else if kf.e < rf.s {
                g += 1;
            } else {
                if equal_ranges(rf, kf, delta) {
                    matched = true;
                }
                // otherwise they overlap without being equal
                g += 1;
            }
        }
        // known introns exhausted: read intron r keeps its value, later ones are not 1
        r == sel.len() || (r == sel.len() - 1 && matched)
    }
}

/// Reference implementation of 4.1 `are_known_introns`: the full profile sweep over every
/// known intron (`long_read_profiles.py:110-142`, read part only). Test oracle.
#[cfg(test)]
pub(crate) fn are_known_introns_reference(known: &[Iv], delta: i32, sel: &[Iv]) -> bool {
    let mut read_profile = vec![0i8; sel.len()];
    let (mut g, mut r) = (0usize, 0usize);
    while g < known.len() && r < sel.len() {
        let (rf, kf) = (sel[r], known[g]);
        if rf.e < kf.s {
            if read_profile[r] == 0 && g > 0 {
                read_profile[r] = -1;
            }
            r += 1;
        } else if kf.e < rf.s {
            g += 1;
        } else if equal_ranges(rf, kf, delta) {
            read_profile[r] = 1;
            g += 1;
        } else {
            g += 1;
        }
    }
    read_profile.iter().all(|&v| v == 1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assigner::AssignParams;
    use iq_types::fixture::schema::{U6AssignIn, U6AssignOut};
    use iq_types::fixture::{read_cases, FixtureCtx};
    use iq_types::enums::match_consts::UNDEFINED_REGION;
    use std::path::PathBuf;

    fn iv(s: Pos, e: Pos) -> Iv {
        Iv::new(s, e)
    }

    fn params() -> MatchingParams {
        MatchingParams::precise()
    }

    fn cmp(jc: &JunctionComparator, rj: &[Iv], rr: Iv, ij: &[Iv], ir: Iv) -> Vec<(EventType, (i64, i64), (i64, i64))> {
        let mut ev = Events::new();
        jc.compare_junctions(rj, rr, ij, ir, &mut ev);
        ev.iter().map(|e| (e.ty, e.iso_region, e.read_region)).collect()
    }

    /// Tiny deterministic generator (no proptest in this crate).
    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self, n: u32) -> i32 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((self.0 >> 33) % n as u64) as i32
        }
    }

    #[test]
    fn known_introns_match_reference_sweep() {
        let mut rng = Lcg(7);
        for _ in 0..20000 {
            // known introns: sorted by (s, e), overlapping, some tiny
            let nk = rng.next(12) as usize;
            let mut known: Vec<Iv> = (0..nk)
                .map(|_| {
                    let s = 100 + rng.next(400);
                    iv(s, s + { let b = if rng.next(3) == 0 { 8 } else { 200 }; rng.next(b) })
                })
                .collect();
            known.sort();
            known.dedup();
            // read introns: disjoint, sorted, small gaps (tiny exons) allowed
            let nr = 1 + rng.next(5) as usize;
            let mut pos = 90 + rng.next(100);
            let mut read = Vec::new();
            for _ in 0..nr {
                let len = { let b = if rng.next(3) == 0 { 6 } else { 150 }; rng.next(b) };
                read.push(iv(pos, pos + len));
                pos += len + 2 + { let b = if rng.next(2) == 0 { 4 } else { 80 }; rng.next(b) };
            }
            // half the time, copy known introns (with jitter) into the read
            if !known.is_empty() && rng.next(2) == 0 {
                for r in read.iter_mut() {
                    if rng.next(2) == 0 {
                        let k = known[rng.next(known.len() as u32) as usize];
                        *r = iv(k.s + rng.next(5) - 2, k.e + rng.next(5) - 2);
                    }
                }
                read.sort();
            }
            let mut p = params();
            p.delta = rng.next(5);
            let jc = JunctionComparator::new(&p, &known);
            for lo in 0..read.len() {
                for hi in lo..read.len() {
                    assert_eq!(
                        jc.are_known_introns(&read, lo, hi),
                        are_known_introns_reference(&known, p.delta, &read[lo..=hi]),
                        "known {known:?} read {read:?} [{lo},{hi}] delta {}",
                        p.delta
                    );
                }
            }
        }
    }

    #[test]
    fn known_introns_tiny_exon_consumes_known() {
        // read intron 0 equals known 0 within delta; known 0 is consumed, and read
        // intron 1 (after a 1-base exon) also equals known 0 but nothing else
        let known = [iv(100, 200)];
        let read = [iv(100, 199), iv(201, 202)];
        let mut p = params();
        p.delta = 3;
        let jc = JunctionComparator::new(&p, &known);
        assert!(jc.are_known_introns(&read, 0, 0));
        assert!(!jc.are_known_introns(&read, 0, 1));
        assert!(!jc.are_known_introns(&read, 1, 1));
        // starting the window late: a long early known intron must not be skipped
        let known = [iv(10, 1000), iv(500, 600)];
        let read = [iv(12, 998)];
        let jc = JunctionComparator::new(&p, &known);
        assert!(jc.are_known_introns(&read, 0, 0));
    }

    #[test]
    fn round_half_even_overlap() {
        let mut p = params();
        p.min_abs_exon_overlap = 100;
        p.min_rel_exon_overlap = 0.5;
        let jc = JunctionComparator::new(&p, &[]);
        assert_eq!(jc.min_exon_overlap(iv(1, 5)), 2); // round(2.5) == 2
        assert_eq!(jc.min_exon_overlap(iv(1, 7)), 4); // round(3.5) == 4
        assert_eq!(jc.min_exon_overlap(iv(1, 1)), 1); // max(1, round(0.5) == 0)
    }

    #[test]
    fn mono_exon_cases() {
        let p = params();
        let jc = JunctionComparator::new(&p, &[]);
        assert_eq!(cmp(&jc, &[], iv(100, 200), &[], iv(50, 300)), vec![(EventType::MonoExonMatch, UNDEFINED_REGION, UNDEFINED_REGION)]);
        assert_eq!(cmp(&jc, &[], iv(100, 200), &[iv(300, 400)], iv(50, 500)), vec![(EventType::MonoExonic, UNDEFINED_REGION, UNDEFINED_REGION)]);
        let ev = cmp(&jc, &[], iv(100, 1000), &[iv(300, 400), iv(950, 1100), iv(1500, 1600)], iv(50, 2000));
        assert_eq!(
            ev,
            vec![
                (EventType::UnsplicedIntronRetention, (0, 0), (ABSENT, 0)),
                (EventType::IncompleteIntronRetentionRight, (1, 1), (ABSENT, 0)),
            ]
        );
        // micro intron well inside
        let ev = cmp(&jc, &[], iv(100, 1000), &[iv(300, 330)], iv(50, 2000));
        assert_eq!(ev, vec![(EventType::FakeMicroIntronRetention, (0, 0), (ABSENT, 0))]);
    }

    #[test]
    fn spliced_basic_events() {
        let p = params(); // delta 4
        let jc = JunctionComparator::new(&p, &[]);
        let ir = iv(1, 1000);
        // identical
        assert_eq!(cmp(&jc, &[iv(100, 200)], iv(50, 300), &[iv(100, 200)], ir)[0].0, EventType::None);
        // intron retention: read lacks iso intron 1
        let ev = cmp(&jc, &[iv(100, 200)], iv(50, 900), &[iv(100, 200), iv(400, 500)], ir);
        assert_eq!(ev, vec![(EventType::IntronRetention, (1, 1), (ABSENT, 1))]);
        // exon skipping: read intron spans iso introns 0..1, exon between is long
        let ev = cmp(&jc, &[iv(100, 500)], iv(50, 900), &[iv(100, 200), iv(400, 500)], ir);
        assert_eq!(ev, vec![(EventType::ExonSkippingNovel, (0, 1), (0, 0))]);
        // intron shift
        let ev = cmp(&jc, &[iv(110, 210)], iv(50, 900), &[iv(100, 200)], ir);
        assert_eq!(ev, vec![(EventType::IntronShift, (0, 0), (0, 0))]);
        // extra flanking intron on the left, long first exon
        let ev = cmp(&jc, &[iv(60, 80), iv(100, 200)], iv(1, 900), &[iv(100, 200)], iv(90, 1000));
        assert_eq!(ev, vec![(EventType::ExtraIntronFlankingLeft, EXTRA_LEFT_REGION, (0, 0))]);
    }

    // ----------------------------------------------------------------------------------
    // Fixture-driven: every U6 case whose expected assignment is inconsistent* comes from
    // `match_inconsistent`, so each match's events = compare_junctions (minus `none`) +
    // exon elongation + polyA verification. Strip the latter two and compare.
    // ----------------------------------------------------------------------------------

    pub(crate) fn fixture_root() -> Option<PathBuf> {
        let p = iq_types::fixture::fixtures_root()
            .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../../fixtures"));
        p.is_dir().then_some(p)
    }

    pub(crate) fn u6_files(root: &std::path::Path) -> Vec<PathBuf> {
        let mut caps = Vec::new();
        iq_types::fixture::find_capture_dirs(root, &mut caps).unwrap();
        caps.into_iter().map(|c| c.join("units/U6.jsonl.gz")).filter(|f| f.is_file()).collect()
    }

    fn is_jc_event(t: EventType, read_region: (i64, i64)) -> bool {
        use EventType::*;
        match t {
            ExonElongationLeft | ExonElongationRight | MajorExonElongationLeft | MajorExonElongationRight
            | TerminalSiteMatchLeft | TerminalSiteMatchRight | TerminalSiteMatchLeftPrecise
            | TerminalSiteMatchRightPrecise => false,
            t if (t.code() >= 200 && t.code() <= 300) => false,
            // polyA verification adds these with an undefined read region
            TerminalExonMisalignmentLeft | TerminalExonMisalignmentRight => read_region != UNDEFINED_REGION,
            None => false,
            _ => true,
        }
    }

    #[test]
    fn fixture_inconsistent_junction_events() {
        let Some(root) = fixture_root() else {
            eprintln!("fixtures not found; skipped");
            return;
        };
        let (mut n_case, mut n_match, mut bad) = (0usize, 0usize, Vec::new());
        for f in u6_files(&root) {
            let ctx = FixtureCtx::open(&f);
            for c in read_cases::<U6AssignIn, U6AssignOut>(&f).unwrap() {
                let c = c.unwrap();
                let Some(want) = &c.out.assignment else { continue };
                if !(30..=32).contains(&want.atype.code()) || c.input.quick_mode {
                    continue;
                }
                let dump = ctx.index_dump(&c.input.index).unwrap();
                let index = ctx.index(&c.input.index).unwrap();
                let mp = AssignParams::from_fixture(&c.input.params, &MatchingParams::precise(), false).matching;
                let jc = JunctionComparator::new(&mp, &index.introns);
                let prof = &c.input.profile;
                let rf = &prof.split_exon.read_features;
                let read_region = iv(rf[0].s, rf[rf.len() - 1].e);
                let read_j = &prof.intron.read_features;
                n_case += 1;
                for m in &want.matches {
                    let Some(tx) = &m.tx else { continue };
                    let k = dump.isoform_exons.iter().position(|(t, _)| t == tx).expect("tx in index");
                    let mut ev = Events::new();
                    jc.compare_junctions(read_j, read_region, index.iso_introns(k), index.iso[k].span, &mut ev);
                    let got: Vec<_> =
                        ev.iter().filter(|e| is_jc_event(e.ty, e.read_region)).map(|e| (e.ty, e.iso_region, e.read_region)).collect();
                    let exp: Vec<_> = m
                        .events
                        .iter()
                        .filter(|e| is_jc_event(e.t, e.read))
                        .map(|e| (e.t, e.iso, e.read))
                        .collect();
                    n_match += 1;
                    if got != exp && bad.len() < 20 {
                        bad.push(format!("{} {tx}\n  got  {got:?}\n  want {exp:?}", c.case));
                    }
                    if got != exp {
                        n_match = n_match.wrapping_add(1usize << 40);
                    }
                }
            }
        }
        let n_bad = n_match >> 40;
        eprintln!("JC fixture check: {n_case} inconsistent cases, {} matches, {n_bad} mismatching", n_match & ((1 << 40) - 1));
        assert!(bad.is_empty(), "{} mismatches, first:\n{}", n_bad, bad.join("\n"));
    }
}
