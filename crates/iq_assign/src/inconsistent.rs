//! U6 inconsistent path: `LongReadAssigner.match_inconsistent` and helpers
//! (`assignment/long_read_assigner.py`), owned by the U6b agent. Uses [`crate::junction`].
//!
//! Python line numbers (`:N`) refer to `long_read_assigner.py`.

use crate::assigner::{
    classify_assignment, inconsistency_classification, mono_exon_classification, Assigner, NucScore,
};
use crate::junction::JunctionComparator;
use iq_types::enums::match_consts::{ABSENT_POSITION, UNDEFINED_REGION};
use iq_types::params::TOP_SCORED_FACTOR;
use iq_types::smallvec::SmallVec;
use iq_types::{
    elongation_cost, Assignment, AssignmentType, CombinedProfile, EventType, Events, IsoformMatch, MatchClass,
    MatchEvent, Matches, Strand,
};

/// `read_matches` of `detect_inconsistensies`: `(isoform, events)` in insertion order
/// (= isoforms sorted by id).
type ReadMatches = SmallVec<[(u32, Events); 4]>;

/// `difference_in_present_features(iso, read, -1, range)` (`common.py:680-710`): Hamming
/// distance over `[lo, hi)`, ignoring positions where either profile is 0.
#[inline]
fn difference_in_present_features(iso: &[i8], read: &[i8], lo: i64, hi: i64) -> u32 {
    let lo = lo.max(0) as usize;
    let hi = (hi.max(0) as usize).min(read.len());
    if lo >= hi {
        return 0;
    }
    iso[lo..hi].iter().zip(&read[lo..hi]).filter(|&(&a, &b)| a != 0 && b != 0 && a != b).count() as u32
}

/// Number of introns an event affects (`select_best_among_inconsistent`, `:767-781`).
#[inline]
fn event_count(e: &MatchEvent) -> i64 {
    use EventType as E;
    let defined = |r: (i64, i64)| r != UNDEFINED_REGION && r.0 != ABSENT_POSITION && r.1 != ABSENT_POSITION;
    if defined(e.iso_region) && matches!(e.ty, E::ExonSkippingKnown | E::ExonSkippingNovel) {
        e.iso_region.1 - e.iso_region.0 + 1
    } else if defined(e.read_region) {
        let c = e.read_region.1 - e.read_region.0 + 1;
        match e.ty {
            E::ExonGainNovel
            | E::ExonGainKnown
            | E::MutuallyExclusiveExonsNovel
            | E::MutuallyExclusiveExonsKnown
            | E::ExonDetachKnown
            | E::ExonDetachNovel => c - 1,
            E::IntronRetention
            | E::UnsplicedIntronRetention
            | E::FakeMicroIntronRetention
            | E::IncompleteIntronRetentionLeft
            | E::IncompleteIntronRetentionRight => 1,
            _ => c,
        }
    } else {
        1
    }
}

impl<'a> Assigner<'a> {
    /// `match_inconsistent(read_id, combined_read_profile)` (`long_read_assigner.py:620-655`).
    /// Called by `assign_to_isoform` when the read has `-1`/`0` features or
    /// `match_consistent` returned `None`. `read_id` is unused in 4.1 beyond logging.
    pub(crate) fn match_inconsistent(&mut self, p: &CombinedProfile) -> Assignment {
        if self.params.quick_mode {
            let m = IsoformMatch::new(MatchClass::Genic, None, None, [], Strand::Unknown, 0.0);
            return Assignment::new(AssignmentType::Inconsistent, Matches::from_elem(m, 1));
        }
        let read_region = self.read_region(p);

        let mut candidates = std::mem::take(&mut self.scratch.inc2);
        candidates.clear();
        self.select_similar_isoforms(p, &mut candidates);
        let read_matches =
            if candidates.is_empty() { ReadMatches::new() } else { self.detect_inconsistencies(p, &mut candidates) };
        self.scratch.inc2 = candidates;
        if read_matches.is_empty() {
            return self.assign_to_overlapping_genes(read_region, MatchClass::Genic);
        }

        let mut best = std::mem::take(&mut self.scratch.inc0);
        let penalty = self.select_best_among_inconsistent(p, &read_matches, &mut best);
        let result = if best.is_empty() {
            self.assign_to_overlapping_genes(read_region, MatchClass::Genic)
        } else {
            self.sort_by_id(&mut best);
            let events_of =
                |iso: u32| -> &Events { &read_matches.iter().find(|x| x.0 == iso).expect("selected isoform").1 };
            let atype = classify_assignment(best.len() > 1, best.iter().flat_map(|&i| events_of(i).iter()));
            let matches: Matches = if p.intron.read_profile.is_empty() {
                // create_monoexon_matches (:678-689)
                best.iter()
                    .map(|&i| {
                        let ev = events_of(i);
                        self.iso_match(mono_exon_classification(ev), i, ev.iter().copied(), 0.0)
                    })
                    .collect()
            } else if atype.is_inconsistent() {
                // create_inconsistent_matches (:701-714)
                best.iter()
                    .map(|&i| {
                        let ev = events_of(i);
                        self.iso_match(inconsistency_classification(ev), i, ev.iter().copied(), penalty)
                    })
                    .collect()
            } else {
                // create_consistent_matches (:691-699)
                best.iter()
                    .map(|&i| {
                        let mut m = self.categorize_correct_splice_match(p, i);
                        for &e in events_of(i) {
                            if e.ty != EventType::None {
                                m.add_event(e);
                            }
                        }
                        m
                    })
                    .collect()
            };
            Assignment::new(atype, matches)
        };
        self.scratch.inc0 = best;
        result
    }

    /// `select_similar_isoforms` (`:163-197`): isoforms sharing a split exon with the read,
    /// with coverage score ≥ -0.5, whose intron-profile Hamming distance (+0.5 per extra
    /// terminal side) is within 3 of the best. Written to `out` in `match_profile` order
    /// (by distance, ties by id); empty when 4.1 returns `None`/`[]`.
    fn select_similar_isoforms(&mut self, p: &CombinedProfile, out: &mut Vec<u32>) {
        let mut overlapping = std::mem::take(&mut self.scratch.inc0);
        let mut significant = std::mem::take(&mut self.scratch.inc1);
        overlapping.clear();
        significant.clear();
        self.find_overlapping_isoforms(p, None, &mut overlapping);
        if !overlapping.is_empty() {
            self.resolve_by_nucleotide_score(p, &overlapping, NucScore::Coverage, 0.0, &mut significant);
        }
        if !significant.is_empty() {
            // match_profile (:108-131): hint order (= sorted ids), stable sort by diff
            let ip = &p.intron;
            let read_region = self.read_region(p);
            let delta = self.m().delta;
            let (lo, hi) = (ip.gene_profile_range.0 as i64, ip.gene_profile_range.1 as i64);
            let mut cands: SmallVec<[(u32, u32, f64); 16]> = significant
                .iter()
                .map(|&iso| {
                    let diff = difference_in_present_features(
                        self.index.intron_prof.row(iso as usize),
                        &ip.gene_profile,
                        lo,
                        hi,
                    );
                    (iso, diff, 0.0)
                })
                .collect();
            cands.sort_by_key(|c| c.1);
            for c in cands.iter_mut() {
                let sp = self.iso_region(c.0);
                let extra_left = if read_region.s as i64 + (delta as i64) < sp.s as i64 { 0.5 } else { 0.0 };
                let extra_right = if read_region.e as i64 - (delta as i64) > sp.e as i64 { 0.5 } else { 0.0 };
                c.2 = c.1 as f64 + extra_right + extra_left;
            }
            let best = cands.iter().map(|c| c.2).fold(f64::INFINITY, f64::min);
            out.extend(cands.iter().filter(|c| c.2 <= best + 3.0).map(|c| c.0));
        }
        self.scratch.inc0 = overlapping;
        self.scratch.inc1 = significant;
    }

    /// `detect_inconsistensies` (`:716-762`) over `sorted(candidates)` (sorts in place).
    fn detect_inconsistencies(&mut self, p: &CombinedProfile, candidates: &mut [u32]) -> ReadMatches {
        self.sort_by_id(candidates);
        let jc = JunctionComparator::new(self.m(), &self.index.introns);
        let read_region = self.read_region(p);
        let read_introns = &p.intron.read_features;
        let mut out = ReadMatches::new();
        for &iso in candidates.iter() {
            let mut ev = Events::new();
            jc.compare_junctions(read_introns, read_region, self.iso_introns(iso), self.iso_region(iso), &mut ev);
            if ev.len() == 1 && ev[0].ty == EventType::Undefined {
                continue;
            }
            self.categorize_exon_elongation_subtype(p, iso, &mut ev);
            self.verify_read_ends(p, iso, &mut ev);
            out.push((iso, ev));
        }
        out
    }

    /// `select_best_among_inconsistent` (`:764-817`): penalty = Σ cost × count in event
    /// order (f64, same summation order as 4.1); keeps the isoforms with the minimal penalty
    /// (exact `==`), resolving ties by coverage score. Writes them to `best` (cleared
    /// first) and returns the minimal penalty.
    fn select_best_among_inconsistent(
        &mut self,
        p: &CombinedProfile,
        read_matches: &ReadMatches,
        best: &mut Vec<u32>,
    ) -> f64 {
        let m = self.m();
        let (minor, major) = (m.minor_exon_extension as i64, m.major_exon_extension as i64);
        let mut scores: SmallVec<[(u32, f64); 8]> = SmallVec::new();
        for (iso, events) in read_matches {
            let mut penalty = 0.0f64;
            for e in events {
                let count = event_count(e);
                let cost = match e.ty {
                    EventType::MajorExonElongationLeft
                    | EventType::MajorExonElongationRight
                    | EventType::ExonElongationRight
                    | EventType::ExonElongationLeft => elongation_cost(minor, major, e.info),
                    t => t.cost().expect("event without a cost (antisense) in inconsistent path"),
                };
                penalty += cost * count as f64;
            }
            scores.push((*iso, penalty));
        }
        // min() keeps the first minimum; NaN cannot occur
        let min_penalty = scores.iter().map(|s| s.1).fold(f64::INFINITY, f64::min);
        best.clear();
        best.extend(scores.iter().filter(|s| s.1 == min_penalty).map(|s| s.0));
        if best.len() > 1 {
            let tied = std::mem::take(best);
            self.resolve_by_nucleotide_score(p, &tied, NucScore::Coverage, TOP_SCORED_FACTOR, best);
        }
        min_penalty
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assigner::{assign_with, AssignParams, AssignScratch};
    use iq_types::fixture::schema::{AssignmentFixture, U6AssignIn, U6AssignOut};
    use iq_types::fixture::{diff_values, find_capture_dirs, read_cases, unit_files, DiffOptions, FixtureCtx};
    use iq_types::params::MatchingParams;
    use std::path::PathBuf;

    #[test]
    fn hamming_ignores_zeros() {
        let iso = [1, -1, 0, 1, -2, 1];
        let read = [1, 1, 1, 0, -1, -2];
        assert_eq!(difference_in_present_features(&iso, &read, 0, 6), 3);
        assert_eq!(difference_in_present_features(&iso, &read, 6, 0), 0);
        assert_eq!(difference_in_present_features(&iso, &read, 1, 2), 1);
    }

    #[test]
    fn event_counts() {
        let skip = MatchEvent::with_regions(EventType::ExonSkippingNovel, (2, 4), (1, 1));
        assert_eq!(event_count(&skip), 3);
        let gain = MatchEvent::with_regions(EventType::ExonGainKnown, (1, 1), (1, 3));
        assert_eq!(event_count(&gain), 2);
        let ir = MatchEvent::with_regions(EventType::IntronRetention, (1, 1), (ABSENT_POSITION, 2));
        assert_eq!(event_count(&ir), 1);
        let extra = MatchEvent::with_regions(EventType::ExtraIntronNovel, (ABSENT_POSITION, 0), (0, 2));
        assert_eq!(event_count(&extra), 3);
        assert_eq!(event_count(&MatchEvent::with_info(EventType::ExonElongationLeft, 20)), 1);
    }

    /// U6 cases routed to `match_inconsistent` by 4.1 can only be told apart by running the
    /// dispatch; this harness runs every case whose expected atype is inconsistent* (or all
    /// cases with `IQRS_U6B_ALL=1`) and reports exact / mismatch with a few diffs.
    /// `cargo test --release -p iq_assign u6b_inconsistent -- --ignored --nocapture`
    #[test]
    #[ignore]
    fn u6b_inconsistent() {
        let root = std::env::var("IQRS_FIXTURES")
            .map(PathBuf::from)
            .unwrap_or_else(|_| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../../fixtures"));
        let all = std::env::var("IQRS_U6B_ALL").is_ok();
        let n_ex: usize = std::env::var("IQRS_U6_EXAMPLES").ok().and_then(|v| v.parse().ok()).unwrap_or(5);
        let mut caps = Vec::new();
        find_capture_dirs(&root, &mut caps).unwrap();
        let (mut exact, mut bad, mut shown) = (0u64, 0u64, 0usize);
        let mut s = AssignScratch::default();
        for cap in caps {
            let ctx = FixtureCtx::new(&cap);
            for f in unit_files(&cap, "U6").unwrap() {
                for c in read_cases::<U6AssignIn, U6AssignOut>(&f).unwrap() {
                    let c = c.unwrap();
                    let Some(want) = &c.out.assignment else { continue };
                    if !all && !(30..=32).contains(&want.atype.code()) {
                        continue;
                    }
                    let i = &c.input;
                    let index = ctx.index(&i.index).unwrap();
                    let params = AssignParams::from_fixture(&i.params, &MatchingParams::default_ont(), i.quick_mode);
                    let got = assign_with(&index, &params, &i.profile, i.profile.polya, &mut s)
                        .map(|a| AssignmentFixture::from_core(&a, &ctx.ids()));
                    let got = U6AssignOut { assignment: got };
                    if got == c.out {
                        exact += 1;
                    } else {
                        bad += 1;
                        if shown < n_ex {
                            shown += 1;
                            let d = diff_values(&got, &c.out, &DiffOptions::default());
                            println!("--- {}\n{}", c.case, d.map(|d| d.to_string()).unwrap_or_default());
                        }
                    }
                }
            }
        }
        println!("u6b: exact {exact}, mismatch {bad}");
        assert_eq!(bad, 0);
    }
}
