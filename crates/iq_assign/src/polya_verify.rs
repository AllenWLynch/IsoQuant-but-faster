//! U6a: `PolyAVerifier` (`terminal_prediction/polya_verification.py:117-390`) and the
//! `shift_polya` / `shift_polyt` helpers (`:18-49`).
//!
//! The verifier edits an isoform match's event list in place, exactly like 4.1 (which
//! mutates the Python list it is given): elongation events on the polyA side are replaced
//! by a polyA-site event, and terminal-exon misalignments may be appended.

use iq_types::params::MatchingParams;
use iq_types::{Events, EventType, IsoformIndex, Iv, MatchEvent, PolyAInfo, Strand};

/// `math.inf` stand-in for absent polyA distances (all real distances are far smaller).
const INF: i64 = i64::MAX;

#[inline]
fn ilen(iv: &Iv) -> i64 {
    iv.e as i64 - iv.s as i64 + 1
}

/// `shift_polya` (`:18-32`): polyA position recomputed when the last `exon_count` read exons
/// are fake.
pub(crate) fn shift_polya(read_exons: &[Iv], exon_count: usize, polya_pos: i64) -> i64 {
    let n = read_exons.len();
    if exon_count == 0 || exon_count >= n || polya_pos == -1 {
        // `exon_count > n` raises in 4.1 (assertion in the caller); unreachable in practice
        return polya_pos;
    }
    let mut dist = 0i64;
    for i in 0..exon_count {
        let exon = read_exons[n - 1 - i];
        if exon.s as i64 > polya_pos {
            continue;
        } else if dist == 0 {
            dist += polya_pos - exon.s as i64;
        } else {
            dist += ilen(&exon);
        }
    }
    read_exons[n - exon_count - 1].e as i64 + dist
}

/// `shift_polyt` (`:35-49`).
pub(crate) fn shift_polyt(read_exons: &[Iv], exon_count: usize, polyt_pos: i64) -> i64 {
    let n = read_exons.len();
    if exon_count == 0 || exon_count >= n || polyt_pos == -1 {
        return polyt_pos;
    }
    let mut dist = 0i64;
    for exon in &read_exons[..exon_count] {
        if (exon.e as i64) < polyt_pos {
            continue;
        } else if dist == 0 {
            dist += exon.e as i64 - polyt_pos;
        } else {
            dist += ilen(exon);
        }
    }
    read_exons[exon_count].s as i64 - dist
}

/// `PolyAVerifier` bound to one index and parameter set.
#[derive(Clone, Copy)]
pub(crate) struct PolyAVerifier<'a> {
    pub index: &'a IsoformIndex,
    pub params: &'a MatchingParams,
}

impl<'a> PolyAVerifier<'a> {
    pub(crate) fn new(index: &'a IsoformIndex, params: &'a MatchingParams) -> Self {
        PolyAVerifier { index, params }
    }

    /// `verify_read_ends` (`:123-148`) for isoform `iso` (index into `index.iso`).
    pub(crate) fn verify_read_ends(&self, read_exons: &[Iv], polya: PolyAInfo, iso: usize, events: &mut Events) {
        let iso_exons = self.index.iso_exons(iso);
        match self.index.iso[iso].strand {
            Strand::Plus => {
                let is_internal = check_internal(polya.int_a, events, EventType::IncompleteIntronRetentionRight, EventType::InternalPolyaRight);
                if !is_internal && (polya.ext_a != -1 || polya.int_a != -1) {
                    self.verify_polya(iso_exons, read_exons, polya, events);
                }
            }
            Strand::Minus => {
                let is_internal = check_internal(polya.int_t, events, EventType::IncompleteIntronRetentionLeft, EventType::InternalPolyaLeft);
                if !is_internal && (polya.ext_t != -1 || polya.int_t != -1) {
                    self.verify_polyt(iso_exons, read_exons, polya, events);
                }
            }
            Strand::Unknown => {}
        }
        if events.is_empty() {
            events.push(MatchEvent::new(EventType::None));
        }
    }

    /// `verify_polya` (`:150-208`).
    fn verify_polya(&self, iso_exons: &[Iv], read_exons: &[Iv], polya: PolyAInfo, events: &mut Events) {
        let isoform_end = iso_exons[iso_exons.len() - 1].e as i64;
        let mut to_remove = None;
        let mut fake_count = 0usize;
        let mut misaligned = 0usize;
        for (i, e) in events.iter().enumerate() {
            match e.ty {
                EventType::MajorExonElongationRight | EventType::ExonElongationRight => to_remove = Some(i),
                EventType::FakeTerminalExonRight => fake_count += 1,
                EventType::TerminalExonMisalignmentRight => misaligned += 1,
                _ => {}
            }
        }
        if let Some(i) = to_remove {
            events.remove(i);
        }
        let (ext0, int0) = (polya.ext_a as i64, polya.int_a as i64);
        if self.check_if_close(isoform_end, ext0, int0, events, EventType::CorrectPolyaSiteRight) {
            return;
        }
        let mut ext = shift_polya(read_exons, fake_count, ext0);
        let mut int = shift_polya(read_exons, fake_count, int0);
        if misaligned > 0 {
            ext = isoform_end;
            int = isoform_end;
        } else {
            (ext, int) = self.detect_reference_exons_beyond_polya(iso_exons, ext, int, events);
        }
        if self.check_if_close(isoform_end, ext, int, events, EventType::CorrectPolyaSiteRight) {
            return;
        }
        let pos = if int == -1 { ext } else { int };
        let ty = if (pos - isoform_end).abs() > self.params.apa_delta as i64 {
            EventType::AlternativePolyaSiteRight
        } else {
            EventType::CorrectPolyaSiteRight
        };
        events.push(MatchEvent::with_info(ty, pos));
    }

    /// `verify_polyt` (`:210-267`).
    fn verify_polyt(&self, iso_exons: &[Iv], read_exons: &[Iv], polya: PolyAInfo, events: &mut Events) {
        let isoform_start = iso_exons[0].s as i64;
        let mut to_remove = None;
        let mut fake_count = 0usize;
        let mut misaligned = 0usize;
        for (i, e) in events.iter().enumerate() {
            match e.ty {
                EventType::MajorExonElongationLeft | EventType::ExonElongationLeft => to_remove = Some(i),
                EventType::FakeTerminalExonLeft => fake_count += 1,
                EventType::TerminalExonMisalignmentLeft => misaligned += 1,
                _ => {}
            }
        }
        if let Some(i) = to_remove {
            events.remove(i);
        }
        let (ext0, int0) = (polya.ext_t as i64, polya.int_t as i64);
        if self.check_if_close(isoform_start, ext0, int0, events, EventType::CorrectPolyaSiteLeft) {
            return;
        }
        let mut ext = shift_polyt(read_exons, fake_count, ext0);
        let mut int = shift_polyt(read_exons, fake_count, int0);
        if misaligned > 0 {
            ext = isoform_start;
            int = isoform_start;
        } else {
            (ext, int) = self.detect_reference_exons_before_polyt(iso_exons, ext, int, events);
        }
        if self.check_if_close(isoform_start, ext, int, events, EventType::CorrectPolyaSiteLeft) {
            return;
        }
        let pos = if int == -1 { ext } else { int };
        let ty = if (pos - isoform_start).abs() > self.params.apa_delta as i64 {
            EventType::AlternativePolyaSiteLeft
        } else {
            EventType::CorrectPolyaSiteLeft
        };
        events.push(MatchEvent::with_info(ty, pos));
    }

    /// `detect_reference_exons_beyond_polya` (`:312-338`).
    fn detect_reference_exons_beyond_polya(&self, iso_exons: &[Iv], ext: i64, int: i64, events: &mut Events) -> (i64, i64) {
        let n = iso_exons.len();
        let polya_pos = if int != -1 { int } else { ext };
        let mut cnt = 0usize;
        while cnt < n && iso_exons[n - cnt - 1].s as i64 >= polya_pos {
            cnt += 1;
        }
        if cnt == n || cnt == 0 {
            return (ext, int);
        }
        let term_len: i64 = iso_exons[n - cnt..].iter().map(ilen).sum();
        let anchor = iso_exons[n - cnt - 1].e as i64;
        let dist = (anchor - ext).abs().min((anchor - int).abs());
        if self.short_terminal_exons(term_len, dist) {
            for i in 0..cnt {
                let r = n as i64 - 2 - i as i64;
                events.push(MatchEvent::with_regions(EventType::TerminalExonMisalignmentRight, (r, r), iq_types::match_consts::UNDEFINED_REGION));
            }
            let end = iso_exons[n - 1].e as i64;
            return (end, end);
        }
        (ext, int)
    }

    /// `detect_reference_exons_before_polyt` (`:341-367`).
    fn detect_reference_exons_before_polyt(&self, iso_exons: &[Iv], ext: i64, int: i64, events: &mut Events) -> (i64, i64) {
        let n = iso_exons.len();
        let polyt_pos = if int != -1 { int } else { ext };
        let mut cnt = 0usize;
        while cnt < n && iso_exons[cnt].e as i64 <= polyt_pos {
            cnt += 1;
        }
        if cnt == 0 || cnt == n {
            return (ext, int);
        }
        let term_len: i64 = iso_exons[..cnt].iter().map(ilen).sum();
        let anchor = iso_exons[cnt].s as i64;
        let dist = (anchor - ext).abs().min((anchor - int).abs());
        if self.short_terminal_exons(term_len, dist) {
            for i in 0..cnt as i64 {
                events.push(MatchEvent::with_regions(EventType::TerminalExonMisalignmentLeft, (i, i), iq_types::match_consts::UNDEFINED_REGION));
            }
            let start = iso_exons[0].s as i64;
            return (start, start);
        }
        (ext, int)
    }

    #[inline]
    fn short_terminal_exons(&self, term_len: i64, dist: i64) -> bool {
        let p = self.params;
        (term_len <= p.max_fake_terminal_exon_len as i64 && dist <= p.max_fake_terminal_exon_len as i64)
            || (term_len <= p.max_missed_exon_len as i64 && (term_len - dist).abs() <= p.delta as i64)
    }

    /// `check_if_close` (`:369-385`): appends `ty` and returns true when a polyA position is
    /// within `apa_delta` of the isoform end.
    fn check_if_close(&self, iso_end: i64, ext: i64, int: i64, events: &mut Events, ty: EventType) -> bool {
        let apa = self.params.apa_delta as i64;
        let d_ext = if ext != -1 { (iso_end - ext).abs() } else { INF };
        let d_int = if int != -1 { (iso_end - int).abs() } else { INF };
        if d_int <= apa && d_int <= d_ext {
            events.push(MatchEvent::with_info(ty, int));
            true
        } else if d_ext <= apa && d_int > d_ext {
            events.push(MatchEvent::with_info(ty, ext));
            true
        } else {
            false
        }
    }
}

/// `check_internal_polya` / `check_internal_polyt` (`:270-296`).
fn check_internal(pos: i32, events: &mut Events, retention: EventType, internal: EventType) -> bool {
    if pos == -1 {
        return false;
    }
    if let Some(e) = events.iter().find(|e| e.ty == retention) {
        let ev = MatchEvent { ty: internal, iso_region: e.iso_region, read_region: iq_types::match_consts::UNDEFINED_REGION, info: pos as i64 };
        events.push(ev);
        return true;
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shifts() {
        let ex = [Iv::new(100, 200), Iv::new(300, 400), Iv::new(500, 510)];
        assert_eq!(shift_polya(&ex, 0, 505), 505);
        // last exon fake: polyA at 505 is 5 bases into it -> 400 + 5
        assert_eq!(shift_polya(&ex, 1, 505), 405);
        assert_eq!(shift_polya(&ex, 3, 505), 505);
        assert_eq!(shift_polya(&ex, 1, -1), -1);
        assert_eq!(shift_polyt(&ex, 1, 150), 300 - 50);
    }
}
