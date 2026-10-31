//! `PolyAFinder` (`terminal_prediction/polya_finder.py:19-231`), `PolyAFixer` and
//! `shift_polya/polyt` (`polya_verification.py:18-115`), and `reconcile_polya` /
//! `inject_external_polya` (`external_polya.py:119-146`). D-P1 / D-P2 sites use
//! `iq_types::quirk`.
//!
//! The finder never materialises the checked window or its reverse complement. It scans
//! `seq` in place through an `is_a(i)` predicate: for the polyA tail, `seq[start + i]` is
//! `A`; for the polyT head, `seq[end - 1 - i]` is `T` (= `A` in the reverse complement).
//! 4.1 upper-cases the window, so lower-case bases count too.

use iq_types::read::cigar_decode;
use iq_types::{quirk, DevId, Iv, PolyAInfo, Pos, RecordLite, Strand};

/// Parameters of the finder/fixer (`MatchingParams.polya_window/polya_fraction`,
/// `max_fake_terminal_exon_len`).
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct PolyAParams {
    pub window: i32,
    pub fraction: f64,
    pub max_fake_terminal_exon_len: i32,
}

impl PolyAParams {
    pub fn from_params(p: &iq_types::Params) -> Self {
        PolyAParams {
            window: p.matching.polya_window,
            fraction: p.matching.polya_fraction,
            max_fake_terminal_exon_len: p.matching.max_fake_terminal_exon_len,
        }
    }
    /// `int(window_size * min_polya_fraction)`.
    #[inline]
    fn threshold(&self) -> usize {
        (self.window as f64 * self.fraction) as usize
    }
}

// ------------------------------------------------------------------------------------------
// find_polya
// ------------------------------------------------------------------------------------------

/// `PolyAFinder.find_polya` (`:212-231`) over a virtual sequence of length `n`.
/// Returns the 0-based start in the virtual sequence.
///
/// D-P1: 4.1's `while i < len(seq) - window` never accepts the final window
/// (start `n - window`); iqrs tests it.
#[inline]
fn find_polya(n: usize, p: &PolyAParams, is_a: impl Fn(usize) -> bool) -> Option<usize> {
    let w = p.window as usize;
    if n < w || w == 0 {
        return None;
    }
    let thr = p.threshold();
    // exclusive bound on accepted window starts
    let limit = if quirk(DevId::P1) { n - w } else { n - w + 1 };
    let mut count = (0..w).filter(|&j| is_a(j)).count();
    let mut found = None;
    for i in 0..limit {
        if count >= thr {
            found = Some(i);
            break;
        }
        if i + w < n {
            count = count + is_a(i + w) as usize - is_a(i) as usize;
        }
    }
    let i = found?;
    // i + max(0, seq[i:].find("AA"))
    let aa = (i..n.saturating_sub(1)).find(|&j| is_a(j) && is_a(j + 1)).map_or(0, |j| j - i);
    Some(i + aa)
}

/// `check_entire_*`: `tail.count('A') < len(tail) * fraction` rejects the call.
#[inline]
fn entire_ok(pos: usize, n: usize, p: &PolyAParams, is_a: impl Fn(usize) -> bool) -> bool {
    let cnt = (pos..n).filter(|&j| is_a(j)).count();
    // operands are finite, so `>=` is the exact negation of 4.1's `<`
    (cnt as f64) >= (n - pos) as f64 * p.fraction
}

#[inline]
fn base_is(b: u8, upper: u8) -> bool {
    b & !0x20 == upper
}

// ------------------------------------------------------------------------------------------
// move_ref_coord_alogn_alignment
// ------------------------------------------------------------------------------------------

/// `move_ref_coord_alogn_alignment(alignment, shift)` (`:19-78`). Always non-negative in
/// practice (`reference_length_consumed - 1`).
fn move_ref_coord(cigar: &[u32], shift: i64) -> i64 {
    if shift == 0 {
        return 0;
    }
    let n = cigar.len() as i64;
    let op = |i: i64| -> (u8, i64) {
        let w = cigar[if i < 0 { (n + i) as usize } else { i as usize }];
        ((w & 0xf) as u8, (w >> 4) as i64)
    };
    let direction: i64 = if shift > 0 { 1 } else { -1 };
    let mut shift = shift.abs();
    let mut cur: i64;
    if direction == 1 {
        cur = 0;
        if n > 1 && op(0).0 == 5 && op(1).0 == 4 {
            cur = 2;
        } else if matches!(op(0).0, 4 | 5) {
            cur = 1;
        }
    } else {
        cur = -1;
        if n > 1 && op(-1).0 == 5 && op(-2).0 == 4 {
            cur -= 2;
        } else if matches!(op(-1).0, 4 | 5) {
            cur -= 1;
        }
    }
    let mut read_consumed: i64 = 0;
    let mut ref_consumed: i64 = 0;
    shift += 1; // 4.1: "a hack to jump to the next event"
    while cur < n && cur >= -n && read_consumed < shift {
        let (ev, len) = op(cur);
        match ev {
            1 => read_consumed += len,
            2 | 3 => ref_consumed += len,
            0 | 7 | 8 => {
                let remaining = shift - read_consumed;
                if len < remaining {
                    ref_consumed += len;
                    read_consumed += len;
                } else {
                    read_consumed += remaining;
                    ref_consumed += remaining;
                }
            }
            4 | 5 => break,
            // 4.1 would raise TypeError in its log call here (P); never seen in practice
            _ => {}
        }
        cur += direction;
    }
    ref_consumed - 1
}

// ------------------------------------------------------------------------------------------
// Tail / head
// ------------------------------------------------------------------------------------------

/// Read-side geometry shared by the four searches.
struct Aln<'a> {
    cigar: &'a [u32],
    seq: &'a [u8],
    /// `reference_start` (0-based).
    ref_start: i64,
    /// `reference_end` (0-based exclusive).
    ref_end: i64,
}

impl Aln<'_> {
    fn op(&self, i: usize) -> u8 {
        (self.cigar[i] & 0xf) as u8
    }
    fn len_at(&self, i: usize) -> usize {
        cigar_decode(self.cigar[i]).1 as usize
    }
    fn soft_tail(&self) -> usize {
        let n = self.cigar.len();
        if n > 1 && self.op(n - 1) == 5 && self.op(n - 2) == 4 {
            self.len_at(n - 2)
        } else if n > 0 && self.op(n - 1) == 4 {
            self.len_at(n - 1)
        } else {
            0
        }
    }
    fn soft_head(&self) -> usize {
        let n = self.cigar.len();
        if n > 1 && self.op(0) == 5 && self.op(1) == 4 {
            self.len_at(1)
        } else if n > 0 && self.op(0) == 4 {
            self.len_at(0)
        } else {
            0
        }
    }

    /// `find_polya_tail(alignment, from_pos, to_pos, check_entire_tail)`.
    fn polya_tail(&self, from_pos: usize, to_pos: usize, check_entire: bool, p: &PolyAParams) -> Pos {
        let l = self.seq.len();
        let sc = self.soft_tail();
        if l == 0 || sc >= l {
            // 4.1: `return -1` for an empty sequence; AssertionError for sc >= len
            return -1;
        }
        let mapped_end = l - sc;
        let cs = mapped_end.saturating_sub(from_pos);
        let ce = l.min(mapped_end + to_pos + 1);
        let n = ce - cs;
        let seq = &self.seq[cs..ce];
        let is_a = |i: usize| base_is(seq[i], b'A');
        let Some(pos) = find_polya(n, p, is_a) else { return -1 };
        if check_entire && !entire_ok(pos, n, p, is_a) {
            return -1;
        }
        let pos = (cs + pos) as i64;
        let shift = pos - mapped_end as i64;
        let r = if shift >= 0 { self.ref_end + shift } else { self.ref_end - move_ref_coord(self.cigar, shift) };
        r as Pos
    }

    /// `find_polyt_head(alignment, from_pos, to_pos, check_entire_head)`.
    fn polyt_head(&self, from_pos: usize, to_pos: usize, check_entire: bool, p: &PolyAParams) -> Pos {
        let l = self.seq.len();
        let sc = self.soft_head();
        if l == 0 || sc >= l {
            return -1;
        }
        let mapped_start = sc;
        let cs = mapped_start.saturating_sub(to_pos);
        let ce = l.min(mapped_start + from_pos + 1);
        let n = ce - cs;
        let seq = &self.seq[cs..ce];
        // reverse complement: rc[i] == 'A' <=> seq[n-1-i] == 'T'
        let is_a = |i: usize| base_is(seq[n - 1 - i], b'T');
        let Some(pos) = find_polya(n, p, is_a) else { return -1 };
        if check_entire && !entire_ok(pos, n, p, is_a) {
            return -1;
        }
        let pos = ce as i64 - pos as i64 - 1;
        let ms = mapped_start as i64;
        let r = if pos <= ms {
            self.ref_start - (ms - pos)
        } else {
            self.ref_start + move_ref_coord(self.cigar, pos - ms)
        };
        r.max(1) as Pos
    }
}

/// `PolyAFinder.detect_polya(alignment)`.
pub fn detect_polya(rec: &RecordLite<'_>, p: &PolyAParams) -> PolyAInfo {
    if rec.cigar.is_empty() {
        return iq_types::NO_POLYA;
    }
    let a = Aln { cigar: rec.cigar, seq: rec.seq, ref_start: rec.pos0, ref_end: rec.end0() };
    let w = p.window.max(0) as usize;
    PolyAInfo {
        ext_a: a.polya_tail(2, 2 * w, false, p),
        ext_t: a.polyt_head(2, 2 * w, false, p),
        int_a: a.polya_tail(4 * w, 2, true, p),
        int_t: a.polyt_head(4 * w, 2, true, p),
    }
}

// ------------------------------------------------------------------------------------------
// PolyAFixer
// ------------------------------------------------------------------------------------------

/// `count_polya_exons` (`:74-95`). D-P2: 4.1 skips a non-qualifying exon and keeps
/// counting older ones; iqrs stops at the first non-qualifying exon. On sorted, disjoint
/// exons (all `get_read_blocks` output) both give the same count: see the
/// `d_p2_unreachable_on_disjoint_exons` test. The gate is kept so §9 stays checkable.
fn count_polya_exons(exons: &[Iv], int_a: Pos, max_fake: i32) -> i32 {
    if int_a == -1 {
        return 0;
    }
    let mut count = 0;
    for ex in exons.iter().rev() {
        if ex.e <= int_a {
            break;
        }
        let len_to = int_a - ex.s;
        if len_to <= 0 || (len_to <= max_fake && ex.e - int_a > 2 * len_to) {
            count += 1;
        } else if !quirk(DevId::P2) {
            break;
        }
    }
    count
}

/// `count_polyt_exons` (`:99-115`), mirror of [`count_polya_exons`] (D-P2).
fn count_polyt_exons(exons: &[Iv], int_t: Pos, max_fake: i32) -> i32 {
    if int_t == -1 {
        return 0;
    }
    let mut count = 0;
    for ex in exons {
        if ex.s >= int_t {
            break;
        }
        let len_to = ex.e - int_t;
        if len_to <= 0 || (len_to <= max_fake && 2 * len_to < int_t - ex.s) {
            count += 1;
        } else if !quirk(DevId::P2) {
            break;
        }
    }
    count
}

/// `PolyAFixer.correct_read_info` (`:57-71`): `(polya_exon_count, polyt_exon_count)`.
///
/// When the two counts together exceed the exon count, 4.1 raises an `IndexError` while
/// trimming. iqrs trims nothing in that case (no fixture can show it).
pub fn correct_read_info(exons: &[Iv], polya: &PolyAInfo, max_fake: i32) -> (i32, i32) {
    let n = exons.len() as i32;
    if n == 1 {
        return (0, 0);
    }
    let mut pa = count_polya_exons(exons, polya.int_a, max_fake);
    let mut pt = count_polyt_exons(exons, polya.int_t, max_fake);
    if pa + pt == n {
        pa -= 1;
        pt -= 1;
    } else if pa + pt > n {
        return (0, 0);
    }
    (pa, pt)
}

/// `shift_polya(read_exons, exon_count, polya_pos)` (`:18-33`).
fn shift_polya(exons: &[Iv], count: usize, pos: Pos) -> Pos {
    if count == 0 || count == exons.len() || pos == -1 {
        return pos;
    }
    let mut dist = 0;
    for ex in exons.iter().rev().take(count) {
        if ex.s > pos {
            continue;
        } else if dist == 0 {
            dist += pos - ex.s;
        } else {
            dist += ex.e - ex.s + 1;
        }
    }
    exons[exons.len() - count - 1].e + dist
}

/// `shift_polyt(read_exons, exon_count, polyt_pos)` (`:37-52`).
fn shift_polyt(exons: &[Iv], count: usize, pos: Pos) -> Pos {
    if count == 0 || count == exons.len() || pos == -1 {
        return pos;
    }
    let mut dist = 0;
    for ex in exons.iter().take(count) {
        if ex.e < pos {
            continue;
        } else if dist == 0 {
            dist += ex.e - pos;
        } else {
            dist += ex.e - ex.s + 1;
        }
    }
    exons[count].s - dist
}

// ------------------------------------------------------------------------------------------
// external polyA (no-ops for --polya_trimmed none: strand None, no_tail false)
// ------------------------------------------------------------------------------------------

/// What an external polyA source says about the read (`resolve_external_strand`,
/// `ExternalPolyAStatus.no_tail`). `Default` = no source (`--polya_trimmed none`).
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct ExternalPolyA {
    pub strand: Option<Strand>,
    pub no_tail: bool,
}

fn reconcile_polya(p: &mut PolyAInfo, ext: ExternalPolyA) {
    if ext.no_tail {
        *p = iq_types::NO_POLYA;
    } else if ext.strand == Some(Strand::Plus) {
        p.ext_t = -1;
        p.int_t = -1;
    } else if ext.strand == Some(Strand::Minus) {
        p.ext_a = -1;
        p.int_a = -1;
    }
}

fn inject_external_polya(p: &mut PolyAInfo, exons: &[Iv], strand: Option<Strand>) {
    match strand {
        Some(Strand::Plus) if p.ext_a == -1 && p.int_a == -1 => p.ext_a = exons[exons.len() - 1].e + 1,
        Some(Strand::Minus) if p.ext_t == -1 && p.int_t == -1 => p.ext_t = exons[0].s - 1,
        _ => {}
    }
}

// ------------------------------------------------------------------------------------------
// AlignmentInfo.add_polya_info
// ------------------------------------------------------------------------------------------

/// `AlignmentInfo.add_polya_info` (`alignment_info.py:121-155`). Trims the three block
/// lists in place (no allocation). `exons` must be non-empty. Returns
/// `(polya_info, exons_changed)`.
pub fn add_polya_info(
    rec: &RecordLite<'_>,
    p: &PolyAParams,
    ext: ExternalPolyA,
    exons: &mut Vec<Iv>,
    read_blocks: &mut Vec<(i32, i32)>,
    cigar_blocks: &mut Vec<(i32, i32)>,
) -> (PolyAInfo, bool) {
    let mut info = detect_polya(rec, p);
    reconcile_polya(&mut info, ext);
    let (pa, pt) = correct_read_info(exons, &info, p.max_fake_terminal_exon_len);
    let mut changed = false;
    if pa > 0 {
        let pa = pa as usize;
        info.int_a = shift_polya(exons, pa, info.int_a);
        info.ext_a = shift_polya(exons, pa, info.ext_a);
        let keep = exons.len() - pa;
        exons.truncate(keep);
        read_blocks.truncate(keep);
        cigar_blocks.truncate(keep);
        changed = true;
    }
    if pt > 0 {
        let pt = pt as usize;
        info.int_t = shift_polyt(exons, pt, info.int_t);
        info.ext_t = shift_polyt(exons, pt, info.ext_t);
        exons.drain(..pt);
        read_blocks.drain(..pt);
        cigar_blocks.drain(..pt);
        changed = true;
    }
    inject_external_polya(&mut info, exons, ext.strand);
    (info, changed)
}

// ------------------------------------------------------------------------------------------
// Fixture runner: U4 process_alignment (the crate entry point, plus `init` blocks)
// ------------------------------------------------------------------------------------------

pub(crate) mod fixture {
    use super::*;
    use iq_types::fixture::schema::{U4BlocksFixture, U4ProcessAlignmentIn, U4ProcessAlignmentOut};
    use iq_types::fixture::{FixError, FixResult, FixtureCtx, UnitFixture};
    use iq_types::AlignedRead;

    pub struct U4ProcessAlignment;

    impl UnitFixture for U4ProcessAlignment {
        const UNIT: &'static str = "U4";
        const FN: &'static str = "process_alignment";
        const DEVIATIONS: &'static [DevId] = &[DevId::P1, DevId::P2];
        type In = U4ProcessAlignmentIn;
        type Out = U4ProcessAlignmentOut;

        fn run(i: &U4ProcessAlignmentIn, _ctx: &FixtureCtx) -> FixResult<U4ProcessAlignmentOut> {
            let rec = i.aln.to_record();
            let rec = rec.as_lite();
            let p = PolyAParams {
                window: i.polya_window,
                fraction: i.polya_fraction,
                max_fake_terminal_exon_len: i.max_fake_terminal_exon_len,
            };
            let strand = match i.external_strand.as_deref() {
                None => None,
                Some(s) => Some(Strand::parse(s).ok_or_else(|| FixError(format!("bad external_strand {s:?}")))?),
            };
            let ext = ExternalPolyA { strand, no_tail: i.external_no_tail };
            let init = crate::cigar::read_blocks(rec.pos0, rec.cigar);
            let init = U4BlocksFixture { read_exons: init.ref_blocks, read_blocks: init.read_blocks, cigar_blocks: init.cigar_blocks };
            let mut ar = AlignedRead {
                exons: vec![],
                read_blocks: vec![],
                cigar_blocks: vec![],
                polya: iq_types::NO_POLYA,
                exons_changed: false,
                deletions: vec![],
            };
            let called = crate::process_alignment_into(&rec, &p, ext, &mut ar);
            Ok(U4ProcessAlignmentOut {
                init,
                polya: called.then_some(ar.polya),
                read_exons: ar.exons,
                read_blocks: ar.read_blocks,
                cigar_blocks: ar.cigar_blocks,
                exons_changed: ar.exons_changed,
                polya_called: called,
            })
        }
    }
}

pub(crate) fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    reg.add::<fixture::U4ProcessAlignment>();
}

#[cfg(test)]
mod tests {
    use super::*;
    use iq_types::read::cigar_encode;
    use iq_types::{with_quirks, CigarOp, QuirkSet};

    const P: PolyAParams = PolyAParams { window: 16, fraction: 0.75, max_fake_terminal_exon_len: 20 };

    fn fp(s: &[u8]) -> Option<usize> {
        find_polya(s.len(), &P, |i| base_is(s[i], b'A'))
    }

    #[test]
    fn find_polya_last_window_d_p1() {
        // C + 16 A: window 0 already has 15 A >= 12, so both builds call it at 1
        let s = [&b"C"[..], &[b'A'; 16]].concat();
        assert_eq!(fp(&s), Some(1));
        let mut t = b"CCCCCCCCCCCCCCCCC".to_vec();
        t.extend_from_slice(b"AAAAAAAAAAAA"); // 17 C + 12 A, len 29; first qualifying window starts 13
        let q41 = with_quirks(QuirkSet::only(DevId::P1), || fp(&t));
        let iq = with_quirks(QuirkSet::NONE, || fp(&t));
        assert_eq!(iq, Some(17)); // window 13 accepted, then "AA" at 17
        if iq_types::QUIRKS41 {
            assert_eq!(q41, None);
        }
        // exact window length: 4.1 never accepts
        let u = [b'A'; 16];
        assert_eq!(with_quirks(QuirkSet::NONE, || fp(&u)), Some(0));
        if iq_types::QUIRKS41 {
            assert_eq!(with_quirks(QuirkSet::only(DevId::P1), || fp(&u)), None);
        }
    }

    /// D-P2 has no reachable effect on sorted, disjoint exons (all `get_read_blocks`
    /// output): an exon that fails the test has `s < int_a`, so every older exon has
    /// `e < s < int_a` and 4.1 breaks there anyway. Mirror for polyT.
    #[test]
    fn d_p2_unreachable_on_disjoint_exons() {
        let mut x: u64 = 0x9e3779b97f4a7c15;
        let mut rnd = |m: u64| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            (x % m) as i32
        };
        for _ in 0..20000 {
            let n = 1 + rnd(6) as usize;
            let mut ex = Vec::with_capacity(n);
            let mut pos = 1 + rnd(50);
            for _ in 0..n {
                let s0 = pos + 1 + rnd(40);
                let e0 = s0 + rnd(60);
                ex.push(Iv::new(s0, e0));
                pos = e0 + 1;
            }
            let ia = rnd(pos as u64 + 20);
            let it = rnd(pos as u64 + 20);
            let cases: [(fn(&[Iv], Pos, i32) -> i32, Pos); 2] = [(count_polya_exons, ia), (count_polyt_exons, it)];
            for (f, pos) in cases {
                let a = with_quirks(QuirkSet::only(DevId::P2), || f(&ex, pos, 20));
                let b = with_quirks(QuirkSet::NONE, || f(&ex, pos, 20));
                assert_eq!(a, b, "{ex:?} pos={pos}");
            }
        }
    }

    /// End to end: a leading 20-bp polyT exon is trimmed; D-P3 moves the deletion anchor.
    #[test]
    fn polyt_trim_and_d_p3_anchor() {
        let cigar: Vec<u32> = [(CigarOp::M, 20), (CigarOp::N, 100), (CigarOp::M, 10), (CigarOp::D, 2), (CigarOp::M, 60)]
            .iter()
            .map(|&(o, l)| cigar_encode(o, l))
            .collect();
        let mut seq = vec![b'T'; 20];
        seq.extend(b"ACGCGTCAGCGCAGCTGACGCGTCAGCGCAGCTGACGCGTCAGCGCAGCTGACGCGTCAGCGCAGCTGACG".iter());
        let rec = iq_types::RecordOwned { qname: "r".into(), flag: 0, pos0: 999, mapq: 60, cigar, seq, tag: None };
        let run = |q| {
            with_quirks(q, || {
                let mut ar = iq_types::AlignedRead {
                    exons: vec![],
                    read_blocks: vec![],
                    cigar_blocks: vec![],
                    polya: iq_types::NO_POLYA,
                    exons_changed: false,
                    deletions: vec![],
                };
                assert!(crate::process_alignment_into(&rec.as_lite(), &P, ExternalPolyA::default(), &mut ar));
                ar
            })
        };
        let iq = run(QuirkSet::NONE);
        assert!(iq.exons_changed);
        assert_eq!(iq.exons, vec![Iv::new(1120, 1191)]);
        assert_eq!(iq.cigar_blocks, vec![(2, 4)]);
        assert_eq!(iq.deletions, vec![(1130, 2)]); // anchored at reference_start + 1
        if iq_types::QUIRKS41 {
            let q = run(QuirkSet::only(DevId::P3));
            // 4.1 walks the full CIGAR from the trimmed exon start: the deletion lands at
            // 1120 + 20 + 100 + 10 = 1250, more than margin from every boundary, so it is lost
            assert!(q.deletions.is_empty());
        }
    }

    #[test]
    fn shift_and_trim() {
        let ex = vec![Iv::new(100, 200), Iv::new(300, 310)];
        assert_eq!(shift_polya(&ex, 1, 305), 205);
        assert_eq!(shift_polyt(&ex, 1, 150), 250);
        assert_eq!(correct_read_info(&ex, &PolyAInfo { int_a: 301, ..iq_types::NO_POLYA }, 20), (1, 0));
    }

    #[test]
    fn move_ref() {
        let c: Vec<u32> = [(CigarOp::S, 3), (CigarOp::M, 10), (CigarOp::D, 2), (CigarOp::M, 10), (CigarOp::S, 4)]
            .iter()
            .map(|&(o, l)| cigar_encode(o, l))
            .collect();
        assert_eq!(move_ref_coord(&c, 0), 0);
        assert_eq!(move_ref_coord(&c, 5), 5);
        assert_eq!(move_ref_coord(&c, 12), 14);
        assert_eq!(move_ref_coord(&c, -5), 5);
    }
}
