//! U7: `ExonCorrector.correct_assigned_read` (`exon_corrector.py:29-213`) with
//! `AlignmentInfo.set_aligned_pairs/get_error_count` (`alignment_info.py:42-97`),
//! `get_assignment_strand` (`alignment_processor.py:538-551`), `StrandDetector`
//! (`gene_info.py:924-967`) with `get_intron_strand` (`common.py:854-864`), and
//! `check_antisense` (`alignment_processor.py:684-695`).
//!
//! Hot-path notes (design §6.4.5):
//! - 4.1 materialises `get_aligned_pairs()` for the whole read. Here the aligned-pair list
//!   stays virtual: the pair index of an exon boundary comes from one monotone CIGAR walk,
//!   and only the ≤24 pairs of the window around a site are generated, on demand.
//! - Reference bases come from [`RefBases`]: a contiguous [`RefSlice`] (FASTA mmap in the
//!   pipeline) or a list of windows (fixtures). There is no per-site FASTA lookup.
//! - `StrandDetector`'s per-intron cache is not ported. It memoises a pure function of the
//!   reference, and four bytes from a slice cost less than a hash lookup, so outputs are
//!   identical.
//!
//! No deviations: the 4.1 behaviour of all three functions is kept as is.

use iq_types::{
    match_consts, AlignedRead, Assignment, AssignmentType, EventType, IsoformIndex, Iv, MatchClass, MatchEvent, Params,
    PolyAInfo, Pos, RecordLite, RefSlice, Strand,
};

// ------------------------------------------------------------------------------------------
// Reference access
// ------------------------------------------------------------------------------------------

/// Source of reference bases (1-based positions, bases verbatim, no case folding).
pub trait RefBases {
    /// Base at 1-based `pos`, or `None` outside the sequence.
    fn base_at(&self, pos: Pos) -> Option<u8>;
    /// No reference at all (4.1 `chr_record` is `None`): mismatches are not counted.
    fn is_absent(&self) -> bool;
}

impl RefBases for RefSlice<'_> {
    #[inline]
    fn base_at(&self, pos: Pos) -> Option<u8> {
        self.base(pos)
    }
    fn is_absent(&self) -> bool {
        self.seq.is_empty()
    }
}

/// Sorted (by `start`) list of windows, possibly overlapping; all windows must agree where
/// they overlap (they are slices of one chromosome).
impl RefBases for [RefSlice<'_>] {
    fn base_at(&self, pos: Pos) -> Option<u8> {
        let k = self.partition_point(|w| w.start <= pos);
        self[..k].iter().rev().find_map(|w| w.base(pos))
    }
    fn is_absent(&self) -> bool {
        self.iter().all(|w| w.seq.is_empty())
    }
}

// ------------------------------------------------------------------------------------------
// Parameters
// ------------------------------------------------------------------------------------------

/// The `ExonCorrector` parameters (`params.delta` and the `SplicSiteCorrectionStrategy`
/// flags, `isoquant.py:1143-1161`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CorrectParams {
    pub delta: i32,
    pub fuzzy_junctions: bool,
    pub intron_shifts: bool,
    pub skipped_exons: bool,
    pub fake_terminal_exons: bool,
    pub terminal_exons: bool,
    pub microintron_retention: bool,
}

impl CorrectParams {
    pub fn from_params(p: &Params) -> Self {
        CorrectParams {
            delta: p.matching.delta,
            fuzzy_junctions: p.splice.fuzzy_junctions,
            intron_shifts: p.splice.intron_shifts,
            skipped_exons: p.splice.skipped_exons,
            fake_terminal_exons: p.splice.fake_terminal_exons,
            terminal_exons: p.splice.terminal_exons,
            microintron_retention: p.splice.microintron_retention,
        }
    }
}

// ------------------------------------------------------------------------------------------
// Virtual aligned pairs (pysam `get_aligned_pairs()`, matches_only=False)
// ------------------------------------------------------------------------------------------

/// `AlignmentInfo.REGION_TO_CHECK_LEN`.
const REGION_TO_CHECK_LEN: usize = 12;

#[derive(Clone, Copy)]
enum PairKind {
    /// M, =, X: `(qpos, rpos)`.
    Match,
    /// I, S: `(qpos, None)`.
    Query,
    /// D, N: `(None, rpos)`.
    Ref,
}

/// Decodes one CIGAR word into its pair kind and length. H and P produce no pairs.
#[inline]
fn pair_op(w: u32) -> Option<(PairKind, u32)> {
    let len = w >> 4;
    match w & 0xf {
        0 | 7 | 8 => Some((PairKind::Match, len)),
        1 | 4 => Some((PairKind::Query, len)),
        2 | 3 => Some((PairKind::Ref, len)),
        _ => None,
    }
}

/// Monotone cursor over the CIGAR that finds the index of the aligned pair with a given
/// 0-based reference position.
struct Seeker<'a> {
    cigar: &'a [u32],
    k: usize,
    pair_base: usize,
    rpos: i64,
}

impl<'a> Seeker<'a> {
    fn new(cigar: &'a [u32], pos0: i64) -> Self {
        Seeker { cigar, k: 0, pair_base: 0, rpos: pos0 }
    }
    /// Pair index of reference position `x` (0-based). Queries must be non-decreasing.
    fn seek(&mut self, x: i64) -> Option<usize> {
        while self.k < self.cigar.len() {
            let Some((kind, len)) = pair_op(self.cigar[self.k]) else {
                self.k += 1;
                continue;
            };
            if !matches!(kind, PairKind::Query) {
                if x < self.rpos {
                    return None;
                }
                if x < self.rpos + len as i64 {
                    return Some(self.pair_base + (x - self.rpos) as usize);
                }
                self.rpos += len as i64;
            }
            self.pair_base += len as usize;
            self.k += 1;
        }
        None
    }
}

/// Total number of aligned pairs.
fn n_pairs(cigar: &[u32]) -> usize {
    cigar.iter().filter_map(|&w| pair_op(w)).map(|(_, l)| l as usize).sum()
}

/// `aligned_pairs_end_index[j]` of `set_aligned_pairs`: pair index of 0-based reference
/// position `exons[j].e` (the base after the exon), provided the ends of all previous exons
/// were found (`exon_index` only advances on an end hit). `None` = 4.1 IndexError.
fn end_index(cigar: &[u32], pos0: i64, exons: &[Iv], j: usize) -> Option<usize> {
    let mut sk = Seeker::new(cigar, pos0);
    let mut idx = None;
    for ex in exons.get(..=j)? {
        idx = Some(sk.seek(ex.e as i64)?);
    }
    idx
}

/// `aligned_pairs_start_index[j]`: the j-th found pair index of 0-based reference position
/// `exons[k].s` (scanning exons while their predecessors' ends are found), followed by two
/// terminal entries `len(aligned_pairs)`. `None` = 4.1 IndexError.
fn start_index(cigar: &[u32], pos0: i64, exons: &[Iv], j: usize) -> Option<usize> {
    let mut sk = Seeker::new(cigar, pos0);
    let mut found = 0usize;
    for ex in exons {
        if let Some(idx) = sk.seek(ex.s as i64) {
            if found == j {
                return Some(idx);
            }
            found += 1;
        }
        if sk.seek(ex.e as i64).is_none() {
            break;
        }
    }
    (j - found < 2).then(|| n_pairs(cigar))
}

/// `get_error_count(ref_start, ref_end, intron_index, left_site, chr_record)`
/// (`alignment_info.py:60-94`): `(indel_count, mismatch_count)` over the aligned pairs
/// `[idx-12, idx+12)` around the site. `None` when 4.1 would raise IndexError.
fn error_count<R: RefBases + ?Sized>(
    rec: &RecordLite<'_>,
    exons: &[Iv],
    reference: &R,
    ref_start: i64,
    ref_end: i64,
    intron_index: usize,
    left_site: bool,
) -> Option<(u32, u32)> {
    let idx = if left_site {
        end_index(rec.cigar, rec.pos0, exons, intron_index)?
    } else {
        start_index(rec.cigar, rec.pos0, exons, intron_index + 1)?
    };
    let total = n_pairs(rec.cigar);
    let lo = idx.saturating_sub(REGION_TO_CHECK_LEN);
    let hi = (idx + REGION_TO_CHECK_LEN).min(total);
    let check_mm = !reference.is_absent() && !rec.seq.is_empty();

    let (mut indels, mut mms) = (0u32, 0u32);
    let mut last_ref_pos: i64 = 0;
    // walk pairs [lo, hi)
    let (mut pair_base, mut qpos, mut rpos) = (0usize, 0i64, rec.pos0);
    for &w in rec.cigar {
        let Some((kind, len)) = pair_op(w) else { continue };
        let len_u = len as usize;
        let op_end = pair_base + len_u;
        if op_end > lo && pair_base < hi {
            let from = lo.max(pair_base) - pair_base;
            let to = hi.min(op_end) - pair_base;
            for off in from..to {
                let off_i = off as i64;
                match kind {
                    PairKind::Query => {
                        // insertion: previous position inside or adjacent to the region
                        if ref_start - 1 <= last_ref_pos && last_ref_pos <= ref_end {
                            indels += 1;
                        }
                    }
                    PairKind::Ref | PairKind::Match => {
                        let r0 = rpos + off_i;
                        let true_ref = r0 + 1;
                        if ref_start <= true_ref && true_ref <= ref_end {
                            if matches!(kind, PairKind::Ref) {
                                indels += 1;
                            } else if check_mm {
                                let q = rec.seq.get((qpos + off_i) as usize).copied();
                                // pyfaidx returns '' outside the chromosome: always unequal
                                let r = reference.base_at(true_ref as Pos);
                                if r.is_none() || q != r {
                                    mms += 1;
                                }
                            }
                        }
                        last_ref_pos = true_ref;
                        if last_ref_pos > ref_end {
                            return Some((indels, mms));
                        }
                    }
                }
            }
        }
        if !matches!(kind, PairKind::Ref) {
            qpos += len as i64;
        }
        if !matches!(kind, PairKind::Query) {
            rpos += len as i64;
        }
        pair_base = op_end;
        if pair_base >= hi {
            break;
        }
    }
    Some((indels, mms))
}

// ------------------------------------------------------------------------------------------
// ExonCorrector
// ------------------------------------------------------------------------------------------

/// Python list indexing: negative indices wrap; out of range = `None` (IndexError).
#[inline]
fn py_get<T: Copy>(v: &[T], i: i64) -> Option<T> {
    let n = v.len() as i64;
    let j = if i < 0 { i + n } else { i };
    if (0..n).contains(&j) {
        Some(v[j as usize])
    } else {
        None
    }
}

/// `v[a..=b]` with Python `[v[i] for i in range(a, b + 1)]` semantics (indices may wrap;
/// any out-of-range index = IndexError → `None`, and nothing is appended).
fn py_extend(out: &mut Vec<Iv>, v: &[Iv], a: i64, b: i64) -> Option<()> {
    let mark = out.len();
    let mut i = a;
    while i <= b {
        match py_get(v, i) {
            Some(x) => out.push(x),
            None => {
                out.truncate(mark);
                return None;
            }
        }
        i += 1;
    }
    Some(())
}

/// `contains_well_inside(bigger, smaller, delta)`.
#[inline]
fn contains_well_inside(bigger: Iv, smaller: Iv, delta: i32) -> bool {
    bigger.e as i64 >= smaller.e as i64 + delta as i64 && bigger.s as i64 <= smaller.s as i64 - delta as i64
}

/// `event_map` of `correct_misalignments`: a dict where the last event with a key wins.
/// Events per match are few, so lookups scan the event list backwards (no allocation).
struct EventMap<'a> {
    events: &'a [MatchEvent],
    microintron: bool,
}

impl<'a> EventMap<'a> {
    #[inline]
    fn key(&self, e: &MatchEvent) -> Option<i64> {
        if e.read_region == match_consts::UNDEFINED_REGION {
            None
        } else if e.read_region.0 == match_consts::ABSENT_POSITION {
            (e.ty == EventType::FakeMicroIntronRetention && self.microintron).then(|| -e.read_region.1 - 1)
        } else {
            Some(e.read_region.0)
        }
    }
    #[inline]
    fn get(&self, k: i64) -> Option<&'a MatchEvent> {
        self.events.iter().rev().find(|e| self.key(e) == Some(k))
    }
}

/// Event types whose read introns are replaced by their fuzzy-corrected versions.
fn is_known_structure(t: EventType) -> bool {
    use EventType::*;
    matches!(
        t,
        ExtraIntronKnown
            | IntronAlternationKnown
            | IntronMigration
            | ExonSkippingKnown
            | ExonMergeKnown
            | TerminalExonShiftKnown
            | MutuallyExclusiveExonsKnown
            | ExonGainKnown
            | ExonDetachKnown
            | AlternativeStructureKnown
            | AlternativeStructureNovel
    )
}

/// `ExonCorrector.correct_assigned_read` with [`CorrectParams`] taken from `params`.
pub fn correct_read<R: RefBases + ?Sized>(
    ar: &AlignedRead,
    rec: &RecordLite<'_>,
    a: &Assignment,
    index: &IsoformIndex,
    reference: &R,
    params: &Params,
) -> Vec<Iv> {
    correct_read_with(&ar.exons, rec, a, index, reference, &CorrectParams::from_params(params))
}

/// `ExonCorrector.correct_assigned_read(alignment_info, read_assignment)`
/// (`exon_corrector.py:29-41`). `exons` are the read exons after polyA/T trimming
/// (`AlignedRead::exons`); `rec` is the alignment they come from (CIGAR and sequence for
/// `get_error_count`). `index` is the block's reference `GeneInfo`.
///
/// Inputs on which 4.1 raises (IndexError / KeyError) do not occur in practice. Here they
/// fall back conservatively: an unknown transcript returns the exons unchanged, an
/// unresolvable error-count window keeps the reference site, and an event with an
/// out-of-range region contributes no introns.
pub fn correct_read_with<R: RefBases + ?Sized>(
    exons: &[Iv],
    rec: &RecordLite<'_>,
    a: &Assignment,
    index: &IsoformIndex,
    reference: &R,
    p: &CorrectParams,
) -> Vec<Iv> {
    if exons.len() == 1 || a.atype == AssignmentType::Noninformative || a.matches.is_empty() {
        return exons.to_vec();
    }
    let m0 = &a.matches[0];
    let Some(iso) = m0.tx.and_then(|t| index.find_iso(t)) else {
        return exons.to_vec();
    };
    let isoform_region = index.iso[iso].span;
    let isoform_introns = index.iso_introns(iso);
    let read_region = Iv::new(exons[0].s, exons[exons.len() - 1].e);
    let read_introns = iq_types::prim::junctions_from_blocks(exons);
    let event_map = EventMap { events: &m0.events, microintron: p.microintron_retention };

    // --- process_events: fuzzy junction correction
    let corrected_introns: Vec<Iv> = if p.fuzzy_junctions {
        let potential = iq_assign::match_genomic_features(&index.introns, &read_introns, p.delta);
        debug_assert_eq!(potential.len(), read_introns.len());
        read_introns
            .iter()
            .zip(&potential)
            .enumerate()
            .map(|(i, (&ri, &pi))| {
                let site_ok = |start: Pos, end: Pos, left: bool| match error_count(
                    rec,
                    exons,
                    reference,
                    start as i64,
                    end as i64,
                    i,
                    left,
                ) {
                    Some((indel, mm)) => indel == 0 && mm <= 1,
                    None => false,
                };
                let left = if ri.s == pi.s || site_ok(ri.s.min(pi.s), ri.s.max(pi.s) - 1, true) { ri.s } else { pi.s };
                let right =
                    if ri.e == pi.e || site_ok(ri.e.min(pi.e) + 1, ri.e.max(pi.e), false) { ri.e } else { pi.e };
                Iv::new(left, right)
            })
            .collect()
    } else {
        read_introns.clone()
    };

    // --- process_events: event-driven corrections
    let mut crr = read_region;
    let mut new_introns: Vec<Iv> = Vec::with_capacity(corrected_introns.len() + 2);
    let n = corrected_introns.len() as i64;
    let mut i: i64 = 0;
    while i < n {
        if let Some(ev) = event_map.get(-i - 1) {
            // special case for fake IR
            if let Some(x) = py_get(isoform_introns, ev.iso_region.0) {
                new_introns.push(x);
            }
        }
        let Some(ev) = event_map.get(i) else {
            new_introns.push(corrected_introns[i as usize]);
            i += 1;
            continue;
        };
        let misalignment = (p.intron_shifts && ev.ty == EventType::IntronShift)
            || (p.skipped_exons && ev.ty == EventType::ExonMisalignment);
        let (r0, r1) = ev.read_region;
        let (i0, i1) = ev.iso_region;
        match ev.ty {
            EventType::FakeTerminalExonLeft if p.fake_terminal_exons => {
                if let Some(x) = py_get(&read_introns, r0) {
                    crr = Iv::new(x.e + 1, crr.e);
                }
            }
            EventType::FakeTerminalExonRight if p.fake_terminal_exons => {
                if let Some(x) = py_get(&read_introns, r0) {
                    crr = Iv::new(crr.s, x.s - 1);
                }
            }
            EventType::TerminalExonMisalignmentLeft if p.terminal_exons => {
                if let Some(x) = py_get(isoform_introns, i0) {
                    new_introns.push(x);
                }
                crr = Iv::new(isoform_region.s, crr.e);
            }
            EventType::TerminalExonMisalignmentRight if p.terminal_exons => {
                if let Some(x) = py_get(isoform_introns, i0) {
                    new_introns.push(x);
                }
                crr = Iv::new(crr.s, isoform_region.e);
            }
            _ if misalignment
                && matches!(
                    (py_get(isoform_introns, i0), py_get(isoform_introns, i1)),
                    (Some(a), Some(b)) if contains_well_inside(read_region, Iv::new(a.s, b.e), p.delta)
                ) =>
            {
                // misalignments but inside read region
                let _ = py_extend(&mut new_introns, isoform_introns, i0, i1);
            }
            t if is_known_structure(t) => {
                let _ = py_extend(&mut new_introns, &corrected_introns, r0, r1);
            }
            _ => {
                let _ = py_extend(&mut new_introns, &read_introns, r0, r1);
            }
        }
        // 4.1 would loop forever on r1 < i; always make progress
        i = (r1 + 1).max(i + 1);
    }

    // --- correct_assigned_read: introns back to exons
    let mut out = Vec::with_capacity(new_introns.len() + 1);
    match (new_introns.first(), new_introns.last()) {
        (Some(f), Some(l)) => {
            out.push(Iv::new(crr.s, f.s - 1));
            for w in new_introns.windows(2) {
                if w[0].e + 1 < w[1].s {
                    out.push(Iv::new(w[0].e + 1, w[1].s - 1));
                }
            }
            out.push(Iv::new(l.e + 1, crr.e));
        }
        _ => out.push(crr),
    }
    out
}

// ------------------------------------------------------------------------------------------
// Strand
// ------------------------------------------------------------------------------------------

/// `get_intron_strand` on the donor `ref[s..=s+1]` and acceptor `ref[e-1..=e]`
/// dinucleotides (upper-cased; a missing base makes the site non-canonical).
pub fn intron_strand_from_sites(donor: [Option<u8>; 2], acceptor: [Option<u8>; 2]) -> Strand {
    let up = |x: [Option<u8>; 2]| match x {
        [Some(a), Some(b)] => Some([a.to_ascii_uppercase(), b.to_ascii_uppercase()]),
        _ => None,
    };
    let (Some(l), Some(r)) = (up(donor), up(acceptor)) else {
        return Strand::Unknown;
    };
    let is_fwd = matches!((&l, &r), (b"GT", b"AG") | (b"GC", b"AG") | (b"AT", b"AC"));
    let is_rev = matches!((&l, &r), (b"CT", b"AC") | (b"CT", b"GC") | (b"GT", b"AT"));
    match (is_fwd, is_rev) {
        (true, false) => Strand::Plus,
        (false, true) => Strand::Minus,
        _ => Strand::Unknown,
    }
}

/// `get_intron_strand(intron, chr_record)` (`common.py:854-864`).
pub fn intron_strand<R: RefBases + ?Sized>(intron: Iv, reference: &R) -> Strand {
    let b = |p: Pos| reference.base_at(p);
    intron_strand_from_sites([b(intron.s), b(intron.s + 1)], [b(intron.e - 1), b(intron.e)])
}

/// `StrandDetector.get_strand(introns, has_polya, has_polyt)` (`gene_info.py:958-966`).
pub fn strand_from_introns<R: RefBases + ?Sized>(
    introns: impl IntoIterator<Item = Iv>,
    reference: &R,
    has_polya: bool,
    has_polyt: bool,
) -> Strand {
    let (mut fwd, mut rev) = (0u32, 0u32);
    for iv in introns {
        match intron_strand(iv, reference) {
            Strand::Plus => fwd += 1,
            Strand::Minus => rev += 1,
            Strand::Unknown => {}
        }
    }
    if fwd == rev {
        return polya_strand(has_polya, has_polyt);
    }
    if rev < fwd {
        Strand::Plus
    } else {
        Strand::Minus
    }
}

#[inline]
fn polya_strand(has_polya: bool, has_polyt: bool) -> Strand {
    match (has_polya, has_polyt) {
        (true, false) => Strand::Plus,
        (false, true) => Strand::Minus,
        _ => Strand::Unknown,
    }
}

/// Introns of sorted exon blocks (`junctions_from_blocks`) as an iterator.
pub fn introns_of(exons: &[Iv]) -> impl Iterator<Item = Iv> + '_ {
    exons.windows(2).filter(|w| w[0].e + 1 < w[1].s).map(|w| Iv::new(w[0].e + 1, w[1].s - 1))
}

/// `AlignmentCollector.get_assignment_strand` (`alignment_processor.py:538-551`).
/// `n_exons` = `len(read_assignment.exons)` (the read exons, not the corrected ones);
/// `corrected_introns` = `junctions_from_blocks(corrected_exons)` (see [`introns_of`]).
/// `reference` must cover the donor/acceptor dinucleotides of the corrected introns.
pub fn assignment_strand<R: RefBases + ?Sized>(
    a: &Assignment,
    n_exons: usize,
    corrected_introns: impl IntoIterator<Item = Iv>,
    polya: PolyAInfo,
    reference: &R,
) -> Strand {
    if !a.matches.is_empty() && matches!(a.atype, AssignmentType::Unique | AssignmentType::UniqueMinorDifference) {
        return a.matches[0].tx_strand;
    }
    let has_polya = polya.ext_a != -1 || polya.int_a != -1;
    let has_polyt = polya.ext_t != -1 || polya.int_t != -1;
    if n_exons == 1 {
        return polya_strand(has_polya, has_polyt);
    }
    strand_from_introns(corrected_introns, reference, has_polya, has_polyt)
}

/// `AlignmentCollector.check_antisense` (`alignment_processor.py:684-695`). The event is
/// appended to the list directly (not via `add_subclassification`).
pub fn check_antisense(a: &mut Assignment, strand: Strand) {
    for m in a.matches.iter_mut() {
        if m.tx.is_none() || strand == m.tx_strand || strand == Strand::Unknown {
            continue;
        }
        if !matches!(m.class, MatchClass::NovelInCatalog | MatchClass::NovelNotInCatalog) {
            m.class = MatchClass::Antisense;
        }
        m.events.push(MatchEvent::new(EventType::Antisense));
    }
}

// ------------------------------------------------------------------------------------------
// Fixtures
// ------------------------------------------------------------------------------------------

pub(crate) mod fixture {
    use super::*;
    use iq_types::fixture::schema::{
        AssignmentFixture, AssignmentOut, U7AntisenseIn, U7CorrectIn, U7CorrectOut, U7StrandIn, U7StrandOut,
    };
    use iq_types::fixture::{FixResult, FixtureCtx, TxSpace, UnitFixture};

    pub struct U7Correct;
    pub struct U7Strand;
    pub struct U7Antisense;

    impl UnitFixture for U7Correct {
        const UNIT: &'static str = "U7";
        const FN: &'static str = "correct_assigned_read";
        type In = U7CorrectIn;
        type Out = U7CorrectOut;

        fn run(i: &U7CorrectIn, ctx: &FixtureCtx) -> FixResult<U7CorrectOut> {
            let index = ctx.index(&i.index)?;
            let a = i.assignment.to_core(&mut ctx.ids(), TxSpace::Ref);
            let rec = i.aln.to_record();
            let ex = &i.read_exons;
            if ex.is_empty() || (i.read_start, i.read_end) != (ex[0].s, ex[ex.len() - 1].e) {
                return Err("read_start/read_end differ from read_exons".into());
            }
            if i.intron_profile.read_features != iq_types::prim::junctions_from_blocks(ex) {
                return Err("intron_profile.read_features differ from junctions(read_exons)".into());
            }
            let wins: Vec<RefSlice<'_>> = i.ref_.iter().map(|w| w.as_slice()).collect();
            let pp = &i.params;
            let p = CorrectParams {
                delta: pp.delta,
                fuzzy_junctions: pp.correct_fuzzy_junctions,
                intron_shifts: pp.correct_intron_shifts,
                skipped_exons: pp.correct_skipped_exons,
                fake_terminal_exons: pp.correct_fake_terminal_exons,
                terminal_exons: pp.correct_terminal_exons,
                microintron_retention: pp.correct_microintron_retention,
            };
            let corrected = correct_read_with(ex, &rec.as_lite(), &a, &index, &wins[..], &p);
            Ok(U7CorrectOut { corrected_exons: corrected })
        }
    }

    impl UnitFixture for U7Strand {
        const UNIT: &'static str = "U7";
        const FN: &'static str = "get_assignment_strand";
        type In = U7StrandIn;
        type Out = U7StrandOut;

        fn run(i: &U7StrandIn, ctx: &FixtureCtx) -> FixResult<U7StrandOut> {
            let a = i.assignment.to_core(&mut ctx.ids(), TxSpace::Ref);
            // the captured dinucleotides as 2-bp windows at their 1-based positions
            let mut wins: Vec<RefSlice<'_>> = Vec::with_capacity(2 * i.intron_sites.len());
            for (iv, donor, acceptor) in &i.intron_sites {
                wins.push(RefSlice { start: iv.s, seq: donor.as_bytes() });
                wins.push(RefSlice { start: iv.e - 1, seq: acceptor.as_bytes() });
            }
            wins.sort_by_key(|w| w.start);
            let strand = assignment_strand(&a, i.exons.len(), i.corrected_introns.iter().copied(), i.polya, &wins[..]);
            Ok(U7StrandOut { strand })
        }
    }

    impl UnitFixture for U7Antisense {
        const UNIT: &'static str = "U7";
        const FN: &'static str = "check_antisense";
        type In = U7AntisenseIn;
        type Out = AssignmentOut;

        fn run(i: &U7AntisenseIn, ctx: &FixtureCtx) -> FixResult<AssignmentOut> {
            let mut a = i.assignment.to_core(&mut ctx.ids(), TxSpace::Ref);
            check_antisense(&mut a, i.strand);
            Ok(AssignmentOut { assignment: AssignmentFixture::from_core(&a, &ctx.ids()) })
        }
    }
}

pub(crate) fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    reg.add::<fixture::U7Correct>();
    reg.add::<fixture::U7Strand>();
    reg.add::<fixture::U7Antisense>();
}

#[cfg(test)]
mod tests {
    use super::*;
    use iq_types::read::cigar_encode;
    use iq_types::{CigarOp, IsoformMatch};
    use CigarOp::*;

    fn cig(ops: &[(CigarOp, u32)]) -> Vec<u32> {
        ops.iter().map(|&(o, l)| cigar_encode(o, l)).collect()
    }

    /// pysam `get_aligned_pairs()` (matches_only=False): (qpos, rpos0).
    fn naive_pairs(pos0: i64, cigar: &[u32]) -> Vec<(Option<i64>, Option<i64>)> {
        let (mut q, mut r) = (0i64, pos0);
        let mut out = vec![];
        for &w in cigar {
            let (op, len) = iq_types::read::cigar_decode(w);
            for _ in 0..len {
                match op.unwrap() {
                    M | Eq | X => {
                        out.push((Some(q), Some(r)));
                        q += 1;
                        r += 1;
                    }
                    I | S => {
                        out.push((Some(q), None));
                        q += 1;
                    }
                    D | N => {
                        out.push((None, Some(r)));
                        r += 1;
                    }
                    H | P => {}
                }
            }
        }
        out
    }

    /// Line-by-line port of `set_aligned_pairs` + `get_error_count`.
    fn naive_errors(
        pos0: i64,
        cigar: &[u32],
        seq: &[u8],
        exons: &[Iv],
        reference: &RefSlice<'_>,
        rs: i64,
        re: i64,
        ii: usize,
        left: bool,
    ) -> Option<(u32, u32)> {
        let pairs = naive_pairs(pos0, cigar);
        let (mut si, mut ei, mut k) = (vec![], vec![], 0usize);
        for (i, &(_, rp)) in pairs.iter().enumerate() {
            if rp == Some(exons[k].s as i64) {
                si.push(i);
            }
            if rp == Some(exons[k].e as i64) {
                ei.push(i);
                k += 1;
            }
            if k == exons.len() {
                break;
            }
        }
        si.push(pairs.len());
        si.push(pairs.len());
        let idx = if left { *ei.get(ii)? } else { *si.get(ii + 1)? };
        let sel = &pairs[idx.saturating_sub(12)..(idx + 12).min(pairs.len())];
        let (mut ind, mut mm, mut last) = (0, 0, 0i64);
        for &(qp, rp) in sel {
            let t = rp.map(|r| r + 1);
            match t {
                None if rs - 1 <= last && last <= re => ind += 1,
                None => {}
                Some(t) if rs <= t && t <= re => match qp {
                    None => ind += 1,
                    Some(q) => {
                        if Some(seq[q as usize]) != reference.base(t as Pos) {
                            mm += 1
                        }
                    }
                },
                _ => {}
            }
            if let Some(t) = t {
                last = t;
            }
            if last > re {
                break;
            }
        }
        Some((ind, mm))
    }

    #[test]
    fn error_count_equals_naive() {
        let pos0 = 99;
        let c = cig(&[(S, 3), (M, 20), (I, 2), (M, 5), (N, 50), (M, 3), (D, 2), (M, 10), (N, 40), (M, 1)]);
        let exons = crate::cigar::read_blocks(pos0, &c).ref_blocks;
        let qlen: u32 =
            c.iter().filter_map(|&w| pair_op(w)).filter(|(k, _)| !matches!(k, PairKind::Ref)).map(|(_, l)| l).sum();
        let seq: Vec<u8> = (0..qlen).map(|i| b"ACGT"[(i * 7 % 4) as usize]).collect();
        let refseq: Vec<u8> = (0..400).map(|i| b"ACGT"[(i * 3 % 4) as usize]).collect();
        let r = RefSlice { start: 1, seq: &refseq };
        let rec = RecordLite { qname: b"r", flag: 0, pos0, mapq: 60, cigar: &c, seq: &seq, tag: None };
        for ii in 0..exons.len() - 1 {
            for left in [true, false] {
                let site = if left { exons[ii].e + 1 } else { exons[ii + 1].s - 1 } as i64;
                for a in -8..=8 {
                    for b in 0..=8 {
                        let (rs, re) = (site + a, site + a + b);
                        assert_eq!(
                            error_count(&rec, &exons, &r, rs, re, ii, left),
                            naive_errors(pos0, &c, &seq, &exons, &r, rs, re, ii, left),
                            "ii={ii} left={left} rs={rs} re={re}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn terminal_start_index_for_one_base_last_exon() {
        // last exon of 1 bp at the alignment end: its 0-based start is past the end, so
        // start_index[last] is the terminal len(aligned_pairs)
        let c = cig(&[(M, 10), (N, 20), (M, 1)]);
        let exons = crate::cigar::read_blocks(0, &c).ref_blocks;
        assert_eq!(start_index(&c, 0, &exons, 1), Some(31));
        assert_eq!(start_index(&c, 0, &exons, 2), Some(31));
        assert_eq!(start_index(&c, 0, &exons, 3), None);
        assert_eq!(end_index(&c, 0, &exons, 0), Some(10));
        assert_eq!(end_index(&c, 0, &exons, 1), None);
    }

    #[test]
    fn intron_strand_sites() {
        // 1-based: intron 3..10 → donor ref[3..=4], acceptor ref[9..=10]
        let r = RefSlice { start: 1, seq: b"AAgtAAAAagAA" };
        assert_eq!(intron_strand(Iv::new(3, 10), &r), Strand::Plus);
        let r = RefSlice { start: 1, seq: b"AACTAAAAACAA" };
        assert_eq!(intron_strand(Iv::new(3, 10), &r), Strand::Minus);
        // GT..AT is reverse only
        let r = RefSlice { start: 1, seq: b"AAGTAAAAATAA" };
        assert_eq!(intron_strand(Iv::new(3, 10), &r), Strand::Minus);
        // window does not cover the acceptor
        let r = RefSlice { start: 1, seq: b"AAGTAAAA" };
        assert_eq!(intron_strand(Iv::new(3, 10), &r), Strand::Unknown);
        let wins = [RefSlice { start: 3, seq: b"GT" }, RefSlice { start: 9, seq: b"AG" }];
        assert_eq!(intron_strand(Iv::new(3, 10), &wins[..]), Strand::Plus);
    }

    #[test]
    fn strand_votes_and_polya() {
        let r = RefSlice { start: 1, seq: b"AAGTAAAAAGAAAAGTAAAAAGAA" };
        let a = Assignment::empty(AssignmentType::Inconsistent);
        let pa = PolyAInfo { ext_a: 100, ..iq_types::NO_POLYA };
        let introns = [Iv::new(3, 10), Iv::new(15, 22)];
        assert_eq!(assignment_strand(&a, 3, introns, iq_types::NO_POLYA, &r), Strand::Plus);
        // no canonical sites: tie → polyA decides
        assert_eq!(assignment_strand(&a, 2, [Iv::new(1, 5)], pa, &r), Strand::Plus);
        assert_eq!(assignment_strand(&a, 1, [], iq_types::NO_POLYA, &r), Strand::Unknown);
        let pt = PolyAInfo { int_t: 5, ..iq_types::NO_POLYA };
        assert_eq!(assignment_strand(&a, 1, [], pt, &r), Strand::Minus);
    }

    #[test]
    fn antisense() {
        let m = |class, strand, tx| IsoformMatch::new(class, None, tx, [], strand, 0.0);
        let mut a = Assignment::new(
            AssignmentType::Unique,
            [
                m(MatchClass::FullSpliceMatch, Strand::Plus, Some(iq_types::IsoRef::Ref(0))),
                m(MatchClass::NovelInCatalog, Strand::Plus, Some(iq_types::IsoRef::Ref(1))),
                m(MatchClass::FullSpliceMatch, Strand::Plus, None),
                m(MatchClass::FullSpliceMatch, Strand::Minus, Some(iq_types::IsoRef::Ref(2))),
            ]
            .into_iter()
            .collect(),
        );
        check_antisense(&mut a, Strand::Minus);
        assert_eq!(a.matches[0].class, MatchClass::Antisense);
        assert_eq!(a.matches[0].events.last().unwrap().ty, EventType::Antisense);
        assert_eq!(a.matches[1].class, MatchClass::NovelInCatalog);
        assert_eq!(a.matches[1].events.len(), 1);
        assert!(a.matches[2].events.is_empty());
        assert!(a.matches[3].events.is_empty());
    }

    /// Runs every U7 fixture under `$IQRS_FIXTURES` (skipped when unset).
    #[test]
    fn fixtures_u7() {
        let Some(root) = iq_types::fixture::fixtures_root() else { return };
        let mut reg = iq_types::fixture::Registry::new();
        register_fixtures(&mut reg);
        let mut sum = iq_types::fixture::Summary::new(3);
        iq_types::fixture::run_path(&reg, "U7", &root, &mut sum).unwrap();
        println!("{sum}");
        assert!(sum.all_ok());
    }
}
