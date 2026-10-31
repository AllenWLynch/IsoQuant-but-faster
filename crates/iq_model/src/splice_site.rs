//! U16: model splice-site correction (`model_filter.py:143-192`,
//! `transcript_splice_site_corrector.py:103-237`).
//!
//! For every exon boundary of a novel model that falls inside a bound read's span, the
//! number of deleted reference bases the read carries in the 8-bp window next to the
//! boundary (on the exon side) is histogrammed. If at least 5 reads span the boundary, the
//! histogram has a unique mode of 3-6 bases, and shifting the boundary by that mode (into
//! the exon) puts a canonical dinucleotide next to it, the boundary is moved.
//!
//! Order: 4.1 keys the per-site accumulators by position in a dict and returns the
//! corrected locations as a list that is immediately turned into a set; the mode is
//! required to be unique (ties → -1, never accepted). So neither dict nor set order can
//! reach the output, and none is emulated here.

use iq_types::smallvec::SmallVec;
use iq_types::{BlockStore, Iv, ModelKind, ModelStore, Pos, RefSlice, Strand, TranscriptModel};

/// `WINDOW_SIZE`.
pub const WINDOW_SIZE: Pos = 8;
/// `MIN_N_OF_ALIGNED_READS`.
pub const MIN_N_OF_ALIGNED_READS: u32 = 5;
/// `ACCEPTED_DEL_CASES` (3, 4, 5, 6).
pub const ACCEPTED_DEL_CASES: std::ops::RangeInclusive<i32> = 3..=6;

/// What U16 needs from each bound read.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SiteRead<'a> {
    /// `exons[0][0]`, `exons[-1][1]` (raw `ReadAssignment.exons`).
    pub start: Pos,
    pub end: Pos,
    /// U4 boundary deletions: `(1-based first deleted base, length)`.
    pub boundary_deletions: &'a [(Pos, u32)],
}

/// `SpliceSiteCase`: histogram `deleted bases → #reads` (insertion order irrelevant, see
/// module doc).
#[derive(Default)]
struct SiteCase {
    location_is_end: bool,
    seen: bool,
    hist: SmallVec<[(u32, u32); 4]>,
}

impl SiteCase {
    fn add(&mut self, count: u32) {
        match self.hist.iter_mut().find(|(k, _)| *k == count) {
            Some((_, n)) => *n += 1,
            None => self.hist.push((count, 1)),
        }
    }

    /// `compute_most_common_case_of_deletions`: signed unique mode, -1 if not unique.
    fn most_common_del(&self) -> i32 {
        let max = self.hist.iter().map(|&(_, n)| n).max().unwrap_or(0);
        let mut modes = self.hist.iter().filter(|&&(_, n)| n == max);
        match (modes.next(), modes.next()) {
            (Some(&(k, _)), None) => {
                let k = k as i32;
                if self.location_is_end {
                    -k
                } else {
                    k
                }
            }
            _ => -1,
        }
    }
}

/// `_deleted_bases_in_window`: deleted bases within the closed window `[low, high]`.
#[inline]
fn deleted_bases_in_window(dels: &[(Pos, u32)], low: Pos, high: Pos) -> u32 {
    let mut total = 0u32;
    for &(s, len) in dels {
        let e = s as i64 + len as i64 - 1;
        let overlap = (high as i64).min(e) - (low as i64).max(s as i64) + 1;
        if overlap > 0 {
            total += overlap as u32;
        }
    }
    total
}

/// `_corrected_site_is_canonical` (upper-cased dinucleotide next to the shifted boundary).
fn corrected_site_is_canonical(
    reference: &RefSlice<'_>,
    boundary: Pos,
    shift: i32,
    is_end: bool,
    strand: Strand,
) -> bool {
    let corrected = boundary + shift;
    // 0-based start of the dinucleotide; 1-based first base = start0 + 1.
    let start0 = if is_end { corrected } else { corrected - 3 };
    if start0 < 0 {
        return false;
    }
    // A slice clipped by the chromosome end is shorter than 2 → not canonical.
    let (Some(a), Some(b)) = (reference.base(start0 + 1), reference.base(start0 + 2)) else {
        return false;
    };
    let site = [a.to_ascii_uppercase(), b.to_ascii_uppercase()];
    let valid: &[&[u8; 2]] = match (strand == Strand::Plus, is_end) {
        (true, true) => &[b"GT", b"GC", b"AT"],   // _FWD_LEFT_SITES
        (true, false) => &[b"AG", b"AC"],         // _FWD_RIGHT_SITES
        (false, true) => &[b"CT", b"GT"],         // _REV_LEFT_SITES
        (false, false) => &[b"AC", b"GC", b"AT"], // _REV_RIGHT_SITES
    };
    valid.iter().any(|v| **v == site)
}

/// `_correct_one_model_splice_sites` on raw exons/strand: the corrected exon list (possibly
/// equal to `exons`) when at least one site qualifies, else `None`.
///
/// `reference` must cover every base 4.1 could read: positions within `[b-9, b+8]` of each
/// boundary `b` (bases outside it count as "not canonical", like a pyfaidx slice clipped at
/// the chromosome end).
pub fn correct_exon_sites<'r>(
    exons: &[Iv],
    strand: Strand,
    reads: impl IntoIterator<Item = SiteRead<'r>>,
    reference: &RefSlice<'_>,
) -> Option<Vec<Iv>> {
    if exons.len() < 2 || !matches!(strand, Strand::Plus | Strand::Minus) {
        return None;
    }
    // Boundary k: exon k/2, start (even) or end (odd). Cases are keyed by position (dict);
    // `canon[k]` is the first boundary index with the same position.
    let nb = exons.len() * 2;
    let pos = |k: usize| if k & 1 == 0 { exons[k / 2].s } else { exons[k / 2].e };
    let canon: Vec<usize> = (0..nb).map(|k| (0..=k).find(|&j| pos(j) == pos(k)).unwrap()).collect();
    let mut cases: Vec<SiteCase> = (0..nb).map(|_| SiteCase::default()).collect();
    let mut any_read = false;
    let mut any_case = false;
    for r in reads {
        any_read = true;
        // `extract_splice_site_locations_within_aligned_read`
        for (i, ex) in exons.iter().enumerate() {
            for (k, loc, is_end) in [(2 * i, ex.s, false), (2 * i + 1, ex.e, true)] {
                if r.start <= loc && loc <= r.end {
                    let (low, high) = if is_end { (loc - WINDOW_SIZE + 1, loc) } else { (loc, loc + WINDOW_SIZE - 1) };
                    let n = deleted_bases_in_window(r.boundary_deletions, low, high);
                    let c = &mut cases[canon[k]];
                    if !c.seen {
                        c.seen = true;
                        c.location_is_end = is_end;
                        any_case = true;
                    }
                    c.add(n);
                }
            }
            if r.end <= ex.e {
                break;
            }
        }
    }
    if !any_read || !any_case {
        return None;
    }
    // `correct_splice_site_errors`: shift per corrected boundary index (canonical ones only).
    let mut shift: Vec<Option<i32>> = vec![None; nb];
    let mut any_corr = false;
    for (k, c) in cases.iter().enumerate() {
        if !c.seen || c.hist.iter().map(|&(_, n)| n).sum::<u32>() < MIN_N_OF_ALIGNED_READS {
            continue;
        }
        let d = c.most_common_del();
        if !ACCEPTED_DEL_CASES.contains(&d.abs()) {
            continue;
        }
        if corrected_site_is_canonical(reference, pos(k), d, c.location_is_end, strand) {
            shift[k] = Some(d);
            any_corr = true;
        }
    }
    if !any_corr {
        return None;
    }
    // `generate_updated_exon_list`
    Some(
        exons
            .iter()
            .enumerate()
            .map(|(i, ex)| {
                let s = ex.s + shift[canon[2 * i]].unwrap_or(0);
                let e = ex.e + shift[canon[2 * i + 1]].unwrap_or(0);
                Iv { s, e }
            })
            .collect(),
    )
}

/// `fn correct_model_sites(&TranscriptModel, reads, &RefSlice) -> Option<Vec<Iv>>`
/// (`_correct_one_model_splice_sites`). Reads with no exons must be left out by the caller.
pub fn correct_model_sites(
    model: &TranscriptModel,
    reads: &[SiteRead<'_>],
    reference: &RefSlice<'_>,
) -> Option<Vec<Iv>> {
    correct_exon_sites(&model.exons, model.strand, reads.iter().cloned(), reference)
}

/// `_valid_exon_chain`: no inverted exon, introns keep at least one base.
pub fn valid_exon_chain(exons: &[Iv]) -> bool {
    exons.iter().enumerate().all(|(i, x)| x.s <= x.e && (i == 0 || x.s > exons[i - 1].e))
}

/// `correct_transcript_splice_sites` (without the `--no_splice_site_correction` check, which
/// the caller does via `params.splice_site_correction`): corrects every non-known
/// model in storage order using its bound reads (`transcript_read_ids`, duplicates count
/// twice as in 4.1). `reference` must cover all model boundaries ±9.
pub fn correct_transcript_splice_sites(ms: &mut ModelStore, reads: &BlockStore, reference: &RefSlice<'_>) {
    let ModelStore { models, transcript_read_ids, .. } = ms;
    for m in models.iter_mut() {
        if matches!(m.kind, ModelKind::Known(_)) {
            continue;
        }
        let ids = transcript_read_ids.get(&m.key).map(|v| v.as_slice()).unwrap_or(&[]);
        let site_reads = ids.iter().filter_map(|&r| {
            let v = reads.view(r as usize);
            let (start, end) = v.exons_span();
            // `if not read.exons: continue` (exons are 1-based, so (0, 0) = none)
            if start == 0 && end == 0 {
                return None;
            }
            Some(SiteRead { start, end, boundary_deletions: v.boundary_deletions() })
        });
        if let Some(c) = correct_exon_sites(&m.exons, m.strand, site_reads, reference) {
            if c != m.exons && valid_exon_chain(&c) {
                m.exons = c;
            }
        }
    }
}

/// Registers the U16 fixture runners.
pub(crate) fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    reg.add::<fixture::U16CorrectOne>();
    reg.add::<fixture::U16CorrectStore>();
}

pub(crate) mod fixture {
    //! U16 runners. `_correct_one_model_splice_sites` carries one reference window per
    //! boundary; they are stitched into one [`RefSlice`] (gaps = `N`, never canonical).
    //! The store-level `correct_transcript_splice_sites` has no windows, so it reads the
    //! capture's genome like U13 (`data/ref/genome.fa`, `$IQRS_REFERENCE` overrides).
    use super::*;
    use crate::fixture_util::{compare_stores, StageInput};
    use iq_types::fixture::schema::{RefWin, StageStoreIn, StoreOut, U16CorrectOneIn, U16CorrectOneOut};
    use iq_types::fixture::{FixError, FixResult, FixtureCtx, UnitFixture, Verdict};
    use std::io::{Read, Seek, SeekFrom};
    use std::path::PathBuf;

    /// Stitches windows into one buffer starting at the smallest window start.
    pub(crate) fn stitch(wins: &[RefWin]) -> (Pos, Vec<u8>) {
        let Some(lo) = wins.iter().map(|w| w.start0).min() else {
            return (1, Vec::new());
        };
        let hi = wins.iter().map(|w| w.start0 + w.seq.len() as i64).max().unwrap();
        let mut buf = vec![b'N'; (hi - lo) as usize];
        for w in wins {
            let o = (w.start0 - lo) as usize;
            buf[o..o + w.seq.len()].copy_from_slice(w.seq.as_bytes());
        }
        ((lo + 1) as Pos, buf)
    }

    pub struct U16CorrectOne;

    impl UnitFixture for U16CorrectOne {
        const UNIT: &'static str = "U16";
        const FN: &'static str = "_correct_one_model_splice_sites";
        type In = U16CorrectOneIn;
        type Out = U16CorrectOneOut;

        fn run(i: &U16CorrectOneIn, _ctx: &FixtureCtx) -> FixResult<U16CorrectOneOut> {
            let (start, seq) = stitch(&i.ref_);
            let reference = RefSlice { start, seq: &seq };
            // The hook records reads with exons only (start/end come from exons[0]/[-1]).
            let reads = i.reads.iter().map(|r| SiteRead {
                start: r.start,
                end: r.end,
                boundary_deletions: &r.boundary_deletions,
            });
            Ok(U16CorrectOneOut { exons: correct_exon_sites(&i.model.exons, i.model.strand, reads, &reference) })
        }
    }

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
        Err(FixError("U16: reference genome not found (set IQRS_REFERENCE)".into()))
    }

    /// Reads `chr:[lo, hi]` (1-based closed, clipped) from an faidx-indexed FASTA.
    /// (Same as the U13 runner's private helper.)
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

    pub struct U16CorrectStore;

    impl UnitFixture for U16CorrectStore {
        const UNIT: &'static str = "U16";
        const FN: &'static str = "correct_transcript_splice_sites";
        type In = StageStoreIn;
        type Out = StoreOut;

        fn run(input: &StageStoreIn, ctx: &FixtureCtx) -> FixResult<StoreOut> {
            let si = StageInput::load(input, ctx)?;
            let mut ms = si.store.clone();
            if si.params.splice_site_correction {
                let (mut lo, mut hi) = (Pos::MAX, Pos::MIN);
                for m in ms.models.iter().filter(|m| !matches!(m.kind, ModelKind::Known(_))) {
                    lo = lo.min(m.start());
                    hi = hi.max(m.end());
                }
                let (lo, seq) = if lo <= hi {
                    let lo = (lo - 16).max(1);
                    (lo, read_region(&genome_path(ctx)?, &input.block.chr, lo, hi + 16)?)
                } else {
                    (1, Vec::new())
                };
                let reference = RefSlice { start: lo, seq: &seq };
                correct_transcript_splice_sites(&mut ms, &si.block.store, &reference);
            }
            Ok(StoreOut { store: si.store_out(&ms, &si.detected_known) })
        }

        fn compare(got: &StoreOut, want: &StoreOut, ctx: &FixtureCtx) -> Verdict {
            compare_stores(&got.store, &want.store, ctx)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rd(start: Pos, end: Pos, dels: &[(Pos, u32)]) -> SiteRead<'_> {
        SiteRead { start, end, boundary_deletions: dels }
    }

    #[test]
    fn window_overlap() {
        assert_eq!(deleted_bases_in_window(&[(10, 3)], 11, 18), 2);
        assert_eq!(deleted_bases_in_window(&[(10, 3), (15, 10)], 11, 18), 6);
        assert_eq!(deleted_bases_in_window(&[(19, 3)], 11, 18), 0);
    }

    #[test]
    fn mode_must_be_unique() {
        let mut c = SiteCase { location_is_end: true, seen: true, ..Default::default() };
        for n in [4, 4, 0, 0] {
            c.add(n);
        }
        assert_eq!(c.most_common_del(), -1);
        c.add(4);
        assert_eq!(c.most_common_del(), -4);
        c.location_is_end = false;
        assert_eq!(c.most_common_del(), 4);
    }

    /// Exon end at 20 with a 4-bp deletion [17, 20] in 5 reads → shift to 16 if the
    /// reference has GT at 17-18.
    #[test]
    fn shifts_exon_end_onto_canonical() {
        let mut seq = vec![b'A'; 100];
        seq[16] = b'g'; // pos 17 (start 1)
        seq[17] = b't'; // pos 18
        let r = RefSlice { start: 1, seq: &seq };
        let exons = [Iv::new(1, 20), Iv::new(50, 80)];
        let dels = [(17, 4)];
        let reads: Vec<_> = (0..5).map(|_| rd(5, 70, &dels)).collect();
        let got = correct_exon_sites(&exons, Strand::Plus, reads.iter().cloned(), &r);
        assert_eq!(got, Some(vec![Iv::new(1, 16), Iv::new(50, 80)]));
        // 4 reads only → not enough support
        assert_eq!(correct_exon_sites(&exons, Strand::Plus, reads[..4].iter().cloned(), &r), None);
        // '-' strand: GT is a valid REV left site too
        assert!(correct_exon_sites(&exons, Strand::Minus, reads.iter().cloned(), &r).is_some());
        assert_eq!(correct_exon_sites(&exons, Strand::Unknown, reads.iter().cloned(), &r), None);
    }

    /// Exon start at 50 with deletion [50, 52] (3 bp) → shift to 53 if AG at 51-52.
    #[test]
    fn shifts_exon_start() {
        let mut seq = vec![b'C'; 100];
        seq[50] = b'A'; // pos 51
        seq[51] = b'G'; // pos 52
        let r = RefSlice { start: 1, seq: &seq };
        let exons = [Iv::new(1, 20), Iv::new(50, 80)];
        let dels = [(50, 3)];
        let reads: Vec<_> = (0..6).map(|_| rd(5, 70, &dels)).collect();
        assert_eq!(
            correct_exon_sites(&exons, Strand::Plus, reads.iter().cloned(), &r),
            Some(vec![Iv::new(1, 20), Iv::new(53, 80)])
        );
        // window outside the reference → not canonical
        let short = RefSlice { start: 60, seq: &seq[59..] };
        assert_eq!(correct_exon_sites(&exons, Strand::Plus, reads.iter().cloned(), &short), None);
    }

    #[test]
    fn exon_chain_validity() {
        assert!(valid_exon_chain(&[Iv::new(1, 5), Iv::new(7, 9)]));
        assert!(!valid_exon_chain(&[Iv::new(1, 5), Iv::new(5, 9)]));
        assert!(!valid_exon_chain(&[Iv::new(6, 5)]));
    }
}
