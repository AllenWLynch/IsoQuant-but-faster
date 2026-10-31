//! U5: `CombinedProfileConstructor` (`assignment/long_read_profiles.py:73-281`) with the
//! windowed sweep `CombinedProfile` (dense) is the output and comparison form.
//!
//! # Windowed sweep (identical outputs to the 4.1 full scan)
//! Known features are sorted by `(start, end)`; ends are not monotone (intron features
//! overlap), so each feature list gets a prefix maximum of ends, `pmax`. Then:
//! - "first feature whose end can be `>= x`" = `pmax.partition_point(|e| e < x)`, and
//!   "features with start `<= y`" = a `partition_point` on starts. A feature overlapping
//!   `[x, y]` lies in that window.
//! - Absence pre-pass (`overlaps_at_least(mapped, f, 20)` / `contains(mapped, f)`): both imply
//!   that `f` overlaps `mapped`, so only the window over `mapped` is tested.
//! - Two-pointer sweeps start at `lo = first end >= read_features[0].s` with `read_pos = 0`.
//!   In 4.1, every feature before `lo` ends before the first read feature, so with
//!   `read_pos == 0` it takes the "feature left of read" branch, which writes nothing.
//!   The 4.1 state on reaching `lo` is therefore exactly `(lo, 0)`. The rest of the sweep
//!   is the 4.1 loop; it stops at the first feature starting after the last read feature,
//!   so it visits only features inside the read's span.
//! - `-2` marks: beyond polyA (`start > a+delta`) is a suffix (binary search). Before polyT
//!   (`end < t-delta`) is a `pmax` prefix plus a scan of features starting before `t-delta`.
//!   Split exons are disjoint, so 4.1's `interval_bin_search(_rev)` reduce to
//!   `partition_point`s.
//! - Every write is non-zero and values never return to 0, so `gene_profile_range` =
//!   (min written index, max written index + 1), or `(K, 0)` when nothing was written.
//!
//! Per read, the cost is O(log K + k) branchy work, where k is the number of features in
//! the windows. The dense output vectors (frozen `CombinedProfile` contract) add an O(K)
//! memset. [`build_profiles`] also validates its per-thread `pmax` cache with an O(K) slice
//! comparison. [`ProfileBuilder`] precomputes `pmax` once per index and skips that check.

use iq_types::fixture::schema::{
    U5ConstructProfilesIn, U5ConstructProfilesOut, U5MatchGenomicFeaturesIn, U5MatchGenomicFeaturesOut,
};
use iq_types::fixture::{FixResult, FixtureCtx, UnitFixture};
use iq_types::params::MatchingParams;
use iq_types::{CombinedProfile, IsoformIndex, Iv, MappedReadProfile, PolyAInfo, Pos};

/// Prefix maximum of ends for one sorted feature list, plus a copy of the list it was
/// built from (the cache key).
#[derive(Debug, Default)]
struct PmaxCache {
    feats: Vec<Iv>,
    pmax: Vec<Pos>,
    valid: bool,
}

impl PmaxCache {
    /// `pmax` for `feats`, recomputed only if `feats` differs from the cached list.
    fn get(&mut self, feats: &[Iv]) -> &[Pos] {
        if !self.valid || self.feats.as_slice() != feats {
            self.feats.clear();
            self.feats.extend_from_slice(feats);
            fill_pmax(feats, &mut self.pmax);
            self.valid = true;
        }
        &self.pmax
    }
}

fn fill_pmax(feats: &[Iv], out: &mut Vec<Pos>) {
    out.clear();
    out.reserve(feats.len());
    let mut m = Pos::MIN;
    for f in feats {
        m = m.max(f.e);
        out.push(m);
    }
}

/// Per-thread scratch buffers reused across reads (and across blocks).
#[derive(Debug, Default)]
pub struct Scratch {
    /// `(read_pos, gene_pos)` matches of the current sweep, in sweep order.
    matches: Vec<(u32, u32)>,
    /// `pmax` caches for intron / exon / split-exon features (`prof_kind` order).
    pmax: [PmaxCache; 3],
    /// U6 assigner buffers (owned by `assigner.rs`).
    pub assign: crate::assigner::AssignScratch,
}

impl Scratch {
    pub fn new() -> Self {
        Self::default()
    }
}

/// Profile parameters used by U5 (the subset of `MatchingParams` that 4.1's
/// `CombinedProfileConstructor` reads).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
struct ProfParams {
    delta: i32,
    minimal_exon_overlap: i32,
    minimal_intron_absence_overlap: i32,
    count_exons: bool,
}

impl ProfParams {
    fn of(p: &MatchingParams) -> Self {
        ProfParams {
            delta: p.delta,
            minimal_exon_overlap: p.minimal_exon_overlap,
            minimal_intron_absence_overlap: p.minimal_intron_absence_overlap,
            count_exons: p.count_exons,
        }
    }
}

/// `fn build_profiles(&IsoformIndex, &Params, &[Iv], PolyAInfo, &mut Scratch) -> CombinedProfile`
/// (`CombinedProfileConstructor(gene_info, params).construct_profiles(blocks, polya, _)`).
pub fn build_profiles(
    index: &IsoformIndex,
    params: &MatchingParams,
    blocks: &[Iv],
    polya: PolyAInfo,
    scratch: &mut Scratch,
) -> CombinedProfile {
    let mut out = CombinedProfile::default();
    build_profiles_into(index, params, blocks, polya, scratch, &mut out);
    out
}

/// [`build_profiles`] writing into `out`, which reuses its vectors' capacity.
pub fn build_profiles_into(
    index: &IsoformIndex,
    params: &MatchingParams,
    blocks: &[Iv],
    polya: PolyAInfo,
    scratch: &mut Scratch,
    out: &mut CombinedProfile,
) {
    let pp = ProfParams::of(params);
    let Scratch { matches, pmax: [pi, pe, ps], .. } = scratch;
    let pmax_intron = pi.get(&index.introns);
    let pmax_split = ps.get(&index.split_exons);
    let pmax_exon: &[Pos] = if pp.count_exons { pe.get(&index.exons) } else { &[] };
    build_with(index, pp, [pmax_intron, pmax_exon, pmax_split], blocks, polya, matches, out);
}

/// A profile constructor bound to one index (4.1 builds one `CombinedProfileConstructor`
/// per `GeneInfo`). It precomputes `pmax` once, so each read costs O(log K + k) plus the
/// dense output fill.
#[derive(Debug)]
pub struct ProfileBuilder<'a> {
    index: &'a IsoformIndex,
    params: ProfParams,
    pmax: [Vec<Pos>; 3],
}

impl<'a> ProfileBuilder<'a> {
    pub fn new(index: &'a IsoformIndex, params: &MatchingParams) -> Self {
        let params = ProfParams::of(params);
        let mut pmax: [Vec<Pos>; 3] = Default::default();
        fill_pmax(&index.introns, &mut pmax[0]);
        if params.count_exons {
            fill_pmax(&index.exons, &mut pmax[1]);
        }
        fill_pmax(&index.split_exons, &mut pmax[2]);
        ProfileBuilder { index, params, pmax }
    }

    /// Same result as [`build_profiles`].
    pub fn build(&self, blocks: &[Iv], polya: PolyAInfo, scratch: &mut Scratch) -> CombinedProfile {
        let mut out = CombinedProfile::default();
        self.build_into(blocks, polya, scratch, &mut out);
        out
    }

    pub fn build_into(&self, blocks: &[Iv], polya: PolyAInfo, scratch: &mut Scratch, out: &mut CombinedProfile) {
        let pm = [self.pmax[0].as_slice(), self.pmax[1].as_slice(), self.pmax[2].as_slice()];
        build_with(self.index, self.params, pm, blocks, polya, &mut scratch.matches, out);
    }
}

fn build_with(
    index: &IsoformIndex,
    pp: ProfParams,
    pmax: [&[Pos]; 3],
    blocks: &[Iv],
    polya: PolyAInfo,
    matches: &mut Vec<(u32, u32)>,
    out: &mut CombinedProfile,
) {
    let gene_region = Iv::new(index.start, index.end);
    let (ext_a, ext_t) = (polya.ext_a, polya.ext_t);
    out.polya = polya;

    // intron profile (`construct_intron_profile`)
    let intron = &mut out.intron;
    junctions_into(blocks, &mut intron.read_features);
    let mapped = match (blocks.first(), blocks.last()) {
        (Some(f), Some(l)) => Iv::new(f.s, l.e),
        _ => Iv::new(0, 0),
    };
    let ofv = OverlapView { feats: &index.introns, pmax: pmax[0] };
    overlapping_profile(
        ofv,
        gene_region,
        Absence::OverlapsAtLeast(pp.minimal_intron_absence_overlap),
        pp.delta,
        mapped,
        ext_a,
        ext_t,
        matches,
        intron,
    );

    // exon profile (`construct_exon_profile`), only with `count_exons`
    if pp.count_exons {
        let ex = out.exon.get_or_insert_with(MappedReadProfile::default);
        ex.read_features.clear();
        ex.read_features.extend_from_slice(blocks);
        let mapped = match (blocks.first(), blocks.last()) {
            (Some(f), Some(l)) => Iv::new(f.e.wrapping_add(pp.delta), l.s.wrapping_sub(pp.delta)),
            _ => Iv::new(0, 0),
        };
        let efv = OverlapView { feats: &index.exons, pmax: pmax[1] };
        overlapping_profile(efv, gene_region, Absence::Contains, pp.delta, mapped, ext_a, ext_t, matches, ex);
    } else {
        out.exon = None;
    }

    // split exon profile (`NonOverlappingFeaturesProfileConstructor.construct_profile`)
    let sp = &mut out.split_exon;
    sp.read_features.clear();
    sp.read_features.extend_from_slice(blocks);
    let sfv = OverlapView { feats: &index.split_exons, pmax: pmax[2] };
    nonoverlapping_profile(sfv, pp.minimal_exon_overlap, pp.delta, ext_a, ext_t, sp);
}

/// `junctions_from_blocks` (`common.py:462-468`) into a reused buffer.
fn junctions_into(blocks: &[Iv], out: &mut Vec<Iv>) {
    out.clear();
    for w in blocks.windows(2) {
        if (w[0].e as i64) + 1 < w[1].s as i64 {
            out.push(Iv::new(w[0].e + 1, w[1].s - 1));
        }
    }
}

// ------------------------------------------------------------------------------------------
// predicates (`common.py`), 64-bit arithmetic like Python ints
// ------------------------------------------------------------------------------------------

#[inline]
fn equal_ranges(a: Iv, b: Iv, delta: i32) -> bool {
    (a.s as i64 - b.s as i64).abs() <= delta as i64 && (a.e as i64 - b.e as i64).abs() <= delta as i64
}

/// `overlaps_at_least(range1, range2, delta)` (`common.py:199-208`).
#[inline]
fn overlaps_at_least(r1: Iv, r2: Iv, delta: i32) -> bool {
    let ovlp1 = r1.e as i64 - r2.s as i64;
    let ovlp2 = r2.e as i64 - r1.s as i64;
    if ovlp1 < 0 || ovlp2 < 0 {
        return false;
    }
    let d = delta as i64 - 1;
    if r1.e < r2.e {
        ovlp1 >= d || r1.s >= r2.s
    } else {
        ovlp2 >= d || r1.s <= r2.s
    }
}

/// `overlaps_at_least_when_overlap(range1, range2, delta)` (`common.py:212-216`).
#[inline]
fn overlaps_at_least_when_overlap(r1: Iv, r2: Iv, delta: i32) -> bool {
    if r1.e < r2.e {
        r1.s >= r2.s || r1.e as i64 - r2.s as i64 + 1 >= delta as i64
    } else {
        r1.s <= r2.s || r2.e as i64 - r1.s as i64 + 1 >= delta as i64
    }
}

#[inline]
fn match_delta(a: Iv, b: Iv) -> i64 {
    (a.s as i64 - b.s as i64).abs() + (a.e as i64 - b.e as i64).abs()
}

#[derive(Copy, Clone, Debug)]
enum Absence {
    /// `partial(overlaps_at_least, delta=..)`.
    OverlapsAtLeast(i32),
    /// `contains`.
    Contains,
}

impl Absence {
    #[inline]
    fn test(self, big: Iv, f: Iv) -> bool {
        match self {
            Absence::OverlapsAtLeast(d) => overlaps_at_least(big, f, d),
            Absence::Contains => big.contains(&f),
        }
    }
}

/// Sorted features with their prefix-max ends.
#[derive(Copy, Clone)]
struct OverlapView<'a> {
    feats: &'a [Iv],
    pmax: &'a [Pos],
}

impl OverlapView<'_> {
    /// First index whose prefix-max end is `>= x` (all features before it end `< x`).
    #[inline]
    fn first_end_ge(&self, x: i64) -> usize {
        self.pmax.partition_point(|&e| (e as i64) < x)
    }
    /// Number of features with start `<= y`.
    #[inline]
    fn start_le(&self, y: i64) -> usize {
        self.feats.partition_point(|f| (f.s as i64) <= y)
    }
    /// Sweep start for a read whose first feature is `r0` (see module docs). Falls back
    /// to 0 for an inverted `r0`, where the skip argument does not hold.
    #[inline]
    fn sweep_start(&self, r0: Option<&Iv>) -> usize {
        match r0 {
            Some(r) if r.s <= r.e => self.first_end_ge(r.s as i64),
            _ => 0,
        }
    }
}

/// Tracks the min/max index written with a non-zero value.
struct Written {
    lo: usize,
    hi: usize,
}

impl Written {
    #[inline]
    fn new() -> Self {
        Written { lo: usize::MAX, hi: 0 }
    }
    #[inline]
    fn at(&mut self, i: usize) {
        self.lo = self.lo.min(i);
        self.hi = self.hi.max(i + 1);
    }
    #[inline]
    fn range(&mut self, a: usize, b: usize) {
        if a < b {
            self.lo = self.lo.min(a);
            self.hi = self.hi.max(b);
        }
    }
    /// `(start_pos, end_pos + 1)` of 4.1.
    fn profile_range(&self, k: usize) -> (i32, i32) {
        if self.lo == usize::MAX {
            (k as i32, 0)
        } else {
            (self.lo as i32, self.hi as i32)
        }
    }
}

/// `OverlappingFeaturesProfileConstructor.construct_profile_for_features` (`:111-184`),
/// windowed. `out.read_features` holds the read features on entry.
fn overlapping_profile(
    fv: OverlapView<'_>,
    gene_region: Iv,
    absence: Absence,
    delta: i32,
    mapped: Iv,
    polya: Pos,
    polyt: Pos,
    matches: &mut Vec<(u32, u32)>,
    out: &mut MappedReadProfile,
) {
    let kf = fv.feats;
    let k = kf.len();
    let rf = &out.read_features;
    let n = rf.len();
    let gp = &mut out.gene_profile;
    gp.clear();
    gp.resize(k, 0);
    let rp = &mut out.read_profile;
    rp.clear();
    rp.extend(rf.iter().map(|&r| if absence.test(gene_region, r) { -1 } else { 0 }));
    let mut w = Written::new();

    // absence pre-pass over the window of features overlapping `mapped`
    let a = fv.first_end_ge(mapped.s as i64);
    let b = fv.start_le(mapped.e as i64);
    for i in a..b.max(a) {
        if absence.test(mapped, kf[i]) {
            gp[i] = -1;
            w.at(i);
        }
    }

    // two-pointer sweep (4.1 loop, entered at the equivalent state)
    matches.clear();
    let mut g = fv.sweep_start(rf.first());
    let mut r = 0usize;
    while g < k && r < n {
        let rfe = rf[r];
        let kfe = kf[g];
        if rfe.e < kfe.s {
            if rp[r] == 0 && g > 0 {
                rp[r] = -1;
            }
            r += 1;
        } else if kfe.e < rfe.s {
            if r > 0 {
                gp[g] = -1;
                w.at(g);
            }
            g += 1;
        } else if equal_ranges(rfe, kfe, delta) {
            gp[g] = 1;
            rp[r] = 1;
            w.at(g);
            matches.push((r as u32, g as u32));
            g += 1;
        } else {
            // overlapping (the two disjoint cases are handled above)
            if absence.test(mapped, kfe) {
                gp[g] = -1;
                w.at(g);
            }
            g += 1;
        }
    }

    // eliminating non-unique matches: keep those with minimal delta
    let mut i = 0;
    while i < matches.len() {
        let r0 = matches[i].0;
        let mut j = i + 1;
        while j < matches.len() && matches[j].0 == r0 {
            j += 1;
        }
        if j - i > 1 {
            let rfe = rf[r0 as usize];
            let best = matches[i..j].iter().map(|&(_, g)| match_delta(rfe, kf[g as usize])).min().unwrap();
            for &(_, g) in &matches[i..j] {
                if match_delta(rfe, kf[g as usize]) > best {
                    gp[g as usize] = -1;
                }
            }
        }
        i = j;
    }
    // (4.1's `corrected_features` is computed but not returned: dead.)

    // beyond polyA: start > polya + delta (a suffix)
    if polya != -1 {
        let x = polya as i64 + delta as i64;
        let s = fv.start_le(x);
        gp[s..].fill(-2);
        w.range(s, k);
    }
    // before polyT: end < polyt - delta (a pmax prefix plus features starting before x)
    if polyt != -1 {
        let x = polyt as i64 - delta as i64;
        let p = fv.first_end_ge(x);
        gp[..p].fill(-2);
        w.range(0, p);
        let q = kf.partition_point(|f| (f.s as i64) < x);
        for i in p..q.max(p) {
            if (kf[i].e as i64) < x {
                gp[i] = -2;
                w.at(i);
            }
        }
    }

    out.gene_profile_range = w.profile_range(k);
}

/// `NonOverlappingFeaturesProfileConstructor.construct_profile` (`:194-244`), windowed.
/// `out.read_features` holds the sorted blocks on entry.
fn nonoverlapping_profile(
    fv: OverlapView<'_>,
    min_overlap: i32,
    delta: i32,
    polya: Pos,
    polyt: Pos,
    out: &mut MappedReadProfile,
) {
    let kf = fv.feats;
    let k = kf.len();
    let rf = &out.read_features;
    let n = rf.len();
    let ep = &mut out.gene_profile;
    ep.clear();
    ep.resize(k, 0);
    let rp = &mut out.read_profile;
    rp.clear();
    rp.resize(n, 0);
    let mut w = Written::new();

    let mut g = fv.sweep_start(rf.first());
    let mut r = 0usize;
    while g < k && r < n {
        let re = rf[r];
        let ge = kf[g];
        if re.e < ge.s {
            if g > 0 && rp[r] == 0 {
                rp[r] = -1;
            }
            r += 1;
        } else if ge.e < re.s {
            if r > 0 && ep[g] == 0 {
                ep[g] = -1;
                w.at(g);
            }
            g += 1;
        } else {
            if overlaps_at_least_when_overlap(re, ge, min_overlap) {
                ep[g] = 1;
                rp[r] = 1;
                w.at(g);
            }
            if re.e < ge.e {
                r += 1;
            } else {
                g += 1;
            }
        }
    }

    if polya != -1 {
        if let Some(ix) = interval_bin_search(kf, polya as i64 + delta as i64) {
            ep[ix + 1..].fill(-2);
            w.range(ix + 1, k);
        }
    }
    if polyt != -1 {
        if let Some(ix) = interval_bin_search_rev(fv, polyt as i64 - delta as i64) {
            ep[..ix].fill(-2);
            w.range(0, ix);
        }
    }

    out.gene_profile_range = w.profile_range(k);
}

/// `interval_bin_search` (`common.py:887-905`): the `ind` with
/// `iv[ind].s <= pos < iv[ind+1].s`, which for sorted starts is unique.
fn interval_bin_search(iv: &[Iv], pos: i64) -> Option<usize> {
    let (first, last) = (iv.first()?, iv.last()?);
    if pos > last.e as i64 || pos < first.s as i64 {
        return None;
    }
    if pos >= last.s as i64 {
        return Some(iv.len() - 1);
    }
    Some(iv.partition_point(|f| f.s as i64 <= pos) - 1)
}

/// `interval_bin_search_rev` (`common.py:908-927`): the `ind` with
/// `iv[ind-1].e < pos <= iv[ind].e`. For disjoint sorted intervals, ends are sorted and
/// `pmax` equals them.
fn interval_bin_search_rev(fv: OverlapView<'_>, pos: i64) -> Option<usize> {
    let iv = fv.feats;
    let (first, last) = (iv.first()?, iv.last()?);
    if pos > last.e as i64 || pos < first.s as i64 {
        return None;
    }
    if pos <= first.e as i64 {
        return Some(0);
    }
    Some(fv.first_end_ge(pos))
}

// ------------------------------------------------------------------------------------------
// match_genomic_features
// ------------------------------------------------------------------------------------------

/// `OverlappingFeaturesProfileConstructor.match_genomic_features` (`:73-109`) with
/// comparator `equal_ranges(delta)`, used by U7. Read features are replaced by their
/// matched known feature (minimal `match_delta`, first in order on ties).
///
/// O(n log K + matches): for each read feature, the 4.1 loop advances over every known
/// feature starting at or before the read feature's end. It then keeps advancing over
/// trailing known features that are still `equal_ranges` (those can start after the read
/// end). Matches can only start within `delta` of the read feature's start.
pub fn match_genomic_features(known_features: &[Iv], read_features: &[Iv], delta: i32) -> Vec<Iv> {
    let kf = known_features;
    let k = kf.len();
    let d = delta as i64;
    let mut out = Vec::with_capacity(read_features.len());
    let mut g = 0usize;
    for (r, &rfe) in read_features.iter().enumerate() {
        if g >= k {
            out.extend_from_slice(&read_features[r..]);
            break;
        }
        let mut best: Option<(i64, usize)> = None;
        let mut consider = |gi: usize| {
            let md = match_delta(rfe, kf[gi]);
            if best.is_none_or(|(b, _)| md < b) {
                best = Some((md, gi));
            }
        };
        if rfe.s > rfe.e {
            // inverted read feature: plain 4.1 loop for this read position
            while g < k {
                let kfe = kf[g];
                if equal_ranges(rfe, kfe, delta) {
                    consider(g);
                } else if !rfe.overlaps(&kfe) && rfe.e < kfe.s {
                    break; // left_of: next read position
                }
                g += 1;
            }
        } else {
            let pe = kf.partition_point(|f| f.s as i64 <= rfe.e as i64);
            if g < pe {
                let lo = g.max(kf.partition_point(|f| (f.s as i64) < rfe.s as i64 - d));
                let hi = pe.min(kf.partition_point(|f| f.s as i64 <= rfe.s as i64 + d));
                for (gi, &f) in kf.iter().enumerate().take(hi).skip(lo) {
                    if equal_ranges(rfe, f, delta) {
                        consider(gi);
                    }
                }
                g = pe;
            }
            while g < k && equal_ranges(rfe, kf[g], delta) {
                consider(g);
                g += 1;
            }
        }
        out.push(match best {
            Some((_, gi)) => kf[gi],
            None => rfe,
        });
    }
    out
}

// ------------------------------------------------------------------------------------------
// Dense reference port (tests, benchmark)
// ------------------------------------------------------------------------------------------

/// Straight line-by-line port of the 4.1 constructors (full scans). Used only to check
/// and benchmark the windowed implementation.
#[doc(hidden)]
pub mod dense {
    use super::*;

    #[allow(clippy::too_many_arguments)]
    fn overlapping(
        known: &[Iv],
        gene_region: Iv,
        absence: Absence,
        delta: i32,
        read_features: &[Iv],
        mapped: Iv,
        polya: Pos,
        polyt: Pos,
    ) -> MappedReadProfile {
        let mut read_profile = vec![0i8; read_features.len()];
        let mut intron_profile = vec![0i8; known.len()];
        let mut matched: Vec<Vec<usize>> = vec![Vec::new(); read_features.len()];
        for i in 0..known.len() {
            if absence.test(mapped, known[i]) {
                intron_profile[i] = -1;
            }
        }
        for i in 0..read_features.len() {
            if absence.test(gene_region, read_features[i]) {
                read_profile[i] = -1;
            }
        }
        let (mut g, mut r) = (0, 0);
        while g < known.len() && r < read_features.len() {
            let rf = read_features[r];
            let kf = known[g];
            if rf.e < kf.s {
                if read_profile[r] == 0 && g > 0 {
                    read_profile[r] = -1;
                }
                r += 1;
            } else if kf.e < rf.s {
                if r > 0 {
                    intron_profile[g] = -1;
                }
                g += 1;
            } else if equal_ranges(rf, kf, delta) {
                intron_profile[g] = 1;
                read_profile[r] = 1;
                matched[r].push(g);
                g += 1;
            } else if rf.overlaps(&kf) {
                if absence.test(mapped, kf) {
                    intron_profile[g] = -1;
                }
                g += 1;
            }
        }
        for (r, m) in matched.iter().enumerate() {
            if m.len() > 1 {
                let ds: Vec<i64> = m.iter().map(|&g| match_delta(read_features[r], known[g])).collect();
                let best = *ds.iter().min().unwrap();
                for (i, &g) in m.iter().enumerate() {
                    if ds[i] > best {
                        intron_profile[g] = -1;
                    }
                }
            }
        }
        if polya != -1 {
            for i in 0..known.len() {
                if known[i].s as i64 > polya as i64 + delta as i64 {
                    intron_profile[i] = -2;
                }
            }
        }
        if polyt != -1 {
            for i in 0..known.len() {
                if (known[i].e as i64) < polyt as i64 - delta as i64 {
                    intron_profile[i] = -2;
                }
            }
        }
        let range = dense_range(&intron_profile);
        MappedReadProfile {
            gene_profile: intron_profile,
            read_profile,
            read_features: read_features.to_vec(),
            gene_profile_range: range,
        }
    }

    fn dense_range(p: &[i8]) -> (i32, i32) {
        let mut s = 0usize;
        while s < p.len() && p[s] == 0 {
            s += 1;
        }
        let mut e = p.len() as i64 - 1;
        while e >= 0 && p[e as usize] == 0 {
            e -= 1;
        }
        (s as i32, (e + 1) as i32)
    }

    /// Python list indexing (negative indices wrap).
    fn py(iv: &[Iv], i: i64) -> Iv {
        iv[if i < 0 { (iv.len() as i64 + i) as usize } else { i as usize }]
    }

    /// 4.1 `interval_bin_search` loop, verbatim.
    fn bin_search(iv: &[Iv], pos: i64) -> i64 {
        if iv.is_empty() || pos > iv[iv.len() - 1].e as i64 || pos < iv[0].s as i64 {
            return -1;
        }
        let s = iv.len() as i64 - 1;
        if pos >= iv[s as usize].s as i64 {
            return s;
        }
        let (mut ind, mut step) = (s / 2, s / 2);
        while !((py(iv, ind).s as i64) <= pos && pos < py(iv, ind + 1).s as i64) {
            step = (step / 2).max(1);
            if pos < py(iv, ind).s as i64 {
                ind -= step;
            } else {
                ind += step;
            }
        }
        ind
    }

    /// 4.1 `interval_bin_search_rev` loop, verbatim.
    fn bin_search_rev(iv: &[Iv], pos: i64) -> i64 {
        if iv.is_empty() || pos > iv[iv.len() - 1].e as i64 || pos < iv[0].s as i64 {
            return -1;
        }
        if pos <= iv[0].e as i64 {
            return 0;
        }
        let s = iv.len() as i64 - 1;
        let (mut ind, mut step) = (s / 2, s / 2);
        while !((py(iv, ind - 1).e as i64) < pos && pos <= py(iv, ind).e as i64) {
            step = (step / 2).max(1);
            if pos > py(iv, ind).e as i64 {
                ind += step;
            } else {
                ind -= step;
            }
        }
        ind
    }

    fn nonoverlapping(
        known: &[Iv],
        min_overlap: i32,
        delta: i32,
        blocks: &[Iv],
        polya: Pos,
        polyt: Pos,
    ) -> MappedReadProfile {
        let mut ep = vec![0i8; known.len()];
        let mut rp = vec![0i8; blocks.len()];
        let (mut g, mut r) = (0, 0);
        while g < known.len() && r < blocks.len() {
            let re = blocks[r];
            let ge = known[g];
            if re.e < ge.s {
                if g > 0 && rp[r] == 0 {
                    rp[r] = -1;
                }
                r += 1;
            } else if ge.e < re.s {
                if r > 0 && ep[g] == 0 {
                    ep[g] = -1;
                }
                g += 1;
            } else if overlaps_at_least_when_overlap(re, ge, min_overlap) {
                ep[g] = 1;
                rp[r] = 1;
                if re.e < ge.e {
                    r += 1;
                } else {
                    g += 1;
                }
            } else if re.e < ge.e {
                r += 1;
            } else {
                g += 1;
            }
        }
        if polya != -1 {
            let ix = bin_search(known, polya as i64 + delta as i64);
            if ix != -1 {
                for v in &mut ep[(ix + 1) as usize..] {
                    *v = -2;
                }
            }
        }
        if polyt != -1 {
            let ix = bin_search_rev(known, polyt as i64 - delta as i64);
            if ix != -1 {
                for v in &mut ep[..ix as usize] {
                    *v = -2;
                }
            }
        }
        let range = dense_range(&ep);
        MappedReadProfile {
            gene_profile: ep,
            read_profile: rp,
            read_features: blocks.to_vec(),
            gene_profile_range: range,
        }
    }

    /// Dense `construct_profiles`.
    pub fn build_profiles(
        index: &IsoformIndex,
        params: &MatchingParams,
        blocks: &[Iv],
        polya: PolyAInfo,
    ) -> CombinedProfile {
        let pp = ProfParams::of(params);
        let gene_region = Iv::new(index.start, index.end);
        let mut introns = Vec::new();
        junctions_into(blocks, &mut introns);
        let mapped = Iv::new(blocks[0].s, blocks[blocks.len() - 1].e);
        let intron = overlapping(
            &index.introns,
            gene_region,
            Absence::OverlapsAtLeast(pp.minimal_intron_absence_overlap),
            pp.delta,
            &introns,
            mapped,
            polya.ext_a,
            polya.ext_t,
        );
        let exon = pp.count_exons.then(|| {
            let mapped = Iv::new(blocks[0].e + pp.delta, blocks[blocks.len() - 1].s - pp.delta);
            overlapping(
                &index.exons,
                gene_region,
                Absence::Contains,
                pp.delta,
                blocks,
                mapped,
                polya.ext_a,
                polya.ext_t,
            )
        });
        let split_exon =
            nonoverlapping(&index.split_exons, pp.minimal_exon_overlap, pp.delta, blocks, polya.ext_a, polya.ext_t);
        CombinedProfile { intron, exon, split_exon, polya }
    }

    /// Dense `match_genomic_features`.
    pub fn match_genomic_features(known: &[Iv], read: &[Iv], delta: i32) -> Vec<Iv> {
        let mut matched: Vec<Vec<usize>> = vec![Vec::new(); read.len()];
        let (mut g, mut r) = (0, 0);
        while g < known.len() && r < read.len() {
            if equal_ranges(read[r], known[g], delta) {
                matched[r].push(g);
                g += 1;
            } else if read[r].overlaps(&known[g]) {
                g += 1;
            } else if read[r].e < known[g].s {
                r += 1;
            } else {
                g += 1;
            }
        }
        read.iter()
            .enumerate()
            .map(|(r, &rf)| {
                let m = &matched[r];
                if m.is_empty() {
                    return rf;
                }
                let best = m.iter().map(|&g| match_delta(rf, known[g])).min().unwrap();
                known[*m.iter().find(|&&g| match_delta(rf, known[g]) == best).unwrap()]
            })
            .collect()
    }
}

// ------------------------------------------------------------------------------------------
// Fixture runners
// ------------------------------------------------------------------------------------------

/// U5 `construct_profiles`.
pub struct U5ConstructProfiles;

impl UnitFixture for U5ConstructProfiles {
    const UNIT: &'static str = "U5";
    const FN: &'static str = "construct_profiles";
    type In = U5ConstructProfilesIn;
    type Out = U5ConstructProfilesOut;

    fn run(i: &U5ConstructProfilesIn, ctx: &FixtureCtx) -> FixResult<U5ConstructProfilesOut> {
        let index = ctx.index(&i.index)?;
        let mut params = MatchingParams::precise();
        params.delta = i.params.delta;
        params.minimal_exon_overlap = i.params.minimal_exon_overlap;
        params.minimal_intron_absence_overlap = i.params.minimal_intron_absence_overlap;
        params.count_exons = i.params.count_exons;
        let mut scratch = Scratch::default();
        let p = build_profiles(&index, &params, &i.blocks, i.polya, &mut scratch);
        debug_assert_eq!(p, dense::build_profiles(&index, &params, &i.blocks, i.polya));
        Ok(U5ConstructProfilesOut { intron: p.intron, exon: p.exon, split_exon: p.split_exon })
    }
}

/// U5 `match_genomic_features`.
pub struct U5MatchGenomicFeatures;

impl UnitFixture for U5MatchGenomicFeatures {
    const UNIT: &'static str = "U5";
    const FN: &'static str = "match_genomic_features";
    type In = U5MatchGenomicFeaturesIn;
    type Out = U5MatchGenomicFeaturesOut;

    fn run(i: &U5MatchGenomicFeaturesIn, ctx: &FixtureCtx) -> FixResult<U5MatchGenomicFeaturesOut> {
        // Schema v1 bug: `in.delta` records the constructor's `self.delta` (default 0 in
        // `exon_corrector.py:24-27`), not the comparator's `equal_ranges` delta
        // (`params.delta`). The exon corrector is built with the run's params, so the
        // comparator delta is `params.json`'s `delta`. Remove once the harness fix
        // is re-captured.
        let delta = ctx.params()?.matching.delta;
        let features = match_genomic_features(&i.known_features, &i.read_features, delta);
        debug_assert_eq!(features, dense::match_genomic_features(&i.known_features, &i.read_features, delta));
        Ok(U5MatchGenomicFeaturesOut { features })
    }
}

/// Registers the U5 fixture runners.
pub fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    reg.add::<U5ConstructProfiles>();
    reg.add::<U5MatchGenomicFeatures>();
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn iv(s: Pos, e: Pos) -> Iv {
        Iv::new(s, e)
    }

    fn params(delta: i32, count_exons: bool) -> MatchingParams {
        let mut p = MatchingParams::precise();
        p.delta = delta;
        p.count_exons = count_exons;
        p
    }

    fn index(start: Pos, end: Pos, introns: Vec<Iv>, exons: Vec<Iv>, split: Vec<Iv>) -> IsoformIndex {
        IsoformIndex { start, end, introns, exons, split_exons: split, ..Default::default() }
    }

    /// Sorted (by `(s, e)`) possibly overlapping features.
    fn arb_feats(max: usize) -> impl Strategy<Value = Vec<Iv>> {
        prop::collection::vec((0i32..400, 0i32..80), 0..max).prop_map(|v| {
            let mut f: Vec<Iv> = v.into_iter().map(|(s, l)| iv(s, s + l)).collect();
            f.sort();
            f.dedup();
            f
        })
    }

    /// Sorted disjoint intervals from (gap, len) pairs.
    fn arb_disjoint(min: usize, max: usize, min_gap: i32) -> impl Strategy<Value = Vec<Iv>> {
        (0i32..60, prop::collection::vec((min_gap..40, 0i32..50), min..max)).prop_map(|(start, v)| {
            let mut out = Vec::new();
            let mut pos = start;
            for (gap, len) in v {
                let s = pos + gap;
                out.push(iv(s, s + len));
                pos = s + len + 1;
            }
            out
        })
    }

    fn arb_pos() -> impl Strategy<Value = Pos> {
        prop_oneof![Just(-1), 0i32..500]
    }

    proptest! {
        #![proptest_config(ProptestConfig { cases: 4000, ..ProptestConfig::default() })]

        #[test]
        fn windowed_equals_dense(
            introns in arb_feats(40),
            exons in arb_feats(40),
            split in arb_disjoint(0, 30, 1),
            blocks in arb_disjoint(1, 12, 0),
            gs in 0i32..200, gl in 0i32..400,
            delta in 0i32..10, min_ov in 0i32..30, min_abs in 0i32..30,
            count_exons: bool,
            a in arb_pos(), t in arb_pos(),
        ) {
            let idx = index(gs, gs + gl, introns, exons, split);
            let mut p = params(delta, count_exons);
            p.minimal_exon_overlap = min_ov;
            p.minimal_intron_absence_overlap = min_abs;
            let polya = PolyAInfo { ext_a: a, ext_t: t, int_a: -1, int_t: -1 };
            let want = dense::build_profiles(&idx, &p, &blocks, polya);
            let mut scratch = Scratch::default();
            prop_assert_eq!(&build_profiles(&idx, &p, &blocks, polya, &mut scratch), &want);
            // cached pmax path and the builder
            prop_assert_eq!(&build_profiles(&idx, &p, &blocks, polya, &mut scratch), &want);
            prop_assert_eq!(&ProfileBuilder::new(&idx, &p).build(&blocks, polya, &mut scratch), &want);
        }

        #[test]
        fn match_features_equals_dense(
            known in arb_feats(40),
            read in prop_oneof![arb_disjoint(0, 12, 0), arb_feats(12)],
            delta in 0i32..12,
        ) {
            prop_assert_eq!(match_genomic_features(&known, &read, delta), dense::match_genomic_features(&known, &read, delta));
        }
    }

    #[test]
    fn scratch_cache_invalidates_across_indices() {
        let p = params(4, false);
        let a = index(0, 1000, vec![iv(100, 200), iv(300, 400)], vec![], vec![iv(50, 99), iv(201, 299), iv(401, 500)]);
        let b = index(0, 1000, vec![iv(100, 250)], vec![], vec![iv(50, 99), iv(251, 500)]);
        let blocks = [iv(60, 99), iv(201, 299), iv(401, 450)];
        let mut s = Scratch::default();
        for idx in [&a, &b, &a] {
            assert_eq!(
                build_profiles(idx, &p, &blocks, NO_POLYA_INFO, &mut s),
                dense::build_profiles(idx, &p, &blocks, NO_POLYA_INFO)
            );
        }
    }

    const NO_POLYA_INFO: PolyAInfo = iq_types::NO_POLYA;

    #[test]
    fn mono_exon_read_and_empty_index() {
        let p = params(4, true);
        let empty = index(0, 0, vec![], vec![], vec![]);
        let r = build_profiles(&empty, &p, &[iv(10, 20)], NO_POLYA_INFO, &mut Scratch::default());
        assert_eq!(r.intron.gene_profile_range, (0, 0));
        assert!(r.intron.read_features.is_empty());
        assert_eq!(r.split_exon.read_profile, vec![0]);
        // no blocks: no panic, empty profiles
        let idx = index(0, 1000, vec![iv(100, 200)], vec![], vec![iv(50, 99)]);
        let r = build_profiles(&idx, &p, &[], NO_POLYA_INFO, &mut Scratch::default());
        assert_eq!(r.intron.gene_profile, vec![0]);
        assert_eq!(r.intron.gene_profile_range, (1, 0));
    }

    #[test]
    fn polyt_before_overlapping_features() {
        // ends are not monotone: (10, 500) ends after polyT, (20, 30) before it
        let idx = index(0, 1000, vec![iv(10, 500), iv(20, 30), iv(40, 600)], vec![], vec![iv(1, 9)]);
        let p = params(0, false);
        let polya = PolyAInfo { ext_a: -1, ext_t: 100, int_a: -1, int_t: -1 };
        let blocks = [iv(100, 200), iv(400, 700)];
        let got = build_profiles(&idx, &p, &blocks, polya, &mut Scratch::default());
        assert_eq!(got, dense::build_profiles(&idx, &p, &blocks, polya));
        assert_eq!(got.intron.gene_profile[1], -2);
        assert_ne!(got.intron.gene_profile[0], -2);
    }

    #[test]
    fn match_after_read_end_within_delta() {
        // equal_ranges can hold for a known feature starting after the read feature's end
        let known = [iv(5, 6), iv(13, 14), iv(20, 30)];
        let read = [iv(10, 11), iv(21, 29)];
        assert_eq!(match_genomic_features(&known, &read, 4), vec![iv(13, 14), iv(20, 30)]);
        assert_eq!(match_genomic_features(&known, &read, 4), dense::match_genomic_features(&known, &read, 4));
    }

    /// Benchmark on a real fixture block: `cargo test --release -p iq_assign -- --ignored bench_profiles --nocapture`.
    #[test]
    #[ignore]
    fn bench_profiles() {
        use iq_types::fixture::{read_cases, FixtureCtx};
        use iq_types::serde_json;
        use std::collections::HashMap;
        use std::time::Instant;
        // `IQ_BENCH_CAP` = capture dir; default: the multigene slice of SRR32588720.
        let cap = std::env::var("IQ_BENCH_CAP").map(std::path::PathBuf::from).unwrap_or_else(|_| {
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../fixtures/SRR32588720/multigene")
        });
        let ctx = FixtureCtx::new(&cap);
        let mut by_index: HashMap<String, Vec<(Vec<Iv>, PolyAInfo, MatchingParams)>> = HashMap::new();
        for c in read_cases::<serde_json::Value, serde_json::Value>(&cap.join("units/U5.jsonl.gz")).unwrap() {
            let c = c.unwrap();
            if c.input["fn"] != "construct_profiles" {
                continue;
            }
            let i: U5ConstructProfilesIn = serde_json::from_value(c.input).unwrap();
            let mut p = MatchingParams::precise();
            p.delta = i.params.delta;
            by_index.entry(i.index.clone()).or_default().push((i.blocks, i.polya, p));
        }
        // the index with the most features among those with >= 100 reads
        let (sha, reads, idx) = by_index
            .into_iter()
            .filter(|(_, v)| v.len() >= 100)
            .map(|(k, v)| {
                let idx = ctx.index(&k).unwrap();
                (k, v, idx)
            })
            .max_by_key(|(_, _, idx)| idx.introns.len() + idx.split_exons.len())
            .expect("no block with enough reads");
        let n = reads.len();
        let reps = (200_000 / n).max(1);
        let mut s = Scratch::default();
        let mut out = CombinedProfile::default();
        let time = |f: &mut dyn FnMut()| {
            let t = Instant::now();
            for _ in 0..reps {
                f();
            }
            t.elapsed().as_secs_f64() * 1e6 / (reps * n) as f64
        };
        let mut sink = 0usize;
        let t_dense = time(&mut || {
            for (b, a, p) in &reads {
                sink += dense::build_profiles(&idx, p, b, *a).intron.gene_profile_range.0 as usize;
            }
        });
        let t_fast = time(&mut || {
            for (b, a, p) in &reads {
                sink += build_profiles(&idx, p, b, *a, &mut s).intron.gene_profile_range.0 as usize;
            }
        });
        let t_into = time(&mut || {
            for (b, a, p) in &reads {
                build_profiles_into(&idx, p, b, *a, &mut s, &mut out);
                sink += out.intron.gene_profile_range.0 as usize;
            }
        });
        let builder = ProfileBuilder::new(&idx, &reads[0].2);
        let t_builder = time(&mut || {
            for (b, a, _) in &reads {
                builder.build_into(b, *a, &mut s, &mut out);
                sink += out.intron.gene_profile_range.0 as usize;
            }
        });
        println!(
            "index {sha}: {n} reads, K = {} introns / {} split exons; µs/read: dense {t_dense:.3}, build_profiles {t_fast:.3}, build_profiles_into {t_into:.3}, ProfileBuilder::build_into {t_builder:.3} (sink {sink})",
            idx.introns.len(),
            idx.split_exons.len()
        );
    }
}
