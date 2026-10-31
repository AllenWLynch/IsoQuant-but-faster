//! U10a: `_histogram_record`, `_rank_peaks`, `_window_count`, `detect_peaks`
//! (`terminal_prediction/terminal_peaks.py:84-214`), numpy pairwise sum, scipy
//! `find_peaks(distance=10)` / `peak_widths(rel_height=0.98)`; D-T1 tie order.
//!
//! numpy 2.4 / scipy 1.17 semantics ported explicitly:
//! - `np.var(data)`: integer mean (exact), `d = x - mean`, `pairwise_sum(d*d) / n`;
//! - `stats.skew(data)` (bias=True): `m3 / m2**1.5` with `m3 = pairwise_sum((d*d)*d) / n`,
//!   NaN when `m2 <= (eps*mean)**2`;
//! - `stats.entropy(padded)`: `pairwise_sum(entr(c / total))`;
//! - `stats.mode`: smallest most frequent value;
//! - `find_peaks(padded, distance=10)`: `_local_maxima_1d` + `_select_by_peak_distance`;
//! - `peak_widths(padded, peaks, rel_height=0.98)`: `_peak_prominences` (wlen=None) +
//!   `_peak_widths`.
//!
//! **D-T1:** scipy's distance filter visits peaks in `np.argsort(heights)` order (from the
//! end) and IsoQuant's ranking uses `np.argsort(-heights)`. iqrs emulates numpy's default
//! (non-SIMD) argsort exactly, see [`np_argsort`]: introsort with an insertion sort for
//! partitions of ≤ 16 elements. On small peak sets that sort is stable, so the distance filter
//! keeps the *rightmost* of equal peaks closer than `distance`, and the ranking keeps equal
//! heights in position order. This equals 4.1 on CPUs without x86-simd-sort; on SIMD CPUs 4.1
//! can still differ for large tied arrays, so the fixture comparator (fixtures were captured on
//! a SIMD CPU) accepts any 4.1 output consistent with some tie order as `Explained(D-T1)`.
//! Deterministic on every CPU.

use crate::xgb::{sigmoid, XgbModel, N_FEATURES};
use iq_types::{Peak, Pos};

/// Position → count, in insertion (dict) order — the order fixes float summation order.
pub type Histogram = iq_types::indexmap::IndexMap<Pos, u32>;

/// `PEAK_DISTANCE` (find_peaks `distance`).
pub const PEAK_DISTANCE: usize = 10;
/// `PEAK_REL_HEIGHT` (peak_widths `rel_height`).
pub const PEAK_REL_HEIGHT: f64 = 0.98;
/// `HISTOGRAM_PAD` zeros on each side of the histogram.
pub const HISTOGRAM_PAD: usize = 10;

/// numpy's `DOUBLE_pairwise_sum` as used by `np.add.reduce` on a contiguous float64 array
/// (`0.0 + pairwise(a)`): n < 8 sequential from -0.0; n ≤ 128 eight strided accumulators
/// combined `((r0+r1)+(r2+r3))+((r4+r5)+(r6+r7))` plus a sequential tail; otherwise split at
/// `n/2 - (n/2 mod 8)`.
pub fn np_sum(a: &[f64]) -> f64 {
    0.0 + pairwise(a)
}

fn pairwise(a: &[f64]) -> f64 {
    let n = a.len();
    if n < 8 {
        let mut r = -0.0;
        for &x in a {
            r += x;
        }
        r
    } else if n <= 128 {
        let mut r = [0.0f64; 8];
        r.copy_from_slice(&a[..8]);
        let mut i = 8;
        while i < n - n % 8 {
            for j in 0..8 {
                r[j] += a[i + j];
            }
            i += 8;
        }
        let mut res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        while i < n {
            res += a[i];
            i += 1;
        }
        res
    } else {
        let mut n2 = n / 2;
        n2 -= n2 % 8;
        pairwise(&a[..n2]) + pairwise(&a[n2..])
    }
}

/// One `find_peaks` peak after `peak_widths` (`record['raw_peaks']` column), start-relative.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RawPeak {
    /// `peak_location` (padded index − 10).
    pub location: i64,
    /// `peak_width` (`widths[0]`).
    pub width: f64,
    /// `peak_left` = `int(left_ips − 10)` (truncation).
    pub left: i64,
    /// `peak_right` = `int(right_ips − 10)`.
    pub right: i64,
    /// `peak_heights` = `padded[location + 10]`.
    pub height: u32,
}

/// `_histogram_record` for one histogram.
#[derive(Clone, Debug, PartialEq)]
pub struct HistRecord {
    pub start: Pos,
    /// Start-relative mode.
    pub mode: i64,
    pub var: f64,
    pub skew: f64,
    pub entropy: f64,
    pub mean_height: f64,
    /// `histogram`: 10 zeros + per-position counts over `[start, max]` + 10 zeros.
    pub padded: Vec<u32>,
    /// `raw_peaks`, in position order (empty = `peak_count == 0`).
    pub peaks: Vec<RawPeak>,
}

/// `_histogram_record(position_counts)`. `None` when no count is positive (4.1 raises on
/// `min([])` there; callers treat it like an empty histogram).
pub fn histogram_record(hist: &Histogram) -> Option<HistRecord> {
    histogram_record_ties(hist, Ties::Numpy)
}

/// Order in which equal heights are visited (distance filter) and ranked. iqrs uses
/// [`Ties::Numpy`] (D-T1); the position orders exist so tests and the D-T1 comparators can
/// build alternative 4.1 (SIMD-CPU) outputs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Ties {
    /// numpy's default argsort ([`np_argsort`]): production.
    Numpy,
    /// Equal heights by ascending position (distance filter keeps the leftmost).
    PositionAsc,
    /// Equal heights by descending position.
    PositionDesc,
}

fn tie_cmp(t: Ties, a: usize, b: usize) -> std::cmp::Ordering {
    match t {
        Ties::PositionAsc | Ties::Numpy => a.cmp(&b),
        Ties::PositionDesc => b.cmp(&a),
    }
}

/// numpy's default `np.argsort(v)` (`kind='quicksort'` = `aquicksort_` in
/// `npysort/quicksort.cpp`, non-SIMD path): median-of-3 introsort, partitions of ≤ 16
/// elements finished by insertion sort (stable there). The depth-limit heapsort fallback
/// (`2 * floor(log2 n)` levels) is ported too, so the result matches for every input.
pub fn np_argsort(v: &[i64]) -> Vec<usize> {
    let n = v.len();
    let mut t: Vec<usize> = (0..n).collect();
    if n < 2 {
        return t;
    }
    const SMALL: isize = 15; // numpy SMALL_QUICKSORT
    let msb = |x: usize| (usize::BITS - 1 - x.leading_zeros()) as i32;
    let mut stack: Vec<(isize, isize, i32)> = Vec::new();
    let (mut pl, mut pr) = (0isize, n as isize - 1);
    let mut cdepth = msb(n) * 2;
    loop {
        if cdepth < 0 {
            np_aheapsort(v, &mut t[pl as usize..=pr as usize]);
        } else {
            while pr - pl > SMALL {
                let pm = pl + ((pr - pl) >> 1);
                let (l, m, r) = (pl as usize, pm as usize, pr as usize);
                if v[t[m]] < v[t[l]] {
                    t.swap(m, l);
                }
                if v[t[r]] < v[t[m]] {
                    t.swap(r, m);
                }
                if v[t[m]] < v[t[l]] {
                    t.swap(m, l);
                }
                let vp = v[t[m]];
                let mut pi = pl;
                let mut pj = pr - 1;
                t.swap(m, pj as usize);
                loop {
                    loop {
                        pi += 1;
                        if v[t[pi as usize]] >= vp {
                            break;
                        }
                    }
                    loop {
                        pj -= 1;
                        if vp >= v[t[pj as usize]] {
                            break;
                        }
                    }
                    if pi >= pj {
                        break;
                    }
                    t.swap(pi as usize, pj as usize);
                }
                t.swap(pi as usize, (pr - 1) as usize);
                cdepth -= 1;
                if pi - pl < pr - pi {
                    stack.push((pi + 1, pr, cdepth));
                    pr = pi - 1;
                } else {
                    stack.push((pl, pi - 1, cdepth));
                    pl = pi + 1;
                }
            }
            // insertion sort (numpy checks the depth limit only after a stack pop)
            let mut pi = pl + 1;
            while pi <= pr {
                let vi = t[pi as usize];
                let vp = v[vi];
                let mut pj = pi;
                while pj > pl && vp < v[t[(pj - 1) as usize]] {
                    t[pj as usize] = t[(pj - 1) as usize];
                    pj -= 1;
                }
                t[pj as usize] = vi;
                pi += 1;
            }
        }
        match stack.pop() {
            Some((a, b, d)) => {
                pl = a;
                pr = b;
                cdepth = d;
            }
            None => break,
        }
    }
    t
}

/// numpy `aheapsort_` on an index slice (1-based heap, as in `npysort/heapsort.cpp`).
fn np_aheapsort(v: &[i64], t: &mut [usize]) {
    let n = t.len();
    if n < 2 {
        return;
    }
    // a[1..=n] view: a(k) = t[k-1]
    let mut l = n >> 1;
    while l > 0 {
        let tmp = t[l - 1];
        let (mut i, mut j) = (l, l << 1);
        while j <= n {
            if j < n && v[t[j - 1]] < v[t[j]] {
                j += 1;
            }
            if v[tmp] < v[t[j - 1]] {
                t[i - 1] = t[j - 1];
                i = j;
                j += j;
            } else {
                break;
            }
        }
        t[i - 1] = tmp;
        l -= 1;
    }
    let mut m = n;
    while m > 1 {
        let tmp = t[m - 1];
        t[m - 1] = t[0];
        m -= 1;
        let (mut i, mut j) = (1usize, 2usize);
        while j <= m {
            if j < m && v[t[j - 1]] < v[t[j]] {
                j += 1;
            }
            if v[tmp] < v[t[j - 1]] {
                t[i - 1] = t[j - 1];
                i = j;
                j += j;
            } else {
                break;
            }
        }
        t[i - 1] = tmp;
    }
}

fn histogram_record_ties(hist: &Histogram, ties: Ties) -> Option<HistRecord> {
    let mut lo = i64::MAX;
    let mut hi = i64::MIN;
    let mut n: u64 = 0;
    let mut sum: i128 = 0;
    for (&p, &c) in hist {
        if c > 0 {
            lo = lo.min(p as i64);
            hi = hi.max(p as i64);
            n += c as u64;
            sum += p as i128 * c as i128;
        }
    }
    if n == 0 {
        return None;
    }
    let span = (hi - lo) as usize + 1;
    let mut padded = vec![0u32; span + 2 * HISTOGRAM_PAD];
    for (&p, &c) in hist {
        if c > 0 {
            padded[(p as i64 - lo) as usize + HISTOGRAM_PAD] += c;
        }
    }
    let nf = n as f64;
    // np.mean of an int64 array: float64 sum of integers (exact below 2**53) / n.
    let mean = sum as f64 / nf;
    let mut sq = Vec::with_capacity(n as usize);
    for (&p, &c) in hist {
        if c > 0 {
            let d = p as f64 - mean;
            sq.extend(std::iter::repeat_n(d * d, c as usize));
        }
    }
    let var = np_sum(&sq) / nf;
    let skew = if n > 3 && var >= 1e-12 {
        // scipy: m2 = mean(d**2) (== var, same ops), m3 = mean((d*d)*d).
        let m2 = var;
        let mut cube = sq;
        let mut k = 0;
        for (&p, &c) in hist {
            if c > 0 {
                let d = p as f64 - mean;
                for v in &mut cube[k..k + c as usize] {
                    *v *= d;
                }
                k += c as usize;
            }
        }
        let m3 = np_sum(&cube) / nf;
        let eps_mean = f64::EPSILON * mean;
        if m2 <= eps_mean * eps_mean {
            f64::NAN
        } else {
            m3 / m2.powf(1.5)
        }
    } else {
        f64::NAN
    };
    let ent: Vec<f64> = padded
        .iter()
        .map(|&c| {
            let p = c as f64 / nf;
            if p > 0.0 {
                -p * p.ln()
            } else {
                0.0
            }
        })
        .collect();
    let entropy = np_sum(&ent);
    let mean_height = nf / padded.len() as f64;
    // stats.mode: smallest most frequent value = first maximum of the histogram.
    let mut mode = 0usize;
    for (i, &c) in padded.iter().enumerate() {
        if c > padded[mode] {
            mode = i;
        }
    }
    let peaks = find_peaks_with_widths(&padded, ties);
    Some(HistRecord {
        start: lo as Pos,
        mode: mode as i64 - HISTOGRAM_PAD as i64,
        var,
        skew,
        entropy,
        mean_height,
        padded,
        peaks,
    })
}

/// scipy `_local_maxima_1d`: plateau midpoints `(left + right) // 2`.
pub fn local_maxima(x: &[u32]) -> Vec<usize> {
    let mut out = Vec::new();
    if x.len() < 3 {
        return out;
    }
    let i_max = x.len() - 1;
    let mut i = 1;
    while i < i_max {
        if x[i - 1] < x[i] {
            let mut ahead = i + 1;
            while ahead < i_max && x[ahead] == x[i] {
                ahead += 1;
            }
            if x[ahead] < x[i] {
                out.push((i + ahead - 1) / 2);
                i = ahead;
            }
        }
        i += 1;
    }
    out
}

/// scipy `_select_by_peak_distance`: peaks (ascending positions) are visited in reverse
/// `np.argsort(heights)` order ([`np_argsort`], D-T1); each kept peak removes its neighbours
/// closer than `distance`.
pub fn select_by_peak_distance(peaks: &[usize], x: &[u32], distance: usize) -> Vec<usize> {
    select_by_peak_distance_ties(peaks, x, distance, Ties::Numpy)
}

fn select_by_peak_distance_ties(peaks: &[usize], x: &[u32], distance: usize, ties: Ties) -> Vec<usize> {
    let order: Vec<usize> = if ties == Ties::Numpy {
        let pri: Vec<i64> = peaks.iter().map(|&p| x[p] as i64).collect();
        let mut o = np_argsort(&pri);
        o.reverse();
        o
    } else {
        let mut o: Vec<usize> = (0..peaks.len()).collect();
        o.sort_by(|&a, &b| x[peaks[b]].cmp(&x[peaks[a]]).then(tie_cmp(ties, a, b)));
        o
    };
    let mut keep = vec![true; peaks.len()];
    for &j in &order {
        if !keep[j] {
            continue;
        }
        let mut k = j;
        while k > 0 && peaks[j] - peaks[k - 1] < distance {
            keep[k - 1] = false;
            k -= 1;
        }
        let mut k = j + 1;
        while k < peaks.len() && peaks[k] - peaks[j] < distance {
            keep[k] = false;
            k += 1;
        }
    }
    peaks.iter().zip(&keep).filter(|(_, &k)| k).map(|(&p, _)| p).collect()
}

/// scipy `_peak_prominences` (wlen=None) for one peak: `(prominence, left_base, right_base)`.
/// `left_min`/`right_min` move only on strictly smaller values.
fn prominence(x: &[f64], peak: usize) -> (f64, usize, usize) {
    let xp = x[peak];
    let mut left_base = peak;
    let mut left_min = xp;
    let mut i = peak as isize;
    while i >= 0 && x[i as usize] <= xp {
        if x[i as usize] < left_min {
            left_min = x[i as usize];
            left_base = i as usize;
        }
        i -= 1;
    }
    let mut right_base = peak;
    let mut right_min = xp;
    let mut i = peak;
    while i < x.len() && x[i] <= xp {
        if x[i] < right_min {
            right_min = x[i];
            right_base = i;
        }
        i += 1;
    }
    (xp - left_min.max(right_min), left_base, right_base)
}

/// scipy `peak_widths(x, [peak], rel_height)` for one peak: `(width, left_ips, right_ips)`.
pub fn peak_width(x: &[f64], peak: usize, rel_height: f64) -> (f64, f64, f64) {
    let (prom, i_min, i_max) = prominence(x, peak);
    let height = x[peak] - prom * rel_height;
    let mut i = peak;
    while i_min < i && height < x[i] {
        i -= 1;
    }
    let mut left_ip = i as f64;
    if x[i] < height {
        left_ip += (height - x[i]) / (x[i + 1] - x[i]);
    }
    let mut i = peak;
    while i < i_max && height < x[i] {
        i += 1;
    }
    let mut right_ip = i as f64;
    if x[i] < height {
        right_ip -= (height - x[i]) / (x[i - 1] - x[i]);
    }
    (right_ip - left_ip, left_ip, right_ip)
}

/// `find_peaks(padded, distance=10)` + `peak_widths(..., rel_height=0.98)` → raw peaks.
fn find_peaks_with_widths(padded: &[u32], ties: Ties) -> Vec<RawPeak> {
    let kept = select_by_peak_distance_ties(&local_maxima(padded), padded, PEAK_DISTANCE, ties);
    if kept.is_empty() {
        return Vec::new();
    }
    let x: Vec<f64> = padded.iter().map(|&c| c as f64).collect();
    kept.iter().map(|&p| raw_peak(&x, padded, p)).collect()
}

/// One raw peak at padded index `p`.
pub(crate) fn raw_peak(x: &[f64], padded: &[u32], p: usize) -> RawPeak {
    let pad = HISTOGRAM_PAD as f64;
    let (width, l, r) = peak_width(x, p, PEAK_REL_HEIGHT);
    RawPeak {
        location: p as i64 - HISTOGRAM_PAD as i64,
        width,
        // Python int() truncates toward zero, as does `as`.
        left: (l - pad) as i64,
        right: (r - pad) as i64,
        height: padded[p],
    }
}

/// `_rank_peaks`: `(index into record.peaks, relative_height)` in `np.argsort(-heights)`
/// order ([`np_argsort`], D-T1).
pub fn rank_peaks(rec: &HistRecord) -> Vec<(usize, f64)> {
    rank_peaks_ties(rec, Ties::Numpy)
}

fn rank_peaks_ties(rec: &HistRecord, ties: Ties) -> Vec<(usize, f64)> {
    let mut order: Vec<usize> = (0..rec.peaks.len()).collect();
    if order.len() <= 1 {
        return order.into_iter().map(|i| (i, 1.0)).collect();
    }
    if ties == Ties::Numpy {
        let neg: Vec<i64> = rec.peaks.iter().map(|p| -(p.height as i64)).collect();
        order = np_argsort(&neg);
    } else {
        order.sort_by(|&a, &b| rec.peaks[b].height.cmp(&rec.peaks[a].height).then(tie_cmp(ties, a, b)));
    }
    let top = rec.peaks[order[0]].height as f64;
    order
        .into_iter()
        .map(|i| {
            let h = rec.peaks[i].height as f64;
            (i, if top != 0.0 { h / top } else { 0.0 })
        })
        .collect()
}

/// Classifier input row (`FEATURE_COLUMNS` order).
pub fn feature_row(rec: &HistRecord, peak: &RawPeak, relative_height: f64) -> [f64; N_FEATURES] {
    [
        rec.var,
        rec.skew,
        rec.peaks.len() as f64,
        peak.width,
        rec.entropy,
        rec.mean_height,
        peak.height as f64,
        relative_height,
    ]
}

/// `_window_count`: reads in the padded window `[left, right]` (start-relative), clamped.
pub fn window_count(padded: &[u32], left: i64, right: i64) -> u32 {
    let low = (left + HISTOGRAM_PAD as i64).max(0);
    let high = (right + HISTOGRAM_PAD as i64 + 1).min(padded.len() as i64);
    if high <= low {
        return 0;
    }
    padded[low as usize..high as usize].iter().sum()
}

/// One scored candidate peak (a classifier row).
#[derive(Clone, Debug, PartialEq)]
pub struct ScoredPeak {
    /// Index into `HistRecord::peaks`.
    pub peak: usize,
    pub features: [f64; N_FEATURES],
    pub margin: f32,
    pub proba: f32,
    pub accepted: bool,
}

/// Everything `detect_peaks_batch` computes for one histogram.
#[derive(Clone, Debug, PartialEq)]
pub struct Detection {
    pub record: HistRecord,
    /// Classifier rows in ranked order (empty for the zero-peak mode fallback).
    pub rows: Vec<ScoredPeak>,
    /// Accepted peaks (genomic coordinates), ranked order.
    pub peaks: Vec<Peak>,
}

/// The peak a raw peak becomes when accepted.
pub fn genomic_peak(rec: &HistRecord, p: &RawPeak) -> Peak {
    let start = rec.start as i64;
    Peak {
        position: (start + p.location) as Pos,
        count: window_count(&rec.padded, p.left, p.right),
        left: (start + p.left) as Pos,
        right: (start + p.right) as Pos,
    }
}

/// `detect_peaks` with every intermediate. `None` for an empty histogram (no positive count).
pub fn detect_peaks_detailed(hist: &Histogram, model: &XgbModel) -> Option<Detection> {
    detect_peaks_ties(hist, model, Ties::Numpy)
}

/// [`detect_peaks_detailed`] with an explicit tie order (tests build alternative 4.1 outputs).
pub fn detect_peaks_ties(hist: &Histogram, model: &XgbModel, ties: Ties) -> Option<Detection> {
    let record = histogram_record_ties(hist, ties)?;
    let mut rows = Vec::new();
    let mut peaks = Vec::new();
    if record.peaks.is_empty() {
        // Unreachable (the zero padding guarantees a local maximum); ported anyway.
        let position = (record.start as i64 + record.mode) as Pos;
        let count = window_count(&record.padded, record.mode, record.mode);
        peaks.push(Peak { position, count, left: position, right: position });
    } else {
        for (i, rel) in rank_peaks_ties(&record, ties) {
            let features = feature_row(&record, &record.peaks[i], rel);
            let margin = model.margin(&features);
            let proba = sigmoid(margin);
            let accepted = proba > 0.5;
            if accepted {
                peaks.push(genomic_peak(&record, &record.peaks[i]));
            }
            rows.push(ScoredPeak { peak: i, features, margin, proba, accepted });
        }
    }
    Some(Detection { record, rows, peaks })
}

/// `fn detect_peaks(&Histogram, &XgbModel) -> Vec<Peak>`: accepted peaks in ranked order
/// (height descending, D-T1 ties by position); empty for an empty histogram.
pub fn detect_peaks(hist: &Histogram, model: &XgbModel) -> Vec<Peak> {
    detect_peaks_detailed(hist, model).map(|d| d.peaks).unwrap_or_default()
}

pub mod fixture {
    use super::*;
    use crate::xgb::bundled_models;
    use crate::xgb::fixture::{judge_pred, verdict};
    use iq_types::fixture::schema::{
        PeakRecordFixture, PeakRowFixture, RawPeaksFixture, U10aDetectPeaksIn, U10aDetectPeaksOut,
    };
    use iq_types::fixture::{DiffEntry, FixResult, FixtureCtx, PyF64, UnitFixture, Verdict};

    pub(crate) fn record_dto(r: &HistRecord) -> PeakRecordFixture {
        PeakRecordFixture {
            start: r.start,
            mode: r.mode as Pos,
            var: r.var,
            skew: PyF64(r.skew),
            entropy: r.entropy,
            mean_height: r.mean_height,
            peak_count: r.peaks.len() as u32,
            histogram: r.padded.clone(),
            raw_peaks: raw_dto(&r.peaks),
        }
    }

    fn raw_dto(peaks: &[RawPeak]) -> Option<RawPeaksFixture> {
        if peaks.is_empty() {
            return None;
        }
        Some(RawPeaksFixture {
            peak_location: peaks.iter().map(|p| p.location).collect(),
            peak_width: peaks.iter().map(|p| p.width).collect(),
            peak_left: peaks.iter().map(|p| p.left).collect(),
            peak_right: peaks.iter().map(|p| p.right).collect(),
            peak_heights: peaks.iter().map(|p| p.height as i64).collect(),
        })
    }

    fn row_dto(hist: usize, r: &ScoredPeak) -> PeakRowFixture {
        PeakRowFixture {
            hist: hist as u32,
            features: r.features.iter().map(|&f| PyF64(f)).collect(),
            proba: r.proba as f64,
            margin: r.margin as f64,
            label: r.accepted as i32,
        }
    }

    fn peak_tuple(p: &Peak) -> (Pos, u32, Pos, Pos) {
        (p.position, p.count, p.left, p.right)
    }

    pub struct U10aDetectPeaks;
    impl UnitFixture for U10aDetectPeaks {
        const UNIT: &'static str = "U10a";
        const FN: &'static str = "detect_peaks_batch";
        type In = U10aDetectPeaksIn;
        type Out = U10aDetectPeaksOut;

        fn run(input: &Self::In, _ctx: &FixtureCtx) -> FixResult<Self::Out> {
            run_ties(input, Ties::Numpy)
        }

        fn compare(got: &Self::Out, want: &Self::Out, _ctx: &FixtureCtx) -> Verdict {
            if got == want {
                return Verdict::Exact;
            }
            compare_batch(got, want)
        }
    }

    /// The runner with an explicit tie order.
    pub fn run_ties(input: &U10aDetectPeaksIn, ties: Ties) -> FixResult<U10aDetectPeaksOut> {
        let model = bundled_models().by_name(&input.model).ok_or_else(|| format!("unknown model {:?}", input.model))?;
        let mut out = U10aDetectPeaksOut { records: Vec::new(), rows: Vec::new(), peaks: Vec::new() };
        for (i, h) in input.histograms.iter().enumerate() {
            let hist: Histogram = h.iter().flatten().copied().collect();
            match detect_peaks_ties(&hist, model, ties) {
                Some(d) => {
                    out.records.push(Some(record_dto(&d.record)));
                    out.rows.extend(d.rows.iter().map(|r| row_dto(i, r)));
                    out.peaks.push(d.peaks.iter().map(peak_tuple).collect());
                }
                None => {
                    out.records.push(None);
                    out.peaks.push(Vec::new());
                }
            }
        }
        Ok(out)
    }

    /// f64 relative difference ≤ 1e-12 (NaN equals NaN).
    pub(crate) fn feq(a: f64, b: f64) -> bool {
        (a.is_nan() && b.is_nan()) || a == b || (a - b).abs() <= 1e-12 * a.abs().max(b.abs())
    }

    fn push(e: &mut Vec<DiffEntry>, path: String, got: impl std::fmt::Debug, want: impl std::fmt::Debug) {
        e.push(DiffEntry { path, got: format!("{got:?}"), want: format!("{want:?}") });
    }

    /// Tie-independent record fields.
    fn cmp_record_base(p: &str, g: &PeakRecordFixture, w: &PeakRecordFixture, e: &mut Vec<DiffEntry>) {
        if (g.start, g.mode, &g.histogram) != (w.start, w.mode, &w.histogram) {
            push(e, format!("{p}.start/mode/histogram"), (g.start, g.mode), (w.start, w.mode));
        }
        for (f, a, b) in [
            ("var", g.var, w.var),
            ("skew", g.skew.0, w.skew.0),
            ("entropy", g.entropy, w.entropy),
            ("mean_height", g.mean_height, w.mean_height),
        ] {
            if !feq(a, b) {
                push(e, format!("{p}.{f}"), a, b);
            }
        }
    }

    fn cmp_raw(p: &str, g: &Option<RawPeaksFixture>, w: &Option<RawPeaksFixture>, e: &mut Vec<DiffEntry>) {
        match (g, w) {
            (None, None) => {}
            (Some(g), Some(w)) => {
                let ints = |r: &RawPeaksFixture| {
                    (r.peak_location.clone(), r.peak_left.clone(), r.peak_right.clone(), r.peak_heights.clone())
                };
                if ints(g) != ints(w) {
                    push(e, format!("{p}.raw_peaks"), ints(g), ints(w));
                } else if !g.peak_width.iter().zip(&w.peak_width).all(|(&a, &b)| feq(a, b)) {
                    push(e, format!("{p}.raw_peaks.peak_width"), &g.peak_width, &w.peak_width);
                }
            }
            _ => push(e, format!("{p}.raw_peaks"), g.is_some(), w.is_some()),
        }
    }

    fn cmp_row(p: &str, g: &PeakRowFixture, w: &PeakRowFixture, e: &mut Vec<DiffEntry>) {
        if g.hist != w.hist || g.features.len() != w.features.len() {
            push(e, format!("{p}.hist"), g.hist, w.hist);
            return;
        }
        for (k, (a, b)) in g.features.iter().zip(&w.features).enumerate() {
            if !feq(a.0, b.0) {
                push(e, format!("{p}.features[{k}]"), a.0, b.0);
            }
        }
        judge_pred(p, (g.margin, g.proba, g.label), (w.margin, w.proba, w.label), e);
    }

    /// Tolerant comparison (§3 U10a pass criteria); histograms that fail it are re-checked
    /// for D-T1 consistency.
    pub fn compare_batch(got: &U10aDetectPeaksOut, want: &U10aDetectPeaksOut) -> Verdict {
        let mut e = Vec::new();
        if got.records.len() != want.records.len() || got.peaks.len() != want.peaks.len() {
            push(&mut e, ".len()".into(), got.records.len(), want.records.len());
            return verdict(e);
        }
        let rows_of = |o: &U10aDetectPeaksOut, i: usize| -> Vec<PeakRowFixture> {
            o.rows.iter().filter(|r| r.hist as usize == i).cloned().collect()
        };
        let mut explained = false;
        // The case's model is not in `out`: it is the bundled model that reproduces iqrs's own
        // first row bit-exactly.
        let model = got.rows.first().and_then(|r| {
            let f: [f64; N_FEATURES] = std::array::from_fn(|i| r.features[i].0);
            let ms = bundled_models();
            [&ms.polya, &ms.tss].into_iter().find(|m| m.margin(&f) as f64 == r.margin)
        });
        for i in 0..want.records.len() {
            let p = format!(".records[{i}]");
            let mut he = Vec::new();
            let (gr, wr) = (&got.records[i], &want.records[i]);
            match (gr, wr) {
                (None, None) => {}
                (Some(g), Some(w)) => {
                    cmp_record_base(&p, g, w, &mut he);
                    if !he.is_empty() {
                        e.extend(he);
                        continue;
                    }
                    let mut te = Vec::new();
                    if g.peak_count != w.peak_count {
                        push(&mut te, format!("{p}.peak_count"), g.peak_count, w.peak_count);
                    }
                    cmp_raw(&p, &g.raw_peaks, &w.raw_peaks, &mut te);
                    let (grows, wrows) = (rows_of(got, i), rows_of(want, i));
                    if grows.len() != wrows.len() {
                        push(&mut te, format!(".rows(hist {i}).len()"), grows.len(), wrows.len());
                    } else {
                        for (k, (a, b)) in grows.iter().zip(&wrows).enumerate() {
                            cmp_row(&format!(".rows(hist {i})[{k}]"), a, b, &mut te);
                        }
                    }
                    if got.peaks[i] != want.peaks[i] {
                        push(&mut te, format!(".peaks[{i}]"), &got.peaks[i], &want.peaks[i]);
                    }
                    if te.is_empty() {
                        continue;
                    }
                    let Some(model) = model else {
                        e.extend(te);
                        continue;
                    };
                    match tie_consistent(w, &wrows, &want.peaks[i], model) {
                        Ok(()) => explained = true,
                        Err(msg) => {
                            e.extend(te);
                            push(&mut e, format!("{p} (D-T1 check)"), "inconsistent", msg);
                        }
                    }
                }
                _ => push(&mut e, p, gr.is_some(), wr.is_some()),
            }
        }
        if e.is_empty() && explained {
            Verdict::Explained("D-T1")
        } else {
            verdict(e)
        }
    }

    /// D-T1: is 4.1's output for one histogram what some order of equal heights (in the
    /// distance filter and in the ranking) produces? The tie-independent fields were already
    /// checked equal. `model` re-evaluates 4.1's rows.
    pub(crate) fn tie_consistent(
        w: &PeakRecordFixture,
        wrows: &[PeakRowFixture],
        wpeaks: &[(Pos, u32, Pos, Pos)],
        model: &XgbModel,
    ) -> Result<(), String> {
        let x = &w.histogram;
        let cand = local_maxima(x);
        let pad = HISTOGRAM_PAD as i64;
        let kept: Vec<usize> = match &w.raw_peaks {
            Some(r) => r.peak_location.iter().map(|&l| (l + pad) as usize).collect(),
            None => Vec::new(),
        };
        // 1. Distance filter: greedy by priority keeps a peak iff no earlier-kept peak lies
        //    within the distance. A kept set K is reachable iff K ⊆ candidates, K is pairwise
        //    ≥ distance apart, and every dropped candidate has a K-neighbour within the
        //    distance that is at least as high (ties can be ordered K-first).
        if kept.windows(2).any(|p| p[1] <= p[0]) || kept.iter().any(|k| !cand.contains(k)) {
            return Err("kept peaks are not sorted local maxima".into());
        }
        if kept.windows(2).any(|p| p[1] - p[0] < PEAK_DISTANCE) {
            return Err("kept peaks closer than the distance".into());
        }
        for &c in cand.iter().filter(|c| !kept.contains(c)) {
            if !kept.iter().any(|&k| k.abs_diff(c) < PEAK_DISTANCE && x[k] >= x[c]) {
                return Err(format!("candidate {c} dropped without a dominating neighbour"));
            }
        }
        // 2. Raw peaks for K are tie-independent.
        let xf: Vec<f64> = x.iter().map(|&c| c as f64).collect();
        let raw: Vec<RawPeak> = kept.iter().map(|&p| raw_peak(&xf, x, p)).collect();
        let mut e = Vec::new();
        cmp_raw("", &raw_dto(&raw), &w.raw_peaks, &mut e);
        if w.peak_count as usize != raw.len() || !e.is_empty() {
            return Err(format!("raw peaks differ: {e:?}"));
        }
        if raw.is_empty() {
            return if wrows.is_empty() && wpeaks.len() == 1 { Ok(()) } else { Err("mode fallback differs".into()) };
        }
        // 3. Ranking: 4.1's rows must be a height-descending order of K. Accepted rows are
        //    pinned to their peak by the reported peak tuple; rejected rows are matched to the
        //    remaining peaks of their height by features.
        let rec = HistRecord {
            start: w.start,
            mode: w.mode as i64,
            var: w.var,
            skew: w.skew.0,
            entropy: w.entropy,
            mean_height: w.mean_height,
            padded: x.clone(),
            peaks: raw,
        };
        let top = rec.peaks.iter().map(|p| p.height).max().unwrap_or(0) as f64;
        let model_rows: Vec<PeakRowFixture> = rec
            .peaks
            .iter()
            .map(|p| {
                let rel = if rec.peaks.len() == 1 {
                    1.0
                } else if top != 0.0 {
                    p.height as f64 / top
                } else {
                    0.0
                };
                let features = feature_row(&rec, p, rel);
                PeakRowFixture {
                    hist: 0,
                    features: features.iter().map(|&f| PyF64(f)).collect(),
                    proba: 0.0,
                    margin: 0.0,
                    label: 0,
                }
            })
            .collect();
        if wrows.len() != rec.peaks.len() {
            return Err("row count differs".into());
        }
        if wrows.windows(2).any(|r| r[1].features[6].0 > r[0].features[6].0) {
            return Err("rows not in height-descending order".into());
        }
        let feats_eq =
            |a: &PeakRowFixture, b: &PeakRowFixture| a.features.iter().zip(&b.features).all(|(x, y)| feq(x.0, y.0));
        let mut used = vec![false; rec.peaks.len()];
        let mut assign = vec![usize::MAX; wrows.len()];
        let mut next_peak = 0;
        for (k, r) in wrows.iter().enumerate() {
            if r.label != 1 {
                continue;
            }
            let Some(&t) = wpeaks.get(next_peak) else {
                return Err("fewer peaks than accepted rows".into());
            };
            next_peak += 1;
            let j = rec
                .peaks
                .iter()
                .position(|p| peak_tuple(&genomic_peak(&rec, p)) == t)
                .ok_or_else(|| format!("accepted peak {t:?} is not a raw peak"))?;
            if used[j] || !feats_eq(&model_rows[j], r) {
                return Err(format!("accepted row {k} does not match peak {t:?}"));
            }
            used[j] = true;
            assign[k] = j;
        }
        if next_peak != wpeaks.len() {
            return Err("more peaks than accepted rows".into());
        }
        for (k, r) in wrows.iter().enumerate() {
            if r.label == 1 {
                continue;
            }
            let j = (0..rec.peaks.len())
                .find(|&j| !used[j] && feats_eq(&model_rows[j], r))
                .ok_or_else(|| format!("rejected row {k} matches no remaining peak"))?;
            used[j] = true;
            assign[k] = j;
        }
        // 4. Predictions: re-evaluate each row with the native model.
        let check_pred = |r: &PeakRowFixture| -> Result<(), String> {
            let f: [f64; N_FEATURES] = std::array::from_fn(|i| r.features[i].0);
            let m = model.margin(&f);
            let mut e = Vec::new();
            judge_pred(
                "",
                (m as f64, sigmoid(m) as f64, (sigmoid(m) > 0.5) as i32),
                (r.margin, r.proba, r.label),
                &mut e,
            );
            if e.is_empty() {
                Ok(())
            } else {
                Err(format!("prediction differs: {e:?}"))
            }
        };
        for r in wrows {
            check_pred(r)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pairwise_small_and_blocked() {
        assert_eq!(np_sum(&[]), 0.0);
        assert_eq!(np_sum(&[1.0, 2.0, 3.0]), 6.0);
        let v: Vec<f64> = (0..300).map(|i| i as f64 * 0.1).collect();
        assert!((np_sum(&v) - 4485.0).abs() < 1e-9);
    }

    #[test]
    fn local_maxima_plateau_midpoint() {
        assert_eq!(local_maxima(&[0, 1, 1, 1, 1, 0]), vec![2]);
        assert_eq!(local_maxima(&[0, 2, 0, 3, 3, 0]), vec![1, 3]);
        assert_eq!(local_maxima(&[0, 1, 1]), Vec::<usize>::new());
    }

    #[test]
    fn distance_ties_by_position() {
        // two equal peaks 5 apart: numpy's stable small-array argsort visits the right one
        // first, so it wins (D-T1); the position-ordered alternatives keep either side.
        let mut x = vec![0u32; 30];
        x[10] = 3;
        x[15] = 3;
        let p = local_maxima(&x);
        assert_eq!(select_by_peak_distance(&p, &x, 10), vec![15]);
        assert_eq!(select_by_peak_distance_ties(&p, &x, 10, Ties::PositionAsc), vec![10]);
        assert_eq!(select_by_peak_distance_ties(&p, &x, 10, Ties::PositionDesc), vec![15]);
    }

    /// `np.argsort` reference vectors from numpy 2.4.6 with
    /// `NPY_DISABLE_CPU_FEATURES="AVX512_SPR AVX512_ICL X86_V4 X86_V3"` (non-SIMD path), random
    /// heights with many ties (`default_rng(7)`); both `argsort(v)` and `argsort(-v)`.
    #[test]
    fn np_argsort_matches_numpy_nosimd() {
        #[rustfmt::skip]
        let cases: &[(&[i64], &[usize], &[usize])] = &[
            (&[1, 1, 1, 1, 1], &[0, 1, 2, 3, 4], &[0, 1, 2, 3, 4]),
            (&[2, 2, 0, 0, 0, 0, 2, 2, 0, 1, 2, 0, 2, 0, 1, 2, 0], &[8, 13, 11, 5, 16, 3, 2, 4, 9, 14, 6, 7, 15, 10, 12, 1, 0], &[0, 1, 12, 6, 7, 15, 10, 14, 9, 13, 11, 8, 4, 3, 2, 5, 16]),
            (&[1, 0, 2, 0, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 2], &[1, 3, 0, 15, 10, 9, 16, 8, 6, 5, 7, 4, 11, 12, 13, 14, 2, 17], &[17, 2, 14, 4, 13, 12, 11, 15, 10, 0, 16, 7, 6, 5, 9, 8, 3, 1]),
            (&[1, 0, 3, 0, 3, 2, 0, 0, 1, 0, 0, 2, 3, 1, 3, 3, 3, 2, 1, 2, 1, 1, 1, 0, 3, 0, 0, 0, 3, 2, 3, 0, 2, 1, 1, 0, 2, 3, 2, 0], &[39, 23, 25, 26, 27, 10, 9, 31, 6, 3, 35, 1, 7, 21, 33, 34, 22, 20, 0, 8, 18, 13, 38, 36, 5, 32, 11, 29, 17, 19, 28, 24, 30, 14, 15, 16, 4, 2, 37, 12], &[37, 2, 24, 4, 16, 15, 14, 30, 28, 12, 19, 17, 29, 38, 32, 5, 36, 11, 33, 34, 22, 0, 20, 18, 13, 21, 8, 1, 3, 35, 6, 7, 9, 10, 27, 26, 25, 23, 31, 39]),
            (&[2, 1, 4, 4, 0, 2, 4, 4, 3, 3, 0, 3, 2, 0, 1, 2, 3, 2, 3, 4, 3, 1, 3, 2, 0, 0, 3, 1, 2, 1, 1, 0, 1, 4, 2, 1, 2, 4, 1, 2, 2, 3, 2, 3, 2, 3, 4, 0, 2, 2, 1, 1, 0, 2, 4, 0, 2, 4, 4, 1, 0, 3, 0, 1, 2, 4, 0, 3, 2, 0, 4, 4, 2, 4, 3, 4, 4, 2, 1, 0, 2, 0, 3, 4, 1, 2, 0, 0, 1, 4, 4, 3, 2, 2, 0, 1, 4, 2, 1, 1], &[69, 52, 55, 60, 62, 31, 66, 25, 24, 79, 81, 13, 47, 86, 10, 87, 4, 94, 50, 51, 1, 59, 38, 95, 35, 63, 32, 29, 30, 27, 14, 84, 78, 21, 98, 99, 88, 72, 64, 93, 92, 56, 97, 53, 80, 68, 85, 77, 0, 49, 23, 28, 17, 34, 36, 48, 39, 12, 15, 40, 42, 44, 5, 9, 82, 8, 16, 91, 18, 20, 11, 22, 74, 26, 45, 67, 61, 43, 41, 73, 2, 96, 3, 6, 7, 90, 89, 54, 57, 58, 37, 83, 33, 65, 19, 76, 75, 70, 71, 46], &[19, 76, 75, 37, 83, 73, 89, 90, 57, 70, 71, 58, 7, 6, 65, 96, 3, 2, 46, 33, 54, 67, 41, 74, 26, 61, 82, 20, 18, 16, 43, 91, 11, 9, 8, 45, 22, 56, 0, 77, 68, 72, 80, 85, 92, 93, 97, 64, 53, 49, 28, 34, 23, 36, 17, 15, 40, 39, 44, 12, 5, 48, 42, 14, 88, 98, 84, 21, 95, 1, 78, 27, 29, 50, 99, 32, 51, 35, 38, 63, 30, 59, 4, 52, 94, 47, 10, 55, 13, 31, 87, 62, 66, 81, 79, 69, 24, 25, 60, 86]),
            (&[4, 0, 5, 5, 3, 2, 4, 3, 2, 1, 2, 4, 0, 0, 4, 2, 4, 0, 5, 0, 2, 5, 0, 3, 4, 2, 0, 3, 2, 5, 4, 2, 0, 3, 1, 4, 5, 2, 5, 3, 3, 4, 4, 5, 5, 0, 2, 5, 1, 0, 3, 4, 4, 4, 3, 0, 5, 2, 5, 4, 5, 0, 5, 3, 2, 4, 5, 3, 0, 4, 1, 1, 1, 1, 5, 2, 0, 1, 2, 2, 1, 5, 5, 3, 5, 2, 0, 1, 0, 5, 1, 2, 1, 5, 0, 3, 5, 3, 5, 5, 4, 4, 3, 3, 4, 2, 5, 5, 2, 2, 4, 5, 0, 0, 1, 2, 0, 3, 3, 5, 4, 1, 4, 4, 4, 4, 1, 4, 0, 3, 5, 5, 1, 1, 3, 2, 4, 1, 1, 0, 2, 1, 1, 5, 4, 5, 3, 0, 4, 3, 5, 2, 1, 3, 2, 3, 2, 1, 0, 5, 4, 2, 1, 3, 2, 3, 0, 1, 4, 0, 4, 2, 0, 4, 4, 4, 5, 4, 5, 0, 5, 5, 4, 4, 4, 5, 5, 3, 4, 5, 1, 0, 1, 0, 2, 0, 1, 1, 4, 1, 5, 1, 5, 3, 1, 0, 1, 3, 0, 0, 0, 4, 0, 0, 2, 1, 3, 5, 5, 3, 3, 4, 2, 3, 1, 3, 1, 1, 4, 3, 2, 0, 1, 4, 2, 5, 5, 5, 0, 0, 3, 2, 5, 1, 0, 0, 1, 3, 0, 4, 1, 1, 3, 5, 2, 4, 1, 0, 4, 4, 4, 5, 2, 1, 3, 3, 2, 3, 1, 3, 4, 0, 2, 3, 0, 3, 5, 4, 4, 4, 2, 1, 0, 4, 1, 5, 2, 5, 0, 0, 2, 0, 2, 3, 3, 1, 0, 3, 3, 5], &[76, 193, 245, 191, 116, 113, 45, 112, 208, 248, 49, 209, 274, 210, 212, 55, 166, 271, 213, 61, 169, 179, 257, 68, 94, 172, 88, 86, 244, 32, 282, 22, 231, 291, 19, 239, 13, 139, 238, 205, 289, 158, 12, 296, 147, 26, 195, 128, 1, 288, 17, 256, 142, 126, 224, 157, 92, 90, 87, 232, 204, 137, 167, 138, 251, 227, 250, 133, 132, 114, 162, 226, 246, 263, 121, 243, 141, 80, 152, 192, 268, 48, 281, 215, 197, 206, 9, 71, 72, 70, 199, 295, 190, 201, 34, 284, 77, 73, 196, 241, 25, 161, 31, 194, 286, 28, 79, 290, 20, 292, 135, 10, 140, 8, 230, 5, 156, 154, 234, 15, 115, 272, 91, 171, 222, 266, 64, 37, 254, 214, 57, 85, 75, 262, 108, 109, 164, 46, 151, 78, 280, 105, 155, 223, 163, 225, 219, 216, 187, 165, 207, 220, 149, 252, 67, 153, 267, 63, 269, 54, 273, 275, 50, 265, 40, 33, 27, 23, 293, 294, 7, 4, 297, 39, 264, 203, 83, 298, 229, 146, 95, 97, 134, 129, 102, 103, 240, 118, 117, 247, 233, 283, 278, 221, 277, 249, 211, 255, 270, 258, 259, 260, 279, 228, 0, 101, 144, 136, 69, 41, 42, 127, 125, 124, 123, 122, 120, 110, 51, 104, 52, 53, 100, 59, 65, 30, 148, 35, 11, 6, 182, 188, 14, 16, 183, 24, 175, 177, 198, 174, 173, 170, 168, 160, 184, 60, 2, 58, 66, 56, 3, 62, 287, 47, 44, 29, 43, 18, 21, 38, 36, 285, 276, 74, 200, 82, 131, 143, 145, 150, 159, 176, 235, 218, 178, 180, 181, 185, 186, 189, 217, 81, 236, 130, 261, 84, 89, 93, 96, 98, 237, 99, 202, 106, 107, 111, 119, 242, 253, 299], &[299, 38, 185, 43, 44, 181, 96, 47, 178, 176, 56, 58, 159, 60, 150, 62, 145, 66, 98, 99, 89, 106, 84, 107, 186, 82, 111, 119, 74, 130, 131, 143, 81, 36, 180, 253, 287, 2, 3, 285, 276, 261, 189, 18, 242, 237, 21, 236, 235, 93, 29, 200, 202, 218, 217, 123, 122, 120, 270, 184, 110, 277, 278, 233, 183, 283, 211, 104, 101, 100, 198, 279, 124, 127, 177, 175, 174, 173, 170, 168, 228, 160, 221, 249, 148, 182, 144, 255, 258, 136, 259, 260, 125, 188, 0, 11, 35, 24, 69, 41, 42, 16, 14, 51, 52, 30, 65, 6, 59, 53, 54, 203, 33, 187, 39, 40, 155, 163, 165, 207, 50, 153, 95, 219, 27, 297, 294, 293, 4, 7, 275, 273, 269, 267, 216, 265, 252, 247, 240, 23, 229, 225, 223, 220, 298, 264, 63, 149, 146, 97, 102, 103, 83, 117, 129, 134, 118, 67, 25, 37, 266, 85, 78, 262, 108, 109, 194, 15, 254, 79, 31, 241, 20, 234, 115, 75, 230, 156, 214, 28, 154, 105, 222, 286, 91, 164, 10, 140, 64, 161, 57, 171, 290, 5, 280, 135, 8, 151, 46, 272, 292, 227, 295, 92, 224, 232, 226, 268, 243, 246, 250, 251, 90, 80, 284, 256, 281, 263, 87, 9, 114, 152, 157, 192, 201, 199, 126, 196, 190, 34, 71, 70, 132, 133, 48, 142, 137, 138, 167, 162, 141, 72, 204, 197, 77, 121, 73, 206, 215, 244, 61, 212, 49, 116, 88, 68, 213, 282, 172, 55, 274, 1, 288, 289, 76, 291, 166, 139, 296, 158, 169, 179, 271, 19, 245, 205, 112, 248, 17, 113, 32, 239, 238, 94, 195, 45, 257, 13, 12, 26, 191, 22, 208, 128, 209, 210, 86, 231, 193, 147]),
        ];
        for (v, asc, desc) in cases {
            assert_eq!(np_argsort(v), *asc, "argsort(v), n={}", v.len());
            let neg: Vec<i64> = v.iter().map(|x| -x).collect();
            assert_eq!(np_argsort(&neg), *desc, "argsort(-v), n={}", v.len());
        }
        // depth-limit heapsort path: adversarial-ish input still yields a sorted permutation
        let v: Vec<i64> = (0..5000).map(|i| ((i * 7919) % 13) as i64).collect();
        let o = np_argsort(&v);
        assert!(o.windows(2).all(|w| v[w[0]] <= v[w[1]]));
    }

    #[test]
    fn single_position() {
        let h: Histogram = [(100, 6)].into_iter().collect();
        let r = histogram_record(&h).unwrap();
        assert_eq!(r.peaks.len(), 1);
        assert_eq!(r.peaks[0].width, 1.9600000000000009);
        assert!(r.skew.is_nan());
        assert_eq!(window_count(&r.padded, 0, 0), 6);
    }

    #[test]
    fn empty_histogram() {
        let m = crate::xgb::bundled_models();
        assert!(detect_peaks(&Histogram::default(), &m.polya).is_empty());
        let h: Histogram = [(5, 0)].into_iter().collect();
        assert!(detect_peaks(&h, &m.polya).is_empty());
    }
}
