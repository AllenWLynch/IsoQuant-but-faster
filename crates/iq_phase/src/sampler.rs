//! Blocked Gibbs sampler for one locus (phased_model/code/gibbs.py; alpha from --phase-alpha, default 5).
//!
//! Model: beta ~ Dir(gamma/T); theta_A, theta_B | beta ~ Dir(alpha beta); rho ~ Beta(a, a);
//! read r: h_r ~ Bern(rho), t_r | h_r ~ Cat(theta_h) restricted to C_r (l_rt uniform), het
//! evidence lambda_r(h).
//!
//! One sweep, in gibbs.py's order:
//!   1. allocate each class multinomially over (h, t in C_r) given rho, theta (+ Rao-Blackwellised
//!      expected counts E[n_ht] and data-only allele evidence per transcript);
//!   2. rho | n ~ Beta(a + n_A, a + n_B);
//!   3. m_ht = sum_{j<n_ht} Bern(alpha beta_t / (alpha beta_t + j)) (CRT draw);
//!   4. beta ~ Dir(gamma/T + m_.t);
//!   5. theta_h ~ Dir(alpha beta + n_h.) (last).
//!
//! Dirichlet draws use the log-space trick of gibbs.py `log_dirichlet` (Gamma(a + 1) U^(1/a)),
//! stable for tiny concentrations. The hot loop allocates nothing: all buffers live in
//! [`Scratch`], reused across loci.

use crate::locus::Locus;
use rand::{Rng, SeedableRng};
use rand_distr::{Binomial, Distribution, StandardNormal};
use rand_xoshiro::Xoshiro256PlusPlus;

pub type PhRng = Xoshiro256PlusPlus;

/// Reads of a class drawn one by one (categorical) below this count, by binomials above.
const CAT_MAX: u64 = 16;

/// Sampler settings.
#[derive(Clone, Copy, Debug)]
pub struct Cfg {
    pub burnin: u32,
    pub sweeps: u32,
    pub gamma: f64,
    pub alpha: f64,
    pub a_rho: f64,
}

impl Default for Cfg {
    fn default() -> Self {
        Cfg { burnin: 200, sweeps: 1000, gamma: 1.0, alpha: 5.0, a_rho: 20.0 }
    }
}

/// Per-transcript summary.
#[derive(Clone, Debug, Default)]
pub struct TxSummary {
    pub a_mean: f64,
    pub a_lo: f64,
    pub a_hi: f64,
    pub b_mean: f64,
    pub b_lo: f64,
    pub b_hi: f64,
    pub phase_uncertainty: f64,
    pub evid_a: f64,
    pub evid_b: f64,
}

/// Per-locus summary.
#[derive(Clone, Debug, Default)]
pub struct LocusSummary {
    pub tx: Vec<TxSummary>,
    pub n_reads: u64,
    pub rho_mean: f64,
    pub rho_lo: f64,
    pub rho_hi: f64,
    pub p_rho_gt_half: f64,
    pub tv_mean: f64,
    pub phase_uncertainty: f64,
    pub evid_a: f64,
    pub evid_b: f64,
}

/// SD of Beta(a, b) relative to SD of Uniform(0, 1).
pub fn phase_uncertainty(a: f64, b: f64) -> f64 {
    (a * b / ((a + b) * (a + b) * (a + b + 1.0))).sqrt() / (1.0f64 / 12.0).sqrt()
}

/// Uniform in the open interval (0, 1).
#[inline]
pub fn unif<R: Rng + ?Sized>(r: &mut R) -> f64 {
    ((r.next_u64() >> 11) as f64 + 0.5) * (1.0 / 9_007_199_254_740_992.0)
}

/// log of a Gamma(shape, 1) draw for shape >= 1 (Marsaglia-Tsang).
#[inline]
pub fn log_gamma_ge1<R: Rng + ?Sized>(r: &mut R, shape: f64) -> f64 {
    debug_assert!(shape >= 1.0);
    let d = shape - 1.0 / 3.0;
    let c = 1.0 / (9.0 * d).sqrt();
    loop {
        let x: f64 = StandardNormal.sample(r);
        let v = 1.0 + c * x;
        if v <= 0.0 {
            continue;
        }
        let v3 = v * v * v;
        let u = unif(r);
        let x2 = x * x;
        if u < 1.0 - 0.0331 * x2 * x2 || u.ln() < 0.5 * x2 + d * (1.0 - v3 + v3.ln()) {
            return d.ln() + v3.ln();
        }
    }
}

/// log of a Gamma(shape, 1) draw for any shape > 0 (Gamma(a + 1) U^(1/a) below 1).
#[inline]
pub fn log_gamma<R: Rng + ?Sized>(r: &mut R, shape: f64) -> f64 {
    if shape >= 1.0 {
        log_gamma_ge1(r, shape)
    } else {
        log_gamma_ge1(r, shape + 1.0) + unif(r).ln() / shape
    }
}

/// Beta(a, b) draw.
#[inline]
pub fn beta<R: Rng + ?Sized>(r: &mut R, a: f64, b: f64) -> f64 {
    let x = log_gamma(r, a);
    let y = log_gamma(r, b);
    // x / (x + y) in log space
    1.0 / (1.0 + (y - x).exp())
}

/// Normalized log Dirichlet draw into `out` (gibbs.py `log_dirichlet`).
pub fn log_dirichlet<R: Rng + ?Sized>(r: &mut R, shape: impl Iterator<Item = f64>, out: &mut [f64]) {
    let mut m = f64::NEG_INFINITY;
    for (o, a) in out.iter_mut().zip(shape) {
        let x = log_gamma_ge1(r, a + 1.0) + unif(r).ln() / a;
        *o = x;
        if x > m {
            m = x;
        }
    }
    let s: f64 = out.iter().map(|&x| (x - m).exp()).sum();
    let z = m + s.ln();
    for o in out.iter_mut() {
        *o -= z;
    }
}

/// ln Gamma(x + a) - ln Gamma(x) for x >= 1, a > 0, accurate for tiny a.
fn lgamma_shift(x: f64, a: f64) -> f64 {
    // recurse x up to >= 16 with exact terms, then the Stirling-series difference
    let mut acc = 0.0;
    let mut x = x;
    while x < 16.0 {
        // lnG(x+a) - lnG(x) = lnG(x+1+a) - lnG(x+1) - ln((x+a)/x)
        acc -= (a / x).ln_1p();
        x += 1.0;
    }
    let y = x + a;
    // (y - 0.5) ln y - (x - 0.5) ln x - a
    let main = (x - 0.5) * (a / x).ln_1p() + a * y.ln() - a;
    let (iy, ix) = (1.0 / y, 1.0 / x);
    // 1/(12z) - 1/(360z^3) + 1/(1260z^5) - 1/(1680 z^7), differences computed directly
    let t1 = -a / (x * y) / 12.0;
    let t3 = -(iy.powi(3) - ix.powi(3)) / 360.0;
    let t5 = (iy.powi(5) - ix.powi(5)) / 1260.0;
    let t7 = -(iy.powi(7) - ix.powi(7)) / 1680.0;
    acc + main + t1 + t3 + t5 + t7
}

/// Number of tables of a Chinese restaurant with n customers and concentration ab:
/// sum_{j=0}^{n-1} Bern(ab / (ab + j)). Exact; skips runs of failures by inverting the
/// survival function S(k) = prod_{i=j}^{k} i / (ab + i) with a binary search, so the cost is
/// O(m log n) instead of O(n).
pub fn crt<R: Rng + ?Sized>(r: &mut R, n: u64, ab: f64) -> u64 {
    if n == 0 {
        return 0;
    }
    const DIRECT: u64 = 24;
    let mut m = 1u64; // j = 0 always opens a table
    let mut j = 1u64;
    while j < n.min(DIRECT) {
        if unif(r) * (ab + j as f64) < ab {
            m += 1;
        }
        j += 1;
    }
    // log S(k) = -(G(k + 1) - G(j)), G(x) = lnGamma(x + ab) - lnGamma(x)
    while j < n {
        let gj = lgamma_shift(j as f64, ab);
        let lu = unif(r).ln();
        let log_s = |k: u64| -(lgamma_shift(k as f64 + 1.0, ab) - gj);
        if log_s(n - 1) >= lu {
            break;
        }
        // smallest k in [j, n-1] with log S(k) < lu. Start from the asymptotic inverse
        // (G(x) ~ ab ln x => k ~ j e^(-lu/ab) - 1), bracket by galloping, then bisect.
        let guess = (j as f64 * (-lu / ab).exp() - 1.0).round();
        let k0 = if guess.is_finite() { (guess.max(j as f64) as u64).min(n - 1) } else { n - 1 };
        let (mut lo, mut hi);
        if log_s(k0) < lu {
            // answer <= k0: find lo with log S(lo - 1) >= lu
            hi = k0;
            let mut step = 1u64;
            loop {
                if hi - j < step {
                    lo = j;
                    break;
                }
                let c = hi - step;
                if log_s(c) >= lu {
                    lo = c + 1;
                    break;
                }
                hi = c;
                step *= 2;
            }
        } else {
            // answer > k0 (and <= n - 1 since log S(n - 1) < lu)
            lo = k0 + 1;
            let mut step = 1u64;
            loop {
                let c = (lo + step - 1).min(n - 1);
                if log_s(c) < lu {
                    hi = c;
                    break;
                }
                lo = c + 1;
                step *= 2;
            }
        }
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            if log_s(mid) < lu {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        m += 1;
        j = lo + 1;
    }
    m
}

/// Reusable buffers.
#[derive(Default)]
pub struct Scratch {
    lt: Vec<f64>,   // [2T] log theta
    w: Vec<f64>,    // [2T] rho_h theta_h,t scaled
    beta: Vec<f64>, // [T]
    lb: Vec<f64>,   // [T]
    n: Vec<u32>,    // [2T] current counts
    m: Vec<f64>,    // [T] table counts per t (summed over h)
    rb: Vec<f64>,   // [2T] sum over kept sweeps of E[n_ht]
    inf: Vec<f64>,  // [2T] sum over kept sweeps of evidence
    hist: Vec<u32>, // [S x 2T]
    rho: Vec<f64>,
    tv: Vec<f64>,
    p: Vec<f64>, // [2K] class cell probabilities
    col: Vec<u32>,
    eh: Vec<(f64, f64)>, // per class exp(+-d/2)
}

/// Per-locus seed: independent of thread scheduling.
pub fn locus_seed(seed: u64, name: &str) -> u64 {
    let mut z = seed ^ xxhash_rust::xxh3::xxh3_64(name.as_bytes());
    // splitmix64 finalizer
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// numpy `quantile` (linear) of sorted data.
fn quantile_sorted(x: &[u32], q: f64) -> f64 {
    let n = x.len();
    if n == 0 {
        return f64::NAN;
    }
    let pos = q * (n - 1) as f64;
    let i = pos.floor() as usize;
    let f = pos - i as f64;
    if i + 1 < n {
        x[i] as f64 + f * (x[i + 1] as f64 - x[i] as f64)
    } else {
        x[n - 1] as f64
    }
}

fn quantile_f64(x: &mut [f64], q: f64) -> f64 {
    x.sort_unstable_by(f64::total_cmp);
    let n = x.len();
    let pos = q * (n - 1) as f64;
    let i = pos.floor() as usize;
    let f = pos - i as f64;
    if i + 1 < n {
        x[i] + f * (x[i + 1] - x[i])
    } else {
        x[n - 1]
    }
}

/// Runs the sampler on one locus.
pub fn run(l: &Locus, cfg: &Cfg, seed: u64, s: &mut Scratch) -> LocusSummary {
    let t = l.tx.len();
    let nc = l.n_classes();
    let mut rng = PhRng::seed_from_u64(seed);
    let kept = cfg.sweeps as usize;
    let total = cfg.burnin as usize + kept;
    let alpha = cfg.alpha;
    let reset = |v: &mut Vec<f64>, n: usize| {
        v.clear();
        v.resize(n, 0.0);
    };
    reset(&mut s.lt, 2 * t);
    reset(&mut s.w, 2 * t);
    reset(&mut s.beta, t);
    reset(&mut s.lb, t);
    reset(&mut s.m, t);
    reset(&mut s.rb, 2 * t);
    reset(&mut s.inf, 2 * t);
    s.n.clear();
    s.n.resize(2 * t, 0);
    s.hist.clear();
    s.hist.reserve(kept * 2 * t);
    s.rho.clear();
    s.tv.clear();
    s.eh.clear();
    s.eh.extend(l.cls_d.iter().map(|&d| ((0.5 * d).exp(), (-0.5 * d).exp())));
    // init: beta from unique-class counts
    for v in s.beta.iter_mut() {
        *v = 1.0;
    }
    for c in 0..nc {
        let tx = l.class_tx(c);
        if tx.len() == 1 {
            s.beta[tx[0] as usize] += l.cls_n[c] as f64;
        }
    }
    let z: f64 = s.beta.iter().sum();
    for i in 0..t {
        s.beta[i] /= z;
        s.lt[i] = s.beta[i].ln();
        s.lt[t + i] = s.beta[i].ln();
    }
    let mut rho = 0.5f64;
    let has_ev: Vec<bool> = (0..nc).map(|c| l.cls_la[c] + l.cls_lb[c] > 0.0).collect();

    for it in 0..total {
        let keep = it >= cfg.burnin as usize;
        // ---- 1. allocation --------------------------------------------------------------
        let (lra, lrb) = (rho.ln(), (-rho).ln_1p());
        let mut mx = f64::NEG_INFINITY;
        for i in 0..t {
            mx = mx.max(s.lt[i] + lra).max(s.lt[t + i] + lrb);
        }
        for i in 0..t {
            s.w[i] = (s.lt[i] + lra - mx).exp();
            s.w[t + i] = (s.lt[t + i] + lrb - mx).exp();
        }
        s.n.iter_mut().for_each(|x| *x = 0);
        let (wa, wb) = s.w.split_at(t);
        let (na, nb) = s.n.split_at_mut(t);
        let (rba, rbb) = s.rb.split_at_mut(t);
        let (infa, infb) = s.inf.split_at_mut(t);
        #[allow(clippy::needless_range_loop)]
        for c in 0..nc {
            let tx = l.class_tx(c);
            let k = tx.len();
            let (ea, eb) = s.eh[c];
            s.p.clear();
            s.p.extend(tx.iter().map(|&ti| wa[ti as usize] * ea));
            s.p.extend(tx.iter().map(|&ti| wb[ti as usize] * eb));
            let mut sum: f64 = s.p.iter().sum();
            if sum.is_nan() || sum <= 1e-250 || !sum.is_finite() {
                // slow path: exact log-space normalisation
                let d = l.cls_d[c];
                let (pa, pb) = s.p.split_at_mut(k);
                for ((x, y), &ti) in pa.iter_mut().zip(pb.iter_mut()).zip(tx) {
                    *x = lra + s.lt[ti as usize] + 0.5 * d;
                    *y = lrb + s.lt[t + ti as usize] - 0.5 * d;
                }
                let m2 = s.p.iter().fold(f64::NEG_INFINITY, |m, &v| m.max(v));
                for v in s.p.iter_mut() {
                    *v = (*v - m2).exp();
                }
                sum = s.p.iter().sum();
            }
            let inv = 1.0 / sum;
            for v in s.p.iter_mut() {
                *v *= inv;
            }
            let nr = l.cls_n[c];
            let (pa, pb) = s.p.split_at(k);
            if keep {
                let nf = nr as f64;
                for ((&x, &y), &ti) in pa.iter().zip(pb).zip(tx) {
                    rba[ti as usize] += nf * x;
                    rbb[ti as usize] += nf * y;
                }
                if has_ev[c] {
                    let (la, lb) = (l.cls_la[c], l.cls_lb[c]);
                    for ((&x, &y), &ti) in pa.iter().zip(pb).zip(tx) {
                        let pt = x + y;
                        infa[ti as usize] += la * pt;
                        infb[ti as usize] += lb * pt;
                    }
                }
            }
            // multinomial(nr, p) over the 2K cells: sequential binomials while many reads are
            // left, then categorical draws (sorted uniforms, one pass over the cdf)
            let mut left = nr as u64;
            let mut mass = 1.0f64;
            let mut j0 = 0usize;
            let mut add = |j: usize, x: u32| {
                if j < k {
                    na[tx[j] as usize] += x;
                } else {
                    nb[tx[j - k] as usize] += x;
                }
            };
            while left > CAT_MAX && j0 + 1 < 2 * k {
                let q = (s.p[j0] / mass).clamp(0.0, 1.0);
                mass -= s.p[j0];
                let x = if q >= 1.0 {
                    left
                } else if q <= 0.0 {
                    0
                } else {
                    Binomial::new(left, q).map(|b| b.sample(&mut rng)).unwrap_or(0)
                };
                add(j0, x as u32);
                left -= x;
                j0 += 1;
            }
            if j0 + 1 == 2 * k {
                add(j0, left as u32);
            } else if left == 1 {
                let u = unif(&mut rng) * mass;
                let mut acc = 0.0;
                let mut j = j0;
                while j + 1 < 2 * k {
                    acc += s.p[j];
                    if u < acc {
                        break;
                    }
                    j += 1;
                }
                add(j, 1);
            } else if left > 1 {
                let mut us = [0.0f64; CAT_MAX as usize];
                let us = &mut us[..left as usize];
                for u in us.iter_mut() {
                    *u = unif(&mut rng) * mass;
                }
                // insertion sort (at most CAT_MAX values)
                for i in 1..us.len() {
                    let v = us[i];
                    let mut q = i;
                    while q > 0 && us[q - 1] > v {
                        us[q] = us[q - 1];
                        q -= 1;
                    }
                    us[q] = v;
                }
                let mut acc = s.p[j0];
                let mut j = j0;
                let mut cnt = 0u32;
                for &u in us.iter() {
                    while u >= acc && j + 1 < 2 * k {
                        if cnt > 0 {
                            add(j, cnt);
                            cnt = 0;
                        }
                        j += 1;
                        acc += s.p[j];
                    }
                    cnt += 1;
                }
                add(j, cnt);
            }
        }
        // ---- 2. rho -------------------------------------------------------------------------
        let na: u64 = s.n[..t].iter().map(|&x| x as u64).sum();
        let nb: u64 = s.n[t..].iter().map(|&x| x as u64).sum();
        rho = beta(&mut rng, cfg.a_rho + na as f64, cfg.a_rho + nb as f64);
        // ---- 3. table counts ----------------------------------------------------------------
        for i in 0..t {
            let ab = alpha * s.beta[i];
            s.m[i] = (crt(&mut rng, s.n[i] as u64, ab) + crt(&mut rng, s.n[t + i] as u64, ab)) as f64;
        }
        // ---- 4. beta ------------------------------------------------------------------------
        let g = cfg.gamma / t as f64;
        log_dirichlet(&mut rng, s.m.iter().map(|&m| g + m), &mut s.lb);
        let mut z = 0.0;
        for i in 0..t {
            s.beta[i] = s.lb[i].exp().max(1e-300);
            z += s.beta[i];
        }
        for b in s.beta.iter_mut() {
            *b /= z;
        }
        // ---- 5. theta (last) ----------------------------------------------------------------
        {
            let (lta, ltb) = s.lt.split_at_mut(t);
            let (beta, n) = (&s.beta, &s.n);
            log_dirichlet(&mut rng, (0..t).map(|i| alpha * beta[i] + n[i] as f64), lta);
            log_dirichlet(&mut rng, (0..t).map(|i| alpha * beta[i] + n[t + i] as f64), ltb);
        }
        if keep {
            s.hist.extend_from_slice(&s.n);
            s.rho.push(rho);
            let tv: f64 = (0..t).map(|i| (s.lt[i].exp() - s.lt[t + i].exp()).abs()).sum::<f64>() * 0.5;
            s.tv.push(tv);
        }
    }

    // ---- summaries -----------------------------------------------------------------------
    let kf = kept.max(1) as f64;
    let mut out = LocusSummary { n_reads: l.n_reads(), ..Default::default() };
    for i in 0..t {
        let mut q = [0.0; 4];
        for (h, qq) in [(0usize, 0usize), (1, 2)] {
            s.col.clear();
            s.col.extend((0..kept).map(|sw| s.hist[sw * 2 * t + h * t + i]));
            s.col.sort_unstable();
            q[qq] = quantile_sorted(&s.col, 0.05);
            q[qq + 1] = quantile_sorted(&s.col, 0.95);
        }
        let (ea, eb) = (s.inf[i] / kf, s.inf[t + i] / kf);
        out.tx.push(TxSummary {
            a_mean: s.rb[i] / kf,
            a_lo: q[0],
            a_hi: q[1],
            b_mean: s.rb[t + i] / kf,
            b_lo: q[2],
            b_hi: q[3],
            phase_uncertainty: phase_uncertainty(1.0 + ea, 1.0 + eb),
            evid_a: ea,
            evid_b: eb,
        });
    }
    if kept > 0 {
        out.rho_mean = s.rho.iter().sum::<f64>() / kf;
        out.p_rho_gt_half = s.rho.iter().filter(|&&r| r > 0.5).count() as f64 / kf;
        out.tv_mean = s.tv.iter().sum::<f64>() / kf;
        out.rho_lo = quantile_f64(&mut s.rho, 0.05);
        out.rho_hi = quantile_f64(&mut s.rho, 0.95);
    }
    out.evid_a = l.cls_la.iter().sum();
    out.evid_b = l.cls_lb.iter().sum();
    out.phase_uncertainty = phase_uncertainty(1.0 + out.evid_a, 1.0 + out.evid_b);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mean_var(x: &[f64]) -> (f64, f64) {
        let n = x.len() as f64;
        let m = x.iter().sum::<f64>() / n;
        (m, x.iter().map(|v| (v - m) * (v - m)).sum::<f64>() / (n - 1.0))
    }

    #[test]
    fn gamma_and_beta_moments() {
        let mut r = PhRng::seed_from_u64(1);
        for &a in &[0.01, 0.3, 1.0, 2.5, 40.0] {
            let x: Vec<f64> = (0..200_000).map(|_| log_gamma(&mut r, a).exp()).collect();
            let (m, v) = mean_var(&x);
            assert!((m - a).abs() < 4.0 * (a / 200_000f64).sqrt() + 1e-3, "gamma({a}) mean {m}");
            assert!((v - a).abs() / a < 0.05, "gamma({a}) var {v}");
        }
        let x: Vec<f64> = (0..200_000).map(|_| beta(&mut r, 20.0, 30.0)).collect();
        let (m, v) = mean_var(&x);
        assert!((m - 0.4).abs() < 1e-3);
        let ev = 20.0 * 30.0 / (50.0f64 * 50.0 * 51.0);
        assert!((v - ev).abs() / ev < 0.03);
    }

    #[test]
    fn dirichlet_moments_tiny_shapes() {
        let mut r = PhRng::seed_from_u64(2);
        let a = [1e-4, 0.5, 2.0, 1e-4];
        let a0: f64 = a.iter().sum();
        let mut out = [0.0; 4];
        let mut sums = [0.0; 4];
        let n = 100_000;
        for _ in 0..n {
            log_dirichlet(&mut r, a.iter().copied(), &mut out);
            let s: f64 = out.iter().map(|x| x.exp()).sum();
            assert!((s - 1.0).abs() < 1e-9);
            for i in 0..4 {
                sums[i] += out[i].exp();
            }
        }
        for i in 0..4 {
            let m = sums[i] / n as f64;
            assert!((m - a[i] / a0).abs() < 0.01, "dir mean {i}: {m} vs {}", a[i] / a0);
        }
    }

    #[test]
    fn lgamma_shift_matches_sum() {
        for &a in &[1e-12, 1e-6, 0.3, 1.0, 7.5] {
            for &x in &[1.0, 3.0, 15.0, 16.0, 100.0, 12345.0] {
                // reference: ln Gamma(x+a) - ln Gamma(x) via sum to x+N plus asymptotic is circular;
                // check the recurrence G(x+1) - G(x) = ln((x+a)/x) instead
                let g0 = lgamma_shift(x, a);
                let g1 = lgamma_shift(x + 1.0, a);
                let want = (a / x).ln_1p();
                assert!((g1 - g0 - want).abs() < 1e-12 + 1e-9 * want.abs(), "a={a} x={x}");
            }
        }
        // integer a = 1: lnGamma(x+1) - lnGamma(x) = ln x
        assert!((lgamma_shift(1000.0, 1.0) - 1000f64.ln()).abs() < 1e-12);
        assert!((lgamma_shift(20.0, 2.0) - (20f64 * 21.0).ln()).abs() < 1e-12);
    }

    #[test]
    fn crt_matches_direct() {
        let mut r = PhRng::seed_from_u64(3);
        for &(n, ab) in &[(1u64, 0.5), (30, 0.01), (500, 0.7), (5000, 1e-4), (3000, 3.0)] {
            let reps = 40_000;
            let fast: Vec<f64> = (0..reps).map(|_| crt(&mut r, n, ab) as f64).collect();
            let direct: Vec<f64> =
                (0..reps).map(|_| (0..n).filter(|&j| unif(&mut r) < ab / (ab + j as f64)).count() as f64).collect();
            let exact_mean: f64 = (0..n).map(|j| ab / (ab + j as f64)).sum();
            let exact_var: f64 = (0..n).map(|j| ab / (ab + j as f64) * (1.0 - ab / (ab + j as f64))).sum();
            let (mf, vf) = mean_var(&fast);
            let (md, _) = mean_var(&direct);
            let se = (exact_var / reps as f64).sqrt();
            assert!((mf - exact_mean).abs() < 5.0 * se + 1e-9, "n={n} ab={ab}: {mf} vs {exact_mean}");
            assert!((md - exact_mean).abs() < 5.0 * se + 1e-9);
            if exact_var * reps as f64 > 2000.0 {
                assert!((vf - exact_var).abs() / exact_var < 0.06, "var n={n} ab={ab}: {vf} vs {exact_var}");
            }
        }
    }

    /// Fully phased toy gene: every read carries strong het evidence, so the A/B split per
    /// transcript is known.
    #[test]
    fn known_answer_fully_phased() {
        // t1: 300 A reads, 100 B reads; t2: 0 A, 200 B
        let mut reads: Vec<(Vec<&str>, f64)> = Vec::new();
        let strong = 20.0; // ~ log-lambda ratio of several Q30 sites
        reads.extend(std::iter::repeat_n((vec!["t1"], strong), 300));
        reads.extend(std::iter::repeat_n((vec!["t1"], -strong), 100));
        reads.extend(std::iter::repeat_n((vec!["t2"], -strong), 200));
        let l = Locus::from_reads("chr1", reads, |_| "g".to_string());
        let cfg = Cfg { burnin: 200, sweeps: 2000, ..Default::default() };
        let mut s = Scratch::default();
        let r = run(&l, &cfg, 11, &mut s);
        let t1 = &r.tx[0];
        let t2 = &r.tx[1];
        assert!((t1.a_mean - 300.0).abs() < 0.5, "{t1:?}");
        assert!((t1.b_mean - 100.0).abs() < 0.5);
        assert!(t2.a_mean < 0.5 && (t2.b_mean - 200.0).abs() < 0.5);
        assert_eq!((t1.a_lo, t1.a_hi), (300.0, 300.0));
        // rho | n ~ Beta(20 + 300, 20 + 300) -> mean 0.5
        assert!((r.rho_mean - 0.5).abs() < 0.01, "{}", r.rho_mean);
        // gene level: Beta(301, 301) -> SD 0.0204 / SD(U) 0.2887
        assert!((r.phase_uncertainty - phase_uncertainty(301.0, 301.0)).abs() < 1e-3);
        assert!((r.phase_uncertainty - 0.0706).abs() < 1e-3);
        assert!(t1.phase_uncertainty < t2.phase_uncertainty + 0.1);
        // evidence sums: l(A) ~ 1 for A reads
        assert!((t1.evid_a - 300.0).abs() < 1e-3 && (t1.evid_b - 100.0).abs() < 1e-3);
        // rho posterior: exact Beta(320, 320) quantiles ~ 0.5 +- 0.0325
        assert!((r.rho_hi - r.rho_lo - 2.0 * 1.645 * (320.0f64 * 320.0 / (640.0 * 640.0 * 641.0)).sqrt()).abs() < 0.01);
    }

    /// Unphased reads split by the prior; with no evidence anywhere, rho stays at its prior and
    /// A/B means are symmetric.
    #[test]
    fn no_evidence_symmetric() {
        let reads: Vec<(Vec<&str>, f64)> = std::iter::repeat_n((vec!["t1"], 0.0), 400)
            .chain(std::iter::repeat_n((vec!["t1", "t2"], 0.0), 100))
            .collect();
        let l = Locus::from_reads("chr1", reads, |_| "g".to_string());
        let mut s = Scratch::default();
        let r = run(&l, &Cfg { burnin: 500, sweeps: 4000, ..Default::default() }, 5, &mut s);
        let a: f64 = r.tx.iter().map(|x| x.a_mean).sum();
        let b: f64 = r.tx.iter().map(|x| x.b_mean).sum();
        assert!((a + b - 500.0).abs() < 1e-6);
        assert!((a / 500.0 - 0.5).abs() < 0.03, "{a} {b}");
        assert!((r.phase_uncertainty - 1.0).abs() < 1e-12);
        assert!(r.tx.iter().all(|x| (x.phase_uncertainty - 1.0).abs() < 1e-12));
    }
}
