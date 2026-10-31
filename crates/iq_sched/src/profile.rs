//! `--profile`: RSS sampler, per-phase wall/CPU, admission estimate-vs-actual log.

use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Arc;
use std::time::{Duration, Instant};

const PAGE: u64 = 4096;
/// `sysconf(_SC_CLK_TCK)`; 100 on every Linux x86-64/aarch64 kernel config we target.
const CLK_TCK: f64 = 100.0;

/// Current resident set size (bytes), from `/proc/self/statm`.
pub fn rss_bytes() -> u64 {
    std::fs::read_to_string("/proc/self/statm")
        .ok()
        .and_then(|s| s.split_whitespace().nth(1).and_then(|v| v.parse::<u64>().ok()))
        .map_or(0, |p| p * PAGE)
}

/// Process CPU time (user + system, all threads), seconds, from `/proc/self/stat`.
pub fn cpu_seconds() -> f64 {
    let Ok(s) = std::fs::read_to_string("/proc/self/stat") else { return 0.0 };
    // comm may contain spaces: fields after the closing paren
    let Some(rest) = s.rfind(')').map(|i| &s[i + 2..]) else { return 0.0 };
    let f: Vec<&str> = rest.split_whitespace().collect();
    // rest[0] is field 3 (state); utime = field 14, stime = field 15
    let ut: f64 = f.get(11).and_then(|v| v.parse().ok()).unwrap_or(0.0);
    let st: f64 = f.get(12).and_then(|v| v.parse().ok()).unwrap_or(0.0);
    (ut + st) / CLK_TCK
}

/// Background thread sampling RSS every `period`.
pub struct RssSampler {
    stop: Arc<AtomicBool>,
    peak: Arc<AtomicU64>,
    handle: Option<std::thread::JoinHandle<Vec<(f32, u64)>>>,
    pub baseline: u64,
}

impl RssSampler {
    pub fn start(period: Duration) -> Self {
        let stop = Arc::new(AtomicBool::new(false));
        let baseline = rss_bytes();
        let peak = Arc::new(AtomicU64::new(baseline));
        let (s2, p2) = (stop.clone(), peak.clone());
        let t0 = Instant::now();
        let handle = std::thread::Builder::new()
            .name("iqrs-rss".into())
            .spawn(move || {
                let mut series = Vec::new();
                while !s2.load(Ordering::Relaxed) {
                    let r = rss_bytes();
                    p2.fetch_max(r, Ordering::Relaxed);
                    series.push((t0.elapsed().as_secs_f32(), r));
                    std::thread::sleep(period);
                }
                series
            })
            .ok();
        RssSampler { stop, peak, handle, baseline }
    }

    pub fn peak(&self) -> u64 {
        self.peak.load(Ordering::Relaxed).max(rss_bytes())
    }

    /// Stops the sampler; returns (peak, time series of (seconds, rss)).
    pub fn stop(mut self) -> (u64, Vec<(f32, u64)>) {
        self.stop.store(true, Ordering::Relaxed);
        let series = self.handle.take().and_then(|h| h.join().ok()).unwrap_or_default();
        (self.peak(), series)
    }
}

impl Drop for RssSampler {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::Relaxed);
    }
}

/// Task kinds, also the phases of the run.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Kind {
    PassA,
    Collect,
    Spill,
    Barrier,
    Construct,
    Merge,
    Finish,
}

pub const KINDS: [Kind; 7] =
    [Kind::PassA, Kind::Collect, Kind::Spill, Kind::Barrier, Kind::Construct, Kind::Merge, Kind::Finish];

/// One finished task.
#[derive(Clone, Debug)]
pub struct TaskRec {
    pub kind: Kind,
    /// Unit id (Collect/Construct/Spill) or chromosome (others).
    pub id: usize,
    pub start: Duration,
    pub end: Duration,
    pub est_bytes: u64,
    /// Measured resident bytes left behind (Collect: Σ store bytes).
    pub actual_bytes: u64,
    pub reads: u64,
}

/// Run report, returned by [`crate::engine::execute`].
#[derive(Clone, Debug, Default)]
pub struct Report {
    pub threads: usize,
    pub budget: u64,
    pub wall: Duration,
    pub cpu_seconds: f64,
    pub peak_rss: u64,
    pub baseline_rss: u64,
    pub peak_accounted: u64,
    pub over_budget_admissions: u64,
    pub n_units: usize,
    pub tasks: Vec<TaskRec>,
    pub rss_series: Vec<(f32, u64)>,
}

impl Report {
    /// Average busy cores over the run (CPU seconds / wall seconds).
    pub fn avg_cores(&self) -> f64 {
        self.cpu_seconds / self.wall.as_secs_f64().max(1e-9)
    }

    /// Fraction of the pool used: avg cores / threads.
    pub fn utilisation(&self) -> f64 {
        self.avg_cores() / self.threads.max(1) as f64
    }

    /// (first start, last end, Σ task wall, n) per kind.
    pub fn phase(&self, k: Kind) -> Option<(Duration, Duration, Duration, usize)> {
        let mut it = self.tasks.iter().filter(|t| t.kind == k).peekable();
        it.peek()?;
        let (mut s, mut e, mut sum, mut n) = (Duration::MAX, Duration::ZERO, Duration::ZERO, 0);
        for t in it {
            s = s.min(t.start);
            e = e.max(t.end);
            sum += t.end - t.start;
            n += 1;
        }
        Some((s, e, sum, n))
    }

    /// Human-readable summary (`--profile`).
    pub fn summary(&self) -> String {
        use std::fmt::Write;
        let gb = |b: u64| b as f64 / (1u64 << 30) as f64;
        let mut s = String::new();
        let _ = writeln!(
            s,
            "[profile] wall {:.2}s  cpu {:.2}s  avg cores {:.2}/{} ({:.0}%)  peak RSS {:.3} GiB (baseline {:.3})  units {}",
            self.wall.as_secs_f64(),
            self.cpu_seconds,
            self.avg_cores(),
            self.threads,
            100.0 * self.utilisation(),
            gb(self.peak_rss),
            gb(self.baseline_rss),
            self.n_units
        );
        let _ = writeln!(
            s,
            "[profile] budget {:.3} GiB  peak accounted {:.3} GiB  over-budget admissions {}",
            gb(self.budget),
            gb(self.peak_accounted),
            self.over_budget_admissions
        );
        for k in KINDS {
            if let Some((st, en, sum, n)) = self.phase(k) {
                let _ = writeln!(
                    s,
                    "[profile] {:<9} tasks {:>6}  span {:>8.2}s..{:>8.2}s  Σtask {:>9.2}s",
                    format!("{k:?}"),
                    n,
                    st.as_secs_f64(),
                    en.as_secs_f64(),
                    sum.as_secs_f64()
                );
            }
        }
        // estimate vs actual for collect (recalibration of the 160 B/read constant)
        let (mut est, mut act, mut reads) = (0u64, 0u64, 0u64);
        for t in self.tasks.iter().filter(|t| t.kind == Kind::Collect) {
            est += t.est_bytes;
            act += t.actual_bytes;
            reads += t.reads;
        }
        if reads > 0 {
            let _ = writeln!(
                s,
                "[profile] collect est {:.1} B/read vs resident {:.1} B/read over {} reads",
                est as f64 / reads as f64,
                act as f64 / reads as f64,
                reads
            );
        }
        let mut big: Vec<&TaskRec> = self.tasks.iter().filter(|t| t.kind == Kind::Collect).collect();
        big.sort_by_key(|t| std::cmp::Reverse(t.reads));
        for t in big.iter().take(5) {
            let _ = writeln!(
                s,
                "[profile]   unit {:>6}: {:>10} reads  collect {:>7.2}s  est {:>8.1} MiB  actual {:>8.1} MiB",
                t.id,
                t.reads,
                (t.end - t.start).as_secs_f64(),
                t.est_bytes as f64 / (1 << 20) as f64,
                t.actual_bytes as f64 / (1 << 20) as f64
            );
        }
        s
    }
}
