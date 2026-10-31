//! Byte-budget admission control (design §5.4) and memory-budget detection.
//!
//! The controller tracks two quantities:
//! - `in_flight`: the *estimated* bytes of running tasks (released when a task ends);
//! - `resident`: the *measured* bytes of finished-but-not-yet-consumed results (collected
//!   `BlockStore`s waiting for Phase M; released as construct consumes them or spill moves
//!   them out of anonymous memory).
//!
//! A task of estimate `e` is admitted iff `resident + in_flight + e <= budget`, or nothing is
//! in flight (always make progress: a single unit larger than the budget still runs, alone).

/// Admission controller state.
#[derive(Debug, Default, Clone)]
pub struct Admission {
    pub budget: u64,
    pub in_flight: u64,
    pub resident: u64,
    /// Number of admitted tasks currently running.
    pub running: usize,
    /// Peak of `in_flight + resident` (accounted, not RSS).
    pub peak_accounted: u64,
    /// Number of admissions that exceeded the budget under the progress rule.
    pub over_budget_admissions: u64,
}

impl Admission {
    pub fn new(budget: u64) -> Self {
        Admission { budget, ..Default::default() }
    }

    /// Would a task of estimated `bytes` be admitted now?
    pub fn fits(&self, bytes: u64) -> bool {
        self.running == 0 || self.resident + self.in_flight + bytes <= self.budget
    }

    /// Admit a unit of estimated `bytes` if it fits or nothing is in flight.
    pub fn try_admit(&mut self, bytes: u64) -> bool {
        if !self.fits(bytes) {
            return false;
        }
        if self.resident + self.in_flight + bytes > self.budget {
            self.over_budget_admissions += 1;
        }
        self.in_flight += bytes;
        self.running += 1;
        self.bump_peak();
        true
    }

    /// A task admitted with estimate `est` finished, leaving `resident_bytes` of results.
    pub fn finish(&mut self, est: u64, resident_bytes: u64) {
        self.in_flight -= est;
        self.running -= 1;
        self.resident += resident_bytes;
        self.bump_peak();
    }

    /// Resident results were consumed / spilled.
    pub fn release(&mut self, bytes: u64) {
        self.resident = self.resident.saturating_sub(bytes);
    }

    /// Resident results grew (outside a task completion).
    pub fn add_resident(&mut self, bytes: u64) {
        self.resident += bytes;
        self.bump_peak();
    }

    fn bump_peak(&mut self) {
        self.peak_accounted = self.peak_accounted.max(self.in_flight + self.resident);
    }
}

/// Cgroup (v2 or v1) memory limit of this process, the minimum over its cgroup ancestors.
pub fn cgroup_limit() -> Option<u64> {
    let cg = std::fs::read_to_string("/proc/self/cgroup").ok()?;
    let mut best: Option<u64> = None;
    let mut take = |v: u64| best = Some(best.map_or(v, |b| b.min(v)));
    for line in cg.lines() {
        let mut it = line.splitn(3, ':');
        let (_, ctrl, path) = (it.next()?, it.next()?, it.next()?);
        let (root, file) = if ctrl.is_empty() {
            ("/sys/fs/cgroup", "memory.max")
        } else if ctrl.split(',').any(|c| c == "memory") {
            ("/sys/fs/cgroup/memory", "memory.limit_in_bytes")
        } else {
            continue;
        };
        let mut p = std::path::PathBuf::from(format!("{root}{path}"));
        loop {
            if let Ok(s) = std::fs::read_to_string(p.join(file)) {
                if let Ok(v) = s.trim().parse::<u64>() {
                    // v1 reports "no limit" as a huge page-aligned number
                    if v < (1u64 << 60) {
                        take(v);
                    }
                }
            }
            if p.as_os_str().len() <= root.len() || !p.pop() {
                break;
            }
        }
    }
    best
}

/// `MemAvailable` from `/proc/meminfo`, bytes.
pub fn mem_available() -> Option<u64> {
    let s = std::fs::read_to_string("/proc/meminfo").ok()?;
    let l = s.lines().find(|l| l.starts_with("MemAvailable:"))?;
    let kb: u64 = l.split_whitespace().nth(1)?.parse().ok()?;
    Some(kb * 1024)
}

/// Default budget: 0.85 × min(cgroup limit, MemAvailable) (design §5.4). Falls back to 8 GiB
/// when neither is readable.
pub fn default_budget() -> u64 {
    let lim = match (cgroup_limit(), mem_available()) {
        (Some(a), Some(b)) => a.min(b),
        (Some(a), None) | (None, Some(a)) => a,
        (None, None) => 8 << 30,
    };
    (lim as f64 * 0.85) as u64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn admission_rules() {
        let mut a = Admission::new(100);
        assert!(a.try_admit(60));
        assert!(!a.try_admit(50));
        assert!(a.try_admit(40));
        a.finish(60, 30); // leaves 30 resident
        assert_eq!((a.in_flight, a.resident), (40, 30));
        assert!(!a.try_admit(31));
        assert!(a.try_admit(30));
        a.finish(40, 0);
        a.finish(30, 0);
        // nothing running: an over-budget unit is still admitted (progress)
        assert!(a.try_admit(500));
        assert_eq!(a.over_budget_admissions, 1);
        assert!(!a.try_admit(0) || a.resident + a.in_flight <= a.budget);
        a.finish(500, 0);
        a.release(30);
        assert_eq!(a.resident, 0);
        assert_eq!(a.peak_accounted, 530);
    }

    #[test]
    fn budget_detection_is_sane() {
        let b = default_budget();
        assert!(b > 64 << 20, "budget {b}");
    }
}
