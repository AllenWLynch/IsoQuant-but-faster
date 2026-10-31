//! Scheduler tests against the mock pipeline (U21 pass criteria a-c).
//!
//! Timing- and RSS-sensitive checks run in a child process (`child_*` tests are `#[ignore]`
//! and re-invoked by their parent test through `current_exe`), so concurrently running unit
//! tests cannot disturb them.

use crate::engine::{execute, EngineOpts};
use crate::mock::{MockBlock, MockChr, MockPipeline, MockResult, MockSpec};
use crate::profile::{Kind, Report};
use crate::units::SplitMode;

fn run(spec: &MockSpec, threads: usize, split: SplitMode, budget: u64) -> (MockResult, Report) {
    let p = MockPipeline::new(spec.clone());
    let opts = EngineOpts { threads, budget, split, sample_rss: false, verbose: false };
    let rep = execute(&p, &opts).expect("run");
    (p.take_result(), rep)
}

/// A genome-like spec with a giant block on chromosome 2 (nested parallelism path).
fn det_spec(seed: u64) -> MockSpec {
    let mut s = MockSpec::synthetic(seed, 6, 250, 3000, 0.35);
    s.dup_rate = 0.01;
    s.par_threshold = 5_000;
    // giant block sharing a gene with both neighbours
    let b = &mut s.chrs[1].blocks;
    let g = b[100].genes.clone();
    b.insert(101, MockBlock { genes: g, reads: 60_000 });
    s.chrs.push(MockChr { name: "MT".into(), blocks: vec![MockBlock { genes: vec![], reads: 10 }] });
    s.chrs.push(MockChr { name: "empty".into(), blocks: vec![] });
    s
}

const BIG: u64 = u64::MAX / 4;

#[test]
fn a_determinism_threads_and_splits() {
    for seed in [1u64, 7] {
        let spec = det_spec(seed);
        let (base, rep) = run(&spec, 1, SplitMode::Packed, BIG);
        assert!(base.gtf.lines().count() > 1000, "mock produced too little output");
        assert!(base.counts.contains("__ambiguous"));
        let units_packed = rep.n_units;
        for threads in [1, 4, 16] {
            for split in [SplitMode::Packed, SplitMode::Max, SplitMode::WholeChr] {
                let (r, rep) = run(&spec, threads, split, BIG);
                assert_eq!(r, base, "outputs differ: seed {seed} threads {threads} split {split:?}");
                if split == SplitMode::WholeChr {
                    assert_eq!(rep.n_units, spec.chrs.iter().filter(|c| !c.blocks.is_empty()).count());
                }
                if split == SplitMode::Max {
                    assert!(rep.n_units > units_packed);
                }
            }
        }
        // a tight budget changes the schedule, never the result
        let (r, rep) = run(&spec, 8, SplitMode::Max, 2 << 20);
        assert_eq!(r, base);
        assert!(rep.over_budget_admissions > 0);
    }
}

#[test]
fn a_determinism_check_has_teeth() {
    // splitting blocks that share a gene breaks detected_known threading
    let spec = det_spec(1);
    let (base, _) = run(&spec, 4, SplitMode::Packed, BIG);
    let (bad, _) = run(&spec, 4, SplitMode::UnsafePerBlock, BIG);
    assert_ne!(bad.gtf, base.gtf);
}

#[test]
fn barrier_polya_decision_and_duplicates() {
    let mut spec = MockSpec::synthetic(3, 3, 50, 2000, 0.2);
    spec.dup_rate = 0.05;
    spec.polya_rate = 0.9;
    let (hi, _) = run(&spec, 4, SplitMode::Packed, BIG);
    assert!(hi.gtf.contains("# requires_polya true"));
    spec.polya_rate = 0.5;
    let (lo, _) = run(&spec, 4, SplitMode::Packed, BIG);
    assert!(lo.gtf.contains("# requires_polya false"));
    // totals exclude suspended duplicate copies
    let total_reads: u64 = spec.chrs.iter().flat_map(|c| &c.blocks).map(|b| b.reads as u64).sum();
    let line = lo.counts.lines().find(|l| l.starts_with("__totals")).unwrap();
    let t: u64 = line.split('\t').nth(1).unwrap().parse().unwrap();
    assert!(t < total_reads && t > total_reads * 9 / 10, "{t} vs {total_reads}");
    // duplicates are resolved per chromosome across units: no duplicates → no ambiguous
    spec.dup_rate = 0.0;
    let (nodup, _) = run(&spec, 4, SplitMode::Max, BIG);
    assert!(nodup.counts.contains("__ambiguous\t0\n"));
}

#[test]
fn errors_propagate() {
    struct Failing(MockPipeline);
    impl crate::Pipeline for Failing {
        type Store = <MockPipeline as crate::Pipeline>::Store;
        type UnitState = <MockPipeline as crate::Pipeline>::UnitState;
        type Global = <MockPipeline as crate::Pipeline>::Global;
        type BlockOut = <MockPipeline as crate::Pipeline>::BlockOut;
        type ChrOut = <MockPipeline as crate::Pipeline>::ChrOut;
        fn n_chromosomes(&self) -> usize {
            self.0.n_chromosomes()
        }
        fn chr_name(&self, c: usize) -> &str {
            self.0.chr_name(c)
        }
        fn chr_weight(&self, c: usize) -> u64 {
            self.0.chr_weight(c)
        }
        fn pass_a(&self, c: usize) -> Result<Vec<iq_types::Block>, String> {
            self.0.pass_a(c)
        }
        fn collect_block(&self, c: usize, b: usize, bl: &iq_types::Block) -> Result<Self::Store, String> {
            if c == 1 && b == 7 {
                return Err("boom".into());
            }
            self.0.collect_block(c, b, bl)
        }
        fn store_bytes(&self, s: &Self::Store) -> u64 {
            self.0.store_bytes(s)
        }
        fn barrier_chr(
            &self,
            c: usize,
            b: &[iq_types::Block],
            s: &mut [Self::Store],
        ) -> Result<crate::ReadTotals, String> {
            self.0.barrier_chr(c, b, s)
        }
        fn decide(&self, t: crate::ReadTotals) -> Self::Global {
            self.0.decide(t)
        }
        fn construct_block(
            &self,
            c: usize,
            bi: usize,
            b: &iq_types::Block,
            s: Self::Store,
            st: &mut Self::UnitState,
            g: &Self::Global,
        ) -> Result<Self::BlockOut, String> {
            self.0.construct_block(c, bi, b, s, st, g)
        }
        fn merge_chr(&self, c: usize, b: &[iq_types::Block], o: Vec<Self::BlockOut>) -> Result<Self::ChrOut, String> {
            self.0.merge_chr(c, b, o)
        }
        fn finish(&self, o: Vec<Self::ChrOut>, t: crate::ReadTotals) -> Result<(), String> {
            self.0.finish(o, t)
        }
    }
    let p = Failing(MockPipeline::new(MockSpec::synthetic(1, 3, 40, 100, 0.2)));
    let e = execute(&p, &EngineOpts { threads: 4, ..Default::default() }).unwrap_err();
    assert!(e.contains("boom"), "{e}");
}

#[test]
fn empty_run() {
    let (r, rep) = run(&MockSpec::default(), 2, SplitMode::Packed, BIG);
    assert_eq!(rep.n_units, 0);
    assert!(r.counts.starts_with("feature_id\tcount\n"));
}

// ------------------------------------------------------------------------------------------
// (b) memory budget under a synthetic giant block
// ------------------------------------------------------------------------------------------

const MIB: u64 = 1 << 20;
const READ_BYTES: u64 = 16; // size_of::<MockRead>()

/// A giant block (2 M reads, 200 MiB store) next to 300 small units.
fn giant_spec(store_b: u64, construct_b: u64) -> MockSpec {
    assert_eq!(std::mem::size_of::<crate::mock::MockRead>() as u64, READ_BYTES);
    let mut s = MockSpec::synthetic(5, 4, 100, 20_000, 0.2);
    s.chrs[0].blocks.insert(50, MockBlock { genes: vec![999_999], reads: 2_000_000 });
    s.store_bytes_per_read = store_b;
    s.construct_bytes_per_read = construct_b;
    s.est_collect_per_read = store_b + READ_BYTES + 8;
    s.est_construct_per_read = construct_b + 8;
    s.par_threshold = 50_000;
    s.collect_ns_per_read = 50;
    s.construct_ns_per_read = 50;
    s
}

fn total_store_bytes(s: &MockSpec) -> u64 {
    s.chrs.iter().flat_map(|c| &c.blocks).map(|b| b.reads as u64 * (s.store_bytes_per_read + READ_BYTES)).sum()
}

#[test]
fn b_budget_accounting() {
    // stores fit (no spill needed): accounted peak never exceeds the budget
    let spec = giant_spec(100, 100);
    let stores = total_store_bytes(&spec);
    let budget = stores + 300 * MIB;
    let (_, rep) = run(&spec, 8, SplitMode::Packed, budget);
    assert_eq!(rep.over_budget_admissions, 0);
    assert!(rep.peak_accounted <= budget, "{} > {budget}", rep.peak_accounted);

    // giant unit larger than the budget: admitted alone (progress rule), nothing overlaps it
    let budget = 150 * MIB;
    let mut spec = giant_spec(100, 0);
    spec.spill_dir = Some(tmpdir("b_giant"));
    let (_, rep) = run(&spec, 8, SplitMode::Packed, budget);
    assert_eq!(rep.over_budget_admissions, 1, "{}", rep.summary());
    let g = rep.tasks.iter().find(|t| t.kind == Kind::Collect && t.reads >= 2_000_000).unwrap();
    for t in rep.tasks.iter().filter(|t| matches!(t.kind, Kind::Collect | Kind::PassA | Kind::Construct)) {
        if std::ptr::eq(t, g) {
            continue;
        }
        let overlap = t.start < g.end && g.start < t.end;
        assert!(!overlap || t.kind != Kind::Collect, "unit {} overlaps the giant collect", t.id);
    }
}

fn tmpdir(name: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("iq_sched_test_{}_{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

static CHILD_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

/// Runs `test` (an `#[ignore]` child test) in a fresh process, one child at a time; returns
/// its output.
fn child(test: &str) -> String {
    let _g = CHILD_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let exe = std::env::current_exe().unwrap();
    let out = std::process::Command::new(exe)
        .args(["--ignored", "--exact", test, "--nocapture", "--test-threads", "1"])
        .output()
        .unwrap();
    let so = String::from_utf8_lossy(&out.stdout).into_owned();
    let se = String::from_utf8_lossy(&out.stderr).into_owned();
    assert!(out.status.success(), "child {test} failed:\n{so}\n{se}");
    so + &se
}

fn field(out: &str, key: &str) -> f64 {
    out.lines()
        .find_map(|l| l.strip_prefix(key))
        .unwrap_or_else(|| panic!("no {key} in child output:\n{out}"))
        .trim()
        .parse()
        .unwrap()
}

#[test]
#[ignore = "child process of b_budget_rss"]
fn child_budget_rss() {
    // the giant unit fits (≈ 350 MiB est) but all stores (≈ 460 MiB) exceed the 400 MiB
    // budget: spill keeps anonymous memory bounded
    let mut spec = giant_spec(150, 60);
    spec.spill_dir = Some(tmpdir("rss"));
    let stores = total_store_bytes(&spec);
    let budget = 400 * MIB;
    let p = MockPipeline::new(spec.clone());
    let opts = EngineOpts { threads: 8, budget, split: SplitMode::Packed, sample_rss: true, verbose: false };
    let rep = execute(&p, &opts).unwrap();
    println!("{}", rep.summary());
    println!("STORES {stores}");
    println!("BUDGET {budget}");
    println!("PEAK_RSS {}", rep.peak_rss);
    println!("BASELINE {}", rep.baseline_rss);
    println!("PEAK_ACC {}", rep.peak_accounted);
    println!("OVER {}", rep.over_budget_admissions);
    let _ = std::fs::remove_dir_all(spec.spill_dir.unwrap());
}

#[test]
fn b_budget_rss() {
    let out = child("tests::child_budget_rss");
    eprintln!("{out}");
    let (stores, budget, peak, base, acc) =
        (field(&out, "STORES"), field(&out, "BUDGET"), field(&out, "PEAK_RSS"), field(&out, "BASELINE"), field(&out, "PEAK_ACC"));
    assert!(stores > budget, "scenario must need spill");
    assert_eq!(field(&out, "OVER"), 0.0);
    assert!(acc <= budget);
    // RSS above baseline stays within the budget (+ allocator/page-cache slack)
    assert!(peak - base <= budget * 1.1 + 32.0 * MIB as f64, "peak RSS {peak} baseline {base} budget {budget}");
}

// ------------------------------------------------------------------------------------------
// (c) core utilisation when one chromosome dominates
// ------------------------------------------------------------------------------------------

/// chr1 holds ~85% of the reads (one 300 k-read block + many small units); 9 small chrs.
fn dominant_spec() -> MockSpec {
    let mut s = MockSpec::synthetic(11, 10, 60, 2_000, 0.3);
    s.chrs[0].blocks = MockSpec::synthetic(12, 1, 1500, 4_000, 0.3).chrs.remove(0).blocks;
    s.chrs[0].blocks.insert(700, MockBlock { genes: vec![888_888], reads: 300_000 });
    s.collect_ns_per_read = 3000;
    s.construct_ns_per_read = 3000;
    s.passa_ns_per_read = 20;
    s.par_threshold = 20_000;
    s
}

#[test]
#[ignore = "child process of c_utilisation"]
fn child_utilisation() {
    let threads = 4;
    for (name, split) in [("PACKED", SplitMode::Packed), ("WHOLECHR", SplitMode::WholeChr)] {
        let p = MockPipeline::new(dominant_spec());
        let opts = EngineOpts { threads, budget: BIG, split, sample_rss: false, verbose: false };
        let rep = execute(&p, &opts).unwrap();
        println!("{}", rep.summary());
        println!("{name}_UTIL {}", rep.utilisation());
        println!("{name}_WALL {}", rep.wall.as_secs_f64());
        println!("{name}_CPU {}", rep.cpu_seconds);
    }
}

#[test]
fn c_utilisation() {
    if std::thread::available_parallelism().map_or(1, |n| n.get()) < 4 {
        eprintln!("fewer than 4 cores; skipped");
        return;
    }
    let out = child("tests::child_utilisation");
    eprintln!("{out}");
    let u = field(&out, "PACKED_UTIL");
    let w = field(&out, "WHOLECHR_UTIL");
    // sub-chromosome units + nested parallelism keep most of 4 cores busy although one
    // chromosome holds most of the work; 4.1's granularity (one task per chromosome, here
    // still with nested parallelism inside the giant block) leaves most cores idle. Other
    // tests may run concurrently, hence the margin.
    assert!(u >= 0.65, "utilisation {u}"); // 0.93 on an idle-ish node; shared nodes add noise
    assert!(u >= 1.8 * w, "utilisation {u} vs whole-chr {w}");
}
