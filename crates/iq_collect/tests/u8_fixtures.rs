//! U8/U9 fixture checks over `$IQRS_FIXTURES` (default `../fixtures`), skipped when absent.
//!
//! - `fixtures_pass`: every U9 case exact.
//! - `parallel_equals_serial`: every U8 case collected serially and on a 4-thread pool.
//! - `fullchr8_throughput` (ignored): single-threaded collect timing by stage.

use iq_collect::fixture::U8Prepared;
use iq_types::fixture::schema::{ReadRecordFixture, U8CollectBlockIn};
use iq_types::fixture::{find_capture_dirs, read_cases, run_path, unit_files, FixtureCtx, Registry, Summary};
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::time::Instant;

fn root() -> Option<PathBuf> {
    let r = iq_types::fixture::fixtures_root()
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../../fixtures"));
    r.is_dir().then_some(r)
}

fn captures() -> Vec<PathBuf> {
    let Some(r) = root() else { return Vec::new() };
    let mut v = Vec::new();
    find_capture_dirs(&r, &mut v).unwrap();
    // finished captures only (a capture in progress has no manifest / a `U8.<pid>` file)
    v.retain(|d| {
        d.join("slice").is_dir() && d.join("fixture_manifest.json").is_file() && d.join("units/U8.jsonl.gz").is_file()
    });
    v.sort();
    v
}

fn u8_inputs(cap: &std::path::Path) -> Vec<(String, U8CollectBlockIn, Vec<ReadRecordFixture>)> {
    let mut out = Vec::new();
    for f in unit_files(cap, "U8").unwrap() {
        for c in read_cases::<U8CollectBlockIn, iq_types::fixture::schema::U8CollectBlockOut>(&f).unwrap() {
            let c = c.unwrap();
            out.push((c.case, c.input, c.out.reads));
        }
    }
    out
}

#[test]
fn fixtures_pass() {
    let mut reg = Registry::new();
    iq_collect::register_fixtures(&mut reg);
    let mut s = Summary::new(3);
    for cap in captures() {
        if cap.join("units/U9.jsonl.gz").is_file() {
            run_path(&reg, "U9", &cap, &mut s).unwrap();
        }
    }
    eprintln!("{s}");
    assert!(s.all_ok(), "{s}");
}

#[test]
fn parallel_equals_serial() {
    let pool = rayon::ThreadPoolBuilder::new().num_threads(4).build().unwrap();
    let mut n_cases = 0;
    let mut n_reads = 0;
    for cap in captures() {
        if cap.ends_with("fullchr8") && std::env::var("IQRS_U8_FULL").is_err() {
            continue; // the slice captures cover the same code; IQRS_U8_FULL=1 runs all
        }
        let ctx = FixtureCtx::new(&cap);
        for (case, input, _) in u8_inputs(&cap) {
            let p = U8Prepared::load(&input, &ctx).unwrap();
            let (ser, sst, _) = p.collect(&ctx, false, false);
            let (par, pst, _) = pool.install(|| p.collect(&ctx, true, false));
            assert_eq!(sst, pst, "{case}: stats differ");
            assert_eq!(ser.to_records(), par.to_records(), "{case}: parallel != serial");
            n_cases += 1;
            n_reads += ser.len();
        }
    }
    eprintln!("parallel == serial on {n_cases} blocks, {n_reads} reads");
}

/// Whole-chromosome collect throughput, single-threaded, from prefetched records (fetch and
/// index construction not timed). `cargo test --release -p iq_collect --test u8_fixtures
/// fullchr8_throughput -- --ignored --nocapture`
#[test]
#[ignore]
fn fullchr8_throughput() {
    for cap in captures().into_iter().filter(|c| c.ends_with("fullchr8")) {
        let ctx = FixtureCtx::new(&cap);
        let inputs = u8_inputs(&cap);
        let t_load = Instant::now();
        let preps: Vec<U8Prepared> = inputs.iter().map(|(_, i, _)| U8Prepared::load(i, &ctx).unwrap()).collect();
        let load_s = t_load.elapsed().as_secs_f64();
        let records: usize = preps.iter().map(|p| p.n_records()).sum();
        // untimed-stage wall clock, best of 3
        let mut best = f64::INFINITY;
        let mut stored = 0u64;
        for _ in 0..3 {
            let t = Instant::now();
            stored = 0;
            for p in &preps {
                let (st, _, _) = p.collect(&ctx, false, false);
                stored += st.len() as u64;
            }
            best = best.min(t.elapsed().as_secs_f64());
        }
        // stage breakdown (timed run)
        let mut times = iq_collect::StageTimes::default();
        let mut stats = iq_collect::CollectStats::default();
        for p in &preps {
            let (_, s, t) = p.collect(&ctx, false, true);
            times.add(&t.unwrap());
            stats.add(&s);
        }
        let us = |ns: u64| ns as f64 / 1e3 / stored as f64;
        eprintln!(
            "{}: {} blocks, {} records, {} stored (prep {:.1}s)\n  wall {:.3}s = {:.2} us/stored read, {:.0} reads/s, {:.0} records/s\n  stages us/stored read: select {:.2} u4 {:.2} u5 {:.2} u6 {:.2} u7 {:.2} store {:.2} (timed total {:.2})\n  stats {:?}",
            cap.display(),
            preps.len(),
            records,
            stored,
            load_s,
            best,
            best * 1e6 / stored as f64,
            stored as f64 / best,
            records as f64 / best,
            us(times.select),
            us(times.u4),
            us(times.u5),
            us(times.u6),
            us(times.u7),
            us(times.store),
            us(times.total()),
            stats
        );
    }
}

/// U8 → U9 on whole captures: collect every block, resolve duplicate names per
/// chromosome, and compare each store with the loader's `block/` file (post-resolution,
/// suspended reads dropped).
#[test]
fn u9_chr_resolution_matches_block_files() {
    use iq_types::fixture::schema::BlockRef;
    let (mut n_blocks, mut n_groups, mut n_susp, mut n_totals) = (0, 0, 0, 0);
    for cap in captures() {
        if cap.ends_with("fullchr8") && std::env::var("IQRS_U8_FULL").is_err() {
            continue;
        }
        let m: iq_types::serde_json::Value =
            iq_types::serde_json::from_str(&std::fs::read_to_string(cap.join("fixture_manifest.json")).unwrap())
                .unwrap();
        if m["skipped_blocks"].as_array().is_some_and(|a| !a.is_empty()) {
            eprintln!("{}: has skipped blocks; not comparable", cap.display());
            continue;
        }
        let sample = m["sample"].as_str().unwrap().to_owned();
        let slice = m["slice"].as_str().unwrap().to_owned();
        let ctx = FixtureCtx::new(&cap);
        // "Total assignments used for analysis: N, polyA tail detected in M (..)"
        let log = std::fs::read_to_string(cap.join("capture.log")).unwrap_or_default();
        let want_totals = log.lines().find_map(|l| {
            let t = l.split("Total assignments used for analysis: ").nth(1)?;
            let (n, rest) = t.split_once(',')?;
            let m = rest.split("detected in ").nth(1)?.split_whitespace().next()?;
            Some((n.trim().parse::<u64>().ok()?, m.parse::<u64>().ok()?))
        });
        let (mut tot, mut pa) = (0u64, 0u64);
        let mut by_chr: BTreeMap<String, Vec<(u32, iq_types::BlockStore)>> = BTreeMap::new();
        for (_, input, _) in u8_inputs(&cap) {
            let p = U8Prepared::load(&input, &ctx).unwrap();
            let (st, _, _) = p.collect(&ctx, true, false);
            by_chr.entry(input.chr.clone()).or_default().push((input.block_no, st));
        }
        for (chr, mut blocks) in by_chr {
            blocks.sort_by_key(|b| b.0);
            let mut stores: Vec<_> = blocks.iter_mut().map(|b| std::mem::take(&mut b.1)).collect();
            let o =
                iq_collect::multimap::resolve_chr_duplicates(&mut stores, None, iq_types::MultimapStrategy::TakeBest);
            n_groups += o.groups;
            tot += o.total;
            pa += o.polya;
            n_susp += o.suspended;
            for ((bno, _), st) in blocks.iter().zip(&stores) {
                let lb = ctx
                    .block(&BlockRef { sample: sample.clone(), slice: slice.clone(), chr: chr.clone(), block_no: *bno })
                    .unwrap();
                // 4.1's dump stores penalties as int(p * 2^20) (serialization.py:16,
                // isoform_assignment.py:523); no construct stage reads them.
                let q = |mut v: Vec<iq_types::ReadRecord>| {
                    for r in &mut v {
                        for m in r.assignment.matches.iter_mut() {
                            m.penalty = (m.penalty * (1u64 << 20) as f64).trunc() / (1u64 << 20) as f64;
                        }
                    }
                    v
                };
                let got = q(st.to_records());
                let want = q(lb.records.clone());
                if got != want {
                    let k = got.iter().zip(&want).position(|(a, b)| a != b).unwrap_or(got.len().min(lb.records.len()));
                    panic!(
                        "{}: block {chr}:{bno}: len {} vs {}; first diff at {k}:\n got {:?}\nwant {:?}",
                        cap.display(),
                        got.len(),
                        want.len(),
                        got.get(k),
                        want.get(k)
                    );
                }
                n_blocks += 1;
            }
        }
        if let Some((n, m)) = want_totals {
            assert_eq!((tot, pa), (n, m), "{}: total/polyA assignments", cap.display());
            n_totals += 1;
        }
    }
    eprintln!("{n_blocks} blocks equal to block/ files; {n_groups} groups resolved, {n_susp} suspended; total/polyA equal to 4.1 log in {n_totals} captures");
}
