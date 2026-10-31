//! Randomized U10a cases generated from the live 4.1 Python (`tests/data/gen_u10a.py`) and
//! replayed through the fixture runners. Python is not needed at test time.
//!
//! - `u10a_random_stable`: generated with numpy's non-SIMD argsort (the order iqrs emulates), so every
//!   case must be bit-exact (no tolerance, no D-T1).
//! - `u10a_random_native`: Python's own (CPU-dependent) tie order; cases may only differ by
//!   D-T1, which the comparator verifies.

use iq_terminal::peaks::fixture::U10aDetectPeaks;
use iq_terminal::xgb::fixture::U10aXgbPredict;
use iq_types::fixture::{read_lines, run_file, Case, FixtureCtx, Registry, Summary, UnitFixture};
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(name)
}

fn summary(name: &str) -> Summary {
    let path = data(name);
    let mut reg = Registry::new();
    iq_terminal::register_fixtures(&mut reg);
    let ctx = FixtureCtx::open(&path);
    let mut s = Summary::new(5);
    run_file(&reg, &path, &ctx, Some("U10a"), &mut s).unwrap();
    for ((u, f), t) in &s.by_key {
        eprintln!("{name}: {u} {f}: exact {} explained {:?} mismatch {}", t.exact, t.explained, t.mismatch);
        for (c, d) in &t.examples {
            eprintln!("  {c}\n{d}");
        }
    }
    s
}

#[test]
fn random_native_tie_order() {
    let s = summary("u10a_random_native.jsonl.gz");
    assert!(s.all_ok());
    let t = &s.by_key[&("U10a".to_string(), "detect_peaks_batch".to_string())];
    assert!(t.exact + t.explained.get("D-T1").copied().unwrap_or(0) >= 600);
}

#[test]
fn random_stable_bit_exact() {
    let s = summary("u10a_random_stable.jsonl.gz");
    assert!(s.all_ok());
    assert!(s.by_key.values().all(|t| t.explained.is_empty()));
    // Bit-exact, beyond the comparator's tolerances. Only on Linux: the reference values came from
    // glibc's libm, and other platforms' powf/ln can differ in the last bit.
    if !cfg!(target_os = "linux") {
        return;
    }
    let ctx = FixtureCtx::open(&data("u10a_random_stable.jsonl.gz"));
    let (mut n_peaks, mut n_xgb) = (0, 0);
    for c in read_lines::<Case>(&data("u10a_random_stable.jsonl.gz")).unwrap() {
        let c = c.unwrap();
        let fn_ = c.input["fn"].as_str().unwrap().to_owned();
        if fn_ == "detect_peaks_batch" {
            let i = iq_types::serde_json::from_value(c.input).unwrap();
            let want: <U10aDetectPeaks as UnitFixture>::Out = iq_types::serde_json::from_value(c.out).unwrap();
            assert_eq!(U10aDetectPeaks::run(&i, &ctx).unwrap(), want, "{}", c.case);
            n_peaks += 1;
        } else {
            let i = iq_types::serde_json::from_value(c.input).unwrap();
            let want: <U10aXgbPredict as UnitFixture>::Out = iq_types::serde_json::from_value(c.out).unwrap();
            assert_eq!(U10aXgbPredict::run(&i, &ctx).unwrap(), want, "{}", c.case);
            n_xgb += 1;
        }
    }
    assert!(n_peaks >= 1500 && n_xgb >= 100);
}

/// The D-T1 checker accepts the reverse tie order and rejects corrupted outputs.
#[test]
fn tie_checker_accepts_other_orders_and_rejects_corruption() {
    use iq_terminal::peaks::fixture::{compare_batch, run_ties};
    use iq_terminal::peaks::Ties;
    use iq_types::fixture::Verdict;
    let (mut alt_differs, mut rejected, mut corrupted) = (0, 0, 0);
    for name in ["u10a_random_stable.jsonl.gz", "u10a_random_native.jsonl.gz"] {
        for c in read_lines::<Case>(&data(name)).unwrap() {
            let c = c.unwrap();
            if c.input["fn"] != "detect_peaks_batch" {
                continue;
            }
            let i = iq_types::serde_json::from_value(c.input).unwrap();
            let got = run_ties(&i, Ties::Numpy).unwrap();
            let alt = run_ties(&i, Ties::PositionDesc).unwrap();
            if alt != got {
                alt_differs += 1;
                let v = compare_batch(&got, &alt);
                assert!(v == Verdict::Exact || v == Verdict::Explained("D-T1"), "{}: {v:?}", c.case);
            }
            // Corruptions that no tie order explains.
            let mut bad = alt.clone();
            if let Some(p) = bad.peaks.iter_mut().flatten().next() {
                p.1 += 1; // window count
                corrupted += 1;
                assert!(!compare_batch(&got, &bad).is_ok(), "{}: count corruption accepted", c.case);
                rejected += 1;
            }
            let mut bad = alt.clone();
            if let Some(r) = bad.records.iter_mut().flatten().find(|r| r.raw_peaks.is_some()) {
                let rp = r.raw_peaks.as_mut().unwrap();
                rp.peak_width[0] += 0.5;
                corrupted += 1;
                assert!(!compare_batch(&got, &bad).is_ok(), "{}: width corruption accepted", c.case);
                rejected += 1;
            }
            // Swapping two rows of different heights breaks the ranking.
            let mut bad = alt.clone();
            if let Some(k) = (1..bad.rows.len()).find(|&k| {
                bad.rows[k].hist == bad.rows[k - 1].hist && bad.rows[k].features[6].0 != bad.rows[k - 1].features[6].0
            }) {
                bad.rows.swap(k, k - 1);
                corrupted += 1;
                assert!(!compare_batch(&got, &bad).is_ok(), "{}: row swap accepted", c.case);
                rejected += 1;
            }
            // Dropping an accepted peak.
            let mut bad = alt.clone();
            if let Some(pl) = bad.peaks.iter_mut().find(|p| !p.is_empty()) {
                pl.pop();
                corrupted += 1;
                assert!(!compare_batch(&got, &bad).is_ok(), "{}: dropped peak accepted", c.case);
                rejected += 1;
            }
        }
    }
    eprintln!("alt order differs in {alt_differs} cases; {rejected}/{corrupted} corruptions rejected");
    assert!(alt_differs > 100);
}
