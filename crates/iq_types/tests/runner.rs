//! End-to-end test of the fixture runner machinery on a synthetic capture dir: registry
//! dispatch by `(unit, fn)`, verdicts, Python exceptions,
//! decode errors, unregistered fns, and panics (`todo!()`).

use iq_types::fixture::schema::{Case, U9ResolveIn, U9ResolveOut};
use iq_types::fixture::{self, FixResult, FixtureCtx, Registry, Summary, UnitFixture};
use std::io::Write;

/// Toy runner: echoes the input assignments; panics on an empty list.
struct ToyU9;
impl UnitFixture for ToyU9 {
    const UNIT: &'static str = "U9";
    const FN: &'static str = "resolve";
    type In = U9ResolveIn;
    type Out = U9ResolveOut;
    fn run(input: &U9ResolveIn, _ctx: &FixtureCtx) -> FixResult<U9ResolveOut> {
        let a = input.assignments.clone();
        assert!(!a.is_empty(), "todo: empty");
        Ok(U9ResolveOut { assignments: a })
    }
}

fn basic(mm: bool) -> String {
    format!(
        r#"{{"assignment_id":1,"read_id":"r0","chr":"8","start":100,"end":400,"region":[99,5000],"multimapper":{mm},"polya_found":false,"atype":10,"gene_atype":10,"penalty":0.0,"genes":["g"],"isoforms":["t"]}}"#
    )
}

fn write_capture(dir: &std::path::Path) {
    std::fs::create_dir_all(dir.join("units")).unwrap();
    let f = std::fs::File::create(dir.join("units/U9.jsonl.gz")).unwrap();
    let mut gz = flate2::write::GzEncoder::new(f, flate2::Compression::fast());
    let line = |case: &str, inp: &str, out: &str| {
        format!(
            r#"{{"unit":"U9","v":1,"case":"S:sl:8:-:{case}","in":{{"fn":"resolve","strategy":2,"assignments":[{inp}]}},"out":{out}}}"#
        )
    };
    // 1: echo matches -> Exact
    writeln!(gz, "{}", line("a", &basic(false), &format!(r#"{{"assignments":[{}]}}"#, basic(false)))).unwrap();
    // 2: differs -> Mismatch
    writeln!(gz, "{}", line("b", &basic(false), &format!(r#"{{"assignments":[{}]}}"#, basic(true)))).unwrap();
    // 3: python exception
    writeln!(gz, "{}", line("c", &basic(false), r#"{"exception":"KeyError: 'x'"}"#)).unwrap();
    // 4: undecodable out
    writeln!(gz, "{}", line("d", &basic(false), r#"{"assignments":7}"#)).unwrap();
    // 5: panic in run
    writeln!(gz, r#"{{"unit":"U9","v":1,"case":"S:sl:8:-:e","in":{{"fn":"resolve","strategy":2,"assignments":[]}},"out":{{"assignments":[]}}}}"#).unwrap();
    // 6: unregistered fn
    writeln!(gz, r#"{{"unit":"U9","v":1,"case":"S:sl:8:-:f","in":{{"fn":"other"}},"out":{{}}}}"#).unwrap();
    gz.finish().unwrap();
}

#[test]
fn registry_and_verdicts() {
    let dir = std::env::temp_dir().join(format!("iq_types_runner_{}", std::process::id()));
    write_capture(&dir);
    let mut reg = Registry::new();
    reg.add::<ToyU9>();
    assert!(reg.get("U9", "resolve").is_some());

    let mut s = Summary::new(5);
    fixture::run_path(&reg, "U9", &dir, &mut s).unwrap();
    let t = &s.by_key[&("U9".to_string(), "resolve".to_string())];
    assert_eq!((t.py_exception, t.decode_error), (1, 1));
    assert_eq!(s.by_key[&("U9".to_string(), "other".to_string())].unregistered, 1);
    assert_eq!(t.exact, 1);
    assert_eq!((t.mismatch, t.errors), (2, 1)); // the diff + the panic

    let mut p = Summary::new(5);
    fixture::run_path(&reg, "all", &dir.join("units/U9.jsonl.gz"), &mut p).unwrap();
    let tp = &p.by_key[&("U9".to_string(), "resolve".to_string())];
    assert_eq!((tp.exact, tp.mismatch, tp.errors), (1, 2, 1));
    let report = s.to_string();
    assert!(report.contains("U9"), "{report}");
    assert!(report.contains(".assignments[0].multimapper"), "{report}");
    assert!(report.contains("panic: todo: empty"), "{report}");

    // typed reading of the same file
    let n =
        fixture::read_cases::<serde_json::Value, serde_json::Value>(&dir.join("units/U9.jsonl.gz")).unwrap().count();
    assert_eq!(n, 6);
    let _: Option<Case> = None;
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn case_keys() {
    let k = fixture::CaseKey::parse("SRR1:ordinary1:8:12:a37.1").unwrap();
    assert_eq!((k.sample, k.slice, k.chr, k.block_no, k.local), ("SRR1", "ordinary1", "8", Some(12), "a37.1"));
    let k = fixture::CaseKey::parse("S:sl:8:-:split_coverage_regions:n0").unwrap();
    assert_eq!((k.block_no, k.local), (None, "split_coverage_regions:n0"));
    assert!(fixture::CaseKey::parse("annotation:8").is_none());
}
