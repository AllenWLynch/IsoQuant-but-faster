//! Runs the U10b runner over `$IQRS_FIXTURES` (or `../fixtures`) when present, and checks the
//! runner's "second flush" shortcut: an empty recorded accumulator with a non-empty block
//! accumulation occurs only for the `.1` case (the `dump()` re-flush of the last block).
//! The authoritative gate is `iqrs-fixture run U10b ../fixtures`.

use iq_terminal::counter::{TerminalAccumulator, TerminalKind};
use iq_types::fixture::schema::{PredictionsOut, U10bFlushIn};
use iq_types::fixture::{read_cases, run_path, unit_files, FixtureCtx, Registry, Summary};
use iq_types::IsoRef;
use std::path::PathBuf;

fn root() -> Option<PathBuf> {
    let root = iq_types::fixture::fixtures_root()
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../../fixtures"));
    root.is_dir().then_some(root)
}

#[test]
fn u10b_fixtures() {
    let Some(root) = root() else {
        eprintln!("no fixtures: skipped");
        return;
    };
    let mut reg = Registry::new();
    iq_terminal::register_fixtures(&mut reg);
    let mut s = Summary::new(5);
    run_path(&reg, "U10b", &root, &mut s).unwrap();
    for ((u, f), t) in &s.by_key {
        eprintln!("{u} {f}: exact {} explained {:?} mismatch {}", t.exact, t.explained, t.mismatch);
        for (c, d) in &t.examples {
            eprintln!("  {c}\n{d}");
        }
    }
    assert!(s.all_ok());
}

#[test]
fn u10b_accumulation_by_case_key() {
    let Some(root) = root() else { return };
    let mut n = [0usize; 2];
    for sample in std::fs::read_dir(&root).unwrap().flatten() {
        for slice in std::fs::read_dir(sample.path()).into_iter().flatten().flatten() {
            let cap = slice.path();
            let Ok(files) = unit_files(&cap, "U10b") else { continue };
            let ctx = FixtureCtx::new(&cap);
            for f in files {
                for case in read_cases::<U10bFlushIn, PredictionsOut>(&f).unwrap() {
                    let case = case.unwrap();
                    let kind = if case.input.kind == "polya" { TerminalKind::PolyA } else { TerminalKind::Tss };
                    let block = ctx.block(&case.input.block).unwrap();
                    let mut acc = TerminalAccumulator::new();
                    for r in block.store.iter() {
                        acc.add_read(kind, &r);
                    }
                    let got: Vec<(String, Vec<i32>)> = {
                        let ids = ctx.ids();
                        acc.transcripts()
                            .into_iter()
                            .map(|(t, d)| (ids.iso_name(IsoRef::Ref(t)).to_owned(), d))
                            .collect()
                    };
                    let want: Vec<(String, Vec<i32>)> =
                        case.input.transcripts.iter().map(|(t, e)| (t.clone(), e.data.clone())).collect();
                    if case.case.ends_with(".0") {
                        assert_eq!(got, want, "{}", case.case);
                        n[0] += 1;
                    } else {
                        assert!(case.case.ends_with(".1"), "{}", case.case);
                        assert!(want.is_empty() && case.out.predictions.is_empty(), "{}", case.case);
                        n[1] += 1;
                    }
                }
            }
        }
    }
    eprintln!("first flushes {}, dump re-flushes {}", n[0], n[1]);
}
