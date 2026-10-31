//! Runs the U10a fixture runners over `$IQRS_FIXTURES` (or `../fixtures`) when present.
//! The authoritative gate is `iqrs-fixture run U10a ../fixtures`; this test keeps the crate
//! self-checking without the binary.

use iq_types::fixture::{run_path, Registry, Summary};
use std::path::PathBuf;

#[test]
fn u10a_fixtures() {
    let root = iq_types::fixture::fixtures_root()
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../../fixtures"));
    if !root.is_dir() {
        eprintln!("no fixtures at {}: skipped", root.display());
        return;
    }
    let mut reg = Registry::new();
    iq_terminal::register_fixtures(&mut reg);
    let mut s = Summary::new(5);
    run_path(&reg, "U10a", &root, &mut s).unwrap();
    for ((u, f), t) in &s.by_key {
        eprintln!("{u} {f}: exact {} explained {:?} mismatch {}", t.exact, t.explained, t.mismatch);
        for (c, d) in &t.examples {
            eprintln!("  {c}\n{d}");
        }
    }
    assert!(s.all_ok());
}
