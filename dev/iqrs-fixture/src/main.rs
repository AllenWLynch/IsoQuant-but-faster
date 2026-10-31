//! `iqrs-fixture`: runs unit fixture runners against golden fixtures (design §3.0, §10.3).
//!
//! ```text
//! iqrs-fixture run <unit|all> <file|dir> [--examples N]
//! iqrs-fixture check <file|dir>          schema conformance (byte-equal re-encoding)
//! iqrs-fixture list                      registered (unit, fn) runners
//! ```
//! Developer tool: fixtures are captured from IsoQuant 4.1 with `dev/harness` (not shipped).

use anyhow::{bail, Result};
use iq_types::fixture::{self, roundtrip, Registry, Summary};
use std::path::PathBuf;

/// Collects every crate's fixture runners. Unit owners: add `reg.add::<YourRunner>()` in
/// your crate's `register_fixtures`, not here.
fn registry() -> Registry {
    let mut reg = Registry::new();
    iq_annot::register_fixtures(&mut reg);
    iq_bam::register_fixtures(&mut reg);
    iq_regions::register_fixtures(&mut reg);
    iq_read::register_fixtures(&mut reg);
    iq_assign::register_fixtures(&mut reg);
    iq_collect::register_fixtures(&mut reg);
    iq_terminal::register_fixtures(&mut reg);
    iq_graph::register_fixtures(&mut reg);
    iq_model::register_fixtures(&mut reg);
    iq_quant::register_fixtures(&mut reg);
    reg
}

const USAGE: &str = "usage: iqrs-fixture <run <unit|all> <file|dir> [--examples N] | check <file|dir> | list>";

fn main() -> Result<()> {
    use lexopt::prelude::*;
    let mut p = lexopt::Parser::from_env();
    let mut pos: Vec<String> = Vec::new();
    let mut examples = 3usize;
    while let Some(arg) = p.next()? {
        match arg {
            Long("examples") => examples = p.value()?.parse()?,
            Short('h') | Long("help") => {
                println!("{USAGE}");
                return Ok(());
            }
            Value(v) => pos.push(v.string()?),
            _ => return Err(arg.unexpected().into()),
        }
    }
    match pos.first().map(String::as_str) {
        Some("run") => {
            let (Some(unit), Some(path)) = (pos.get(1), pos.get(2)) else { bail!("{USAGE}") };
            let reg = registry();
            let mut summary = Summary::new(examples);
            fixture::run_path(&reg, unit, &PathBuf::from(path), &mut summary).map_err(|e| anyhow::anyhow!("{e}"))?;
            print!("{summary}");
            if !summary.all_ok() {
                std::process::exit(1);
            }
            Ok(())
        }
        Some("check") => {
            let Some(path) = pos.get(1) else { bail!("{USAGE}") };
            let path = PathBuf::from(path);
            let files = if path.is_dir() { roundtrip::fixture_files(&path) } else { vec![path] };
            let (mut lines, mut bad, mut unreadable) = (0u64, 0u64, 0u64);
            let mut by_kind: std::collections::BTreeMap<String, (u64, u64)> = Default::default();
            for f in &files {
                let fc = match roundtrip::check_file(f, 5) {
                    Ok(fc) => fc,
                    Err(e) => {
                        // e.g. a truncated gzip of a capture still being written
                        println!("   error              {}: {e}", f.display());
                        unreadable += 1;
                        continue;
                    }
                };
                lines += fc.lines;
                bad += fc.lines - fc.ok;
                let e = by_kind.entry(file_kind(f)).or_default();
                e.0 += fc.lines;
                e.1 += fc.ok;
                println!("{:>8} lines {:>6} bad  {}", fc.lines, fc.lines - fc.ok, f.display());
                for (no, m) in fc.failures {
                    println!("    line {no}: {m}");
                }
            }
            println!("{:<12} {:>10} {:>10}", "kind", "lines", "roundtrip");
            for (k, (n, ok)) in &by_kind {
                println!("{k:<12} {n:>10} {ok:>10}");
            }
            println!("{} files ({unreadable} unreadable), {lines} lines, {bad} non-conforming", files.len());
            if bad > 0 || unreadable > 0 {
                std::process::exit(1);
            }
            Ok(())
        }
        Some("list") => {
            for (u, f) in registry().keys() {
                println!("{u}\t{f}");
            }
            Ok(())
        }
        _ => bail!("{USAGE}"),
    }
}

/// Tally key for `check`: the unit of a unit file, or `block` / `index` / `params`.
fn file_kind(f: &std::path::Path) -> String {
    let name = f.file_name().and_then(|n| n.to_str()).unwrap_or("");
    if f.components().any(|c| c.as_os_str() == "block") {
        "block".into()
    } else if name == "params.json" {
        "params".into()
    } else if name.ends_with(".json") {
        "index".into()
    } else {
        name.split('.').next().unwrap_or(name).to_owned()
    }
}
