//! Loads every real fixture file and checks it re-encodes to byte-equal JSON through the
//! typed DTOs (U0 pass criterion). Fixture root: `$IQRS_FIXTURES`, else
//! `<workspace>/../fixtures`. Skips (passes) when no fixtures exist. Lines per file are capped
//! by `$IQRS_FIXTURE_MAX_LINES` (default 2000) to bound debug-build runtime; use
//! `iqrs-fixture check <dir>` (release) for the exhaustive check.

use iq_types::fixture::roundtrip::{check_file_limited, fixture_files};
use std::path::PathBuf;

fn root() -> Option<PathBuf> {
    let r = iq_types::fixture::fixtures_root()
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../../fixtures"));
    r.is_dir().then_some(r)
}

#[test]
fn every_fixture_file_roundtrips_byte_equal() {
    let Some(root) = root() else {
        eprintln!("no fixtures found; skipping");
        return;
    };
    let files = fixture_files(&root);
    eprintln!("checking {} fixture files under {}", files.len(), root.display());
    let todo: Vec<&PathBuf> = files.iter().filter(|f| !in_progress(f)).collect();
    let skipped = files.len() - todo.len();
    // files are independent: check them on all cores
    let next = std::sync::atomic::AtomicUsize::new(0);
    let results = std::sync::Mutex::new((0u64, Vec::<String>::new()));
    let n_threads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4).min(16);
    std::thread::scope(|sc| {
        for _ in 0..n_threads {
            sc.spawn(|| loop {
                let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                let Some(f) = todo.get(i) else { break };
                let fc = check_file_limited(f, 3, max_lines()).unwrap_or_else(|e| panic!("{}: {e}", f.display()));
                let mut r = results.lock().unwrap();
                r.0 += fc.lines;
                for (no, msg) in fc.failures {
                    r.1.push(format!("{}:{no}: {msg}", f.display()));
                }
            });
        }
    });
    let (lines, bad) = results.into_inner().unwrap();
    eprintln!("{lines} lines checked; {skipped} files of in-progress captures skipped");
    assert!(bad.is_empty(), "{} non-conforming lines:\n{}", bad.len(), bad.join("\n"));
}

/// A capture still being written: per-PID unit files (renamed on completion, schema §0.1)
/// or a capture dir without `fixture_manifest.json` (written last).
fn in_progress(f: &std::path::Path) -> bool {
    let name = f.file_name().and_then(|n| n.to_str()).unwrap_or("");
    let pid_file = name.split('.').nth(1).is_some_and(|p| !p.is_empty() && p.bytes().all(|b| b.is_ascii_digit()));
    let cap = f.ancestors().find(|a| a.join("units").is_dir() || a.join("index").is_dir());
    pid_file || cap.is_some_and(|c| !c.join("fixture_manifest.json").exists())
}

/// Side tables of every complete capture load into core types: `params.json` → `Params`,
/// every `index/*.json` → `IsoformIndex`, every block → `BlockStore` whose records convert
/// back to the identical DTOs.
#[test]
fn side_tables_load_into_core_types() {
    use iq_types::fixture::{find_capture_dirs, load_block, FixtureCtx, ReadRecordFixture};
    let Some(root) = root() else { return };
    let mut caps = Vec::new();
    find_capture_dirs(&root, &mut caps).unwrap();
    let mut n_idx = 0;
    let mut n_reads = 0;
    for cap in caps.iter().filter(|c| c.join("fixture_manifest.json").exists()) {
        let ctx = FixtureCtx::new(cap);
        if cap.join("params.json").exists() {
            ctx.params().unwrap_or_else(|e| panic!("{}: params: {e}", cap.display()));
        }
        if let Ok(rd) = std::fs::read_dir(cap.join("index")) {
            for e in rd.flatten().take(100) {
                // bounded: up to 100 indices per capture
                let p = e.path();
                let sha = p.file_stem().unwrap().to_str().unwrap().to_owned();
                let idx = ctx.index(&sha).unwrap_or_else(|e| panic!("{}: {e}", p.display()));
                assert_eq!(idx.id_rank.len(), idx.n_iso());
                n_idx += 1;
            }
        }
        // bounded: first 5 block files per capture
        for f in iq_types::fixture::roundtrip::fixture_files(&cap.join("block")).into_iter().take(5) {
            let mut ids = ctx.ids();
            let b = load_block(&f, &mut ids).unwrap_or_else(|e| panic!("{}: {e}", f.display()));
            assert_eq!(b.records.len() as u32, b.header.n_reads, "{}", f.display());
            for (i, (fx, ex)) in b.fixtures.iter().zip(&b.extras).enumerate() {
                let back = ReadRecordFixture::from_core(&b.store.record(i), &ids, ex);
                assert_eq!(&back, fx, "{} read {i}", f.display());
            }
            n_reads += b.records.len();
        }
    }
    eprintln!("{n_idx} indices, {n_reads} block reads loaded");
}

fn max_lines() -> u64 {
    std::env::var("IQRS_FIXTURE_MAX_LINES").ok().and_then(|v| v.parse().ok()).unwrap_or(2_000)
}
