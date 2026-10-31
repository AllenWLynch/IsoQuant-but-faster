# 0.1.0 release preparation

Baseline: commit a67223f. That commit also holds the internal coordination notes deleted here:
STATUS.md, UNITS.md, CONTRACT_CHANGES.md, docs/progress_*.md, PILOT_U4.md and vtg_proto.md.

## What changed

- **One behaviour only.** The `quirks41` cargo feature is removed, along with `iq_types::dev`
  (`DevId`, `QuirkSet`, `quirk()`, `with_quirks()`), the `IQRS_QUIRKS` env var, the fixture
  `RunMode` and `UnitFixture::DEVIATIONS`. At each of the 17 former `quirk()` sites only the 4.1
  branch is left, marked "matches IsoQuant 4.1". `Verdict::Explained` now carries a
  `&'static str`. Only `"D-T1"` (peak tie order) is still used. See `deviations.md`.
- `--gene-types` now defaults to `isoquant`, and `GeneTypes::default()` is `Isoquant`.
- iqrs-compare: `--mode parity|default` is replaced by `--exclude-nongene-loci`.
- CLI: full `--help` and `--version`. Usage errors exit with code 2, run errors with code 1.
  Clear messages for a missing `.bai` (both locations listed), `@HD SO:queryname|unsorted`,
  out-of-order records (checked in pass A), a missing `.fai`, a compressed FASTA and missing
  input files. The output directory is created up front.
- Layout: `harness/` moved to `dev/harness/` and `bins/iqrs-fixture` to `dev/iqrs-fixture`.
  Hard-coded cluster paths were replaced by `$IQRS_PROJECT` or script-relative paths. Capture
  logs were removed.
- Small unit-test data (`crates/*/tests/data/**/*.jsonl.gz`, about 3.8 MB) is now tracked. It
  was gitignored before, so a fresh clone could not run those tests.
- `rustfmt.toml` sets `max_width = 120` and `use_small_heuristics = "Max"`. The tree had never
  been fmt-clean, so this release includes a one-time whitespace-only reformat.
- `[profile.release]`: fat LTO, `codegen-units = 1`, `strip = "debuginfo"`. The new
  `[profile.profiling]` keeps line tables. A clean release build takes about 1.5 min.
- musl builds use mimalloc as the global allocator. This is `cfg(target_env = "musl")` only,
  approved by the maintainer, and credited in NOTICE.
- Added README, LICENSE (GPL-2.0 text from IsoQuant's GPL2.txt), NOTICE, CHANGELOG,
  CONTRIBUTING, `.github/workflows/{ci,release}.yml` and `docs/dev/deviations.md`.

## Validation

Byte-identity ignores only the GTF `# Command line:` header.

- **chr8, 5 samples** (2 ONT, Iso-Seq, 2 Revio; `validation/run_chr8_dr2.sh` invocation):
  - The glibc and musl builds are both byte-identical to the current parity outputs
    (`transcript_models.gtf`, `discovered_transcript_counts.tsv`).
  - Checked on the pre-reformat build at 8 and 16 threads, and on the final binaries at
    16 threads: 30/30 identical.
- **Fixture gate** (`iqrs-fixture run all`): every unit exact, apart from the usual D-T1 explained
  cases, with exit code 0. The per-function counts are identical to the previous `quirks41`
  parity run.
- **Whole genome, SRR32588726**, 16 threads, `--read-group tag:HP`, same node for all runs.
  Every run below is byte-identical to `bench_full/scale3_parity_SRR32588726` (GTF, counts and
  grouped HP counts). Peak RSS is `/usr/bin/time` max RSS ÷ 1e6, as in the 19.7 GB baseline.

| build | wall | peak RSS | user + sys CPU |
|---|---|---|---|
| glibc (final) | 1:39, 1:35 | 20.0, 20.5 GB | 1,016 s, 1,004 s |
| musl + mimalloc (final) | 1:59, 1:56 | 19.3, 19.5 GB | 1,124 s, 1,093 s |
| musl, system malloc | 11:44 | 14.4 GB | 8,328 s (6,071 s sys) |
| glibc + mimalloc (test only) | 1:47 | 19.8 GB | 1,029 s |
| previous baseline (glibc) | 1:48 | 19.7 GB | 1,016 s |

- glibc stays within the target: about 10% faster than the baseline and 2-4% more peak memory,
  inside node-to-node noise.
- musl + mimalloc is about 22% slower than glibc on the same node, almost all of it in the
  Collect phase (Σtask 962-1242 s vs 822-896 s). glibc + mimalloc does not show that slowdown,
  so the cause is musl's libc routines (likely memcpy), not the allocator. mimalloc does not help
  glibc, so glibc keeps the system allocator.
- The release workflow therefore also builds an `x86_64-unknown-linux-gnu` binary against
  glibc 2.17 with cargo-zigbuild, and keeps musl as the fully static fallback.

`../validation/iqrs_release` is the final glibc binary (built on the validation host), and
`../validation/iqrs_release_musl` is the final static build. Locally, the musl build compiles
mimalloc's C code with the host gcc (`CC_x86_64_unknown_linux_musl=gcc`). CI uses `musl-gcc`
from `musl-tools`.

## Tests

`cargo test --workspace --release` passes in the tree with fixtures present. It also passes in a
fresh-clone copy (tracked plus new files only) with no fixtures: fixture tests skip.
`cargo clippy --workspace --all-targets -- -D warnings` is clean, also for the musl target.
`cargo fmt --all --check` is clean.

## Not tested

The GitHub workflows have not been run. Untested parts: the aarch64 musl build on the
`ubuntu-24.04-arm` runner, cargo-zigbuild for the glibc 2.17 build, and the macOS build.
