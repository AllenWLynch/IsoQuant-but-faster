# U10a progress (terminal peaks + XGBoost)

## Status (2026-10-29)
- Implemented `iq_terminal::{peaks, xgb}`: histogram features (numpy pairwise sum, var/skew/entropy/mode),
  scipy `_local_maxima_1d`, `_select_by_peak_distance` (D-T1 ties by position), `_peak_prominences`/`_peak_widths`,
  `_rank_peaks`, `_window_count`, native xgboost f32 tree eval (hard-coded base-margin bits), `Sigmoid` port.
- Fixture runners `U10a/detect_peaks_batch`, `U10a/xgb_predict`; comparator with tolerances + D-T1 consistency check.
- Fixtures (via crate test `tests/fixtures_u10a.rs`, since iqrs-fixture was broken by another unit's file):
  detect_peaks_batch exact 43 / T1 38 / mismatch 0; xgb_predict exact 143.
- Random suite `tests/random_u10a.rs` with `tests/data/gen_u10a.py` (live Python):
  - stable file (numpy argsort monkeypatched to D-T1 order): 1500 batch + 100x25 xgb rows, all **bit-exact**;
  - native file: 600 batch cases, 212 exact / 388 T1;
  - tie checker accepts reverse tie order, rejects 5204/5204 corruptions.

## Findings
- np.add.reduce on contiguous f64 = `0.0 + pairwise(all)` (no 8192 chunking) — verified on 3000 random arrays.
- Int sums (mean) exact as integers.
- Base margin bits 0x3ea048e8 / 0xbf283a5a confirmed by bit-exact margins (log-based recompute gives ...e5/e6).
- D-T1 is frequent on real data: 38/81 fixture batch cases differ by tie order (AVX-512 capture node).

## Done
- iqrs-fixture (quirks41 parity, quirks41 --mode default, plain default build): detect_peaks_batch 43 exact + 38 D-T1,
  xgb_predict 143 exact, 0 mismatch in all three. clippy clean with/without quirks41; tests pass in both builds.
- Note: `cargo fmt -p iq_terminal` reformats counter.rs (U10b); restored the stub. The repo has no rustfmt.toml.
