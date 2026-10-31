# U10b progress (per-block polyA/TSS accumulation -> predictions)

## Status (2026-10-29): done
- `crates/iq_terminal/src/counter.rs`: `TerminalKind`, `TerminalAccumulator` (add_read / add / flush_into /
  transcripts / clear), `read_order_moments`, `block_predictions` (+ `block_predictions_with` reusing buffers).
- Runner `U10b/flush` registered from `lib.rs` (`counter::register_fixtures`); crate-level
  `#![allow(unused_variables)]` removed.
- Fixtures: 140 cases, exact 138 + D-T1 2, mismatch 0 (quirks41 parity, quirks41 --mode default, plain default build).
- Tests: `tests/fixtures_u10b.rs` (runner + case-key check of the accumulation incl. `.1` re-flush), unit tests in
  counter.rs (moments vs U10a bitwise, accumulator order/reuse, flush vs detect_peaks, D-T1 checker rejects shifts).

## Findings
- 4.1 `_predict_rows` does NOT call `terminal_peaks.detect_peaks`; it has its own pandas pipeline, identical except
  that `np.var`/`stats.skew` run over `data` in **read (add) order**, while `detect_peaks` expands the histogram in
  dict order. Pairwise-sum bits differ by order (random test: 119/500 var values differ in the last bits; no
  prediction differed). So the accumulator keeps a flat `(slot, pos)` buffer in read order (8 B per accepted read,
  no per-read allocation, reused across blocks) rather than a histogram; histogram/peaks/entropy/mode are order-free
  and reuse U10a (`histogram_record` + overridden var/skew, `rank_peaks`, `feature_row`, `XgbModel::predict`).
- Last block of each chromosome has a second case `.1`: `dump()` calls `flush()` again on the empty accumulator ->
  always `[]` (25 cases). The runner detects it as "recorded accumulator empty but block accumulation non-empty";
  the crate test verifies that only happens for `.1` keys.
- Rank is not an XGBoost feature, so acceptance depends only on the kept peak set; D-T1 check enumerates
  distance-filter tie orders (exhaustive up to 5040, else seeded random) and requires 4.1's per-transcript rows be
  an accepted set of some order, height-descending. U10a's `tie_consistent` needs 4.1 records/rows (not in U10b out).
- Zero-peak fallback unreachable; ported with 4.1's frame order (zero-peak rows first).
- `IsoRef::Model` matches are skipped (`Predictions` is TxIdx-keyed; never occurs on reference assignments).
