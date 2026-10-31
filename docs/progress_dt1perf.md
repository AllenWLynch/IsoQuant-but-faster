# D-T1 slowdown on SRR32588726: root cause + fix (2026-10-29)

Brief: full genome on SRR32588726 went 1:43 → 2:39 with the D-T1 numpy-argsort change (A/B 57988884,
`../bench_full/ab_*.log`). Find the cause, restore runtime without changing outputs.
Build dirs: `../.cargo-target/dt1perf` (quirks41), `../.cargo-target/dt1perf-def` (default).

## Diagnosis
- Temporary per-stage timing in `construct_block` (env-gated) + a temporary env switch back to the pre-D-T1 tie order
  (`Ties::PositionAsc`). Both are removed again. Full-genome run, SLURM 57989755 (`../bench_full/dt1perf/diag_*.log`):
  numpy ties 2:36, asc ties 1:37. The construct tail was a single block, chr1 ref window from 24,868,015 (fbxw7 locus,
  1,057,793 reads): **filter = 62.3 s** with numpy ties (7 models) vs 0.02 s with asc ties (6 models).
- Reproduced on a 1.08 M-read chr1:24-27 Mb BAM (`../bench_full/dt1perf/vtg2/blk.bam`): 64 s in `filter`, 4 threads.
- Cause: with 4.1's real tie order, `_add_known_alternative_end_models` creates one more alternative-end NIC on the deep
  known fbxw7 model, and that NIC claims ~10^5 reads. Each claimed read went through `ModelStore::move_read_to_model`,
  which does `iter().position()` + `Vec::remove` on the source list (a port of Python `list.remove`). That is O(reads moved ×
  list length), about 10^11 element ops here. So the extra work comes from the new (correct, 4.1) output, and the quadratic
  read move is what makes it expensive. `np_argsort` itself is not the cost.

## Fix
- `crates/iq_model/src/ends.rs`: `move_reads_batched` replaces the per-read move loop in `add_known_alt_end_models`.
  The loop moves `min(len(alt.reads), len0 - 1)` reads (each move shrinks the source by one, so the `len <= 1` break lands
  at that point). The first occurrence of each moved read is removed in one `retain` pass, which keeps the remaining order.
  Alt reads are distinct (`claimed`). `internal_counter[src]` is decremented (and inserted if missing) before the target
  entries, as on the first sequential move. Then `save_assigned_read` runs per read in order. If a read is missing, it
  panics like `list.remove` does.
- Unit test `batched_moves_match_sequential_list_remove` compares against the sequential `move_read_to_model` loop:
  duplicates in the source, break at len 1, counter insertion order, and IndexMap order.
- `iq_types` is frozen and was not touched. `ModelStore::move_read_to_model` is still quadratic for any future bulk caller.
- Block result: 64 s → 0.14 s in `filter`. Outputs are byte-identical (only the command-line header line differs).

## Validation (SLURM 57991269, `../bench_full/dt1perf/validate.*.out`)
- chr8, 5 samples × {parity, default}: GTF + discovered counts byte-identical to `../validation/runs/` (header line excluded).
- Fixtures (quirks41 build, parity and `--mode default`): U10a, U10b, U13, U14, U15a, U15b, U16, U17, U19 give the same counts
  as `../validation/fixtures_dr2/q41_*.txt` (U10a 215+66 D-T1, U10b 361+7 D-T1, U15b all exact, U19 default 163+22 D-A4).
- iq_model tests pass and clippy is clean, with and without quirks41.
- Full-genome outputs from scale3 are byte-identical to scale2 for both samples and both builds.

## Full genome (SLURM 57991270, `../bench_full/scale3_*`), 16 threads, tag:HP; RSS = max RSS in GiB
| sample | parity scale2 → scale3 | default scale2 → scale3 |
|---|---|---|
| SRR32588726 | 2:26 / 19.4 → **1:48 / 18.8** | 2:25 / 19.1 → **1:42 / 19.3** |
| SRR32588716 | 2:48 / 14.1 → 2:51 / 14.1 | 2:46 / 14.4 → 2:47 / 14.3 |
- 726 construct span 63 s → 16 s, Σtask 361 → 260 s (parity). 716 was not affected. Its wall time is collect-bound:
  construct starts at about 155 s.
- `../validation/iqrs_parity` / `iqrs_default` are updated to these builds. The previous binaries are kept as
  `../validation/attrib/iqrs_{parity,default}_dr2`.

## Open observations (not changed)
- In the per-stage traces, big blocks show 9-18 s in `reassign1` / `final` (e.g. 4 M-read block final ≈ 9-18 s). Some of
  that is probably rayon work-stealing inside nested `par_iter` (a waiting task runs other blocks' work, so its wall time
  grows). This was the same with both tie orders, so it is the next construct-tail candidate.
