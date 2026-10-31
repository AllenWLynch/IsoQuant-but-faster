# D-R2 retirement + D-M6 "explained" reporting (2026-10-29)

Brief: revert D-R2 to 4.1 double processing in both builds (user-approved after `../validation/deviation_effects.md`),
mark D-M6's tx/gene effect as explained in iqrs-compare, re-validate (fixtures, chr8 E2E, full-genome timing).
Build dirs: `../.cargo-target/d-r2` (default), `../.cargo-target/d-r2-q41` (quirks41).

## Code changes
- `crates/iq_collect/src/lib.rs`: `ClusterOwner::owns` = 4.1 fetch-window test, unconditional (no `quirk(R2)`).
  New `Owner::shared` (default false; ClusterOwner: the record also overlaps an adjacent block region) and
  `BlockCollector::finish_shared()` → name hashes of selected shared reads.
- `crates/iq_collect/src/multimap.rs`: `resolve_chr_duplicates(stores, candidates: Option<&[u64]>, strategy)`.
  `Some` = pass-A dup candidates ∪ shared names (only those reads are grouped); `None` = scan all. Quirk gate removed.
- `crates/iq_collect/src/fixture.rs`, `tests/u8_fixtures.rs`: R2 removed from deviation lists; test uses `None`.
- `crates/iq_sched/src/real.rs`: `RealStore::shared`; barrier merges it into the U9 candidates.
- `crates/iq_regions/src/pass_a.rs`: `owned_reads` estimate adds each boundary bin's coverage to both neighbours
  (upper bound on straddlers), so scheduler cost/memory estimates (`pipeline.rs` est_* use owned_reads) include the copies.
  `owner.rs`/`lib.rs` docs: owner_of is now only the estimate partition.
- `crates/iq_types/src/dev.rs`: D-R2 verdict Drop → Keep (variant kept so IQRS_QUIRKS specs parse); `read.rs` Block doc.
- `bins/iqrs-compare`: `--explained METRIC=DEV` (repeatable / comma list; metric-name prefix). A failing gated metric
  with a rule is EXPLAINED (DEV) in summary.md, `"verdict": "explained"`, `"explained_by"` in compare.json, plus
  `gated_explained`, `gated_fail`, `all_pass_or_explained`, `explained_rules`; TSV `EXPLAINED:DEV`. `--fail-on-gate` ignores explained.
- Docs: design §5.1/§5.2/§6.2/§9 (D-R2 KEEP, D-C1, D-M5 rows), CONTRACT_CHANGES.md section.

## Validation
- Unit tests (iq_collect, iq_regions, iq_sched, iq_types, iqrs-compare) pass; clippy clean both builds.
- Note: the interactive job has 8 GB; `iqrs-fixture run all` was OOM-killed at 4.6 GB there, so fixtures run via
  SLURM: `../validation/fixtures_dr2/run.sbatch` → `q41_parity.txt`, `q41_default.txt`, `default_build.txt`, `u9_chr_q41.txt`.
- chr8 E2E (`../validation/run_chr8_dr2.sh`, log `../validation/run_chr8_dr2.log`; previous runs kept in
  `../validation/runs_pre_dr2/`, `compare_pre_dr2/`):
  - parity: 15/15 on all 5 (also vs 4.1 no-SIMD for 720/726: 15/15). Parity GTF+counts byte-identical to
    `attrib/runs/t1np_parity` (old full-scan U9 under quirk(R2)) on all 5 → the candidate-list U9 path is exact.
  - default (`--explained tx_per_gene=D-M6`): ONT 249/254 15/15; 720, 716, 726 14/15 + tx/gene EXPLAINED (D-M6), 0 fail.
    ONT 254 fate assigned −0.019 pp (was −1.24); Revio tx/locus +0.30 / +0.32 % (were +1.53 / +5.06 %);
    tx/gene +1.63 % (720), +4.10 % (716), +1.77 % (726); ONT 249 Spearman 0.9973.
- Fixtures (SLURM 57987498, `../validation/fixtures_dr2/`):
  - quirks41 parity: every unit 100 % exact or D-T1-explained (U8 1801/1801, U9 164/164, U3 all exact).
  - quirks41 `--mode default`: U3 and U9 all exact. U8 1596 exact + 194 explained (D-P1 129, D-A4 65) + 11 mismatch, with **no D-R2**
    (the pre-change binary gave 1514 + 245 (D-R2 82) + 42). The 11 are blocks with a D-P1 read and a separate D-A4 read, which no
    single quirk restores. Per-read attribution (`u8_per_read_attrib.txt`): 757/739,047 reads differ, P1 608 / A4 140 / P3 6 / two-dev 3,
    **0 unexplained**. Other units as before (D-R1, D-G1, D-M6, D-A4, D-O1, D-P1, D-T1).
  - default build `--mode default`: U8 has the same 1596 exact; U3/U9 exact.
  - U8→U9 whole-capture (`u9_chr_q41.txt`, includes fullchr8): 1801 blocks equal to 4.1 block files; total/polyA equal to the 4.1 log in 39 captures.
- Full genome (SLURM 57987499, `../bench_full/scale2_*`), 16 threads, tag:HP:
  | sample | parity before → after | default before → after |
  |---|---|---|
  | SRR32588726 | 1:52 / 22.2 GB → 2:26 / 20.3 GB | 1:40 / 16.5 GB → 2:25 / 20.0 GB |
  | SRR32588716 | 2:50 / 17.1 GB → 2:48 / 14.8 GB | 2:30 / 12.6 GB → 2:46 / 15.1 GB |
  - D-R2 itself costs about 10 % of default time on 716 (barrier Σ 1.3 → 15.7 s, since `apply_resolution` rebuilds every store that has a
    suspended copy; collect +6 %). Default memory now matches parity, as expected since both do the same work.
  - **The 726 slowdown (+35-45 s) is not D-R2.** A same-node A/B (SLURM 57988884, `../bench_full/ab_*.log`, parity builds) gave
    old 20:25 binary 1:43, `attrib/iqrs_parity_t1np` (D-T1 numpy-argsort fix, before D-R2) 2:39, new 2:36. Construct span went 18 s → 64 s with the
    D-T1 change (one long construct tail). `np_argsort` alone is fast (1M elements in 108 ms, ties fine), so the cost is likely downstream
    (different TSS peaks → more alt-end work in one giant block). Open item for the D-T1 owner. The candidate-list U9 path cut the parity barrier Σ from 41-44 s to 20 s.
