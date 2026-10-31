# iqrs project status (coordinator handoff file)

Last updated: 2026-10-29 (coordinator session "isoquant-rust").

## Where we are
- Specs: `../specs/00_scope.md` (scope, decisions log), `01-05` (component specs), `06` (go verdict), `07_rust_design.md` (design; §3 = units, §12 = milestones, §13.1 = resolved questions).
- **M0 done**: `iq_types` v1 + fixture schema v1 **frozen**; `iqrs-fixture check ../fixtures` → 694,483 lines, 0 non-conforming. Harness: `harness/iqcap/` (README, SCHEMA.md). Fixtures: `../fixtures/` (3 samples × 5 slices + fullchr8 for ONT/Iso-Seq).
- Reference IsoQuant 4.1 outputs (chr8): `../reference_outputs/<sample>/iq/OUT/` (5 samples + Iso-Seq no-SIMD noise-floor run: 1/1375 transcripts differ).
- Rules for unit agents: `UNITS.md`. Pending contract changes: `CONTRACT_CHANGES.md`.

## Running (wave 1, max 4 concurrent agents per user)
| Unit | Files owned | Progress notes |
|---|---|---|

## Done
- D-T1 perf regression fixed (iq_model/src/ends.rs `move_reads_batched`; the quadratic per-read list.remove port hit ~1e5 reads on one fbxw7 alt-end model). Byte-identical outputs (chr8 x5 x2 builds, fixtures, full genome). scale3 (57991270): 726 parity 1:48/19.7 GB, default 1:42/20.2 GB; 716 2:51 / 2:47 (unaffected). Graphic: 65x faster than 4.1. Next candidate: 9-18 s in U17/U19 on the largest blocks. docs/progress_dt1perf.md.
- 4.1 full-genome SRR32588726 (job 57982465): 1:55:52, 81.4 GB SLURM MaxRSS (40.9 GB largest process), 34,052 CPU-s. vs iqrs scale2: parity 17/17 gates, chain 1.0000/1.0000, Spearman 1.0000, 73,023 vs 73,032 models; default 17/17, 0.9989/0.9990/0.9997. `../validation/compare_full/41_vs_*`. Graphic updated (4.1 only, no 3.10).
- D-R2 retired (4.1 double processing in both builds; U9 candidates = pass-A dups ∪ collector-flagged shared reads) + iqrs-compare `--explained tx_per_gene=D-M6`. chr8: parity 15/15 ×5 (byte-identical to previous parity); default 15/15 on ONT, 14/15 + tx/gene EXPLAINED (D-M6) on 720/716/726. Fixtures: U3/U9 exact both builds, U8 default has no D-R2 left (11 two-deviation P1+A4 blocks, 0 unexplained reads). Full genome (57987499): 726 2:26/20.3 GB parity, 2:25/20.0 GB default; 716 2:48/14.8 GB, 2:46/15.1 GB. The 726 slowdown comes from the D-T1 change (A/B 57988884), not D-R2. See docs/progress_dr2.md.
- Full-scale (job 57980044): SRR32588726 parity 1:52/22.2 GB, default 1:40/16.5 GB; SRR32588716 parity 2:50/17.1 GB, default 2:30/12.6 GB. vs 3.10 outputs large diffs = 3.10↔4.1 version gap (chr8 3.10 vs 4.1 Spearman 0.93-0.94). Graphic: https://claude.ai/artifact/7FZZR6Evb6MTkUGdKBWQGX
- D-T1 revised (user): emulate numpy argsort for peak ties; attribution agent applying it.
- Full-genome benchmark (job 57977452): SRR32588726 1:34 / 16.3 GB (3.10: 2:54h / 58 GB PSS); SRR32588716 2:30 / 12.6 GB (3.10: 3:53h / 132 GB). `../bench_full/`.
- M3 (U8+U9): U8 1801/1801 blocks exact (quirks41); U9 164/164; collect 1.5-2.1 µs/read single-thread (4.1: 180-300). Fixed rayon thread-local quirk propagation bug.
- Capture batch 2: all units 0 mismatches on 22 new dirs (see CONTRACT_CHANGES).
- iqrs-compare: done; noise floor saved `../validation/noise_floor/`. chr8 E2E (`../validation/`): **parity mode passes 15/15 gates on all 5 samples**; default mode fails 1-2 gates/sample (tx/gene +1.6-4.2%, ONT fates ±1.3 pp, one Spearman 0.9949) — attribution running. iqrs chr8 runtimes 1.8-9.9 s @8 thr.
- U21 (scheduler + iqrs CLI): end-to-end wired. chr8 SRR32588726: iqrs 5.97 s/1.05 GB @16 thr (coordinator repro 9.4 s/0.79 GB @8 thr) vs 4.1 17:17/2.45 GB. vtg.bam 1:58/4.26 GB. Deterministic across threads/splits. ONT chr8 parity-identical; PacBio small residuals (triage in M5). Open: real BlockStore spill, streaming pass A, U8 collector serial section (~31k reads/s on 2.5M-read block; core use 32-39%).
- Full-genome iqrs benchmark submitted: SLURM 57977452 → ../bench_full/.
- U1 (annotation): 25/25 + 67/67 + 913/913 exact both builds; +300 generated D-R1 cases; exon_ids 320,082/320,082; parse 1.3 s / 42.5 MB.
- U4 (CIGAR/polyA/deletions): 55,683/55,683 exact (quirks41); default: 55 explained D-P1. Playbook: `docs/PILOT_U4.md`.
- U6 (assigner, U6a+U6b): 68,945/68,945 exact (quirks41); default: 63 explained D-A4.
- U5 (profiles): construct_profiles 32,635/32,635 and match_genomic_features 12,802/12,802 exact, both builds. ProfileBuilder 0.1-0.3 µs/read.
- U7 (exon correction/strand/antisense): 78,157/78,157 exact, both builds, no deviations. U8 integration notes in CONTRACT_CHANGES.md.
- U10a (peaks + XGBoost): detect_peaks_batch 43 exact + 38 D-T1; xgb_predict 143/143; random suites bit-exact. API: detect_peaks, XgbModels::bundled().
- U2 (BAM/FASTA IO): 2,059/2,059 exact + 19k edge cases (IUPAC, all HP types). noodles beats htslib 1.1-1.9x (1 thr 269k rec/s; 4 thr 1.05M rec/s no-seq). noodles stays default.
- U3 (regions/pass A): 193/193 exact both builds. 8.8M rec/s, 3.6 B/record peak.
- U12 (path threading): 134/134 exact both builds; per-graph precomputed terminal lookups, one thread per interned chain.
- U11 (intron graph): 67/67 exact (quirks41); default 66 + 1 D-G1. Python set-order emulation removed; adjacency compared unordered (no output effect).
- vtg memory prototype: giant block 10.44M reads → 62.8 B/read, 811 reads/chain, 1.07 GB peak RSS; projected full SRR32588716 @16 thr 3.5-7.5 GB vs 4.1's 138 GB. See docs/vtg_proto.md.
- U14 (assignment-based models): 67/67 exact both builds; shared `iq_model/src/fixture_util.rs` (StageInput, store_out, compare_stores, next_novel_gene). Novel/known mono-exon paths unit-tested only (fixture gap).
- U13 (FL-path models): 67/67 exact both builds; parallel per-path assign (rayon ≥32 paths). `fl::LiveAssigner` for live assignment.
- U10b (per-block polyA/TSS predictions): 138 exact + 2 D-T1 of 140, all builds.
- U16 (model splice-site correction): 172/172 + 67/67 exact both builds; 40k randomized cases vs live Python exact.
- U17 (read reassignment + ownership): 134/134 + 67/67 exact all builds; parallel == serial on every case.
- U15a (model filters): 134+67+67 exact all builds; pair-cache cuts assignments 6728→2456.
- U18 (gene joiner): 67/67 exact (quirks41); default 64 + 3 D-M6 (all novel+novel merges).
- U15b (end processing): 556/556 exact both builds; D-E1 added; U15a filter with live Ends 67/67.
- U19 (final per-read assignment): 67/67 exact (quirks41), default 58 + 9 D-A4; 92k reads; parallel == serial.
- U20 (counters + writers): 156,656 + 68 + 1,675 fixture cases exact (quirks41; default 740 D-O1). Final files from fullchr8 captures byte-identical to real 4.1 (modulo novel ids).

## Queue (start as slots free, ≤4 at once)
M5 chr8 E2E (`iqrs-compare`) → M6 full-scale benchmarks.

## How to resume after an interruption
1. `cd iqrs && cargo build --workspace --release && cargo test --workspace` to see what compiles.
2. For each unit: read `docs/progress_<unit>.md`, then `cargo run --release --features quirks41 -p iqrs-fixture -- run <UNIT> ../fixtures`.
3. Relaunch an agent per unfinished unit with: its row above + "read UNITS.md and docs/progress_<unit>.md, continue".

## User preferences (also in coordinator memory)
- ≤4 concurrent agents; ping user with genuine questions as they arise.
- Parity: "as close as possible without doing weird things"; drop quirks unless history marks them deliberate.
- Don't run full IsoQuant end-to-end; chr8/slice runs only. SLURM: `-p park,short`, default account.
- No git commits / email config unless asked; no GitHub remote.
