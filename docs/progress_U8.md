# U8 (collect driver) + U9 (multimapper resolution) — progress

Build: `export CARGO_TARGET_DIR=$HOME/.cache/iqrs-target/U8`.
Deps added to iq_collect (approved by main): iq_annot, iq_bam, iq_regions (+ their quirks41 features).

## Status (2026-10-29)
- U8 collect_block: 1675/1675 exact (quirks41). Default mode: 1533 exact + 142 explained (D-R2 92, D-P1 30, D-A4 20), 0 mismatch.
- U9 resolve: 81/81 exact, both modes. D-C1 gate is in `BasicRead::from_view` (U9 fixture inputs carry penalty already).
- Runner: real BAM records from `<capdir>/slice/*.bam`, live index from the GFF3 (iq_annot, D-R1), whole-chr reference
  from data/ref/genome.fa, ownership among all U8 regions of the chromosome (ClusterOwner, D-R2 gate).

## API (crates/iq_collect/src/lib.rs)
- `BlockCollector::new(&Block, &IsoformIndex, &RefSlice, &Params, &dyn Owner, &dyn GroupIds)`, `.push(&[RecordLite])`,
  `.finish() -> (BlockStore, CollectStats)`; `.parallel(bool)`, `.timed()`.
- `collect_block(...) -> BlockStore` (design signature), `collect_block_with(..., groups) -> (BlockStore, CollectStats)`.
- Owners: `AllRecords`, `ClusterOwner{regions, me}`. Groups: `NoGroups`, `SharedGroups`, `FnGroups`.
- multimap: `duplicate_groups(&[&BlockStore])`, `BasicRead::from_view`, `resolve(&mut [BasicRead], strategy)`, `apply_resolution`.

## Done since
- Rayon tasks inherit the thread-local quirk set (`with_quirks(current_quirks())`); before that, default-mode runs of
  ≥64-read batches silently ran with all quirks on.
- tests/u8_fixtures.rs: parallel == serial (1701 blocks, 292,954 reads); per-read attribution (ignored test):
  873/277,170 reads differ in default build: D-P1 499, D-R2 250 (0.09%), D-A4 114, D-P3 6, 4 reads need two deviations
  (R2+A4, P1+A4), 0 unexplained.
- `multimap::resolve_chr_duplicates(stores, dup_candidates, strategy) -> DupOutcome` for U21; suspended reads are
  removed from stores (dense re-index = 4.1 loader). Collect+resolve == `block/` files on 1701 blocks / 22 captures;
  total/polyA == 4.1 log line in all 22. (Penalties quantized int(p*2^20) in the comparison: 4.1 dump artifact.)
- Told main: U21 RefSlice must cover index [start,end] too (terminal-misalignment isoform introns).

## Latest (all finished captures incl. batch 2, 1794 blocks / 734,395 reads)
- parallel == serial on all 1794 blocks; collect+resolve == block/ files on all 1794; total/polyA == 4.1 log in 38 captures.
- per-read default attribution: 2015 differing reads: D-R2 1258 (0.17%), D-P1 608, D-A4 137, D-P3 6, 6 two-deviation; 0 unexplained.
- clippy clean (default + quirks41).

## Throughput (default build, single thread, prefetched records; fetch/index build excluded; node load avg ~83 on 8 cores)
| capture | blocks | records | stored | µs/stored read | reads/s | select | U4 | U5 | U6 | U7 | store |
|---|---|---|---|---|---|---|---|---|---|---|---|
| SRR32588720 Iso-Seq fullchr8 | 708 | 122,456 | 113,670 | 2.08 | 480k | 0.07 | 0.59 | 0.26 | 0.55 | 0.42 | 0.30 |
| SRR33046254 ONT fullchr8 | 900 | 74,427 | 70,956 | 1.53 | 655k | 0.08 | 0.54 | 0.12 | 0.37 | 0.22 | 0.29 |
Run: `RAYON_NUM_THREADS=1 cargo test --release -p iq_collect --test u8_fixtures fullchr8_throughput -- --ignored --nocapture`.

## Open
- Penalty quantization int(p*2^20) in 4.1 dump not emulated (no construct reader).
- U21 must pass a RefSlice covering index [start,end] ∪ cluster ±1000 (terminal-misalignment isoform introns).
