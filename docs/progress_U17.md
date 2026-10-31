# U17 progress (reassign.rs) — DONE

- `assign_reads_to_models(ms, ctx)` (all block reads) → `assign_reads(ms, ctx, &reads)` → `assign_reads_with_index(.., index, par_threshold)`.
  Chunks of 32768 reads; plan per read (bound / gate / profile) before the serial 4.1 loop; profile+assign (quick_mode) in
  parallel when ≥256 reads to profile in a chunk. Exact because the loop only touches the current read's slot.
- `resolve_ownership(&mut ModelStore)` (stub had an unused ctx arg; removed). Defaultdict-read inserts mirrored.
- `reference_unique_known_isoform(&ReadView)`.
- Fixtures: quirks41 134/134 + 67/67 exact; `--mode default` and plain default build identical (0 mismatches).
  Runner runs serial (threshold MAX) and forced-parallel (threshold 0) and errors if the snapshots differ: equal on all 134.
  Production path (provisional-id `index_from_models`) also 134/134 (checked by temporary edit).
- D-M5 not exercised (no duplicate read names in any U17 case); not emulated under quirks41 (as U13/U14).
- Resolve fixture coverage thin (7 cases with nonfl knowns, 1 changes, 0 drops) → 5 unit tests.
- clippy clean (both builds) for reassign.rs.
