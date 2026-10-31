# U19 progress (final_assign.rs) — DONE

- `final_assignments(ms, ctx, sink: impl FnMut(u32, FinalAssign))` → `for_each_final(ms, ctx, reads, full, single, par_threshold, map, emit)`.
  Route per read (construction / known-reuse / single re-profile / discarded / full / noninformative). Each needed model's
  single-model index + ProfileBuilder is built once per call (in parallel). Reads are evaluated in parallel in chunks of 16384
  and emitted in read order, so no per-read assignment is kept.
- Single path: a `LiveAssigner` with `delta = Params::final_single_delta()`, because the ReadAssigner trait can't carry delta.
  Leftovers use ctx.assigner.
- The runner compares full assignments, which is stronger than the compact contract. It runs serial vs forced-parallel
  and must get the same result. It also runs the production entry with provisional ids, compared compactly, ignoring tx order.
- quirks41: 67/67 exact. --mode default: 58 + 9 D-A4 (U6). Plain default: 9 mismatches (= D-A4 count).
- Route coverage over 92311 fixture reads: reuse 70179, single 7760 (1996 forced), full 11877, discarded+gene 1248,
  discarded 174, noninformative 1073.
- 4 unit tests, clippy clean in both builds.
