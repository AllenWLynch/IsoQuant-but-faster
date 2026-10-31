# U6a progress (consistent path + PolyAVerifier + runner)

Status: DONE. With U6b's inconsistent.rs, U6 parity (quirks41) = 68945/68945 exact.
Default build (quirks41 attribution) = 68882 exact + 63 explained by D-A4, 0 mismatches.

## Files
- `crates/iq_assign/src/assigner.rs`: entry `assign` / `assign_with`, `Assigner<'a>` (index, params, polya, scratch),
  `route()` dispatch, the consistent path, shared helpers, the `U6Assign` fixture runner, and the ignored test `u6_breakdown`.
- `crates/iq_assign/src/polya_verify.rs`: `PolyAVerifier` and `shift_polya`/`shift_polyt`.
- Scratch: `profile::Scratch.assign: AssignScratch` (c0..c4 for the consistent path, inc0..inc2 for U6b, scores).

## Breakdown (parity, combined, all exact)
collect 41828 · final_full 7199 · final_single 5100 · fl 619 · reassign 10784 · similar 3415.
Route as computed here: consistent 37450, early 6542, inconsistent 24953.
D-A4 (default build) hits only the inconsistent route: tie resolution in select_best_among_inconsistent
(collect 4, final_full 58, fl 1). D-A3 never fires in these fixtures.

## Breakdown tool
`cargo test --release -p iq_assign --features quirks41 u6_breakdown -- --ignored --nocapture`
Env vars: IQRS_U6_CALL_SITE, IQRS_U6_ROUTE (early|consistent|inconsistent), IQRS_U6_ATYPE (expected code),
IQRS_U6_DEFAULT=1 (default mode with attribution), IQRS_U6_EXAMPLES=N, IQRS_FIXTURES.

## Workarounds
- Model/region indices have no gene regions (4.1 `get_gene_regions()` iterates the empty `gene_db_list`), so
  `assign_to_overlapping_genes` ignores `index.genes` when the isoform keys are `IsoRef::Model`. A flag request is in CONTRACT_CHANGES.md.
