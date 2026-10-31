# U15b progress (terminal end processing)

## State: done
- `crates/iq_model/src/ends.rs`: `Ends` (implements `filter::EndProcessor`, `ties` field for D-T1),
  inherent `derive_alternative_end_models` / `alternative_ends` / `nic_model`; 5 fixture runners with a
  D-T1 fallback (reverse tie order → Explained(T1)) that also honours parity/default quirk sets.
  lib.rs: `ends::register_fixtures(reg);`. The stub's own unused `ends::EndProcessor` trait was removed.
- D-E1 (approved by main, DevId added by main): mono-exon end edit overwrote the start edit in 4.1;
  gated at one site in `correct_ends`; unit test covers both branches. 0 fixture cases exercise it.
- Fixtures: 556/556 exact, parity and default (0 explained). No D-T1 fallback needed.
- Integration: ignored test `ends::tests::filter_with_live_ends` (U15a filter with live Ends): 67/67.
  Ignored `ends::tests::u15b_fixtures` runs U15b via the crate registry (iqrs-fixture binary was
  blocked by an unrelated iq_quant compile error at the time).
- Perf: add_known derives all eligible models in parallel (rayon), then sequential dedup/moves (exact);
  mark confirmation does no peak calls; peak calls are direct tree walks (no per-call setup).
- Coverage gaps: only 2 non-empty derive cases, 1 changing add_known case, novel_apa never on.
