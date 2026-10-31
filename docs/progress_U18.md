# U18 progress (gene joiner) — DONE

- `crates/iq_model/src/joiner.rs`: `join_genes` / `join_genes_by(models, index, novel_cmp)`; registered in lib.rs.
- quirks41 parity: 67/67 exact. Default (attributed): 64 exact + 3 explained by D-M6. Plain default: 3 mismatches (= D-M6 count).
- D-M6 cases (all novel-novel; none novel-annotated): SRR32588720:nogene:8:1 (novel_gene_8_15 joins novel_gene_8_7, score 1.286),
  SRR32588720:nogene:8:2 (_18 + _23, score 1.333; 4.1 loop never runs with exactly one pair), SRR32588726:nogene:8:0 (_46 joins _6, score 1.0).
- Novel-novel merge direction in 4.1 = larger id string survives. Fixture runner passes the 4.1 string order; production
  `join_genes` uses block-local key (creation) order. 3 fixture cases differ in surviving key between the two orders,
  partition identical in all 67 (both builds) — labels only, renumbered in Phase W.
- Tests: 2 unit tests (D-M6 single annotated+novel pair; novel direction + weak overlap). Clippy clean both builds.
- Build with CARGO_TARGET_DIR=$HOME/.cache/iqrs-target/U18.
