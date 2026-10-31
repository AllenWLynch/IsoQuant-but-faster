# U6b progress (inconsistent path + JunctionComparator): DONE

## Status (2026-10-29)
- `crates/iq_assign/src/junction.rs`: `JunctionComparator` (compare_junctions, contradiction classification,
  mono-exon subtype, extra flanking events, are_suspicious_introns, are_known_introns). `are_known_introns` is an
  exact windowed version of 4.1's profile sweep: it starts at the first known intron the sweep can stop at
  (binary search on a lazily built prefix max of known-intron ends). The full reference sweep is kept as a test
  oracle and checked against random inputs, including tiny-exon cases.
- `crates/iq_assign/src/inconsistent.rs`: match_inconsistent, select_similar_isoforms (+ match_profile Hamming),
  detect_inconsistencies, select_best_among_inconsistent (f64 penalty in event order), create_*_matches.
- Results: parity (quirks41) U6 68945/68945 exact. Default with quirks41: 68882 exact + 63 explained by D-A4.
  D-A3 never fires in the fixtures. JC fixture test: 14608 matches across 11830 inconsistent cases, 0 mismatches.
- Tests: `cargo test --release -p iq_assign` (junction + inconsistent unit tests); harness
  `u6b_inconsistent` (ignored; `IQRS_U6B_ALL=1` runs all U6 cases).

## Perf follow-ups (not blocking)
- A JunctionComparator is built per read, so the prefix max (O(K)) is rebuilt per inconsistent read that needs
  a known-intron check. Cache it per index (in AssignScratch keyed by index, or in Assigner) when convenient.
- ReadMatches allocates per read (SmallVec<[(u32, Events);4]>); could move to scratch.
