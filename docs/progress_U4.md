# U4 progress (pilot): DONE

- Files: crates/iq_read/src/{cigar.rs, polya.rs, deletions.rs}, lib.rs (process_alignment[_into], register_fixtures);
  docs/PILOT_U4.md.
- Parity (quirks41): extract_relevant_deletions 27831/27831, process_alignment 27852/27852 exact.
- Default (quirks41 build, --mode default): deletions 27831 exact; process_alignment 27797 exact + 55 D-P1.
  Plain default build: the same 55 show up as mismatches.
- D-P2 unreachable (randomized test; §9 now NONE); inert gate kept. D-P3: 0 fixture diffs, covered by the
  polyt_trim_and_d_p3_anchor unit test.
- Workaround: deletions runner infer_pos0 (Err on ambiguity) until `pos0` is in the schema (CONTRACT_CHANGES.md).
- Tests 9 pass both builds; clippy clean both builds.
