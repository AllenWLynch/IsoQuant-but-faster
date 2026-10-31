# U5 progress (read profiles): crates/iq_assign/src/profile.rs. DONE

- Windowed sweep (§6.4.2): build_profiles, build_profiles_into, ProfileBuilder (per-index pmax), match_genomic_features
  (O(n log K + matches)), plus the `dense` reference port (tests/bench only).
- Fixtures (shared tree): construct_profiles 32635/32635, match_genomic_features 12802/12802 exact in the parity, default
  and debug builds (the debug build also debug-asserts windowed == dense on every case).
- Tests: proptest windowed == dense (4000 cases, incl. count_exons and polyA/T), mgf == dense, plus edge cases.
  Clippy is clean with and without quirks41.
- Bench: `IQ_BENCH_CAP=<capdir> cargo test --release -p iq_assign -- --ignored bench_profiles --nocapture`.
- Scratch has `pub assign: AssignScratch` for U6a. Cargo.toml has the proptest dev-dep (approved).
- Workaround: the mgf runner uses ctx.params().matching.delta (CONTRACT_CHANGES item 6). Remove after re-capture.
