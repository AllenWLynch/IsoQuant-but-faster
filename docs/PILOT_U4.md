# Pilot playbook (U4): fixture → Rust → pass

## Structure that worked
- Keep the core fns free of fixture types. U4 has `cigar::read_blocks_into`, `polya::add_polya_info`,
  `deletions::extract_relevant_deletions_into`, and `process_alignment_into(rec, &PolyAParams, ext, &mut AlignedRead)`.
  `process_alignment(rec, &Params)` is a thin wrapper around them.
- Each module has `pub(crate) mod fixture { pub struct U4Xxx; impl UnitFixture for U4Xxx {...} }` and a
  `pub(crate) fn register_fixtures(reg)`. Crate `lib.rs::register_fixtures` only calls the module fns.
- Runners should call the **real entry point**, not internals, and convert at the edge. Build parameters from the case
  `in` (e.g. `polya_window`), not from `ctx.params()`: the case is the ground truth for that call. If the fixture needs
  an intermediate state the entry fn doesn't return (U4 `init` = blocks before trimming), recompute it in the runner.
- Set `DEVIATIONS` per runner (U4: `[P1, P2]` / `[P3]`) so an unrelated ID can't "explain" a mismatch.

## Commands I actually used
```
export PATH=$HOME/.cargo/bin:$PATH
cargo build --release --features quirks41 -p iqrs-fixture
./target/release/iqrs-fixture run U4 ../fixtures                      # parity: must be 100% exact
./target/release/iqrs-fixture run U4 ../fixtures --mode default       # same binary: attribution per DevId
cargo run --release -p iqrs-fixture -- run U4 ../fixtures --mode default --examples 0   # plain build
cargo test -p iq_read [--features quirks41]; cargo clippy -p iq_read --all-targets [--features quirks41]
```
- Calling the binary directly avoids rebuild noise (~56k cases in about 1 s). `--examples N` prints field diffs.
- The plain default build can't attribute, so its mismatch count must **equal** the explained count from the quirks41
  `--mode default` run (U4: 55 = D-P1 55).
- Survey a unit's fixtures with a 15-line python/zcat script first: fn counts, exception counts, input values that never
  vary. U4: no H/P ops, no IUPAC bases, `external_strand` always null, and only 2 reads with `exons_changed`.

## Quirk gating and verifying attribution
- One site per DevId: `if quirk(DevId::P1) { 4.1 } else { iqrs }`. Make the switch a *bound*, not duplicated loops
  (U4 find_polya: `limit = if quirk(P1) { n-w } else { n-w+1 }`).
- Fixtures may never exercise a DevId (U4: D-P2 and D-P3 differ in 0 cases). Then write a unit test that runs both
  sides via `with_quirks(QuirkSet::only(id), ..)` vs `QuirkSet::NONE`. Guard the 4.1-side asserts with
  `if iq_types::QUIRKS41`, because without the feature `with_quirks` is a no-op.
- Check that each quirk is reachable at all. D-P2 turned out to be unreachable on sorted, disjoint exons, which a
  randomized test proved. §9 now lists it as NONE, and the inert gate stays.
- If 4.1 *crashes* on an input, choose sane behaviour, document it at the site, and tell main. No DevId is needed.

## Fixture / ID pitfalls
- Case keys look like `<sample>:<slice>:<chr>:<block>:<fn>:a<aln_no>.<n>`. Two fns of one unit share an alignment but
  live in separate cases, and the runner can't reach the sibling case.
- A DevId's site can sit in the *caller* while the fixture captures the *callee*. Here D-P3 (deletion anchor) is
  chosen in process_alignment, but the deletions fixture `in` carries only the 4.1 anchor (`ref_start`), not `pos0`.
  Workaround: `deletions::fixture::infer_pos0` aligns the trimmed exons to `get_read_blocks(0, cigar)`. It runs only when
  `quirk(P3)` is off and **returns Err (= mismatch) when the reconstruction is ambiguous or fails**, so it can't hide
  errors. Queued fix: add `pos0` to that `in` (CONTRACT_CHANGES.md); then delete `infer_pos0`.
- `AlnFixture::to_record()` gives `RecordOwned`. Call `.as_lite()` on it. Cigar in fixtures is `[[op,len]]`, and
  `(len << 4) | op` turns it into raw words.
- 4.1 `pysam.reference_end` = `RecordLite::end0()`. `reference_start` = `pos0`, 0-based. Watch the ±1 at every boundary.
- Check arithmetic in hand-written test expectations twice. Both of my early test failures were wrong expectations,
  not wrong code. The fixtures are the real oracle.

## Awkward bits in iq_types / runner (suggested fixes)
1. No sibling lookup between cases: add `FixtureCtx::sibling(case_key, fn)` (queued).
2. `iq_read/lib.rs` didn't delegate registration per module (iq_assign did). U7 must add
   `correct::register_fixtures(reg);` there.
3. `AlignedRead` has no `Default`, so every call site spells out six empty fields. A derive would help (contract change).
4. `quirk()` only works with the `quirks41` feature on the *dependent* crate. Add `quirks41 = ["iq_types/quirks41"]` to
   your crate's features (iq_read already has it).
5. A runner `Err` and a panic both count as `error` inside `mismatch`. Use `Err` deliberately for "can't judge".
