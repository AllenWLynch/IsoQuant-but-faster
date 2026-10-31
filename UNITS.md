# Rules for unit agents (read before touching code)

You are one of up to 4 agents working **concurrently in this single working tree** (no git commits, no branches).
Project background: `../specs/00_scope.md` (scope, parity policy), `../specs/07_rust_design.md` (§3 is your unit's contract;
§9 deviation log), `harness/iqcap/SCHEMA.md` (fixture JSON), `crates/iq_types/src/lib.rs` docs (worked `UnitFixture` example).

## Contract freeze
- `iq_types` v1 and fixture schema v1 are **frozen** (2026-10-29). Do not edit `crates/iq_types/`, `harness/`,
  `Cargo.toml` (workspace root), or another unit's files.
- Need a contract change (new field, type fix, schema bug)? SendMessage to `main` with the concrete diff and why. Work around it
  locally (private helper types in your crate) until the coordinator applies it.
- Do not add crate dependencies without asking `main` first (Cargo.lock is shared). Dependencies already declared are fine.

## File ownership
Only edit the files assigned to you in your brief. Shared `lib.rs` files already delegate `register_fixtures` to each module, so
register your runners inside your own module.

## The loop
1. Read your unit in design §3 and the 4.1 Python it ports (`../IsoQuant/isoquant_lib/...`). The Python is the spec; the design
   doc and specs 01-05 are guides, and may be wrong — trust the code, and report doc errors to `main`.
2. Implement a `UnitFixture` per hooked `fn`; register it.
3. Parity build: `cargo run --release --features quirks41 -p iqrs-fixture -- run <UNIT> ../fixtures` must reach **100% exact**.
4. Default build: `cargo run --release -p iqrs-fixture -- run <UNIT> ../fixtures --mode default` — every mismatch must be
   attributed to a deviation ID (gate dropped quirks with `iq_types::quirk(DevId::...)` at the single site where they live).
5. Write focused unit tests for edge cases you discover. `cargo clippy -p <crate> --all-targets` clean, with and without `quirks41`.

## Parity policy (from the user)
"As close as possible without doing weird things." Do not port things that look like bugs/quirks unless a code comment, changelog
or commit history marks them deliberate (`git -C ../IsoQuant log -S ...`). Every dropped quirk is a `DevId` behind `quirk()`. If
you find a quirk not in §9, ask `main` before adding a DevId.
- Whole-run attribution (quirks41 builds only): env `IQRS_QUIRKS=all,-D-M6` (parity minus one) or `IQRS_QUIRKS=none,+D-M6`
  (default plus one) sets the starting quirk set on every thread; unset = all. See `iq_types::dev` docs.

## Performance
This rewrite exists to be ≥5x faster and leaner. Prefer flat arrays, no per-read allocation in hot paths (use `Scratch`),
asymptotically better algorithms when they give identical outputs. Correctness first, but don't write throwaway slow code.

## Communication
- Questions that are genuinely the user's (scope, parity tradeoffs): SendMessage `main` with your recommendation; keep working.
- Engineering decisions: make them yourself.
- Final report (≤20 lines): exact/explained/mismatch counts per fn (both builds), deviations used/added, contract change requests,
  open issues, files touched.

## Build isolation (added after repeated link collisions)
Concurrent agents sharing `target/` collide (spurious link errors). Always build/test with your own target dir:
`export CARGO_TARGET_DIR=$HOME/.cache/iqrs-target/<your-unit>` (e.g. `.../U15b`). Final verification must still be in
this shared source tree (just not the shared target dir).
