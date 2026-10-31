# Contributing to iqrs

iqrs must produce the same output as IsoQuant 4.1 for the configuration in the README. When a
change could affect output, check it against IsoQuant fixtures and end-to-end outputs before you
send it.

## Layout

```
bins/iqrs            the CLI (the only shipped binary)
bins/iqrs-compare    scores two IsoQuant-format output directories (models matched by structure)
crates/iq_types      shared types, BlockStore, IsoformIndex, Params presets, fixture IO and runner
crates/iq_annot      GFF3/GTF parsing, gene index, per-block isoform index
crates/iq_bam        BAM fetch (noodles; optional rust-htslib backend), FASTA mmap
crates/iq_regions    pass A: coverage clusters, block splitting, read ownership
crates/iq_read       CIGAR blocks, polyA/T, deletions, exon correction, strand, read groups
crates/iq_assign     read profiles, isoform assignment
crates/iq_collect    per-block read collection, multimapper resolution
crates/iq_terminal   polyA/TSS peak calling and the bundled XGBoost models
crates/iq_graph      intron graph and path threading
crates/iq_model      transcript model construction, filtering, read reassignment, gene joining
crates/iq_quant      counters, GTF/TSV writers, novel ID numbering
crates/iq_sched      work units, memory admission, rayon orchestration
dev/iqrs-fixture     runs the unit fixture runners against captured IsoQuant fixtures
dev/harness          iqcap, the Python harness that captures fixtures from a live IsoQuant 4.1 run
docs/dev             developer notes (deviation log, release notes)
```

Code comments name the IsoQuant 4.1 source each function ports (`file.py:lines`, relative to
`isoquant_lib/`). Where iqrs keeps a surprising 4.1 behaviour on purpose, the comment says "matches IsoQuant 4.1".

## Build and test

```bash
cargo build --release
cargo test --workspace --release
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all --check
```

CI runs these steps. The unit tests are self-contained. Tests that need the large fixture
captures skip with a message when the captures are absent. They look in `$IQRS_FIXTURES` or
in `../fixtures` next to the checkout.

## Fixtures (IsoQuant parity at unit level)

`dev/harness` runs IsoQuant 4.1.0 unmodified and records the inputs and outputs of each ported
function as gzip JSON Lines. See `dev/harness/README.md`, and `dev/harness/iqcap/SCHEMA.md` for
the format. A capture needs an IsoQuant 4.1.0 checkout, a Python environment with IsoQuant's
requirements, and the input BAM, genome and gene DB. The captures are too large for the
repository. To replay them:

```bash
cargo run --release -p iqrs-fixture -- run all ../fixtures        # every unit
cargo run --release -p iqrs-fixture -- run U6 ../fixtures         # one unit
```

The U8 runner also needs the genome: `$IQRS_REFERENCE`, or the first `data/ref/genome.fa` found
in a directory above the fixtures. Every case must be `exact`. The only accepted `explained` cases are `D-T1` peak tie-order cases,
listed in `docs/dev/deviations.md`.

## End-to-end check

Run the old and new binaries on the same inputs and compare. `transcript_models.gtf` and the
counts files must be byte-identical apart from the `# Command line:` header. Against an IsoQuant
4.1 output directory, use:

```bash
iqrs-compare --ref isoquant_out/OUT --test iqrs_out --annotation genes.gff3 -o compare_dir
```

## Performance

Prefer flat arrays and per-thread scratch buffers, and avoid allocating per read in hot loops.
Results must never depend on the thread count or on how work is split. Set
`IQRS_DEBUG_UNIT_SPLIT=whole-chr|max` to check that. `cargo build --profile profiling` keeps line
tables for profilers, and `iqrs --profile` prints time and memory per phase.

## License

By contributing you agree that your contributions are licensed under GPL-2.0-only.
