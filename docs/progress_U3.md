# U3 progress (pass A, split_coverage_regions, ownership)

## Status
- [x] split.rs: split_coverage_regions verbatim, skip_read_fraction, coverage_cutoff/is_small_chr
- [x] pass_a.rs: streaming PassA (diff-array bins, gap check via GeneOracle, emit clusters/blocks)
- [x] owner.rs: owner_block/owner_of (D-R2), overlapping_blocks (4.1 side for U8), OwnerCounter (midpoint histogram + explicit long reads)
- [x] dups.rs: scalable blocked Bloom for duplicate-name candidates
- [x] fixture.rs: 4 runners; cluster/process run real PassA over <capdir>/slice/*.bam (mini BGZF reader) + gene double from fixtures/annotation
- [x] quirks41 parity: 193/193 exact (cluster 33, genes_in_region 60, process 67, split 33)
- [x] default build: 193/193 exact (D-R1 never fires on fixtures; unit test fixture::tests::d_r1_gap_bridge_both_sides covers both sides)
- [x] clippy clean both builds; 15 unit/property tests
- [x] bench: examples/bench_pass_a.rs, 10.0M records -> 8.8 M rec/s (1.13 s, incl. xxh3 name hashing), peak heap 35.9 MB = 3.6 B/record (mostly the dup Bloom filter; ~1.5 B/rec with expected_primary hint)
DONE.

## Facts
- 4.1 crashes on records with reference_end None (unmapped flag / no CIGAR); iqrs skips them (stats.skipped).
- Short reads (span <= 32768) straddle at most one block boundary (inner blocks >= 32768 bp), so owner = block of midpoint key ((s+e)/2 - 1) div 256; long/empty spans resolved explicitly.
- D-R1 site is iq_annot::AnnotationIndex::genes_for_region0; GeneOracle::genes_for_region0(Iv0) should delegate to it (no extra shift in U3).
