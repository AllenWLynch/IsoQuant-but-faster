# U20 progress (counters + writers + IDs)

Build: `export CARGO_TARGET_DIR=$HOME/.cache/iqrs-target/U20`.
Deps added to iq_quant (approved by main): iq_annot, rustc-hash, flate2; quirks41 forwards to iq_annot.

## Done
- counter.rs: FeatureCounter/TxCounter (dense ModelIdx, groups first-touch SmallVec), Strategy (UniqueOnly; USC only for gene fixtures).
- ids.rs: IdDistributor (forbidden ids), NovelNamer (tx number then new-gene number, storage order), ExonIds (D-O1 gate).
- gtf.rs: GtfChrWriter::dump_block + BlockAnnot trait; AnnotBlock (iq_annot adapter), print_models.
- output.rs: ChrCounts (merge by final id), ChrPhaseW, natural fragment order (per file suffix), assemble/write_outputs, fmt_fixed2.
- fixture.rs: 3 runners. Prescan of the unit file rebuilds the exon-id dict (from recorded bare ids) and case-key lookup for confirm check.
- Fixtures: parity 156656/68/1675 exact. Default: gff_dump 935 exact + 740 D-O1; plain default 740 mismatch.

## Findings
- 4.1 gene strand in GTF = *last* model's strand (record rebuilt per model), not first (design §8.1 says first).
- Gene-model counters use unique_splicing_consistent (non_intronic single feature counts 1.0).
- counter_dump .usable cannot be reproduced (reads_for_tpm not recorded) -> excluded.

- E2E (fixture.rs `e2e`): for both fullchr8 captures, sequential printing over all blocks plus counts from counter_dump plus `assemble`
  gives byte-identical OUT.transcript_models.gtf (incl. header), discovered_transcript_counts.tsv and the grouped linear tsv
  (only D-O1 exon ids normalized). The fixture outputs equal reference_outputs/*/iq/OUT (minus header), so the same holds against real 4.1.
  With D-M1 renumbering applied, the output equals 4.1's under a consistent bijective id map.
- Production path test (output.rs phase_w_tests): AnnotBlock/print_models/forbidden ids/counts; write_outputs plain and gzip.
- clippy clean in both builds; 13 tests.

## Status: DONE (pending final report)
