# U7 progress (done)

Status: complete. Parity and default runs are both 100% exact. Fixtures: 25163 `correct_assigned_read`, 27831
`get_assignment_strand` and 25163 `check_antisense` cases.

## Files
- `crates/iq_read/src/correct.rs`: core functions, runners and tests.
- `crates/iq_read/src/groups.rs`: `group_string`.
- `crates/iq_read/src/lib.rs`: one line, `correct::register_fixtures(reg)`.
- `crates/iq_read/Cargo.toml`: `iq_assign` dependency (approved by main) for `match_genomic_features`.

## API (differs from the design §3 sketch)
- `correct_read<R: RefBases + ?Sized>(&AlignedRead, &RecordLite, &Assignment, &IsoformIndex, &R, &Params)`
  and `correct_read_with(exons, rec, a, index, ref, &CorrectParams)`.
- `assignment_strand(a, n_exons, corrected_introns: impl IntoIterator<Item=Iv>, polya, &R)`. 4.1 uses
  `len(read_assignment.exons)`, i.e. the read exons, not the corrected ones. `introns_of(&exons)` gives the
  introns without allocating.
- `RefBases` trait: implemented for `RefSlice` (pipeline: FASTA mmap) and `[RefSlice]` (fixture windows).

## Findings
- The aligned-pairs list is virtual: a monotone CIGAR seek gives the pair index, and only the ≤24-pair window
  is generated. A randomized-window test checks it against a naive port of pysam pairs plus `get_error_count`.
- pysam pairs include N (ref-skip) positions. 4.1 compares 0-based `ref_pos` with 1-based exon coords
  (start index = the exon's 2nd base, end index = the base after the exon). Ported verbatim.
- Instrumented run: 2212 error-count windows are exercised, some with mismatches. 0 hits on the
  IndexError/KeyError fallbacks and 0 terminal start-index cases (only reachable with a 1-bp last exon at the
  alignment end; covered by a unit test).
- The StrandDetector cache is not ported (pure memo; outputs identical).
- `fixtures_u7` test runs all U7 fixtures when `IQRS_FIXTURES` is set (useful when another crate breaks
  the iqrs-fixture build).
