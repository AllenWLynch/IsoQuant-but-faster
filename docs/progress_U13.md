# U13 progress (FL-path model construction)

## Status
- [x] Read contract, Python, fixtures survey (67 cases, all in.store has 0 models; nonfl empty; fl assign quick_mode=false,
      params = run args).
- [x] fl.rs implementation + runner (uses shared fixture_util from U14)
- [x] parity 67/67 exact; default 67/67 exact (both binaries)
- [x] equivalence checks: count pre-filter off and all-parallel / all-serial give the same 67/67; tests (4) pass; clippy clean both builds
- Reference sensitivity: with an empty reference 11/67 cases mismatch, so the genome workaround is load-bearing.
- Plumbing switched to the shared `fixture_util` (StageInput / store_out / compare_stores / next_novel_gene).

## Findings
- Strand of novel introns needs the reference (`get_intron_strand`); U13 `in` carries no RefWin. Runner workaround:
  read `data/ref/genome.fa` (+.fai; same file the capture staged, capture_slice.sh:21), env `IQRS_REFERENCE` overrides.
  Contract request: add per-FL-path-intron `intron_sites` (like U7) to U13 `in`.
- `select_reference_gene` ties need gene-id string order; IsoformIndex has no gene id rank. Workaround: `FlAux.gene_rank`
  slice aligned with `index.genes`. Contract request: `IsoformIndex::gene_id_rank`.
- IDs: compared by creation order (novel models and novel genes mapped k-th to k-th of the expected store).
