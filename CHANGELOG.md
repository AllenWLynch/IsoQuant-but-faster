# Changelog

## 0.1.0 (unreleased)

First public release.

- Optional phased quantification (`--phase-vcf`): allele-specific transcript counts from a
  phased VCF via a per-gene Bayesian model, with a per-transcript `phase_uncertainty`.

- Transcript discovery and quantification matching IsoQuant 4.1.0 for one bulk,
  coordinate-sorted, indexed BAM. PacBio is `sensitive_pacbio` with `--fl_data`; ONT is
  `sensitive_ont`.
- Options: `--report-novel-unspliced`, per-tag grouped counts (`--read-group tag:TAG`) and
  `--gene-types isoquant|all` (default `isoquant`; `all` is not validated against 4.1).
- `--threads` defaults to all CPUs available to the process (respects SLURM/cgroup limits).
- Outputs: `transcript_models.gtf`, `discovered_transcript_counts.tsv` and the grouped linear
  counts table, in IsoQuant 4.1's formats.
- IsoQuant 4.1's quirks are reproduced (see `docs/dev/deviations.md`). There is no separate
  "default" behaviour.
- Clear errors for missing BAM or FASTA indexes, name-sorted or out-of-order BAMs, compressed
  references and missing input files.
- Benchmark: whole-genome PacBio Revio sample, 16 threads. iqrs takes 1:48 and 19.7 GB.
  IsoQuant 4.1 takes 1:55:52 and 81.4 GB.
