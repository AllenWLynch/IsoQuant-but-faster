# iqrs-compare progress notes

## Findings (session 1)
- Reference dirs: `../reference_outputs/<s>/iq/OUT/OUT.*`. No grouped files there (unphased runs) -> grouped path tested synthetically.
- 4.1 per-read file vs models: `OUT.transcript_model_reads.tsv.gz` (written by default): 2 `#` lines, header
  `read_id chr strand gene_id gene_assignment_type isoform_id isoform_assignment_type ...`, one row per (read, model);
  isoform_id `.` when none. `OUT.read_info.tsv.gz` is the same layout but vs the reference annotation (not what we want).
- Fates seen: unique, unique_minor_difference, ambiguous, inconsistent, inconsistent_ambiguous, inconsistent_non_intronic,
  noninformative, discarded (+ intergenic in read_info). `__no_feature` = noninformative+intergenic+discarded (checked 2152+2123=4275).
- Noise floor pair: 1374 vs 1375 transcripts (one extra .nnic in nosimd).
- Novel IDs: `transcriptN.<chr>.(nic|nnic)` with source IsoQuant; known: `transcript:ENSDART...`, annotation source.

## Plan / status
- [ ] modules: io (gz-aware readers), model (keys + union-find matching), metrics, annot, report
- [x] modules written: io.rs, matching.rs, annot.rs, metrics.rs, report.rs, main.rs (builds clean)
- [x] noise floor run (SRR32588720 vs nosimd): 15/15 gated pass, 0 / 1 unmatched (transcript2092.8.nnic, 2 reads,
      1 intron shifted 5 bp vs ENSDART00000005989); chain t->r 0.99998; spearman 1.0; per-read: 2 reads inconsistent->unique.
      1.4 s incl. full GFF3 parse.
- [x] self-comparisons (6 refs): every deficit exactly 0, 15/15 gated pass.
- [x] perturbed (SRR32588716: drop 114, shift first-exon end +3bp on 108, x2 counts on 217, novel IDs renamed, 2% reads
      flipped to noninformative, synthetic HP grouped A/B/NA, test reads in proposed iqrs format): 4/24 pass;
      chain r->t 0.9435, t->r 0.9829, spearman 0.984, tx/gene -3.6%, fate unique -1.72pp / no_feature +1.93pp,
      grouped per-group spearman 0.977-0.988; residuals explain shifts ("1 intron differs, max shift 3 bp").
      Generator: scratchpad scripts/perturb.py (not kept in repo).
- [x] per-read storage: sorted Vec<(hash, fine<<28|model)> (16 B/read). 25x replicated revio726 (40M reads/side): 24.5 s, 1.3 GB.
      Reads at >1 locus (36 in 716 chr8): keep the most-assigned record (smallest fine index), order-independent.
- [x] unit tests (6) pass; clippy clean.
- Deps added to bins/iqrs-compare/Cargo.toml (all already workspace-declared, no Cargo.lock change): iq_annot, flate2, rustc-hash.

## Usage
iqrs-compare --ref DIR --test DIR [--annotation ../data/ref/annotation.gff3] [--mode parity|default]
             [--noise-floor floor/compare.json] [--ref-reads F --test-reads F | --no-reads] [-o OUT] [--fail-on-gate]
Per-read files auto-detected: *transcript_model_reads.tsv[.gz] (4.1) or *read_fates.tsv[.gz] (iqrs).

## Proposed iqrs --emit-read-fates format (for U21)
`OUT.read_fates.tsv.gz`: line 1 `#iqrs read fates v1`, header `read_id\tchr\tmodel_ids\tassignment_type\tgroup`,
one row per (read, locus) after U19; model_ids = final printed model IDs comma-joined (`*` if none);
assignment_type = 4.1 names (unique, unique_minor_difference, ambiguous, inconsistent, inconsistent_ambiguous,
inconsistent_non_intronic, noninformative, intergenic, discarded); group = read-group value or NA.
