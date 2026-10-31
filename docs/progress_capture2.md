# Fixture capture batch 2 — progress notes

## Harness changes (hooks.py) — done
- #1 U4 extract_relevant_deletions `in.pos0` (caller's `alignment.reference_start`; fallback CTX.cur_pos0)
- #6 U5 match_genomic_features `delta` = `self.comparator.keywords["delta"]` (new captures only)
- #14 U13 `in.intron_sites`, U16 store-level `in.ref`
- TODO: SCHEMA.md changelog (+ #3 entries), capture_slice.sh env overrides

## Captures
(see below)
- SCHEMA.md: U4 pos0, U5 delta note, U13 intron_sites, U16 ref, manifest slice_def.keep_names, changelog (#3 + batch 2) — done
- slice.py: `--keep-names FILE`; capture_slice.sh env: NOVEL_UNSPLICED, SRC_BAM, KEEP_NAMES, EXTRA_IQARGS, BLOCK_RATE/CAP
- Stage refs: scratchpad/cap2/stage. Diagnostics (full chr8, units U3,U13,U14,U15b,U17, novel_unspliced true):
  scratchpad/cap2/diag.sh SRR33046254 nanopore / SRR32588720 pacbio_ccs -> scratchpad/cap2/diag/<S>/fx
- (c) trigger: polya_fraction = polyA_found assignments / total >= 0.7 (args.polya_percentage_threshold, auto).
  Built b2_polyarich (ONT, valley+ordinary+multigene regions, keep all-polyA reads (from M0a U4) + 1/4 as many non-polyA
  -> ~80%); names in fixtures/SRR33046254/b2_polyarich/slice/keep_names.txt. Launched locally.
- (f) SLURM 57975407: SRR32588716 (vtg.bam = SRR32588716 PacBio) b2_vtg22sub 22:24462404-25048722 sub 0.02, U3,U8,U11, BLOCK_CAP=1e6
- b2_valley (ONT, all units) launched locally to validate harness changes
- b2_valley ONT: validate 0 errors; pos0 present; U5 exon_corrector delta=6 (was 0); U13 intron_sites, U16 ref OK.
  Rust check: non-roundtrip only on U4 erd (pos0), U13 (intron_sites), U16 store-level (ref) -> DTO fields needed.
- b2_polyarich: polyA 80.0%, requires_polya_for_construction=True, require_mono*_polya=True. validate 0 errors. (c) CLOSED.
- Diagnostics done (scratchpad/cap2/analyze.py): ONT chr8: 21 novel mono, 3 known mono, 0 U17 changes, 0 derive.
  Iso-Seq chr8: 1 novel mono, 11 known mono, 8 U17 changes (4 drops), 11 derive, 10 add-known.
- capture_b2.sh holds the full b2 table. Launched: 14 re-captures (b2_<M0a slice>) + b2_unspliced (ONT, PB), b2_ownership,
  b2_altends, b2_altends_napa (PB, --novel_apa) locally.
- NOTE: this interactive job has mem=8G cgroup; >5 concurrent local captures -> SIGKILL. PB b2_unspliced/b2_ownership were
  killed, cleaned and resubmitted to SLURM (57976939/40). Keep local N<=3.
- b2_vtg22sub done (SLURM, 7 min, 241 MB): giant block 208k reads, blocks 0-7, U3=10 U8=8 U11=8.
- b2_unspliced ONT: 11 novel mono, 2 known mono. b2_altends(+_napa): 7 derive, 7 add-known.
- PB b2_unspliced re-captured with region 8:28286000-28304534 (novel mono model lies left of its block region): 1 novel mono, 6 known mono.
- b2_ownership (SLURM): 6 U17 resolve changes, 4 drops.
- ALL DONE: 22 b2 dirs, iqcap.validate 0 errors each, total 363 MB.
- Rust check ../fixtures: 0 unreadable; all M0a dirs conform; 43142 non-conforming lines, all in b2 dirs, all due to
  U4 extract_relevant_deletions in.pos0, U13 in.intron_sites, U16 correct_transcript_splice_sites in.ref (DTO fields needed).
