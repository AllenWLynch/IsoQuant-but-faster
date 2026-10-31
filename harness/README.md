# iqrs harness: `iqcap`, the golden-fixture capture harness for IsoQuant 4.1.0

`iqcap` records the inputs and outputs of IsoQuant 4.1.0 functions, at the hook points of design doc
07 §10.2, while IsoQuant runs unmodified. Monkeypatches are installed only at runtime; the IsoQuant
repo is never edited. The resulting fixtures are the contract for the Rust units U1-U20. Their exact
format is in [`iqcap/SCHEMA.md`](iqcap/SCHEMA.md) (schema `v = 1`).

```
harness/
  iqcap/
    SCHEMA.md            the fixture contract (read this first)
    run.py               python -m iqcap.run ... -- <isoquant args>      (capture)
    hooks.py             all recording wrappers, one _install_uN() per unit
    enc.py               object -> JSON encoders (ReadRecord, Assignment, Model, Store, index dump)
    sink.py              shared dumps(), deterministic sampling, per-PID gzip JSONL sinks, index side table
    slice.py             python -m iqcap.slice     BAM region slice + deterministic HP:Z synthesis (+ subsample)
    dump_bam.py          python -m iqcap.dump_bam  U2: pysam walk of a region
    dump_annotation.py   python -m iqcap.dump_annotation  U1(a): gffutils DB view, no IsoQuant run
    find_slices.py       python -m iqcap.find_slices  proposes slice regions (approximate pass A)
    verify.py            python -m iqcap.verify A B   compares IsoQuant outputs (no-perturbation check)
    validate.py          python -m iqcap.validate DIR structural/round-trip checks of a fixture dir
    report.py            python -m iqcap.report       coverage table (cases per unit per sample)
  capture_slice.sh       slice + capture driver (local or sbatch)
  capture_all.sh         the M0a fixture set (all samples/slices)
```

All commands run from `harness/` with the project env: `../../env/bin/python -m iqcap.<tool>`.

## Capturing

```bash
# 1. slice a BAM (adds HP:Z = A|B|unphased from sha1(qname), so --read_group tag:HP is exercised)
python -m iqcap.slice --bam ../../data/chr8/SRR33046254.chr8.bam --regions 8:21104496-21142297 \
    --out /path/slice.bam [--subsample 0.15]
# 2. capture (always --threads 1; stage genome.fa + the gene DB on local disk first)
python -m iqcap.run --units all --out ../../fixtures/SRR33046254/valley --sample SRR33046254 --slice valley \
    --sample-rate "U4=0.05,U5=0.02,default=1" --block-rate 1.0 --block-cap 200000 \
    --slice-def /path/slice.slice.json -- \
    --reference $STAGE/genome.fa --genedb $STAGE/isoquant.annotation.db --bam /path/slice.bam \
    --data_type nanopore --model_construction_strategy sensitive_ont --report_novel_unspliced false \
    --read_group tag:HP --threads 1 --force -o /tmp/iqout
# or both steps:
STAGE=/local/stage ./capture_slice.sh SAMPLE SLICE REGIONS nanopore|pacbio_ccs [RATES] [SUBSAMPLE] [UNITS]
sbatch -p park,short -t 4:00:00 -c 2 --mem=16G capture_slice.sh ...   # long ones; stages to $TMPDIR
# run-independent dumps
python -m iqcap.dump_annotation --db .../isoquant.annotation.db --gff3 .../annotation.gff3 --chr 8 \
    --out ../../fixtures/annotation/U1.annotation.jsonl.gz
python -m iqcap.dump_bam --bam slice.bam --region 8:1000-2000 --out ../../fixtures/S/ordinary --sample S
python -m iqcap.report            # coverage table
```

PacBio runs use `--data_type pacbio_ccs --model_construction_strategy sensitive_pacbio --fl_data`.
ONT runs use `--data_type nanopore --model_construction_strategy sensitive_ont`.

Sampling options:
- `--sample-rate` sets the per-unit rate for per-read or per-call cases (U3, U4-U7, U9, U10a, U20 per-read, and
  the nested U5/U6 calls of construction).
- `--block-rate` / `--block-cap` select whole blocks for the construct-side units U8 and U10b-U19, so every
  cross reference within a selected block resolves.
- All sampling is a deterministic hash of the case key.

### Capture does not change IsoQuant's behaviour
- Wrappers encode their inputs (a deep copy into JSON-ready structures) before calling the original.
  They never alter arguments or return values. Their only side effects are private `_iq_*` attributes (read index, cached index hash).
- Side computations (the XGBoost `predict_proba`/margin and the U10a feature recomputation) call pure
  functions under a reentrancy guard.
- Check: `NO_CAPTURE=1 VERIFY_DIR=... ./capture_slice.sh ...` and a capture run with the same `VERIFY_DIR`, then
  `python -m iqcap.verify $VERIFY_DIR/<slice>/plain $VERIFY_DIR/<slice>/capture`. Every `.gtf`/`.tsv`/`.tsv.gz` output
  is compared, ignoring only the `# Command line:` header. Results are in "Verification" below.

## Adding a unit hook

1. Find the function and the state you need. Read the IsoQuant source, and specs 01-05 for context.
2. In `hooks.py`, add (or extend) an `_install_uN()` and call it from `install()`. The pattern is:
   ```python
   def f(orig):
       def w(self, *args):
           if not CTX.selected("UN", CTX.bkey()):        # or a per-read key; see Ctx.selected
               return orig(self, *args)
           inp = {"fn": "name", ...encode inputs with enc.* BEFORE the call...}
           res = orig(self, *args)
           _write("UN", CTX.local("UN", "name", "n"), inp, {...encode outputs...})
           return res
       return w
   patch(SomeClass, "name", f)
   ```
   - For module-level functions, rebind the attribute in every module that imported it by name
     (see `extract_relevant_deletions`).
   - For static methods, re-wrap with `staticmethod`.
   - Read identity: `CTX.read_ident()` (collect: `{qname, aln_no}`) or `enc.ridx(read)` (construct: the
     `_iq_idx` set in `ReadAssignmentLoader.get_next`).
3. Encode only with `enc.*` / `sink.fval`, and write only through `Sinks.write` (the one `dumps`).
   Pair lists keep dict order. Sets are sorted. Enums become `.value`.
4. Document the payload in `SCHEMA.md` **before** shipping: the Rust round-trip test fails on undocumented
   fields. A breaking change bumps `SCHEMA_VERSION` in `iqcap/__init__.py` and in SCHEMA.md. Tell the coordinator.
5. Re-run the no-perturbation check on one slice.

## Fixture set (M0a)

Location: `/n/data1/hms/dbmi/park/allen_l/lariat/isoquant_r/fixtures/<sample>/<slice>/`. The slice table, with
regions, rates and subsampling, is in `capture_all.sh`. Every slice dir holds its slice BAM (`slice/`),
`fixture_manifest.json`, `params.json`, `units/`, `index/`, `block/` and the IsoQuant final outputs (`outputs/`).

| slice | design 10.2 | ONT SRR33046254 | Iso-Seq SRR32588720 | Revio SRR32588726 |
|---|---|---|---|---|
| ordinary | (a) 5 ordinary loci (multi-isoform genes) | 5 genes, 1.3k reads | 5 genes, 1.5k reads | 5 genes, 6.8k reads |
| deepest | (b) deepest chr8 cluster/block | 8:30.70-30.83M, 6.0k reads | 8:47.10-47.21M, 8.1k reads | deepest block 8:13.75-13.97M (285k reads), subsampled to 0.15 -> 42k-read block |
| multigene | (c) multi-gene both-strand cluster | 15 genes, 3.9k reads | 24 genes, 3.5k reads | 14 genes, 9.7k reads |
| nogene | (d) annotation-free regions | 3 intergenic clusters | 3 | 3 |
| valley | (e) coverage-valley split | 2 blocks from 1 cluster | 2 | 2 |
| fullchr8 | (f) full chr8, U8 + U20 only | 900 blocks, 71k reads | 708 blocks, 114k reads | not run (1.7M reads; not cheap) |

All captures add `--read_group tag:HP` (HP synthesized by `iqcap.slice`) and run `--threads 1`. Sample
rates: small slices capture 100% of everything. Revio ordinary/multigene/valley and ONT deepest use
`U4=0.25,U5=0.1,U6=0.5,U7=0.25,U20=0.25`. Revio deepest uses `U4=0.05,U5=0.02,U6=0.25,U7=0.05,U20=0.1`. Construct-side
units (U8, U10b-U19, U20 gff_dump, `block/`) always cover every block of every slice. U2: one region of each
sample's ordinary slice in `fixtures/<sample>/u2/`. U1(a): `fixtures/annotation/U1.annotation.jsonl.gz`
(all 25 chromosomes of the gene DB).

### Coverage report (cases per unit per sample; `python -m iqcap.report [--by-slice]`)

| unit | SRR32588720 | SRR32588726 | SRR33046254 |
|---|---|---|---|
| U1 | 306 | 549 | 125 |
| U2 | 484 | 1161 | 414 |
| U3 | 53 | 69 | 71 |
| U4 | 19747 | 13615 | 22321 |
| U5 | 20009 | 7489 | 17939 |
| U6 | 17992 | 33285 | 17668 |
| U7 | 26895 | 20187 | 31075 |
| U8 | 725 | 31 | 919 |
| U9 | 37 | 4 | 40 |
| U10a | 99 | 80 | 45 |
| U10b | 44 | 72 | 24 |
| U11 | 17 | 31 | 19 |
| U12 | 34 | 62 | 38 |
| U13 | 17 | 31 | 19 |
| U14 | 17 | 31 | 19 |
| U15a | 68 | 124 | 76 |
| U15b | 159 | 328 | 69 |
| U16 | 63 | 153 | 23 |
| U17 | 51 | 93 | 57 |
| U18 | 17 | 31 | 19 |
| U19 | 17 | 31 | 19 |
| U20 | 62661 | 35885 | 59853 |

- SRR32588720: slices deepest, fullchr8, multigene, nogene, ordinary, valley; size 131.6 MB
- SRR32588726: slices deepest, multigene, nogene, ordinary, valley; size 57.3 MB
- SRR33046254: slices deepest, fullchr8, multigene, nogene, ordinary, valley; size 62.5 MB
- U1 annotation dump: 25 chromosome case(s), /n/data1/hms/dbmi/park/allen_l/lariat/isoquant_r/fixtures/annotation/U1.annotation.jsonl.gz (6.9 MB)
- total fixtures size: 258.4 MB

U1 counts are the (b)/(c) cases. The annotation dump is reported separately.

### Verification
- **No perturbation.** Capture runs were checked against plain IsoQuant runs on the same slice:
  - `SRR33046254/valley`: all 26 output files identical (`iqcap.verify`, ignoring the `# Command line` header).
  - `SRR32588720/valley` (PacBio, `--fl_data`): the 9 copied outputs identical.
  - The full-chr8 captures (`--threads 1`, HP tags, `--read_group tag:HP`) reproduce the reference runs in
    `reference_outputs/{ont_minion_SRR33046254,pb_isoseq_SRR32588720}` (16 threads, original BAMs) byte for byte
    for `transcript_models.gtf` and `discovered_transcript_counts.tsv`.
- **Structure.** `python -m iqcap.validate <slice dir>` checks every slice. All are 0 errors. It checks the envelope, `v`,
  unique keys, the byte round trip `dumps(loads(line)) == line`, index hashes, block files and graph/paths cross references.
- **Annotation consistency.** Transcript exon lists in `index/` dumps equal those in the U1(a) annotation dump.


## Fixture set batch 2 (`b2_*` dirs)

Driver: `capture_b2.sh local|sbatch [regex]` (table inside; per-entry env overrides for `capture_slice.sh`:
`NOVEL_UNSPLICED`, `EXTRA_IQARGS`, `SRC_BAM`, `KEEP_NAMES`, `BLOCK_CAP`). The M0a dirs are never touched.
Harness changes (additive, SCHEMA.md changelog): U4 `extract_relevant_deletions.in.pos0`; U5 `match_genomic_features.in.delta`
= comparator's effective delta (bug fix, b2 dirs only); U13 `in.intron_sites`; U16 store-level `in.ref`;
`iqcap.slice --keep-names`.

| dir | sample | purpose |
|---|---|---|
| `b2_{ordinary,deepest,multigene,nogene,valley}` | all 3 | M0a slice definitions/rates re-captured (all units) so the new/fixed U4/U5/U13/U16 fields have coverage |
| `b2_unspliced` | ONT 254, Iso-Seq 720 | `--report_novel_unspliced true`; loci with novel mono-exon models and known mono-exon models (gaps a, b) |
| `b2_polyarich` | ONT 254 | read-name-restricted slice (all polyA reads of M0a valley/ordinary/multigene + 1/4 as many non-polyA) -> sample polyA fraction 80% >= 0.7 -> `requires_polya_for_construction = true` (gap c); `--report_novel_unspliced true` |
| `b2_ownership` | Iso-Seq 720 | blocks where `resolve_read_ownership_conflicts` changes/drops models (gap d) |
| `b2_altends`, `b2_altends_napa` | Iso-Seq 720 | blocks where known-model alt-end derivation emits models (gap e); `_napa` adds `--novel_apa` |
| `b2_vtg22sub` | SRR32588716 (vtg.bam) | chr22 vtg locus subsampled to 2%, U3/U8/U11 only, `--block-cap 1000000` (gap f) |

Loci for a, b, d, e were chosen from full-chr8 diagnostic captures (`--units U3,U13,U14,U15b,U17`,
`--report_novel_unspliced true`) by comparing store snapshots before/after each stage.
`requires_polya_for_construction` is set in `dataset_processor.py:290` from `polya_found / total_assignments >=
args.polya_percentage_threshold (0.7)` under `--polya_requirement auto` (or forced by `always`); none of the real
samples reaches it (ONT 8-23%, PacBio ~0%).
