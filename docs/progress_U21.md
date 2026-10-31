# U21 progress (scheduler / orchestration)

Build: `export CARGO_TARGET_DIR=$HOME/.cache/iqrs-target/U21`.
While `iq_collect` (U8, concurrent) did not compile, stage 1 was checked in a scratch crate that symlinks
`crates/iq_sched/src` and depends only on iq_types + iq_quant + rayon + crossbeam-channel
(`$SCRATCH/schedcheck`, target `~/.cache/iqrs-target/U21-scratch`): `cargo test --release`.

## Stage 1 (mock units): DONE
Files: `crates/iq_sched/src/{lib,units,admission,pipeline,engine,profile,mock,tests}.rs`.
- `units.rs`: legal units = maximal runs of consecutive blocks sharing ≥1 gene (§5.3). Default `Packed` mode packs
  consecutive small legal units up to 20 k reads / 512 blocks (never changes outputs). `IQRS_DEBUG_UNIT_SPLIT`
  = `whole-chr` | `max` (finest legal split, no packing) | `unsafe-per-block` (illegal; test teeth only).
- `pipeline.rs`: `Pipeline` trait = the stage interface (pass_a, collect_block, store_bytes, spill, barrier_chr,
  decide, construct_block, merge_chr, finish + estimates). Determinism contract documented on the trait.
- `engine.rs`: coordinator on the calling thread, tasks on one rayon pool (`in_place_scope`), results over a
  crossbeam channel; no shared mutable state. Queue classes: now (barrier/merge/spill, no memory check) >
  pass A (largest chr first) > work (collect/construct, LPT by cost, strict head-of-line byte admission).
  Running cap = threads; nested rayon inside tasks uses idle workers. Pass A of chr X overlaps collect of
  others. Barrier per chr when all its units are collected; global `decide` after every chr; then Phase M;
  merge per chr as soon as its units are constructed; `finish` after all, inside the pool.
  Spill: after a unit's collect, if `resident > budget/2` and `can_spill()`, a Spill task runs.
- `admission.rs`: `fits = running==0 || resident + in_flight + est <= budget`; budget default
  0.85 × min(cgroup memory.max over ancestors (v2/v1), MemAvailable).
- `profile.rs`: RSS sampler (/proc/self/statm, 100 ms), CPU from /proc/self/stat, per-kind task spans,
  collect est-vs-actual B/read, top-5 units. `Report::summary()` = `--profile` output.
- `mock.rs`: mock stages with all real couplings (detected_known per gene, chr-local duplicate names →
  suspended/ambiguous, sample-wide polyA decision, block-local novel ids renumbered per chr), busy-spin
  costs, touched ballast memory, nested rayon above `par_threshold`, optional spill to files.

Tests (13 pass, ~8 s release):
- (a) byte-identical outputs for threads 1/4/16 × Packed/Max/WholeChr, 2 seeds, incl. a 60 k-read nested block,
  an empty chr and a tight 2 MiB budget; `unsafe-per-block` differs (teeth).
- barrier: polyA 0.9 → requires_polya true, 0.5 → false; totals exclude suspended dups.
- (b) accounting: peak accounted ≤ budget, 0 over-budget admissions; a giant unit > budget runs alone.
  RSS (child process): 2 M-read giant + 300 units, stores 433 MiB > budget 400 MiB → spill; peak RSS 446 MiB
  incl. 4 MiB baseline (+6% allocator slack), peak accounted ≤ budget.
- (c) utilisation (child process, 4 threads, chr1 = 88% of reads incl. a 300 k-read block):
  packed 93% of 4 cores (wall 1.94 s) vs whole-chr (4.1 granularity) 36% (wall 5.03 s).

## Stage 2 (real units): wired end to end
Files: `crates/iq_sched/src/{real,construct}.rs`, `bins/iqrs/src/main.rs`. Cargo: added `iq_assign` (internal crate,
for `AssignParams` of the live assigner) to `crates/iq_sched/Cargo.toml` (+ its quirks41 feature).
- `real.rs` `RealPipeline: Pipeline`:
  - setup: `run_params` (preset or hidden `--params-json`, + `--report-novel-unspliced`), GFF3/GTF by extension
    with `GeneTypes::{All,Isoquant}`, FASTA mmap, 3 BAM handles (pass A: ≤4 inflate threads, no tag; collect: inline
    inflate; blocks ≥ 20 k reads: ≤4 inflate threads). Chromosomes = BAM ∩ FASTA ∩ annotation (`ann.chrom`), longest
    FASTA record first (4.1 `get_chromosome_ids`).
  - pass A: `try_fetch(seq=false)` over the chromosome → `iq_regions::PassA` with an `AnnotationIndex::genes_for_region0`
    oracle and `PassAConfig::for_chr` (coverage caps from params). Keeps per chr: block regions, cluster of each block,
    dup candidates.
  - collect: block index `index_from_genes(block.gene_ids)` / `index_from_region` → reference = (cluster ∪ index span)
    ±1000 bp, clipped (U8 note: U7 strand detection reads index introns outside the read span) →
    `iq_collect::BlockCollector` with `ClusterOwner{chromosome-wide regions, me}` and `SharedGroups` when tagged.
  - barrier: U9 `iq_collect::multimap::resolve_chr_duplicates` (removes suspended copies, re-indexes rows) → totals;
    `decide` = 4.1
    `set_polya_requirement_strategy` for requires_polya / monointronic / monoexonic (`dataset_processor.py:283-300`).
  - construct: `construct.rs` block driver in 4.1 order (U10b preds, U11 graph, U12 paths, U13 FL, U14, U15a pre_filter,
    U17 reassign, U15a filter with `Ends::default()`, U16 if `splice_site_correction`, U17 reassign + ownership,
    U18 join, U19 final → `TxCounter::add(fa, view.group())`). `UnitState` = `DetectedKnown`.
  - merge: `ChrPhaseW` per chr; finish: `write_outputs` (`__not_aligned` from BAI).
  - `store_bytes` = `BlockStore::len_bytes()`; `fixed_bytes` = annotation heap + 16 MiB/thread.
- CLI (`bins/iqrs`): design §7.0 flags + `--read-group tag:X` (and `--read-group-tag X`), `--profile`, `--verbose`,
  `IQRS_DEBUG_UNIT_SPLIT`, GTF header command line = argv.
- Packing target scales: clamp(chr_reads / (32 × threads), 200, 5 k) reads (`units::pack_target`); on Revio chr8 the
  finest split was no slower than coarse packing, so packing is kept light.
- Read-group ids (SharedGroups) are first-seen order (scheduling-dependent); U20 grouped rows are written in
  first-touch read order with names resolved through the pool, so the output does not depend on ids.

### SLURM benchmark (job 57976530, `-c 16 --mem=32G`, inputs on node-local tmp, default build, 16 threads)
Script: `../iqrs_bench_U21/bench.sbatch`, log `../iqrs_bench_U21/bench.57976530.log`.
| input | iqrs wall / max RSS | 4.1 wall / max RSS |
|---|---|---|
| SRR32588726 chr8 (1.6 M reads) | 5.97 s / 1.05 GB | 17:17 / 2.45 GB |
| SRR32588716 chr8 (0.72 M) | 3.27 s / 0.60 GB | (see reference_outputs) |
| vtg.bam (chr22 slice, 17.7 M reads; units: 15.26 M-read unit + 2.47 M-read unit) | 1:58 / 4.26 GB | 4.1: ~138 GB (design §5.7) |
- vtg critical path: pass A of chr22 is serial 30 s (vtg_proto predicted ~33 s), then collect of the two giant units
  runs concurrently (15.26 M reads in 59 s; the 2.47 M-read unit oddly needed 79.5 s, i.e. ~31 k reads/s even after
  the big unit finished; likely the serial part of `BlockCollector::push` (fetch/decode/select/append between
  parallel batches) — an Amdahl limit inside U8, not the scheduler). Construct of both: 7 s.
- Whole-run utilisation 32-39%: pass A serial per chromosome + serial parts of the collector. Next steps:
  stream closed blocks out of pass A into the queue (needs an iq_regions callback API) and pipeline fetch/decode with
  per-read work inside `BlockCollector` (U8).
- `--max-mem 4` on vtg: same wall (1:58), peak RSS 4.07 GiB, budget after fixed bytes 3.71 GiB, peak accounted
  3.76 GiB, 1 over-budget admission (progress rule). The real pipeline has no BlockStore spill yet, so the floor is
  Σ resident stores (2.85 + 0.53 GB here); spill (`Pipeline::spill`, columnar file + mmap, design §5.4) is the
  remaining memory lever (mock path implemented and tested).
- Admission: vtg unit 1 est 2.33 GB vs actual store 2.85 GB (≈187 B/read); recalibrate `EST_COLLECT_BYTES_PER_READ`
  to ~220. Done: `EST_COLLECT_BYTES_PER_READ = 224`.

## Open items / next steps
1. BlockStore spill for the real pipeline (needs an iq_types columnar write/mmap-read of BlockStore, or a U21-private
   format over `ReadView` → `ReadRecord` round-trip).
2. Stream closed blocks out of pass A (iq_regions callback) so collect of chr X starts before pass A of X ends
   (vtg: 30 s serial).
3. U8 collector Amdahl limit on giant blocks (2.47 M-read vtg block at ~31 k reads/s).
4. PacBio quirks41 residual diffs (TSS/terminal) — M5 triage.
5. Pass criterion (c) end-to-end compare: `iqrs-compare` is still a skeleton (M5); used a scratch script.

### E2E parity, quirks41 build (vs real 4.1 chr8 outputs; IDs normalised: transcript*/novel_gene*/exon_id)
| sample | GTF | counts |
|---|---|---|
| SRR33046254 ONT | identical | identical (incl. __ambiguous/__no_feature) |
| SRR33046249 ONT | identical | identical |
| SRR32588716 Revio | 70 of 28,104 lines differ (a few models' terminal positions, e.g. start 2456876 vs 2456872) | identical |
| SRR32588720 Iso-Seq | 126 of 21,593 lines | 4 lines |
| SRR32588726 Revio | 152 of 32,693 lines | 13 lines (__ambiguous 45981 vs 45982, __no_feature 68214 vs 68233) |
PacBio-only residuals point at TSS/terminal handling (U10b/U15b, `use_tss_model` is pacbio-only) or nondeterminism
floor (Iso-Seq noise floor 1/1375 transcripts); for M5 to triage. ONT is exact.
Barrier totals equal 4.1's log line (e.g. ONT 254: 70,975 / 7,909).

### Real-data results (default build, `--gene-types isoquant --report-novel-unspliced false`, login node, load ~90)
| sample (chr8) | iqrs wall / RSS | 4.1 wall / RSS (16 thr) | models iqrs/4.1 | read-weighted chain share ref→test / test→ref |
|---|---|---|---|---|
| SRR33046254 ONT (71 k) | 2.4 s / 0.16 GB (8 thr) | 5:33 / 0.44 GB | 449 / 452 | 0.9985 / 1.0000 |
| SRR32588726 Revio (1.6 M) | 12.3 s / 0.76 GB (8 thr) | 17:17 / 2.45 GB | 3128 / 3085 | 0.9999 / 0.9953 |
- Barrier totals: 1,596,312 assignments / 1,392 polyA vs 4.1 1,596,326 / 1,334 (requires_polya false in both).
- Determinism on real data: SRR33046254 outputs byte-identical (minus header line 2) for threads 1/4/8, whole-chr, max.
- Default build `__no_feature` differs (ONT 14442 vs 13471): D-R2 (no valley double processing) — under quirks41 it
  matches exactly.
- Real-data determinism also on SRR32588726 (threads 1/4/8, whole-chr, max: byte-identical).
- Timings on a node with load ~90 are noisy; proper benchmarks (SLURM, vtg.bam, full BAMs) still to do.
- Profile notes: collect is CPU-bound (~31 µs/read on Revio chr8 under load; U8 per-read cost), giant block
  (284 k reads) is not the critical path; pass A of the single data chromosome (~4 s) is serial before collect.
  Collect resident ≈ 217-249 B/read (est constant 160 is low; recalibrate after U9/real assignments).

## Interface for U8/U9 (sent to main 2026-10-29; both now implemented by the collect agent and wired)
**Collect (exists in iq_collect; keep stable):** `BlockCollector::new(&Block, &IsoformIndex, &RefSlice, &Params,
&dyn Owner, &dyn GroupIds)`, `.push(&[RecordLite])` per fetched batch (BAM order), `.finish() -> (BlockStore,
CollectStats)`; `ClusterOwner{regions, me}` (U21 passes all block regions of the chromosome, sorted, disjoint);
`SharedGroups` (pub `.0` RwLock; U21 clones the pool at finish). U21 supplies the block index, a RefSlice over the
block's cluster ±1000 bp (clipped), records from `try_fetch(block.region, seq=true)`. Blocks may be collected
concurrently; the collector may use nested rayon.

**U9 (needed):**
```rust
pub fn resolve_chr_duplicates(stores: &mut [BlockStore], dup_candidates: &[u64], strategy: MultimapStrategy) -> DupOutcome;
pub struct DupOutcome { pub total: u64, pub polya: u64, pub groups: u64, pub suspended: u64 }
```
- `stores` = all blocks of ONE chromosome in genomic order (4.1 NoSecondary is chromosome-local);
  `dup_candidates` = `PassAResult.dup_candidates` (sorted superset of name hashes).
- Confirm groups by name; resolve groups of ≥2 copies with `multimap::resolve`; patch via `set_types`/`set_multimapper`.
- `total`/`polya` = 4.1 `total_assignments`/`polya_assignments` of the chromosome after resolution
  (`processed_read_manager.py:20-87`: suspended copies excluded).
- Groups and copies in (block,row) order (never hash order).
- Document how suspended reads look afterwards and that U10b-U19 skip them as 4.1's loader does
  (`assignment_loader.py:50-70`).
- Call site: `real.rs::RealPipeline::barrier_chr` → replace the body of `real::resolve_duplicates` (signature
  there takes `&mut [&mut BlockStore]`; adapting is trivial). Until then: totals = Σstored/Σpolya_found, confirmed
  unresolved duplicate names are counted and warned about.
