# M1 vtg memory prototype: report

Status: **done** (2026-10-29). Tool: `crates/iq_sched/examples/vtg_proto.rs`.
Build: `cargo build --release -p iq_sched --example vtg_proto`. This adds `[dev-dependencies] iq_read` to
`crates/iq_sched/Cargo.toml`, which main approved.
Run: `vtg_proto passa|store <bam> <chr> [inflate_threads]`.
SLURM job 57965975 (`-p park,short -c 8 --mem=32G`, BAMs copied to node-local `/tmp`) ran on compute-p-17-175.
Its log is `target/vtg_proto_57965975.log`. That file is under target/, so it may be cleaned; the key numbers are copied below.

## What the tool does
- **Pass A.** `iq_bam::fetch_with(seq=false)` over the whole chromosome, fed into `iq_regions::PassA` with `NoGenes` and the 4.1 coverage defaults (`--max_coverage_normal_chr -1`, so no cap).
- **Store.** Pass A, then one more streaming fetch with sequences. Each record goes through these steps:
  - 4.1 collect filter: drop unmapped, supplementary and secondary reads (`use_secondary` is off).
  - `iq_read::process_alignment_into` (CIGAR → exons, polyA/T, boundary deletions; `Params::pacbio()`).
  - One reused `ReadRecord` with a stub assignment: `Assignment::empty(Noninformative)`, `corrected_exons = exons`, `strand = mapped_strand`.
  - `owner_block` (D-R2), then `BlockStore::push` into that block's store.
- **Resident set.** A block's store is finalized and dropped as soon as the coordinate-sorted stream passes its end. Only open blocks are resident, and all `ReadRecord`s are never held at once.
- **RSS** is sampled from `/proc/self/statm` every 100 ms in a thread. VmHWM and `/usr/bin/time -v` are reported as cross-checks.
- **Accounting per block:**
  - **Hot columns.** Length-based bytes of the v1 `BlockStore` columns. That is 57 B fixed: name_hash 8, name_off 4, chain 4, start/end 8, raw 4, polya 16, flags 2, mapq 1, group 4, atype ×2 2, match_off 4, dels_off 4. Add 8 B per boundary deletion.
  - **ChainTable.** 8 B per stored intron, plus about 34 B per chain (offs 4, pure 1, FxHashMap entry about 29).
  - **Cold names.** QNAME bytes.
  - **`heap_bytes()`.** Capacity-based, as implemented in iq_types.
  - **Interning.** Distinct chains, singleton chains, and the share of reads covered by the 100 most frequent chains.

## 1. Pass A over vtg.bam (chr22:24.46-25.05 Mb, 17,775,376 records)
| inflate threads | wall | records/s | peak RSS |
|---|---|---|---|
| 4 | 32.5 s | 547 k | 85 MB (pass A internal 69.6 MB, mostly the dup filter) |
| 1 | 114.8 s | 155 k | 82 MB |

- The slice has 17.73 M primary, 42 k secondary and 5 k supplementary records.
- It forms **1 cluster with 7 blocks** after valley splitting:

| block (0-based) | len | max_cov | owned primary |
|---|---|---|---|
| 24462403-24495104 | 32.7 kb | 1 | 0 |
| 24495105-24527872 | 32.8 kb | 1 | 0 |
| **24527873-24738304** | 210 kb | 10,481,086 | **10,438,643** |
| 24738305-24858112 | 120 kb | 1,658,964 | 4,823,111 |
| 24858113-24890880 | 32.8 kb | 2,470,932 | 2,465,021 |
| 24890881-24923648 | 32.8 kb | 963 | 1,087 |
| 24923649-25048722 | 125 kb | 8 | 4 |

- **The 10.4 M-read block stays one block.** Valley splitting does not break it, which matches 4.1's 10.4 M.
- Two more multi-million-read blocks sit next to it in the same cluster.
- Gene bridging can only merge clusters. Here there is already a single cluster, so the `NoGenes` stub does not change the blocks.
- With annotation, the vtg genes almost certainly make this one work unit of about 17.7 M reads (§5.3).
- **Throughput is inflate-bound.** The slice is 11.4 GB for 17.8 M records, about 640 B compressed per record. U3's clustering alone runs at 8.8 M rec/s.
- chr22 is a serial pass-A critical path of about 33 s (4 threads) before collect of the giant unit can start, and pass B inflates the same 11 GB again.

## 2. BlockStore for the vtg blocks (stub assignment)
Store build over the whole slice: 48.8 s wall, 363 k stored reads/s. That uses 4 inflate threads and 1 consumer thread (U4 only, not parallelized).

| | giant block | block 2 | block 3 | all vtg blocks |
|---|---|---|---|---|
| reads | 10,438,643 | 4,823,111 | 2,465,021 | 17,727,866 |
| introns/read | 15.65 | 15.85 | 16.97 | 15.89 |
| distinct chains | **12,868** | 18,628 | 4,963 | 36,729 |
| reads/chain (interning ratio) | **811** | 259 | 497 | 483 |
| singleton chains (% of reads) | 0.06% | 0.17% | 0.09% | 0.10% |
| top-100 chains cover | 99.6% | 97.5% | 99.4% | 99.0% |
| hot columns B/read | 62.6 | 62.7 | 62.7 | 62.7 |
| ChainTable B/read | **0.2** (2.25 MB total) | 0.7 | 0.4 | 0.4 |
| cold names B/read (+4 B offset in hot) | 19.8 | 19.8 | 19.8 | 19.8 |
| `heap_bytes()` B/read (capacity) | 133.7 | 143.2 | 139.9 | 137.2 |
| chains without interning | 129 B/read | 131 | 140 | 131 |

- **Peak RSS of the whole store run: 1,068 MB.** VmHWM was 1,065 MB and `time -v` reported 1,040 MB; baseline RSS before the build was 36 MB.
- The peak falls just before the giant block is retired. That gives (1068 − 36) MB / 10.44 M = **≈ 99 B/read actually resident**:
  - 63 B hot columns + chains;
  - 20 B QNAMEs;
  - the rest is touched Vec-growth slack and allocator overhead.
- `heap_bytes()` (134 B/read) overstates resident memory. Doubled Vec capacity is mostly untouched pages, and large reallocations use mremap.
- **Interning works far better than §5.7 assumed.** There are 12.9 k chains against the ≤10⁵ assumption, and the chain table costs 2.25 MB.
- Without interning, the chains would cost 129 B/read, which would double the store.

## 3. Typical cases (chr8, whole chromosome)
| sample | records | pass A | stored reads | blocks | distinct chains | reads/chain | singletons | hot + chains B/read | names B/read | peak RSS (store run) |
|---|---|---|---|---|---|---|---|---|---|---|
| **SRR32588726** Revio | 1.69 M | 1.8 s, 939 k rec/s | 1,599,885 | 373 (319 with reads) | 33,478 | 47.8 | 1.2% | 62.0 + 1.7 = **63.7** | 19.8 | 63 MB |
| SRR32588716 Revio (vtg sample) | 0.76 M | 0.76 s, 1.0 M rec/s | 724,616 | 853 | 12,080 | 60.0 | 0.5% | 62.0 + 1.5 = 63.4 | 19.8 | 40 MB |
| SRR32588720 Iso-Seq | 0.12 M | 0.18 s | 113,870 | 789 | 9,559 | 11.9 | 5.6% | 62.3 + 8.7 = 71.0 | 18.7 | 28 MB |
| SRR33046249 ONT MinION | 0.21 M | 0.43 s | 187,436 | 1306 | 11,608 | 16.1 | 4.5% | 70.2 + 4.0 = 74.2 | 27.6 | 25 MB |
| SRR33046254 ONT MinION | 0.07 M | 0.08 s | 71,010 | 1423 | 7,872 | 9.0 | 7.3% | 67.2 + 6.4 = 73.6 | 27.1 | 12 MB |

- **SRR32588726's largest block** has 283,897 reads, 912 chains, and 62.5 B/read.
- **ONT and Iso-Seq rows** were run locally, not under SLURM, because they are small. They use PacBio polyA parameters, which does not affect the memory numbers.
- **ONT hot bytes are higher** (70 vs 62 B) because of 1.1-1.5 boundary deletions per read, against 0.12-0.2 for PacBio, at 8 B each.

## 4. Extrapolation: full SRR32588716 at 16 threads
**Inputs:**
- idxstats gives 41.97 M mapped records, about 41.5 M primary.
- vtg unit: 17.7 M reads in 3 large blocks plus 4 small ones. The giant block has 10.44 M.
- Measured resident cost with the stub assignment is about 99 B/read including names. Adding the design's ~40 B/read for real matches/events gives **≈ 140 B/read** (≈ 115 B/read once names are spilled).

**Projected components:**
- **Annotation and index:** about 0.1-0.3 GB (estimate; U1 not measured).
- **The vtg unit's stores**, built in Phase C and resident through Phase M: 17.7 M × 140 B ≈ **2.5 GB**. With the stub assignment the measured figure would be 1.75 GB; the tool held only one block at a time, at 1.07 GB.
- **Construct-side per-read arrays** for the giant block: 10.44 M × 20-48 B ≈ 0.2-0.5 GB. §5.7 assumes 20 B; `est_M` uses 48 B.
- **The other ~24 M reads' stores:**
  - kept resident: about 3.4 GB;
  - spilled under §5.4: mmapped and reclaimable, so close to 0 GB of anonymous memory.
  - On chr8 the largest non-vtg block is 146 k reads, so the 15 other in-flight units cost < 0.5 GB.
- **Batch buffers:** 16 threads × 4096 records at about 2-3 KB decoded (seq at 1 B/base, no quals) ≈ 0.2-0.4 GB. §5.7 assumed 1.3 GB.

**Projected peak:**

| scenario | peak | reduction vs 4.1's 138 GB |
|---|---|---|
| no spill | **≈ 6.5-7.5 GB** | ≈ 19-21× |
| with spill | **≈ 3.5-4.5 GB** | ≈ 30-40× |

The design's 3-7 GB claim holds, and only the no-spill case touches its upper end.

**Typical full SRR32588726** (53.9 M mapped records): about 51 M × 140 B ≈ 7 GB of stores unspilled, against 5.5 GB in §5.7, plus small in-flight units. It is still about 8× under 4.1's 58 GB, and spill brings it to ≈ 2-3 GB.

## Design assumptions checked
1. **Bytes/read holds.**
   - Measured fixed cost is 62-63 B/read (PacBio) and 67-73 B/read (ONT), against the 60-70 B claim.
   - Real assignments add about 40 B/read, so 100-115 B/read hot.
   - **Gap:** the §2.3 table omits QNAMEs. That is 20-28 B/read plus a 4 B offset, and BlockStore v1 stores them inline (`names: String`) instead of in a spillable cold arena (§2.3.2). Until U21 splits them out they add about 20% to every resident store.
2. **Interning holds, and better than assumed.**
   - vtg giant block: 811 reads/chain, 12.9 k chains (assumed ≤ 10⁵), 0.2 B/read.
   - Revio chr8: 48-60 reads/chain overall, with singleton chains ≤ 1.2% of reads.
   - **Noisy chains do not defeat interning on these data:**
     - ONT MinION still gets 9-16 reads/chain, with singletons at 4.5-7.3% of reads and the chain table at 4-6 B/read.
     - Iso-Seq gets 12 reads/chain.
   - **Caveat:** these are uncorrected exons. With annotation, U7 exon correction makes `corrected_exons` intern better, but every corrected read whose raw exons differ adds a raw chain (`raw_side`, 12 B, plus an interned raw chain). For ONT that raw-chain table could approach one chain per corrected read (about 8 B/intron); measure in M3/M5.
3. **"~34 exons" for vtg is wrong.** Measured 16.9 exons per read (15.9 introns). This is irrelevant to memory because of interning.
4. **"Boundary deletions rare" is not true.**
   - PacBio: 0.12-0.2 per read (1-2 B/read).
   - ONT: 1.1-1.5 per read (9-12 B/read).
   - Small, but include it in `est_C`.
5. **`BlockStore::heap_bytes()` reports capacity**, about 1.35-2× resident here. For admission calibration use length-based sizes, or `shrink_to_fit` the columns after collect. The `est_C` value of 160 B/read is conservative against the measured ≈ 140 B/read and is fine.
6. **Batch buffers** are probably overestimated in §5.7 (1.3 GB): they should be about 0.3 GB. Measure in M6.
7. **Pass A is inflate-bound, not clustering-bound.** It runs at 155 k rec/s on 1 thread and 547 k rec/s on 4 inflate threads for long HiFi reads. chr22 of SRR32588716 alone is about 33 s of serial pass A before the giant unit's collect can start, and the giant unit is the critical path. Consider giving pass A of the heaviest chromosome more inflate threads, or fusing pass A with BGZF reading for pass B (a later optimization).
8. **Collect throughput.** With a single consumer thread, U4 plus store append runs at 363 k reads/s, so about 29 s for the giant block. §5.7's rayon per-read pipeline is still needed once U5-U7 are added.

## Not measured / caveats
- Real U5/U6/U7 per-read work and assignment rows; the +40 B/read is the design's figure, not measured.
- Construct-side (Phase M) memory.
- The annotation index.
- Full-genome runs (M6).
- `cargo clippy` on the example was blocked at the time of writing because another agent's in-progress `iq_model` edits did not compile. The release build of the example compiled cleanly before that.

## Reproduce
`/tmp/.../scratchpad/vtg_proto.sbatch` (scratch copy), in short:
```bash
cp target/release/examples/vtg_proto $TMPDIR/
for s in SRR32588726 SRR32588716; do cp ../data/chr8/$s.chr8.bam* $TMPDIR/; /usr/bin/time -v $TMPDIR/vtg_proto store $TMPDIR/$s.chr8.bam 8 4; done
cp -L ../data/vtg.bam ../data/vtg.bam.bai $TMPDIR/
/usr/bin/time -v $TMPDIR/vtg_proto passa $TMPDIR/vtg.bam 22 4
/usr/bin/time -v $TMPDIR/vtg_proto passa $TMPDIR/vtg.bam 22 1
/usr/bin/time -v $TMPDIR/vtg_proto store $TMPDIR/vtg.bam 22 4
```
Reference names in these BAMs are `1..25`, without a `chr` prefix.
