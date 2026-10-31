# U2 progress (iq_bam)

## Status (2026-10-29)
- bam.rs: own BAM record decoder over noodles-bgzf (single-thread Reader or MultithreadedReader) +
  noodles BAI; htslib `hts_itr_next` region semantics (reg2bins + linear-index min_off,
  optimize_chunks, stop at other tid / pos >= end, keep iff bam_endpos > beg). CG:B:I long CIGAR
  handled. Flat RecordBatch (names/cigars/seqs/meta arrays, tag values interned per batch).
- fasta.rs: mmap + fai arithmetic, Cow borrowed within a line. Identical to pyfaidx on 450 random
  slices (chr 8/1/22, incl. clipped past-end).
- htslib.rs: `HtsBam` (rust-htslib IndexedReader) filling the same RecordBatch; test proves identical
  batches to noodles on the edge BAM.
- Fixtures: U2/fetch 2059/2059 exact, parity + default + plain default build.
- tests/gen_edge_bam.py (pysam) → tests/data/edge.bam + a U2 fixture file (19092 cases): IUPAC,
  every HP type (i,C,c,s,S,I,Z,A,f,B), '*' seq, missing qual, placed unmapped, zero-rlen CIGAR,
  70k-op CIGAR (CG), 150 random regions. All exact; record counts per region equal.
- Real data cross-check (examples/bam_bench check vs bam_check.py, CRC over qname/pos/flag/end/seq/HP):
  identical on SRR33046254 chr8 (300 regions) and SRR32588726 slice; big BAM in the sbatch job.

## rust-htslib build recipe (works on O2)
  LIBCLANG_PATH=<repo>/env/lib  BINDGEN_EXTRA_CLANG_ARGS="-I/usr/lib/gcc/x86_64-redhat-linux/11/include"
  cargo build --release -p iq_bam --features htslib

## Benchmark (DONE, job 57958639, Xeon Platinum 8468, 8 cpus, BAM in node-local $TMPDIR, page cache warm)
Script: target/u2_bench/job.sh; raw output target/u2_bench/iqrs_u2_bench.57958639.out.
SRR32588726.chr8.bam, 631 MB, 1,688,357 records. Full-chr fetch to RecordBatch (qname, flag, pos, mapq,
cigar, HP tag; end0 computed); MB/s = compressed bytes.

| backend | threads | no seq: rec/s | MB/s | with seq: rec/s | MB/s |
|---|---|---|---|---|---|
| noodles | 1 | 269k | 100 | 222k | 83 |
| htslib  | 1 | 142k | 53  | 127k | 48 |
| noodles | 2 | 546k | 204 | 549k | 205 |
| htslib  | 2 | 287k | 107 | 285k | 106 |
| noodles | 4 | 1.05M | 391 | 805k | 301 |
| htslib  | 4 | 576k | 215 | 565k | 211 |
| noodles | 8 | 1.38M | 517 | 746k | 279 |
| htslib  | 8 | 864k | 323 | 696k | 260 |
| pysam (4.1 loop) | 1 | 197k | - | (seq-only loop 10.0 s) | |

Random regions (2000 regions of 10-200 kb, 209 Mb in total, 6.98M records incl. overlaps, seq decoded, 1 thread):
noodles 38.6 s = 51.8 regions/s, 181k rec/s; htslib 67.9 s = 29.5 regions/s, 103k rec/s.
Real-data check vs pysam: 300 random regions, 219,509 records, identical CRC (qname/pos/flag/end/seq/HP).
htslib here = hts-sys 2.2.1 with bundled zlib (no libdeflate, default-features=false).
Verdict: noodles is 1.1-1.9x faster than htslib in every config, so keep noodles as the default; the htslib feature stays a fallback only.
With seq at 8 threads the consumer thread is the bottleneck (seq decode + 8 inflate workers on 8 cores).
