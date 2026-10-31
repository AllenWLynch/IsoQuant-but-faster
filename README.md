# iqrs

iqrs is a fast Rust reimplementation of [IsoQuant](https://github.com/ablab/IsoQuant) 4.1
transcript discovery and quantification, for bulk long-read RNA-seq. It reads one
coordinate-sorted, indexed BAM (PacBio or Oxford Nanopore), a reference genome and a gene
annotation. It writes IsoQuant 4.1's transcript models GTF and transcript count tables.

It aims to give the same results as IsoQuant 4.1 run with the matching options, in about 1/60
of the wall time and 1/4 of the memory. That includes IsoQuant's quirks: where 4.1 does
something surprising, iqrs does the same.

## Credit and citation

The algorithms, parameters, models and output formats are IsoQuant's. IsoQuant is developed at
the University of Helsinki and was previously developed at Saint Petersburg State University.
The polyA and TSS XGBoost models in `crates/iq_terminal/data/` are IsoQuant 4.1's own, embedded
unchanged. iqrs is an independent port and is not affiliated with the IsoQuant authors. If you
use iqrs, please cite IsoQuant:

> Prjibelski, A.D., Mikheenko, A., Joglekar, A. et al. Accurate isoform discovery with IsoQuant
> using long reads. *Nat Biotechnol* 41, 915–918 (2023). https://doi.org/10.1038/s41587-022-01565-y

## Scope

iqrs covers this IsoQuant 4.1 configuration and nothing else:

| iqrs | IsoQuant 4.1 equivalent |
|---|---|
| `--data-type pacbio` | `--data_type pacbio_ccs --model_construction_strategy sensitive_pacbio --fl_data` |
| `--data-type ont` | `--data_type nanopore --model_construction_strategy sensitive_ont` |
| `--report-novel-unspliced true\|false` | `--report_novel_unspliced` |
| `--read-group tag:TAG` | `--read_group tag:TAG` |
| `--gene-types isoquant` (default) | IsoQuant's gene loading (GFF3 type `gene`) |

Validation used GFF3 annotations. GTF input is supported but has been tested less.

**Not supported:** FASTQ input or read alignment (bring your own BAM), multiple BAMs or
samples, single-cell and spatial modes and barcode calling, other data types and model
construction strategies, read groups other than BAM tags, unsorted or unindexed BAMs, `.csi`
indexes, compressed FASTA, and the extra IsoQuant outputs (per-read assignments, gene,
exon and intron counts, TPM tables, `.usable` files). IsoQuant's `--genedb` is replaced by
reading the GFF3/GTF directly.

`--gene-types all` also loads `ncRNA_gene`, `pseudogene` and similar features as genes. It is
opt-in and has not been validated against IsoQuant 4.1.

## Install

Each [GitHub release](https://github.com/AllenWLynch/iqrs/releases) has prebuilt binaries:

- `x86_64-unknown-linux-gnu`: needs glibc 2.17 or newer, the fastest option on Linux.
- `x86_64-unknown-linux-musl` and `aarch64-unknown-linux-musl`: fully static, run on any Linux,
  about 20% slower at 16 threads.
- `aarch64-apple-darwin`: macOS on Apple silicon.

```bash
tar xzf iqrs-v0.1.0-x86_64-unknown-linux-musl.tar.gz
sha256sum -c iqrs-v0.1.0-x86_64-unknown-linux-musl.tar.gz.sha256   # optional
./iqrs-v0.1.0-x86_64-unknown-linux-musl/iqrs --version
```

### Build from source

You need Rust 1.95 or newer (`rustup`). There are no C dependencies.

```bash
git clone https://github.com/AllenWLynch/iqrs && cd iqrs
cargo build --release          # binary: target/release/iqrs
```

## Usage

```bash
samtools index sample.bam      # iqrs needs sample.bam.bai
samtools faidx genome.fa       # and genome.fa.fai
iqrs --bam sample.bam --reference genome.fa --annotation genes.gff3 \
     --data-type pacbio --read-group tag:HP --threads 16 -o iqrs_out
```

`iqrs --help` lists every option. Output does not depend on `--threads`. The validated runs
used `--report-novel-unspliced false`. The default is `true`, which is IsoQuant 4.1's default
for these strategies.

## Outputs

With the default `--prefix OUT`, the output directory gets:

| file | content |
|---|---|
| `OUT.transcript_models.gtf` | discovered transcript models, known and novel, in IsoQuant 4.1's GTF format |
| `OUT.discovered_transcript_counts.tsv` | read count per discovered transcript |
| `OUT.discovered_transcript_grouped_tag_TAG_counts.linear.tsv` | per-group counts (with `--read-group`) |
| `OUT.phased_transcript_counts.tsv` | allele-specific counts per transcript, with 90% intervals (with `--phase-vcf`) |
| `OUT.phased_gene_summary.tsv` | per-gene allelic balance and allele-specific isoform usage (with `--phase-vcf`) |

Novel transcripts and genes use IsoQuant's naming (`transcriptN.CHR.nic|nnic`,
`novel_gene_CHR_N`). The numbers `N` can differ from an IsoQuant run, so compare models by
structure, not by ID. The GTF header records the iqrs version (`4.1.0-iqrs-0.1.0`).

## Phased quantification (optional)

`--phase-vcf donor.vcf.gz` adds allele-specific transcript quantification. iqrs reads each read's
bases at the VCF's phased heterozygous SNVs during its normal BAM pass, then a Bayesian model
(Gibbs sampling, per gene) splits every transcript's reads between haplotype A (the first GT
allele) and haplotype B. Unlike counting reads by an HP tag, it uses every read: reads with no
het site and reads ambiguous between transcripts are shared out by the model.

Each transcript gets a `phase_uncertainty` between 0 and 1. Near 0, its allele split is pinned by
reads that carry het evidence. At 1, no read of that transcript covers a het site, so its split
comes from the prior and should not be interpreted. On whole GTEx genomes the step adds about 6%
CPU time (25-45 CPU-s per sample) and under 1 GB of memory. Phase blocks are not used: GT is taken
as one genome-wide phase. Other outputs are unchanged.

## Benchmark

Test data: PacBio Revio sample SRR32588726, whole genome, 16 threads, `--read-group tag:HP`,
one compute node.

| | wall time | peak memory | CPU time |
|---|---|---|---|
| IsoQuant 4.1 | 1:55:52 | 81.4 GB job total (40.9 GB largest worker) | 34,052 s |
| iqrs 0.1.0 | 1:48 | 19.7 GB | 1,016 s |

Agreement with IsoQuant 4.1 on the same sample, scored with the bundled `iqrs-compare` tool:

- 17/17 gated metrics pass.
- Read-weighted intron-chain match is 1.0000 in both directions.
- Spearman correlation of counts over matched models is 1.0000.
- iqrs reports 73,023 transcript models and IsoQuant 73,032. 4 models on each side have no
  structural match.

On chr8 of five samples (2 ONT, 1 Iso-Seq, 2 Revio), all 15 gated metrics pass for every sample.

## Development

See [CONTRIBUTING.md](CONTRIBUTING.md) for the workspace layout, tests and the IsoQuant
fixture harness.

## License

GPL-2.0-only, the same license as IsoQuant. See [LICENSE](LICENSE) and [NOTICE](NOTICE).
