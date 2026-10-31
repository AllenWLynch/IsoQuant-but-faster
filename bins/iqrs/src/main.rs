//! `iqrs`: a fast Rust reimplementation of IsoQuant 4.1 transcript discovery and quantification
//! for bulk long-read RNA-seq (one coordinate-sorted, indexed BAM; PacBio or ONT).
//!
//! Run `iqrs --help` for usage. Hidden options (tests only): `--params-json FILE` (full IsoQuant
//! parameter dump instead of the built-in preset); env `IQRS_DEBUG_UNIT_SPLIT=whole-chr|max`
//! forces the work-unit split (outputs never depend on it).

use iq_sched::{Config, SplitMode};
use iq_types::DataType;
use std::ffi::OsString;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

#[cfg(target_env = "musl")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

const VERSION: &str = env!("CARGO_PKG_VERSION");

fn help() -> String {
    format!(
        "\
iqrs {VERSION}: transcript discovery and quantification from long RNA reads,
a Rust reimplementation of IsoQuant 4.1 (bulk mode, one BAM, PacBio or ONT).

USAGE:
    iqrs --bam FILE --reference FILE --annotation FILE --data-type pacbio|ont -o DIR [OPTIONS]

REQUIRED:
    --bam FILE                 Coordinate-sorted BAM of reads aligned to the genome, with a .bai
                               index next to it (FILE.bai or FILE with .bam replaced by .bai)
    --reference FILE           Genome FASTA (uncompressed) with a samtools .fai index (FILE.fai)
    --annotation FILE          Gene annotation, GFF3 (.gff3/.gff) or GTF (.gtf)
    --data-type TYPE           pacbio (alias pacbio_ccs): IsoQuant sensitive_pacbio + --fl_data
                               ont (alias nanopore):      IsoQuant sensitive_ont
    -o, --output DIR           Output directory (created if missing)

OPTIONS:
    -t, --threads N            Worker threads [default: 1]; outputs do not depend on N
    --report-novel-unspliced true|false
                               Report novel mono-exonic transcripts [default: true, as IsoQuant
                               4.1 for the sensitive strategies]
    --read-group tag:TAG       Also count transcripts per value of BAM tag TAG (e.g. tag:HP);
                               reads without the tag are counted under NA
    --gene-types isoquant|all  Genes to load from the annotation [default: isoquant]
                               isoquant: features of type `gene` (IsoQuant 4.1 behaviour)
                               all: every gene-like feature with transcripts (ncRNA_gene,
                               pseudogene, ...); not validated against IsoQuant 4.1
    --prefix NAME              Output file name prefix [default: OUT]
    --grouped-out FILE         Path of the grouped counts file [default: in DIR, see below]
    --max-mem GIB              Soft memory budget for scheduling work [default: 85% of the
                               smaller of the cgroup limit and available memory]
    --profile                  Print per-phase wall/CPU time and peak memory to stderr
    -v, --verbose              Log progress and scheduling decisions to stderr
    -h, --help                 Print this help
    -V, --version              Print the version

OUTPUTS (in DIR, with the default prefix OUT):
    OUT.transcript_models.gtf                     discovered transcript models (known and novel)
    OUT.discovered_transcript_counts.tsv          read counts per discovered transcript
    OUT.discovered_transcript_grouped_tag_TAG_counts.linear.tsv
                                                  per-group counts (only with --read-group)
These files have IsoQuant 4.1's formats; novel transcript and gene ids follow IsoQuant's naming
(transcriptN.CHR.nic|nnic, novel_gene_CHR_N).

EXAMPLE:
    iqrs --bam sample.bam --reference genome.fa --annotation genes.gff3 \\
         --data-type pacbio --read-group tag:HP --threads 16 -o iqrs_out
"
    )
}

/// A command-line error (exit code 2) as opposed to a run failure (exit code 1).
struct Usage(String);

fn usage<T>(msg: impl Into<String>) -> Result<T, Usage> {
    Err(Usage(msg.into()))
}

fn value(p: &mut lexopt::Parser, flag: &str) -> Result<OsString, Usage> {
    p.value().or_else(|_| usage(format!("{flag} needs a value")))
}

fn string(p: &mut lexopt::Parser, flag: &str) -> Result<String, Usage> {
    value(p, flag)?.into_string().or_else(|_| usage(format!("{flag}: value is not valid UTF-8")))
}

fn parse_bool(flag: &str, v: &str) -> Result<bool, Usage> {
    match v {
        "true" | "yes" | "1" => Ok(true),
        "false" | "no" | "0" => Ok(false),
        _ => usage(format!("{flag} expects true or false, got {v:?}")),
    }
}

/// `Ok(None)` = help or version was printed.
fn parse() -> Result<Option<Config>, Usage> {
    use lexopt::prelude::*;
    let mut p = lexopt::Parser::from_env();
    let (mut bam, mut reference, mut annotation, mut out_dir) = (None, None, None, None);
    let mut data_type = None;
    let mut threads = 1usize;
    let mut max_mem = None;
    let mut report_novel_unspliced = true;
    let mut read_group_tag = None;
    let mut gene_types_isoquant = true;
    let mut prefix = "OUT".to_string();
    let mut grouped_out = None;
    let mut params_json = None;
    let (mut profile, mut verbose) = (false, false);
    while let Some(arg) = p.next().map_err(|e| Usage(e.to_string()))? {
        match arg {
            Long("bam") => bam = Some(PathBuf::from(value(&mut p, "--bam")?)),
            Long("reference") => reference = Some(PathBuf::from(value(&mut p, "--reference")?)),
            Long("annotation") => annotation = Some(PathBuf::from(value(&mut p, "--annotation")?)),
            Short('o') | Long("output") => out_dir = Some(PathBuf::from(value(&mut p, "--output")?)),
            Long("data-type") => {
                data_type = Some(match string(&mut p, "--data-type")?.as_str() {
                    "pacbio" | "pacbio_ccs" => DataType::Pacbio,
                    "ont" | "nanopore" => DataType::Ont,
                    d => return usage(format!("--data-type must be pacbio or ont, got {d:?}")),
                })
            }
            Short('t') | Long("threads") => {
                let v = string(&mut p, "--threads")?;
                threads = match v.parse::<usize>() {
                    Ok(n) if n >= 1 => n,
                    _ => return usage(format!("--threads expects a positive integer, got {v:?}")),
                };
            }
            Long("max-mem") => {
                let v = string(&mut p, "--max-mem")?;
                match v.parse::<f64>() {
                    Ok(g) if g > 0.0 && g.is_finite() => max_mem = Some((g * (1u64 << 30) as f64) as u64),
                    _ => return usage(format!("--max-mem expects a positive number of GiB, got {v:?}")),
                }
            }
            Long("report-novel-unspliced") => {
                let v = string(&mut p, "--report-novel-unspliced")?;
                report_novel_unspliced = parse_bool("--report-novel-unspliced", &v)?;
            }
            Long("read-group") => {
                let v = string(&mut p, "--read-group")?;
                match v.strip_prefix("tag:") {
                    Some(t) if t.len() == 2 => read_group_tag = Some(t.to_string()),
                    Some(t) => return usage(format!("--read-group: a BAM tag has 2 characters, got {t:?}")),
                    None => return usage(format!("--read-group: only tag:TAG is supported, got {v:?}")),
                }
            }
            // older spelling of --read-group tag:TAG
            Long("read-group-tag") => {
                let t = string(&mut p, "--read-group-tag")?;
                if t.len() != 2 {
                    return usage(format!("--read-group-tag: a BAM tag has 2 characters, got {t:?}"));
                }
                read_group_tag = Some(t);
            }
            Long("gene-types") => {
                gene_types_isoquant = match string(&mut p, "--gene-types")?.as_str() {
                    "isoquant" => true,
                    "all" => false,
                    v => return usage(format!("--gene-types must be isoquant or all, got {v:?}")),
                }
            }
            Long("prefix") => prefix = string(&mut p, "--prefix")?,
            Long("grouped-out") => grouped_out = Some(PathBuf::from(value(&mut p, "--grouped-out")?)),
            Long("params-json") => params_json = Some(PathBuf::from(value(&mut p, "--params-json")?)),
            Long("profile") => profile = true,
            Short('v') | Long("verbose") => verbose = true,
            Short('h') | Long("help") => {
                print!("{}", help());
                return Ok(None);
            }
            Short('V') | Long("version") => {
                println!("iqrs {VERSION} (IsoQuant 4.1.0 compatible)");
                return Ok(None);
            }
            Value(v) => return usage(format!("unexpected argument {:?}", v.to_string_lossy())),
            _ => return Err(Usage(arg.unexpected().to_string())),
        }
    }
    let mut missing = Vec::new();
    for (v, n) in
        [(&bam, "--bam"), (&reference, "--reference"), (&annotation, "--annotation"), (&out_dir, "-o/--output")]
    {
        if v.is_none() {
            missing.push(n);
        }
    }
    if data_type.is_none() {
        missing.push("--data-type");
    }
    if !missing.is_empty() {
        return usage(format!("missing required option(s): {}", missing.join(", ")));
    }
    let split = SplitMode::from_env().map_err(|e| Usage(format!("IQRS_DEBUG_UNIT_SPLIT: {e}")))?;
    Ok(Some(Config {
        bam: bam.unwrap(),
        reference: reference.unwrap(),
        annotation: annotation.unwrap(),
        data_type: data_type.unwrap(),
        out_dir: out_dir.unwrap(),
        threads,
        max_mem,
        report_novel_unspliced,
        read_group_tag,
        gene_types_isoquant,
        prefix,
        grouped_out,
        params_json,
        profile,
        verbose,
        split,
        cmd_line: std::env::args().collect::<Vec<_>>().join(" "),
    }))
}

/// Checks the inputs before any work, so mistakes are reported in one clear line.
fn preflight(cfg: &Config) -> Result<(), String> {
    let file = |p: &Path, what: &str| -> Result<(), String> {
        match std::fs::metadata(p) {
            Ok(m) if m.is_file() => Ok(()),
            Ok(_) => Err(format!("{what} {} is not a file", p.display())),
            Err(e) => Err(format!("{what} {}: {e}", p.display())),
        }
    };
    file(&cfg.bam, "BAM")?;
    file(&cfg.reference, "reference")?;
    file(&cfg.annotation, "annotation")?;
    if let Some(p) = &cfg.params_json {
        file(p, "params JSON")?;
    }
    std::fs::create_dir_all(&cfg.out_dir)
        .map_err(|e| format!("cannot create output directory {}: {e}", cfg.out_dir.display()))?;
    Ok(())
}

fn main() -> ExitCode {
    if std::env::args_os().len() == 1 {
        eprint!("{}", help());
        return ExitCode::from(2);
    }
    let cfg = match parse() {
        Ok(Some(cfg)) => cfg,
        Ok(None) => return ExitCode::SUCCESS,
        Err(Usage(msg)) => {
            eprintln!("iqrs: error: {msg}\nRun `iqrs --help` for usage.");
            return ExitCode::from(2);
        }
    };
    if let Err(e) = preflight(&cfg).and_then(|_| iq_sched::run(&cfg).map(|_| ())) {
        eprintln!("iqrs: error: {e}");
        return ExitCode::FAILURE;
    }
    ExitCode::SUCCESS
}
