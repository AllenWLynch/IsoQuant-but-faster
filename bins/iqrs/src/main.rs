//! `iqrs`: streamlined pure-Rust reimplementation of IsoQuant 4.1.0 (design §7.0 CLI).
//!
//! ```text
//! iqrs --bam X --reference fa --annotation gff3 --data-type {pacbio|ont} -o DIR
//!      [--threads N] [--max-mem G] [--report-novel-unspliced true|false]
//!      [--read-group tag:HP | --read-group-tag HP] [--gene-types all|isoquant]
//!      [--prefix OUT] [--grouped-out PATH] [--profile] [--verbose]
//! ```
//! `IQRS_DEBUG_UNIT_SPLIT=whole-chr|max` forces the work-unit split (determinism tests).

use anyhow::{bail, Result};
use iq_sched::{Config, SplitMode};
use iq_types::DataType;
use std::path::PathBuf;

fn parse() -> Result<Config> {
    use lexopt::prelude::*;
    let mut p = lexopt::Parser::from_env();
    let (mut bam, mut reference, mut annotation, mut out_dir) = (None, None, None, None);
    let mut data_type = None;
    let mut cfg_threads = 1usize;
    let mut max_mem = None;
    let mut report_novel_unspliced = true;
    let mut read_group_tag = None;
    let mut gene_types_isoquant = false;
    let mut prefix = "OUT".to_string();
    let mut grouped_out = None;
    let mut params_json = None;
    let (mut profile, mut verbose) = (false, false);
    while let Some(arg) = p.next()? {
        match arg {
            Long("bam") => bam = Some(PathBuf::from(p.value()?)),
            Long("reference") => reference = Some(PathBuf::from(p.value()?)),
            Long("annotation") => annotation = Some(PathBuf::from(p.value()?)),
            Short('o') | Long("output") => out_dir = Some(PathBuf::from(p.value()?)),
            Long("data-type") => {
                data_type = Some(match p.value()?.string()?.as_str() {
                    "pacbio" | "pacbio_ccs" => DataType::Pacbio,
                    "ont" | "nanopore" => DataType::Ont,
                    d => bail!("unknown --data-type {d}"),
                })
            }
            Long("threads") | Short('t') => cfg_threads = p.value()?.parse()?,
            Long("max-mem") => {
                let g: f64 = p.value()?.parse()?;
                max_mem = Some((g * (1u64 << 30) as f64) as u64);
            }
            Long("report-novel-unspliced") => {
                report_novel_unspliced = match p.value()?.string()?.as_str() {
                    "true" => true,
                    "false" => false,
                    v => bail!("--report-novel-unspliced expects true|false, got {v}"),
                }
            }
            Long("read-group-tag") => read_group_tag = Some(p.value()?.string()?),
            Long("read-group") => {
                let v = p.value()?.string()?;
                match v.strip_prefix("tag:") {
                    Some(t) => read_group_tag = Some(t.to_string()),
                    None => bail!("--read-group: only tag:<TAG> is supported, got {v}"),
                }
            }
            Long("profile") => profile = true,
            Long("verbose") | Short('v') => verbose = true,
            Long("gene-types") => {
                gene_types_isoquant = match p.value()?.string()?.as_str() {
                    "all" => false,
                    "isoquant" => true,
                    v => bail!("--gene-types expects all|isoquant, got {v}"),
                }
            }
            Long("prefix") => prefix = p.value()?.string()?,
            Long("grouped-out") => grouped_out = Some(PathBuf::from(p.value()?)),
            Long("params-json") => params_json = Some(PathBuf::from(p.value()?)),
            Short('h') | Long("help") => {
                println!(
                    "usage: iqrs --bam X --reference fa --annotation gff3|gtf --data-type pacbio|ont -o DIR\n  \
                     [--threads N] [--max-mem GiB] [--report-novel-unspliced true|false]\n  \
                     [--read-group tag:HP] [--gene-types all|isoquant] [--prefix OUT] [--grouped-out PATH]\n  \
                     [--profile] [--verbose]"
                );
                std::process::exit(0);
            }
            _ => return Err(arg.unexpected().into()),
        }
    }
    let need = |o: Option<PathBuf>, n: &str| o.ok_or_else(|| anyhow::anyhow!("missing --{n}"));
    Ok(Config {
        bam: need(bam, "bam")?,
        reference: need(reference, "reference")?,
        annotation: need(annotation, "annotation")?,
        data_type: data_type.ok_or_else(|| anyhow::anyhow!("missing --data-type"))?,
        out_dir: need(out_dir, "output")?,
        threads: cfg_threads.max(1),
        max_mem,
        report_novel_unspliced,
        read_group_tag,
        gene_types_isoquant,
        prefix,
        grouped_out,
        params_json,
        profile,
        verbose,
        split: SplitMode::from_env().map_err(|e| anyhow::anyhow!(e))?,
        cmd_line: std::env::args().collect::<Vec<_>>().join(" "),
    })
}

fn main() -> Result<()> {
    let cfg = parse()?;
    iq_sched::run(&cfg).map_err(|e| anyhow::anyhow!(e))?;
    Ok(())
}
