//! # iq_sched — U21: work units, admission control, rayon orchestration, spill,
//! deterministic merge
//!
//! Replaces 4.1's `dataset_processor.py` / `parallel_workers.py`.
//!
//! - [`units`]: work units and the `IQRS_DEBUG_UNIT_SPLIT` modes;
//! - [`admission`]: byte-budget admission and budget detection;
//! - [`pipeline`]: the [`Pipeline`] trait the scheduler drives (two-phase run);
//! - [`engine`]: the coordinator ([`execute`]);
//! - [`profile`]: `--profile` (RSS sampler, per-phase wall/CPU, est-vs-actual);
//! - [`mock`]: mock stages with configurable cost/memory (tests, benchmarks).

pub mod admission;
pub mod construct;
pub mod engine;
pub mod mock;
pub mod pipeline;
pub mod profile;
pub mod real;
pub mod units;

pub use engine::{execute, EngineOpts, UnitSpec};
pub use pipeline::{Pipeline, ReadTotals};
pub use profile::Report;
pub use units::SplitMode;

use iq_types::DataType;
use std::path::PathBuf;

/// Run configuration (the `iqrs` CLI).
#[derive(Clone, Debug)]
pub struct Config {
    pub bam: PathBuf,
    pub reference: PathBuf,
    pub annotation: PathBuf,
    pub data_type: DataType,
    pub out_dir: PathBuf,
    pub threads: usize,
    /// Bytes; `None` = 0.85 × min(cgroup limit, MemAvailable).
    pub max_mem: Option<u64>,
    pub report_novel_unspliced: bool,
    pub read_group_tag: Option<String>,
    pub gene_types_isoquant: bool,
    pub prefix: String,
    pub grouped_out: Option<PathBuf>,
    /// Hidden, tests only.
    pub params_json: Option<PathBuf>,
    /// `--profile`: print per-phase wall/CPU, peak RSS, est-vs-actual bytes.
    pub profile: bool,
    /// Log scheduling decisions.
    pub verbose: bool,
    /// Unit split mode (`IQRS_DEBUG_UNIT_SPLIT`).
    pub split: SplitMode,
    /// Command line for the GTF header.
    pub cmd_line: String,
}

/// Runs the whole pipeline.
pub fn run(cfg: &Config) -> Result<Report, String> {
    let t0 = std::time::Instant::now();
    let mut p = real::RealPipeline::open(cfg)?;
    p.cmd_line = cfg.cmd_line.clone();
    if cfg.verbose || cfg.profile {
        eprintln!("[iqrs] setup {:.2}s: {} chromosomes", t0.elapsed().as_secs_f64(), p.n_chromosomes());
    }
    let opts = EngineOpts {
        threads: cfg.threads.max(1),
        budget: cfg.max_mem.unwrap_or_else(admission::default_budget),
        split: cfg.split,
        sample_rss: cfg.profile,
        verbose: cfg.verbose,
    };
    let rep = execute(&p, &opts)?;
    if cfg.profile {
        eprint!("{}", rep.summary());
    }
    Ok(rep)
}

#[cfg(test)]
mod tests;
