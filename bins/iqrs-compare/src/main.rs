//! `iqrs-compare`: end-to-end scoring of two IsoQuant-format output sets (design §10.1).
//!
//! Any two output sets (IsoQuant 4.1, iqrs, or two iqrs runs) can be compared. Models are matched by structure only (transcript and gene IDs are never assumed
//! to correspond): multi-exon by `(chr, strand, exact intron chain)`, mono-exon by same-strand
//! reciprocal overlap >= 0.8.
//!
//! ```text
//! iqrs-compare --ref DIR --test DIR [options]
//!   --ref-reads F / --test-reads F   per-read fate files (default: auto-detect
//!                                    *transcript_model_reads.tsv[.gz] or *read_fates.tsv[.gz])
//!   --no-reads                       skip per-read files
//!   --annotation F                   GFF3/GTF: structural known/NIC/NNIC classes and locus categories
//!   --exclude-nongene-loci           exclude loci whose GFF3 gene type is not `gene` from gated metrics
//!                                    (for a `--gene-types all` run vs IsoQuant 4.1; needs --annotation)
//!   --noise-floor F                  compare.json of a noise-floor comparison; thresholds widened by its deficits
//!   --ref-label S / --test-label S   names used in the report
//!   --top N                          residual rows in the markdown (default 40)
//!   -o DIR                           write summary.md, compare.json, compare.tsv, mismatches.tsv
//!   --explained METRIC=DEV           a failing gated metric whose name starts with METRIC is reported
//!                                    as "explained by DEV" (e.g. tx_per_gene=D-M6); repeatable
//!   --fail-on-gate                   exit 2 if a gated metric fails (explained ones do not count)
//! ```

#![allow(clippy::needless_range_loop)] // index loops over (side, class) tables read clearer

mod annot;
mod io;
mod matching;
mod metrics;
mod report;

use anyhow::{bail, Context, Result};
use io::Side;
use report::Report;
use rustc_hash::FxHashMap;
use std::path::{Path, PathBuf};
use std::time::Instant;

const USAGE: &str = "usage: iqrs-compare --ref DIR --test DIR [--ref-reads F] [--test-reads F] [--no-reads] \
[--annotation GFF3|GTF] [--exclude-nongene-loci] [--noise-floor compare.json] [--ref-label S] [--test-label S] \
[--top N] [-o OUTDIR] [--explained METRIC=DEV]... [--fail-on-gate]";

struct Args {
    r: PathBuf,
    t: PathBuf,
    rr: Option<PathBuf>,
    tr: Option<PathBuf>,
    no_reads: bool,
    annotation: Option<PathBuf>,
    exclude_nongene: bool,
    floor: Option<PathBuf>,
    rl: String,
    tl: String,
    top: usize,
    out: Option<PathBuf>,
    fail_on_gate: bool,
    explained: Vec<(String, String)>,
}

fn parse_args() -> Result<Args> {
    use lexopt::prelude::*;
    let mut p = lexopt::Parser::from_env();
    let (mut r, mut t) = (None, None);
    let mut a = Args {
        r: PathBuf::new(),
        t: PathBuf::new(),
        rr: None,
        tr: None,
        no_reads: false,
        annotation: None,
        exclude_nongene: false,
        floor: None,
        rl: String::new(),
        tl: String::new(),
        top: 40,
        out: None,
        fail_on_gate: false,
        explained: Vec::new(),
    };
    while let Some(arg) = p.next()? {
        match arg {
            Long("ref") => r = Some(PathBuf::from(p.value()?)),
            Long("test") => t = Some(PathBuf::from(p.value()?)),
            Long("ref-reads") => a.rr = Some(p.value()?.into()),
            Long("test-reads") => a.tr = Some(p.value()?.into()),
            Long("no-reads") => a.no_reads = true,
            Long("annotation") => a.annotation = Some(p.value()?.into()),
            Long("exclude-nongene-loci") => a.exclude_nongene = true,
            Long("noise-floor") => a.floor = Some(p.value()?.into()),
            Long("ref-label") => a.rl = p.value()?.string()?,
            Long("test-label") => a.tl = p.value()?.string()?,
            Long("top") => a.top = p.value()?.parse()?,
            Short('o') | Long("out") => a.out = Some(p.value()?.into()),
            Long("fail-on-gate") => a.fail_on_gate = true,
            Long("explained") => {
                let v = p.value()?.string()?;
                for item in v.split(',') {
                    let Some((m, d)) = item.split_once('=') else { bail!("--explained wants METRIC=DEV, got {item}") };
                    if m.is_empty() || d.is_empty() {
                        bail!("--explained wants METRIC=DEV, got {item}");
                    }
                    a.explained.push((m.to_string(), d.to_string()));
                }
            }
            Short('h') | Long("help") => {
                println!("{USAGE}");
                std::process::exit(0);
            }
            _ => return Err(arg.unexpected().into()),
        }
    }
    let (Some(r), Some(t)) = (r, t) else { bail!("{USAGE}") };
    let label = |p: &Path| p.display().to_string();
    if a.rl.is_empty() {
        a.rl = label(&r);
    }
    if a.tl.is_empty() {
        a.tl = label(&t);
    }
    a.r = r;
    a.t = t;
    Ok(a)
}

fn reads_path(explicit: &Option<PathBuf>, dir: &Path, no_reads: bool) -> Result<Option<PathBuf>> {
    if no_reads {
        return Ok(None);
    }
    if explicit.is_some() {
        return Ok(explicit.clone());
    }
    io::find_one(dir, "per-read", io::is_reads)
}

fn main() -> Result<()> {
    let a = parse_args()?;
    let t0 = Instant::now();
    let rr = reads_path(&a.rr, &a.r, a.no_reads)?;
    let tr = reads_path(&a.tr, &a.t, a.no_reads)?;

    // parse both sides and the annotation concurrently
    let (rs, ts, ann) = std::thread::scope(|sc| {
        let h1 = sc.spawn(|| Side::load(&a.rl, &a.r, rr.as_deref()));
        let h2 = sc.spawn(|| Side::load(&a.tl, &a.t, tr.as_deref()));
        let h3 = sc.spawn(|| a.annotation.as_deref().map(annot::load).transpose());
        (h1.join().unwrap(), h2.join().unwrap(), h3.join().unwrap())
    });
    let (mut rs, mut ts, ann) = (rs.context("ref side")?, ts.context("test side")?, ann?);
    let t_load = t0.elapsed();

    let mut chroms = Vec::new();
    let mut cmap = FxHashMap::default();
    rs.remap_chroms(&mut chroms, &mut cmap);
    ts.remap_chroms(&mut chroms, &mut cmap);

    let mt = matching::match_models(&rs.models, &ts.models);
    let nr = rs.models.len();
    let n = nr + ts.models.len();
    let model = |g: usize| if g < nr { &rs.models[g] } else { &ts.models[g - nr] };

    let (mut cat, mut cats, mut sclass) = (None, Vec::new(), None);
    if let Some(ann) = &ann {
        let sets = annot::build_sets(ann, &chroms);
        let mut cidx: FxHashMap<String, u8> = FxHashMap::default();
        for k in ["gene", "intergenic"] {
            cidx.insert(k.into(), cats.len() as u8);
            cats.push(k.to_string());
        }
        let mut cv = Vec::with_capacity(n);
        let mut sv = Vec::with_capacity(n);
        for g in 0..n {
            let m = model(g);
            let name = annot::locus_category(ann, &sets, m);
            let ci = match cidx.get(name) {
                Some(&i) => i,
                None => {
                    let i = cats.len() as u8;
                    cats.push(name.to_string());
                    cidx.insert(name.to_string(), i);
                    i
                }
            };
            cv.push(ci);
            sv.push(annot::struct_class(&sets, m));
        }
        cat = Some(cv);
        sclass = Some(sv);
    }

    let ctx = metrics::Ctx {
        r: &rs,
        t: &ts,
        mt: &mt,
        chroms: &chroms,
        cat,
        cats,
        sclass,
        exclude_nongene: a.exclude_nongene,
        top: a.top,
    };
    let mut rep = Report {
        title: format!("iqrs-compare: {} vs {}", a.rl, a.tl),
        explained: a.explained.clone(),
        ..Default::default()
    };
    rep.header_lines.push(format!("ref: `{}`", rs.dir.display()));
    rep.header_lines.push(format!("test: `{}`", ts.dir.display()));
    rep.header_lines.push(format!(
        "gated loci: {}{}",
        if a.exclude_nongene { "gene type `gene` and intergenic only" } else { "all" },
        a.annotation.as_ref().map(|p| format!("; annotation `{}`", p.display())).unwrap_or_default()
    ));
    if a.exclude_nongene && ctx.cat.is_none() {
        rep.warnings.push("--exclude-nongene-loci without --annotation: locus types unknown, all loci gated".into());
    }
    for (s, side) in [("ref", &rs), ("test", &ts)] {
        for w in &side.warnings {
            rep.warnings.push(format!("{s}: {w}"));
        }
    }
    metrics::overview(&ctx, &mut rep);
    metrics::chain(&ctx, &mut rep);
    metrics::counts(&ctx, &mut rep);
    metrics::fates(&ctx, &mut rep);
    metrics::tpg(&ctx, &mut rep);
    metrics::novelty(&ctx, &mut rep);
    metrics::grouped(&ctx, &mut rep);
    metrics::breakdowns(&ctx, &mut rep);
    let mismatches = metrics::residuals(&ctx, &mut rep);

    let floor = match &a.floor {
        Some(p) => {
            let s = std::fs::read_to_string(p).with_context(|| format!("reading noise floor {}", p.display()))?;
            rep.floor_source = Some(p.display().to_string());
            Some(iq_types::serde_json::from_str(&s).with_context(|| format!("parsing {}", p.display()))?)
        }
        None => None,
    };
    rep.judge(floor.as_ref());
    for (m, d) in &a.explained {
        if !rep.metrics.iter().any(|x| x.threshold.is_some() && x.name.starts_with(m.as_str())) {
            rep.warnings.push(format!("--explained {m}={d}: no gated metric matches"));
        }
    }
    rep.header_lines.push(format!("time: load {:.2}s, total {:.2}s", t_load.as_secs_f64(), t0.elapsed().as_secs_f64()));
    let md = rep.markdown();
    match &a.out {
        Some(o) => {
            std::fs::create_dir_all(o)?;
            std::fs::write(o.join("summary.md"), &md)?;
            std::fs::write(o.join("compare.json"), report::to_pretty(&rep.json()))?;
            std::fs::write(o.join("compare.tsv"), rep.tsv())?;
            std::fs::write(o.join("mismatches.tsv"), mismatches)?;
            eprintln!(
                "iqrs-compare: gated metrics {}; wrote {}/{{summary.md,compare.json,compare.tsv,mismatches.tsv}}",
                rep.gate_line(),
                o.display()
            );
        }
        None => print!("{md}"),
    }
    let (p, g) = rep.n_gated();
    if a.fail_on_gate && p + rep.n_explained() < g {
        std::process::exit(2);
    }
    Ok(())
}
