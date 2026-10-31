//! # iq_phase — optional phased (allele-specific) transcript quantification
//!
//! Turned on by `iqrs --phase-vcf`. Pieces:
//! - [`vcf`]: phased het SNVs of one sample (GT; A = first allele);
//! - [`evidence`]: per-read A/B calls with base qualities, collected inline in the read pass;
//! - [`locus`]: model reads -> loci (gene components) -> read equivalence classes;
//! - [`sampler`]: blocked Gibbs sampler per locus (phased_model/code/gibbs.py, alpha = 1);
//! - [`run_loci`] / [`write_outputs`]: rayon over loci and the two TSV outputs.

pub mod evidence;
pub mod locus;
pub mod sampler;
pub mod vcf;

pub use evidence::{BlockEvidence, EvidenceBuf};
pub use locus::{ChrAcc, Locus};
pub use sampler::{Cfg, LocusSummary};
pub use vcf::{ChrSites, HetSites};

use rayon::prelude::*;
use std::io::Write;
use std::path::Path;
use std::time::Instant;

/// Model reads of one block (construct output): read-name hash, compatible final-model keys
/// (block-local [`iq_types::ModelIdx`]) and log lambda(A) - log lambda(B).
#[derive(Clone, Debug, Default)]
pub struct PhaseReads {
    pub hash: Vec<u64>,
    /// Start of each read's models in `models`; one extra entry at the end (empty = no reads).
    pub off: Vec<u32>,
    pub models: Vec<u32>,
    pub d: Vec<f64>,
}

impl PhaseReads {
    /// True for the assignment types the model uses (unique, unique_minor_difference,
    /// ambiguous).
    pub fn is_model_type(a: iq_types::AssignmentType) -> bool {
        use iq_types::AssignmentType as A;
        matches!(a, A::Unique | A::UniqueMinorDifference | A::Ambiguous)
    }
    pub fn push(&mut self, hash: u64, models: &[u32], d: f64) {
        if self.off.is_empty() {
            self.off.push(0);
        }
        self.hash.push(hash);
        self.models.extend_from_slice(models);
        self.off.push(self.models.len() as u32);
        self.d.push(d);
    }
    pub fn len(&self) -> usize {
        self.hash.len()
    }
    pub fn is_empty(&self) -> bool {
        self.hash.is_empty()
    }
    pub fn models_of(&self, i: usize) -> &[u32] {
        &self.models[self.off[i] as usize..self.off[i + 1] as usize]
    }
}

/// Timing and size of a phasing run.
#[derive(Clone, Debug, Default)]
pub struct RunStats {
    pub loci: usize,
    pub genes: usize,
    pub transcripts: usize,
    pub reads: u64,
    pub classes: usize,
    pub wall_s: f64,
    /// Sum of per-locus sampler times (CPU seconds, single-threaded work).
    pub cpu_s: f64,
    /// The slowest loci: (seconds, locus index).
    pub slowest: Vec<(f64, usize)>,
}

/// Samples every locus (rayon, largest first) and returns the summaries in input order.
pub fn run_loci(loci: &[Locus], cfg: &Cfg, seed: u64) -> (Vec<LocusSummary>, RunStats) {
    let t0 = Instant::now();
    // largest loci first for load balance
    let mut order: Vec<usize> = (0..loci.len()).collect();
    let cost = |l: &Locus| l.cls_tx.len() as u64 * 4 + l.n_reads() / 4 + l.tx.len() as u64 * 8 + 16;
    order.sort_by_key(|&i| std::cmp::Reverse(cost(&loci[i])));
    // dynamic scheduling: every worker takes the next-largest locus (LPT order)
    let next = std::sync::atomic::AtomicUsize::new(0);
    let nw = rayon::current_num_threads().max(1);
    let res: Vec<(usize, LocusSummary, f64)> = (0..nw)
        .into_par_iter()
        .flat_map_iter(|_| {
            let mut s = sampler::Scratch::default();
            let mut mine = Vec::new();
            loop {
                let k = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                let Some(&i) = order.get(k) else { break };
                let t = Instant::now();
                let l = &loci[i];
                let r = sampler::run(l, cfg, sampler::locus_seed(seed, &l.name()), &mut s);
                mine.push((i, r, t.elapsed().as_secs_f64()));
            }
            mine
        })
        .collect();
    let mut out: Vec<Option<LocusSummary>> = vec![None; loci.len()];
    let mut cpu = 0.0;
    let mut slowest = Vec::with_capacity(res.len());
    for (i, r, t) in res {
        out[i] = Some(r);
        cpu += t;
        slowest.push((t, i));
    }
    slowest.sort_by(|a, b| b.0.total_cmp(&a.0));
    slowest.truncate(5);
    let st = RunStats {
        loci: loci.len(),
        genes: loci.iter().map(|l| l.genes.len()).sum(),
        transcripts: loci.iter().map(|l| l.tx.len()).sum(),
        reads: loci.iter().map(|l| l.n_reads()).sum(),
        classes: loci.iter().map(|l| l.n_classes()).sum(),
        wall_s: t0.elapsed().as_secs_f64(),
        cpu_s: cpu,
        slowest,
    };
    (out.into_iter().map(|x| x.unwrap()).collect(), st)
}

/// `%.6g`-like formatting (shortest of fixed / scientific, trailing zeros trimmed).
pub fn fmt_g(x: f64) -> String {
    if x == 0.0 {
        return "0".into();
    }
    if !x.is_finite() {
        return format!("{x}");
    }
    let e = x.abs().log10().floor() as i32;
    let s = if (-5..6).contains(&e) {
        let dec = (5 - e).max(0) as usize;
        format!("{x:.dec$}")
    } else {
        format!("{x:.5e}")
    };
    if s.contains('e') {
        let (m, ex) = s.split_once('e').unwrap();
        let m = if m.contains('.') { m.trim_end_matches('0').trim_end_matches('.') } else { m };
        format!("{m}e{ex}")
    } else if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        s
    }
}

/// Posterior share of a transcript's reads on haplotype A; `NA` when it has no reads.
fn phase_prob(a: f64, b: f64) -> String {
    if a + b > 0.0 {
        fmt_g(a / (a + b))
    } else {
        "NA".to_string()
    }
}

/// Writes `PREFIX.phased_transcript_counts.tsv` and `PREFIX.phased_gene_summary.tsv`.
pub fn write_outputs(dir: &Path, prefix: &str, loci: &[Locus], res: &[LocusSummary]) -> std::io::Result<()> {
    let tp = dir.join(format!("{prefix}.phased_transcript_counts.tsv"));
    let gp = dir.join(format!("{prefix}.phased_gene_summary.tsv"));
    let mut t = std::io::BufWriter::new(std::fs::File::create(tp)?);
    let mut g = std::io::BufWriter::new(std::fs::File::create(gp)?);
    writeln!(t, "transcript\tgene\tA_mean\tA_lo90\tA_hi90\tB_mean\tB_lo90\tB_hi90\tphase_prob_A\tphase_uncertainty\tevid_A\tevid_B")?;
    writeln!(
        g,
        "gene\tn_reads\trho_mean\trho_lo90\trho_hi90\tP_rho_gt_0.5\ttv_theta\tphase_uncertainty\tevid_A\tevid_B"
    )?;
    for (l, r) in loci.iter().zip(res) {
        for (i, x) in r.tx.iter().enumerate() {
            writeln!(
                t,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                l.tx[i],
                l.tx_gene[i],
                fmt_g(x.a_mean),
                fmt_g(x.a_lo),
                fmt_g(x.a_hi),
                fmt_g(x.b_mean),
                fmt_g(x.b_lo),
                fmt_g(x.b_hi),
                phase_prob(x.a_mean, x.b_mean),
                fmt_g(x.phase_uncertainty),
                fmt_g(x.evid_a),
                fmt_g(x.evid_b)
            )?;
        }
        writeln!(
            g,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            l.name(),
            r.n_reads,
            fmt_g(r.rho_mean),
            fmt_g(r.rho_lo),
            fmt_g(r.rho_hi),
            fmt_g(r.p_rho_gt_half),
            fmt_g(r.tv_mean),
            fmt_g(r.phase_uncertainty),
            fmt_g(r.evid_a),
            fmt_g(r.evid_b)
        )?;
    }
    t.flush()?;
    g.flush()
}

#[cfg(test)]
mod tests {
    #[test]
    fn fmt_g_like_printf() {
        use super::fmt_g;
        assert_eq!(fmt_g(0.5), "0.5");
        assert_eq!(fmt_g(1234.5678), "1234.57");
        assert_eq!(fmt_g(3.0), "3");
        assert_eq!(fmt_g(1.5e-7), "1.5e-7");
        assert_eq!(fmt_g(0.000123456789), "0.000123457");
        assert_eq!(fmt_g(-2.25), "-2.25");
    }
}
