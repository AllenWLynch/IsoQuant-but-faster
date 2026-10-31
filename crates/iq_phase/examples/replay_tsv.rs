//! Sampler replay on the Python prototype's per-gene inputs (unit-pinned check of the sampler
//! alone): `replay_tsv IQRS_READS.tsv HET_READS.tsv GENE_ID OUT_DIR PREFIX [BURNIN SWEEPS SEED]`.
//! Unit construction follows phased_model/code/analyze.py `load_gene` + model-read filter.
use iq_phase::{locus::Locus, run_loci, sampler::Cfg, write_outputs};
use std::collections::{HashMap, HashSet};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let (reads, het, gene, out, prefix) = (&a[1], &a[2], &a[3], &a[4], &a[5]);
    let num = |i: usize, d: u64| a.get(i).map_or(d, |s| s.parse().unwrap());
    let cfg = Cfg { burnin: num(6, 1000) as u32, sweeps: num(7, 5000) as u32, ..Default::default() };
    let seed = num(8, 7);
    // het: read_id -> logLA - logLB (first row per read)
    let mut d: HashMap<String, f64> = HashMap::new();
    let txt = std::fs::read_to_string(het).unwrap();
    let mut lines = txt.lines();
    let h: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let (ia, ib) = (h.iter().position(|&x| x == "logLA").unwrap(), h.iter().position(|&x| x == "logLB").unwrap());
    for l in lines {
        let f: Vec<&str> = l.split('\t').collect();
        d.entry(f[0].to_string()).or_insert(f[ia].parse::<f64>().unwrap() - f[ib].parse::<f64>().unwrap());
    }
    let txt = std::fs::read_to_string(reads).unwrap();
    let mut lines = txt.lines();
    let h: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let col = |n: &str| h.iter().position(|&x| x == n).unwrap();
    let (ir, ig, it, ic) = (col("read_id"), col("gene_ids"), col("assignment_type"), col("compatible_transcripts"));
    let mut seen = HashSet::new();
    let mut rows: Vec<(Vec<String>, f64)> = Vec::new();
    for l in lines {
        let f: Vec<&str> = l.split('\t').collect();
        if f[ig] != gene || !seen.insert(f[ir].to_string()) {
            continue;
        }
        if !matches!(f[it], "unique" | "unique_minor_difference" | "ambiguous") || f[ic] == "." {
            continue;
        }
        rows.push((f[ic].split(',').map(String::from).collect(), *d.get(f[ir]).unwrap_or(&0.0)));
    }
    let l = Locus::from_reads("chr", rows.iter().map(|(c, d)| (c.iter().map(|s| s.as_str()).collect(), *d)), |_| {
        gene.clone()
    });
    eprintln!(
        "{gene}: {} reads, {} transcripts, {} classes, sum K {}, classes n>16 {}",
        l.n_reads(),
        l.tx.len(),
        l.n_classes(),
        l.cls_tx.len(),
        l.cls_n.iter().filter(|&&n| n > 16).count()
    );
    let loci = vec![l];
    let (res, st) = run_loci(&loci, &cfg, seed);
    eprintln!("sampling {:.3}s ({:.2} us/sweep)", st.cpu_s, st.cpu_s * 1e6 / (cfg.burnin + cfg.sweeps) as f64);
    std::fs::create_dir_all(out).unwrap();
    write_outputs(std::path::Path::new(out), prefix, &loci, &res).unwrap();
}
