//! Parse time / footprint of a GFF3, plus optional checks.
//!
//! `cargo run --release -p iq_annot --example annot_bench -- <gff3> [isoquant|all] [exon_ids.tsv]`
//!
//! `exon_ids.tsv` (chr, start, end, strand, id; from 4.1 `FeatureIdStorage` over the gffutils
//! DB) is compared entry by entry with `AnnotationIndex::exon_id`.

use iq_annot::{parse_gff3, GeneTypes};
use iq_types::{Iv, Strand};
use std::collections::BTreeMap;
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args.get(1).expect("usage: annot_bench <gff3> [isoquant|all] [exon_ids.tsv]");
    let gt = match args.get(2).map(String::as_str) {
        Some("all") => GeneTypes::All,
        _ => GeneTypes::Isoquant,
    };
    let t0 = Instant::now();
    let a = parse_gff3(path.as_ref(), gt).expect("parse");
    let dt = t0.elapsed();
    let n_ex: usize = (0..a.n_txs() as u32).map(|t| a.tx_exons(t).len()).sum();
    println!(
        "mode={gt:?} parse={:.2}s chroms={} genes={} txs={} exons={} exon_ids={} heap≈{:.1} MB",
        dt.as_secs_f64(),
        a.n_chroms(),
        a.n_genes(),
        a.n_txs(),
        n_ex,
        a.n_exon_ids(),
        a.heap_bytes() as f64 / 1e6
    );
    let mut by_type: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    for g in 0..a.n_genes() as u32 {
        let e = by_type.entry(a.gene_type(g).to_owned()).or_default();
        e.0 += 1;
        e.1 += a.gene_txs(g).len();
    }
    for (t, (g, tx)) in &by_type {
        println!("  gene type {t}: {g} genes, {tx} transcripts");
    }
    let mut tx_types: BTreeMap<String, usize> = BTreeMap::new();
    for t in 0..a.n_txs() as u32 {
        *tx_types.entry(a.tx_type(t).to_owned()).or_default() += 1;
    }
    println!("  transcript types: {tx_types:?}");
    if let Some(tsv) = args.get(3) {
        let txt = std::fs::read_to_string(tsv).expect("tsv");
        let (mut ok, mut bad, mut n) = (0usize, 0usize, 0usize);
        for l in txt.lines() {
            let f: Vec<&str> = l.split('\t').collect();
            let chr = a.chrom(f[0]).expect("chr");
            let iv = Iv::new(f[1].parse().unwrap(), f[2].parse().unwrap());
            let st = match f[3] {
                "+" => Strand::Plus,
                "-" => Strand::Minus,
                _ => Strand::Unknown,
            };
            n += 1;
            if a.exon_id(chr, iv, st) == Some(f[4]) {
                ok += 1;
            } else {
                bad += 1;
                if bad <= 5 {
                    println!("  exon_id mismatch {l}: got {:?}", a.exon_id(chr, iv, st));
                }
            }
        }
        println!("exon_ids: {ok}/{n} equal, {bad} differ; rust table has {} entries", a.n_exon_ids());
    }
}
