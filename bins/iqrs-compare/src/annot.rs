//! Annotation-derived model classes (ID-agnostic).
//!
//! - Structural novelty: `known` (intron chain equals an annotated transcript's chain; mono-exon:
//!   reciprocal overlap >= 0.8 with an annotated mono-exon transcript on the same strand), `nic`
//!   (every intron annotated, or every donor and acceptor annotated), `nnic` (otherwise),
//!   `novel_mono` (unmatched mono-exon).
//! - Locus category: the GFF3 type of the model's gene when its `gene_id` is an annotated gene,
//!   else of overlapping same-strand genes (preferring type `gene`), else `intergenic`. Types other
//!   than `gene` (ncRNA_gene, pseudogene, ...) are the D-A1 categories.
//!
//! The annotation is parsed with `GeneTypes::All` so every gene type is visible.

use crate::io::Model;
use crate::matching::{mono_hits, mono_list};
use anyhow::{Context, Result};
use iq_annot::{AnnotationIndex, GeneTypes};
use iq_types::{Iv, Strand};
use rustc_hash::{FxHashMap, FxHashSet};
use std::path::Path;

pub const STRUCT_CLASSES: [&str; 4] = ["known", "nic", "nnic", "novel_mono"];

pub fn load(p: &Path) -> Result<AnnotationIndex> {
    let name = p.to_string_lossy();
    let gtf = name.ends_with(".gtf") || name.ends_with(".gtf.gz");
    let r = crate::io::open(p)?;
    let idx = if gtf { iq_annot::parse_gtf_reader(r, GeneTypes::All) } else { iq_annot::parse_gff3_reader(r, GeneTypes::All) };
    idx.map_err(|e| anyhow::anyhow!("{e:?}")).with_context(|| format!("parsing annotation {}", p.display()))
}

fn strand_byte(s: Strand) -> u8 {
    s.as_char() as u8
}

pub struct AnnotSets {
    chains: FxHashSet<Vec<i32>>,
    introns: FxHashSet<(u32, u8, i32, i32)>,
    donors: FxHashSet<(u32, u8, i32)>,
    acceptors: FxHashSet<(u32, u8, i32)>,
    monos: Vec<(u32, u8, i32, i32, u32)>,
    /// global chr idx -> annotation chr idx
    chr_to_ann: FxHashMap<u32, u16>,
    gene_by_id: FxHashMap<String, u32>,
}

/// Builds lookup sets for the chromosomes present in `chroms` (global table).
pub fn build_sets(ann: &AnnotationIndex, chroms: &[String]) -> AnnotSets {
    let mut s = AnnotSets {
        chains: FxHashSet::default(),
        introns: FxHashSet::default(),
        donors: FxHashSet::default(),
        acceptors: FxHashSet::default(),
        monos: Vec::new(),
        chr_to_ann: FxHashMap::default(),
        gene_by_id: FxHashMap::default(),
    };
    let mut mono_models = Vec::new();
    for (gc, name) in chroms.iter().enumerate() {
        let Some(ac) = ann.chrom(name) else { continue };
        s.chr_to_ann.insert(gc as u32, ac);
        for g in ann.chrom_genes(ac) {
            s.gene_by_id.insert(ann.gene_id(g).to_string(), g);
            for t in ann.gene_txs(g) {
                let mut ex: Vec<(i32, i32)> = ann.tx_exons(t).iter().map(|iv| (iv.s, iv.e)).collect();
                ex.sort_unstable();
                let st = strand_byte(ann.tx_strand(t));
                if ex.len() == 1 {
                    mono_models.push(Model {
                        id: String::new(),
                        gene: String::new(),
                        chr: gc as u32,
                        strand: st,
                        exons: ex,
                        source: String::new(),
                        count: 0.0,
                    });
                    continue;
                }
                let mut key = vec![gc as i32, st as i32];
                for w in ex.windows(2) {
                    let (l, r) = (w[0].1, w[1].0);
                    key.push(l);
                    key.push(r);
                    s.introns.insert((gc as u32, st, l, r));
                    s.donors.insert((gc as u32, st, l));
                    s.acceptors.insert((gc as u32, st, r));
                }
                s.chains.insert(key);
            }
        }
    }
    s.monos = mono_list(mono_models.iter().enumerate().map(|(i, m)| (i as u32, m)));
    s
}

/// Structural class index into [`STRUCT_CLASSES`].
pub fn struct_class(s: &AnnotSets, m: &Model) -> u8 {
    if m.exons.is_empty() {
        return 2;
    }
    if m.is_mono() {
        let mut hit = false;
        mono_hits(&s.monos, m.chr, m.strand, m.exons[0].0, m.exons[0].1, |_, _| hit = true);
        return if hit { 0 } else { 3 };
    }
    if s.chains.contains(&crate::matching::chain_key(m)) {
        return 0;
    }
    let all_introns = m.introns().all(|(l, r)| s.introns.contains(&(m.chr, m.strand, l, r)));
    let all_sites = m
        .introns()
        .all(|(l, r)| s.donors.contains(&(m.chr, m.strand, l)) && s.acceptors.contains(&(m.chr, m.strand, r)));
    if all_introns || all_sites {
        1
    } else {
        2
    }
}

/// Locus category name (GFF3 gene type, or `intergenic`).
pub fn locus_category<'a>(ann: &'a AnnotationIndex, s: &AnnotSets, m: &Model) -> &'a str {
    if let Some(&g) = s.gene_by_id.get(&m.gene) {
        return ann.gene_type(g);
    }
    let Some(&ac) = s.chr_to_ann.get(&m.chr) else { return "intergenic" };
    let (a, b) = m.span();
    let mut best: Option<&str> = None;
    for g in ann.genes_overlapping(ac, Iv::new(a, b)) {
        let (_, _, gs) = ann.gene_span(g);
        if strand_byte(gs) != m.strand {
            continue;
        }
        let t = ann.gene_type(g);
        if t == "gene" {
            return t;
        }
        best.get_or_insert(t);
    }
    best.unwrap_or("intergenic")
}

/// Label-based class from IDs/source as 4.1 writes them: `known` (annotation source, i.e. not
/// `IsoQuant`), `nic`/`nnic` by the `.nic`/`.nnic` ID suffix, `novel_other` otherwise.
pub const LABEL_CLASSES: [&str; 4] = ["known", "nic", "nnic", "novel_other"];
pub fn label_class(m: &Model) -> u8 {
    if m.source != "IsoQuant" {
        0
    } else if m.id.ends_with(".nic") {
        1
    } else if m.id.ends_with(".nnic") {
        2
    } else {
        3
    }
}
