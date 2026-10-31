//! `GFFPrinter.dump` (`model_construction/transcript_printer.py:19-138`), byte-exact
//! modulo IDs (design §8.1).
//!
//! [`GtfChrWriter`] holds the per-chromosome printer state (`printed_gene_ids`, the exon
//! id storage). Blocks are fed in genomic order through [`GtfChrWriter::dump_block`] with
//! their models already named ([`PrintModel`]) and a [`BlockAnnot`] answering the
//! annotation questions the printer asks of the block's `gene_info`. The production adapter
//! over `iq_annot` is [`AnnotBlock`]; the fixture runner has its own adapter over the
//! captured `GeneInfo` dump.

use crate::ids::ExonIds;
use iq_annot::AnnotationIndex;
use iq_types::{ChrIdx, GeneKey, IsoformIndex, Iv, ModelKind, Pos, Strand, TranscriptModel};
use rustc_hash::{FxHashMap, FxHashSet};
use std::borrow::Cow;
use std::fmt::Write;

/// One model as the printer sees it (ids already final).
#[derive(Clone, Debug)]
pub struct PrintModel<'a> {
    pub tid: Cow<'a, str>,
    pub gid: Cow<'a, str>,
    pub strand: Strand,
    pub exons: Cow<'a, [Iv]>,
    /// GFF column 2 (`IsoQuant` for novel models, annotation source for known ones).
    pub source: Cow<'a, str>,
    /// `other_features` (CDS/codons/UTR) of known models: `(start, end, type)`.
    pub other: Vec<(Pos, Pos, Cow<'a, str>)>,
    /// `additional_info` (ordered); `exons` is appended when missing.
    pub additional: Vec<(Cow<'a, str>, String)>,
}

/// What `GFFPrinter.dump` reads from the block's `gene_info`. Methods take the position of
/// a model in the block's model slice (for gene questions: any model of that gene).
pub trait BlockAnnot {
    /// `gene_info.get_gene_regions()[gene_id]` (empty when `gene_info.empty()`).
    fn gene_region(&self, model: usize) -> Option<Iv>;
    /// `gene_info.sources[gene_id]`.
    fn gene_source(&self, model: usize) -> Option<Cow<'_, str>>;
    /// `gene_info.feature_attributes[gene_id]` (present = some attribute qualified).
    fn gene_attrs(&self, model: usize) -> Option<Cow<'_, str>>;
    /// `gene_info.feature_attributes[transcript_id]`.
    fn tx_attrs(&self, model: usize) -> Option<Cow<'_, str>>;
    /// `"<tid>_<s>_<e>_<strand>" in gene_info.feature_attributes`.
    fn exon_key(&self, model: usize, s: Pos, e: Pos, strand: Strand) -> bool;
    /// Annotation `exon_id` of `(chr, s, e, strand)` (`FeatureIdStorage` initial dict).
    fn annotated_exon_id(&self, s: Pos, e: Pos, strand: Strand) -> Option<Cow<'_, str>>;
}

/// `validate_exons` (`transcript_printer.py:19-20`).
pub fn validate_exons(exons: &[Iv]) -> bool {
    exons.windows(2).all(|w| (w[0].s, w[0].e) <= (w[1].s, w[1].e)) && exons.iter().all(|x| 0 < x.s && x.s <= x.e)
}

/// Per-chromosome GTF printer state.
#[derive(Clone, Debug)]
pub struct GtfChrWriter {
    pub chr: String,
    pub printed_gene_ids: FxHashSet<String>,
    pub exon_ids: ExonIds,
}

impl GtfChrWriter {
    pub fn new(chr: &str) -> Self {
        GtfChrWriter { chr: chr.to_owned(), printed_gene_ids: FxHashSet::default(), exon_ids: ExonIds::new() }
    }

    /// `GFFPrinter.dump(gene_info, storage)`: appends the block's lines to `out`.
    pub fn dump_block(&mut self, annot: &impl BlockAnnot, models: &[PrintModel<'_>], out: &mut String) {
        if models.is_empty() {
            return;
        }
        // gene -> (strand, region, model positions), dict insertion order
        struct G {
            first: usize,
            strand: Strand,
            region: (Pos, Pos),
            models: Vec<usize>,
        }
        let mut genes: Vec<G> = Vec::new();
        let mut gpos: FxHashMap<&str, usize> = FxHashMap::default();
        for (i, m) in models.iter().enumerate() {
            if !validate_exons(&m.exons) {
                continue;
            }
            let tr = (m.exons[0].s, m.exons[m.exons.len() - 1].e);
            match gpos.get(m.gid.as_ref()) {
                None => {
                    let region = match annot.gene_region(i) {
                        Some(g) => (g.s.min(tr.0), g.e.max(tr.1)),
                        None => tr,
                    };
                    gpos.insert(m.gid.as_ref(), genes.len());
                    genes.push(G { first: i, strand: m.strand, region, models: vec![i] });
                }
                Some(&k) => {
                    let g = &mut genes[k];
                    g.models.push(i);
                    // 4.1 rebuilds the record with the *current* model's strand
                    g.strand = m.strand;
                    g.region = (g.region.0.min(tr.0), g.region.1.max(tr.1));
                }
            }
        }
        let mut order: Vec<usize> = (0..genes.len()).collect();
        order.sort_by_key(|&k| genes[k].region); // stable: ties keep insertion order

        let chr = self.chr.as_str();
        for k in order {
            let g = &genes[k];
            let gid = models[g.first].gid.as_ref();
            if !self.printed_gene_ids.contains(gid) {
                let source = annot.gene_source(g.first);
                let attrs = annot.gene_attrs(g.first);
                let _ = writeln!(
                    out,
                    "{chr}\t{}\tgene\t{}\t{}\t.\t{}\t.\tgene_id \"{gid}\"; transcripts \"{}\"; {}",
                    source.as_deref().unwrap_or("IsoQuant"),
                    g.region.0,
                    g.region.1,
                    g.strand,
                    g.models.len(),
                    attrs.as_deref().unwrap_or("")
                );
                self.printed_gene_ids.insert(gid.to_owned());
            }
            for &mi in &g.models {
                write_model(chr, &mut self.exon_ids, annot, mi, &models[mi], out);
            }
        }
    }
}

fn write_model(
    chr: &str,
    exon_ids: &mut ExonIds,
    annot: &impl BlockAnnot,
    mi: usize,
    m: &PrintModel<'_>,
    out: &mut String,
) {
    let (tid, gid) = (m.tid.as_ref(), m.gid.as_ref());
    let n_ex = m.exons.len();
    let _ = write!(
        out,
        "{chr}\t{}\ttranscript\t{}\t{}\t.\t{}\t.\tgene_id \"{gid}\"; transcript_id \"{tid}\"; ",
        m.source,
        m.exons[0].s,
        m.exons[n_ex - 1].e,
        m.strand
    );
    // additional_attributes_str(): '%s "%s";' joined by ' ', with `exons` added if missing
    let mut first = true;
    let mut has_exons = false;
    for (k, v) in &m.additional {
        if !first {
            out.push(' ');
        }
        first = false;
        has_exons |= k == "exons";
        let _ = write!(out, "{k} \"{v}\";");
    }
    if !has_exons {
        if !first {
            out.push(' ');
        }
        let _ = write!(out, "exons \"{n_ex}\";");
    }
    let tx_attrs = annot.tx_attrs(mi);
    if let Some(a) = &tx_attrs {
        out.push(' ');
        out.push_str(a);
    }
    out.push('\n');

    // exon-level lines: other_features + exons, sorted by (s, e, type); reversed on '-'
    let mut feats: Vec<(Pos, Pos, &str)> = Vec::with_capacity(m.other.len() + n_ex);
    feats.extend(m.other.iter().map(|(s, e, t)| (*s, *e, t.as_ref())));
    feats.extend(m.exons.iter().map(|x| (x.s, x.e, "exon")));
    feats.sort_unstable();
    if m.strand == Strand::Minus {
        feats.reverse();
    }
    for (i, &(s, e, ty)) in feats.iter().enumerate() {
        let _ = write!(
                out,
                "{chr}\t{}\t{ty}\t{s}\t{e}\t.\t{}\t.\tgene_id \"{gid}\"; transcript_id \"{tid}\"; exon_number \"{}\"; exon_id \"",
                m.source,
                m.strand,
                i + 1
            );
        match annot.annotated_exon_id(s, e, m.strand) {
            Some(x) => out.push_str(&x),
            None => exon_ids.write_id(chr, s, e, m.strand, out),
        }
        out.push_str("\"; ");
        if annot.exon_key(mi, s, e, m.strand) {
            // 4.1 prints the *transcript* attributes here (" " + feature_attributes[tid];
            // a missing key reads as "" through the defaultdict)
            out.push(' ');
            if let Some(a) = &tx_attrs {
                out.push_str(a);
            }
        }
        out.push('\n');
    }
}

// ------------------------------------------------------------------------------------------
// Production adapter
// ------------------------------------------------------------------------------------------

/// [`BlockAnnot`] over the run's [`AnnotationIndex`] and the block's [`IsoformIndex`]
/// (whose `genes` are the block's `gene_db_list`).
pub struct AnnotBlock<'a> {
    pub ann: &'a AnnotationIndex,
    /// Annotation chromosome of the block (`None` if the annotation lacks it).
    pub chr: Option<ChrIdx>,
    pub index: &'a IsoformIndex,
    pub models: &'a [TranscriptModel],
}

impl AnnotBlock<'_> {
    fn ref_gene(&self, i: usize) -> Option<iq_types::GeneIdx> {
        match self.models[i].gene {
            GeneKey::Ref(g) => Some(g),
            GeneKey::Novel(_) => None,
        }
    }
    fn known(&self, i: usize) -> Option<iq_types::TxIdx> {
        match self.models[i].kind {
            ModelKind::Known(t) => Some(t),
            _ => None,
        }
    }
}

fn non_empty(s: String) -> Option<Cow<'static, str>> {
    if s.is_empty() {
        None
    } else {
        Some(Cow::Owned(s))
    }
}

impl BlockAnnot for AnnotBlock<'_> {
    fn gene_region(&self, i: usize) -> Option<Iv> {
        if self.index.is_empty() {
            return None;
        }
        let key = self.models[i].gene;
        self.index.genes.iter().find(|g| g.key == key).map(|g| g.region)
    }
    fn gene_source(&self, i: usize) -> Option<Cow<'_, str>> {
        self.ref_gene(i).map(|g| Cow::Borrowed(self.ann.gene_source(g)))
    }
    fn gene_attrs(&self, i: usize) -> Option<Cow<'_, str>> {
        self.ref_gene(i).and_then(|g| non_empty(self.ann.gene_attrs(g)))
    }
    fn tx_attrs(&self, i: usize) -> Option<Cow<'_, str>> {
        self.known(i).and_then(|t| non_empty(self.ann.tx_attrs(t)))
    }
    fn exon_key(&self, i: usize, s: Pos, e: Pos, strand: Strand) -> bool {
        // set_gene_attributes: key exists when the transcript has an attribute outside the
        // exon skip list and (s, e) is an exon of its gene (any transcript of the gene).
        let Some(t) = self.known(i) else { return false };
        if self.ann.tx_exon_attrs(t).is_empty() {
            return false;
        }
        let g = self.ann.tx_gene(t);
        self.ann
            .gene_txs(g)
            .any(|u| self.ann.tx_strand(u) == strand && self.ann.tx_exons(u).iter().any(|x| x.s == s && x.e == e))
    }
    fn annotated_exon_id(&self, s: Pos, e: Pos, strand: Strand) -> Option<Cow<'_, str>> {
        self.chr.and_then(|c| self.ann.exon_id(c, Iv::new(s, e), strand)).map(Cow::Borrowed)
    }
}

/// Builds the printer view of a block's models (production): names from `names`, known
/// models' source and `other_features` from the annotation.
pub fn print_models<'a>(
    ann: &'a AnnotationIndex,
    models: &'a [TranscriptModel],
    names: &'a crate::ids::BlockNames,
) -> Vec<PrintModel<'a>> {
    models
        .iter()
        .enumerate()
        .map(|(i, m)| {
            let (source, other) = match m.kind {
                ModelKind::Known(t) => (
                    Cow::Borrowed(ann.tx_source(t)),
                    ann.tx_other_features(t).map(|(iv, ty)| (iv.s, iv.e, Cow::Borrowed(ty))).collect(),
                ),
                _ => (Cow::Borrowed("IsoQuant"), Vec::new()),
            };
            PrintModel {
                tid: Cow::Borrowed(names.tx[i].as_str()),
                gid: Cow::Borrowed(names.gene[i].as_str()),
                strand: m.strand,
                exons: Cow::Borrowed(m.exons.as_slice()),
                source,
                other,
                additional: Vec::new(),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    struct NoAnnot;
    impl BlockAnnot for NoAnnot {
        fn gene_region(&self, _: usize) -> Option<Iv> {
            None
        }
        fn gene_source(&self, _: usize) -> Option<Cow<'_, str>> {
            None
        }
        fn gene_attrs(&self, _: usize) -> Option<Cow<'_, str>> {
            None
        }
        fn tx_attrs(&self, _: usize) -> Option<Cow<'_, str>> {
            None
        }
        fn exon_key(&self, _: usize, _: Pos, _: Pos, _: Strand) -> bool {
            false
        }
        fn annotated_exon_id(&self, _: Pos, _: Pos, _: Strand) -> Option<Cow<'_, str>> {
            None
        }
    }

    fn pm<'a>(tid: &'a str, gid: &'a str, strand: Strand, ex: &[(Pos, Pos)]) -> PrintModel<'a> {
        PrintModel {
            tid: tid.into(),
            gid: gid.into(),
            strand,
            exons: Cow::Owned(ex.iter().map(|&(s, e)| Iv::new(s, e)).collect()),
            source: "IsoQuant".into(),
            other: vec![],
            additional: vec![],
        }
    }

    #[test]
    fn novel_block() {
        let mut w = GtfChrWriter::new("8");
        let ms = vec![
            pm("t2", "g2", Strand::Minus, &[(300, 400), (500, 600)]),
            pm("t1", "g1", Strand::Plus, &[(100, 200)]),
            pm("bad", "g1", Strand::Plus, &[(50, 20)]),
        ];
        let mut out = String::new();
        w.dump_block(&NoAnnot, &ms, &mut out);
        let want = "8\tIsoQuant\tgene\t100\t200\t.\t+\t.\tgene_id \"g1\"; transcripts \"1\"; \n\
8\tIsoQuant\ttranscript\t100\t200\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\"; exons \"1\";\n\
8\tIsoQuant\texon\t100\t200\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\"; exon_number \"1\"; exon_id \"1\"; \n\
8\tIsoQuant\tgene\t300\t600\t.\t-\t.\tgene_id \"g2\"; transcripts \"1\"; \n\
8\tIsoQuant\ttranscript\t300\t600\t.\t-\t.\tgene_id \"g2\"; transcript_id \"t2\"; exons \"2\";\n\
8\tIsoQuant\texon\t500\t600\t.\t-\t.\tgene_id \"g2\"; transcript_id \"t2\"; exon_number \"1\"; exon_id \"2\"; \n\
8\tIsoQuant\texon\t300\t400\t.\t-\t.\tgene_id \"g2\"; transcript_id \"t2\"; exon_number \"2\"; exon_id \"3\"; \n";
        assert_eq!(out, want);
        // a second block: genes already printed are not repeated; exon ids are reused
        let mut out2 = String::new();
        let ms2 = vec![pm("t3", "g1", Strand::Plus, &[(100, 200)])];
        w.dump_block(&NoAnnot, &ms2, &mut out2);
        assert!(!out2.contains("\tgene\t"));
        assert!(out2.contains("exon_id \"8.1\""));
    }
}
