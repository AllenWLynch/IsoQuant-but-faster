//! Phase-W ID assignment (design §5.6): `transcript<N>.<chr>.nic|nnic` and
//! `novel_gene_<chr>_<N>` from one per-chromosome counter that skips `forbidden_ids`
//! (`ExcludingIdDistributor`, `utils/id_policy.py:20-51`; KEEP, #219), and the novel
//! `exon_id` storage (`FeatureIdStorage`, `:55-75`; a new id is printed bare, as in 4.1).
//!
//! D-M1: numbers are consumed only by emitted models and genes, in genomic block order and
//! storage order within a block, so they depend only on the blocks (not on scheduling or
//! work-unit splitting). Within a model the transcript number comes first, then the gene
//! number when the model introduces a new novel gene, as in `fl_graph_constructor.py:102,185`.

use iq_types::{GeneKey, ModelKind, Pos, Strand, TranscriptModel};
use rustc_hash::FxHashMap;

/// `ExcludingIdDistributor`: increments, skipping forbidden numbers.
#[derive(Clone, Debug, Default)]
pub struct IdDistributor {
    pub value: u64,
    /// Sorted forbidden numbers.
    forbidden: Vec<u64>,
}

impl IdDistributor {
    /// `forbidden` must be sorted (as `AnnotationIndex::forbidden_ids` returns it).
    pub fn new(forbidden: &[u64]) -> Self {
        debug_assert!(forbidden.windows(2).all(|w| w[0] <= w[1]));
        IdDistributor { value: 0, forbidden: forbidden.to_vec() }
    }
    pub fn increment(&mut self) -> u64 {
        self.value += 1;
        while self.forbidden.binary_search(&self.value).is_ok() {
            self.value += 1;
        }
        self.value
    }
}

/// Per-chromosome namer for novel transcripts and genes.
#[derive(Clone, Debug)]
pub struct NovelNamer {
    chr: String,
    dist: IdDistributor,
}

/// Final names of one block's models (storage order) and their genes.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BlockNames {
    /// `transcript_id` per model, storage order.
    pub tx: Vec<String>,
    /// `gene_id` per model, storage order.
    pub gene: Vec<String>,
}

impl NovelNamer {
    pub fn new(chr: &str, forbidden: &[u64]) -> Self {
        NovelNamer { chr: chr.to_owned(), dist: IdDistributor::new(forbidden) }
    }

    /// Next transcript id with the kind's suffix.
    pub fn next_transcript(&mut self, kind: ModelKind) -> String {
        let n = self.dist.increment();
        let suffix = if kind == ModelKind::Nic { "nic" } else { "nnic" };
        format!("transcript{n}.{}.{suffix}", self.chr)
    }

    /// Next novel gene id.
    pub fn next_gene(&mut self) -> String {
        let n = self.dist.increment();
        format!("novel_gene_{}_{n}", self.chr)
    }

    /// Names one block's models in storage order. `tx_name(t)` / `gene_name(g)` give the
    /// annotation ids of known transcripts and annotated genes.
    pub fn name_block<'a>(
        &mut self,
        models: &[TranscriptModel],
        tx_name: impl Fn(iq_types::TxIdx) -> &'a str,
        gene_name: impl Fn(iq_types::GeneIdx) -> &'a str,
    ) -> BlockNames {
        let mut out = BlockNames { tx: Vec::with_capacity(models.len()), gene: Vec::with_capacity(models.len()) };
        let mut novel: FxHashMap<u32, usize> = FxHashMap::default();
        for m in models {
            let t = match m.kind {
                ModelKind::Known(t) => tx_name(t).to_owned(),
                k => self.next_transcript(k),
            };
            let g = match m.gene {
                GeneKey::Ref(g) => gene_name(g).to_owned(),
                GeneKey::Novel(k) => match novel.get(&k) {
                    Some(&i) => out.gene[i].clone(),
                    None => {
                        novel.insert(k, out.gene.len());
                        self.next_gene()
                    }
                },
            };
            out.tx.push(t);
            out.gene.push(g);
        }
        out
    }

    pub fn last_value(&self) -> u64 {
        self.dist.value
    }
}

/// `FeatureIdStorage` for novel exon ids of one chromosome (annotation ids are looked up
/// by the caller first). `SimpleIDDistributor`: no exclusions.
#[derive(Clone, Debug, Default)]
pub struct ExonIds {
    pub counter: u64,
    ids: FxHashMap<(Pos, Pos, Strand), u64>,
}

impl ExonIds {
    pub fn new() -> Self {
        Self::default()
    }
    /// Registers an id generated earlier (fixture seeding).
    pub fn seed(&mut self, s: Pos, e: Pos, strand: Strand, n: u64) {
        self.ids.insert((s, e, strand), n);
    }
    /// Appends the exon id of a novel exon to `out`: `<chr>.<N>` (stored or new).
    /// Matches IsoQuant 4.1: a bare `N` is printed when the id is created (`id_policy.py:73-80`).
    pub fn write_id(&mut self, chr: &str, s: Pos, e: Pos, strand: Strand, out: &mut String) {
        use std::fmt::Write;
        match self.ids.get(&(s, e, strand)) {
            Some(&n) => {
                let _ = write!(out, "{chr}.{n}");
            }
            None => {
                self.counter += 1;
                let n = self.counter;
                self.ids.insert((s, e, strand), n);
                let _ = write!(out, "{n}");
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use iq_types::Iv;

    fn model(key: u32, kind: ModelKind, gene: GeneKey) -> TranscriptModel {
        TranscriptModel {
            key,
            strand: Strand::Plus,
            gene,
            kind,
            exons: vec![Iv::new(1, 10)],
            intron_path: Box::new([]),
            polya_confirmed: false,
            tss_confirmed: false,
            source_known: None,
        }
    }

    #[test]
    fn excluding_distributor() {
        let mut d = IdDistributor::new(&[2, 3, 5]);
        let v: Vec<u64> = (0..4).map(|_| d.increment()).collect();
        assert_eq!(v, [1, 4, 6, 7]);
    }

    #[test]
    fn naming_order() {
        let mut n = NovelNamer::new("8", &[2]);
        let models = vec![
            model(0, ModelKind::Nnic, GeneKey::Novel(7)),
            model(1, ModelKind::Known(0), GeneKey::Ref(0)),
            model(2, ModelKind::Nic, GeneKey::Novel(7)),
            model(3, ModelKind::Nic, GeneKey::Ref(0)),
            model(4, ModelKind::Nnic, GeneKey::Novel(1)),
        ];
        let b = n.name_block(&models, |_| "transcript:K", |_| "gene:G");
        assert_eq!(
            b.tx,
            ["transcript1.8.nnic", "transcript:K", "transcript4.8.nic", "transcript5.8.nic", "transcript6.8.nnic"]
        );
        assert_eq!(b.gene, ["novel_gene_8_3", "gene:G", "novel_gene_8_3", "gene:G", "novel_gene_8_7"]);
        // counter continues into the next block
        assert_eq!(n.next_transcript(ModelKind::Nic), "transcript8.8.nic");
    }

    #[test]
    fn exon_ids_bare_on_creation() {
        let run = || {
            let mut x = ExonIds::new();
            let mut s = String::new();
            x.write_id("8", 1, 5, Strand::Plus, &mut s);
            s.push('|');
            x.write_id("8", 1, 5, Strand::Plus, &mut s);
            s.push('|');
            x.write_id("8", 1, 5, Strand::Minus, &mut s);
            s
        };
        assert_eq!(run(), "1|8.1|2");
    }
}
