//! `IsoformIndex` builders (design §4.3): `GeneInfo.__init__` (`gene_info.py:233-280`),
//! `GeneInfo.from_models` (`:283-356`) and `GeneInfo.from_region`.
//!
//! 4.1 dict semantics are reproduced where ids can repeat: an isoform keeps the position of
//! its first insertion and the value of its last.

use crate::index::AnnotationIndex;
use iq_types::{
    junctions_from_blocks, ChrIdx, GeneIdx, GeneKey, GeneRec, IsoRec, IsoRef, IsoformIndex, Iv, Pos, ProfMatrix,
    TranscriptModel,
};
use std::collections::HashMap;

/// `IsoformIndex::from_genes(&AnnotationIndex, &[GeneIdx], delta)` (U1): 4.1
/// `GeneInfo(gene_db_list, db, delta, prepare_profiles=True)`. `genes` is the gene list in
/// the order 4.1 passes it (sorted by start, `AnnotationIndex::genes_overlapping` order).
///
/// Panics on an empty gene list (4.1 asserts it; use [`index_from_region`]).
pub fn index_from_genes(ann: &AnnotationIndex, genes: &[GeneIdx], delta: i32) -> IsoformIndex {
    assert!(!genes.is_empty(), "index_from_genes: empty gene list (use index_from_region)");
    // get_gene_region (gene_info.py:519-530)
    let (chr, iv0, _) = ann.gene_span(genes[0]);
    let (mut start, mut end) = (iv0.s, iv0.e);
    for &g in genes {
        let (_, iv, _) = ann.gene_span(g);
        start = start.min(iv.s);
        end = end.max(iv.e);
    }
    let mut idx = IsoformIndex { chr, start, end, delta, ..Default::default() };

    // gene_id_map: transcript id → last gene listing it (set_gene_ids, :569-573)
    let mut gene_of: HashMap<&str, GeneIdx> = HashMap::new();
    for &g in genes {
        for t in ann.gene_txs(g) {
            gene_of.insert(ann.tx_id(t), g);
        }
    }
    // all_isoforms_exons in insertion order (set_introns_and_exons, :600-630)
    let mut pos_of: HashMap<&str, usize> = HashMap::new();
    let mut entries: Vec<(u32, &[Iv])> = Vec::new(); // (TxIdx of first insertion, exons of last)
    for &g in genes {
        for t in ann.gene_txs(g) {
            let ex = ann.tx_exons(t);
            if ex.is_empty() {
                continue;
            }
            match pos_of.get(ann.tx_id(t)) {
                Some(&p) => entries[p].1 = ex,
                None => {
                    pos_of.insert(ann.tx_id(t), entries.len());
                    entries.push((t, ex));
                }
            }
        }
    }
    let mut names = Vec::with_capacity(entries.len());
    for &(t, ex) in &entries {
        let id = ann.tx_id(t);
        push_iso(&mut idx, IsoRef::Ref(t), GeneKey::Ref(gene_of[id]), ann.tx_strand(t), ex);
        names.push(id);
    }
    for &g in genes {
        let (_, iv, strand) = ann.gene_span(g);
        idx.genes.push(GeneRec { key: GeneKey::Ref(g), strand, region: iv });
    }
    finish_profiles(&mut idx, &names);
    idx
}

/// `GeneInfo.from_region` (no genes): an empty index over `[start, end]`.
pub fn index_from_region(chr: ChrIdx, start: Pos, end: Pos, delta: i32) -> IsoformIndex {
    let mut idx = IsoformIndex { chr, start, end, delta, ..Default::default() };
    idx.finish_derived(&[]);
    idx
}

/// `IsoformIndex::from_models(&[TranscriptModel], delta)` (U1; used by U15/U17/U19):
/// 4.1 `GeneInfo.from_models`. `chr` is the models' chromosome (`0` for an empty list).
///
/// `id_rank` follows `ModelIdx` (creation) order, the stand-in for transcript-id order
/// while ids are provisional (design §5.6); [`index_from_models_named`] ranks real ids.
pub fn index_from_models(chr: ChrIdx, models: &[TranscriptModel], delta: i32) -> IsoformIndex {
    let keys: Vec<String> = models.iter().map(|m| format!("{:010}", m.key)).collect();
    let names: Vec<&str> = keys.iter().map(String::as_str).collect();
    index_from_models_named(chr, models, delta, &names)
}

/// [`index_from_models`] with `names[i]` = the transcript id of `models[i]` (for `id_rank`).
pub fn index_from_models_named(chr: ChrIdx, models: &[TranscriptModel], delta: i32, names: &[&str]) -> IsoformIndex {
    assert_eq!(models.len(), names.len());
    let (mut start, mut end) = match models.first() {
        Some(m) => (m.start(), m.end()),
        None => (0, 0),
    };
    for m in models {
        start = start.min(m.start());
        end = end.max(m.end());
    }
    let mut idx = IsoformIndex { chr, start, end, delta, ..Default::default() };
    // dict semantics keyed by transcript id (= ModelIdx)
    let mut pos_of: HashMap<u32, usize> = HashMap::new();
    let mut last: Vec<usize> = Vec::new(); // per unique model: index of its last occurrence
    let mut first_name: Vec<&str> = Vec::new();
    for (i, m) in models.iter().enumerate() {
        match pos_of.get(&m.key) {
            Some(&p) => last[p] = i,
            None => {
                pos_of.insert(m.key, last.len());
                last.push(i);
                first_name.push(names[i]);
            }
        }
    }
    for &i in &last {
        let m = &models[i];
        push_iso(&mut idx, IsoRef::Model(m.key), m.gene, m.strand, &m.exons);
    }
    // gene_strands: first insertion order, last strand
    let mut gpos: HashMap<GeneKey, usize> = HashMap::new();
    for m in models {
        match gpos.get(&m.gene) {
            Some(&p) => idx.genes[p].strand = m.strand,
            None => {
                gpos.insert(m.gene, idx.genes.len());
                idx.genes.push(GeneRec { key: m.gene, strand: m.strand, region: Iv::default() });
            }
        }
    }
    // gene region = span of its isoforms (4.1 has no gene regions here)
    let mut seen = vec![false; idx.genes.len()];
    for r in &idx.iso {
        let p = gpos[&r.gene];
        let g = &mut idx.genes[p].region;
        *g = if seen[p] { Iv::new(g.s.min(r.span.s), g.e.max(r.span.e)) } else { r.span };
        seen[p] = true;
    }
    finish_profiles(&mut idx, &first_name);
    idx
}

fn push_iso(idx: &mut IsoformIndex, key: IsoRef, gene: GeneKey, strand: iq_types::Strand, exons: &[Iv]) {
    let introns = junctions_from_blocks(exons);
    let span = match (exons.first(), exons.last()) {
        (Some(a), Some(b)) => Iv::new(a.s, b.e),
        _ => Iv::default(),
    };
    idx.iso.push(IsoRec {
        key,
        gene,
        strand,
        span,
        exon_off: idx.iso_exons.len() as u32,
        n_exon: exons.len() as u32,
        intron_off: idx.iso_introns.len() as u32,
        n_intron: introns.len() as u32,
    });
    idx.iso_exons.extend_from_slice(exons);
    idx.iso_introns.extend_from_slice(&introns);
}

/// Feature unions, `split_exons`, the three profile matrices with their ranges, and the
/// derived fields (`set_introns_and_exons` tail, `set_junction_profiles`; D-X1: once).
fn finish_profiles(idx: &mut IsoformIndex, names: &[&str]) {
    let mut introns = idx.iso_introns.clone();
    introns.sort_unstable();
    introns.dedup();
    let mut exons = idx.iso_exons.clone();
    exons.sort_unstable();
    exons.dedup();
    idx.split_exons = split_exons(&exons);
    idx.introns = introns;
    idx.exons = exons;
    let n = idx.iso.len();
    idx.intron_prof = ProfMatrix::new(n, idx.introns.len());
    idx.exon_prof = ProfMatrix::new(n, idx.exons.len());
    idx.split_prof = ProfMatrix::new(n, idx.split_exons.len());
    idx.prof_range = Vec::with_capacity(n);
    for i in 0..n {
        let span = idx.iso[i].span;
        let (ex, intr) = (idx.iso_exons(i).to_vec(), idx.iso_introns(i).to_vec());
        let ri = set_profile(idx.intron_prof.row_mut(i), &idx.introns, &intr, span, |a, b| a == b);
        let re = set_profile(idx.exon_prof.row_mut(i), &idx.exons, &ex, span, |a, b| a == b);
        let rs = set_profile(idx.split_prof.row_mut(i), &idx.split_exons, &ex, span, |a, b| a.contains(b));
        idx.prof_range.push([ri, re, rs]);
    }
    idx.finish_derived(names);
}

/// `FeatureProfiles.set_profiles` (`gene_info.py:96-118`) into a row pre-filled with `-1`.
/// `cmp(transcript_feature, known_feature)`. Returns `profile_ranges` (`(F, 0)` when the
/// row has no `1`).
fn set_profile(row: &mut [i8], feats: &[Iv], tfeats: &[Iv], region: Iv, cmp: impl Fn(&Iv, &Iv) -> bool) -> (u32, u32) {
    for (v, f) in row.iter_mut().zip(feats) {
        if !f.overlaps(&region) {
            *v = -2;
        }
    }
    let mut pos = 0;
    for tf in tfeats {
        while pos < feats.len() && !cmp(tf, &feats[pos]) {
            pos += 1;
        }
        while pos < feats.len() && cmp(tf, &feats[pos]) {
            row[pos] = 1;
            pos += 1;
        }
    }
    let first = row.iter().position(|&v| v >= 1).unwrap_or(row.len());
    let last = row.iter().rposition(|&v| v >= 1).map_or(0, |p| p + 1);
    (first as u32, last as u32)
}

/// `GeneInfo.split_exons` (`gene_info.py:830-875`), ported literally: non-overlapping
/// covering blocks of the exon set.
pub fn split_exons(exons: &[Iv]) -> Vec<Iv> {
    let mut starts: Vec<Pos> = exons.iter().map(|e| e.s).collect();
    let mut ends: Vec<Pos> = exons.iter().map(|e| e.e).collect();
    starts.sort_unstable();
    ends.sort_unstable();
    let mut state = 0i64;
    let (mut sp, mut ep) = (0usize, 0usize);
    let mut last_border: Pos = -1;
    let mut blocks = Vec::new();
    while sp < starts.len() {
        if starts[sp] <= ends[ep] {
            if sp == 0 || starts[sp] > starts[sp - 1] {
                let cur = starts[sp];
                if last_border != -1 && state > 0 {
                    blocks.push(Iv::new(last_border, cur - 1));
                }
                last_border = cur;
            }
            state += 1;
            sp += 1;
        } else {
            if ep == 0 || ends[ep] > ends[ep - 1] {
                let cur = ends[ep];
                blocks.push(Iv::new(last_border, cur));
                last_border = cur + 1;
            }
            state -= 1;
            ep += 1;
        }
    }
    while ep < ends.len() {
        if ep == 0 || ends[ep] > ends[ep - 1] {
            let cur = ends[ep];
            blocks.push(Iv::new(last_border, cur));
            last_border = cur + 1;
        }
        ep += 1;
    }
    blocks
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(s: Pos, e: Pos) -> Iv {
        Iv::new(s, e)
    }

    #[test]
    fn split() {
        // overlapping exons [1,10] [5,20] [5,8] [30,40]
        let ex = [iv(1, 10), iv(5, 8), iv(5, 20), iv(30, 40)];
        assert_eq!(split_exons(&ex), vec![iv(1, 4), iv(5, 8), iv(9, 10), iv(11, 20), iv(30, 40)]);
        assert!(split_exons(&[]).is_empty());
    }

    #[test]
    fn profile() {
        let feats = [iv(1, 10), iv(20, 30), iv(40, 50), iv(60, 70)];
        let mut row = vec![-1i8; 4];
        let r = set_profile(&mut row, &feats, &[iv(20, 30), iv(60, 70)], iv(15, 65), |a, b| a == b);
        assert_eq!(row, vec![-2, 1, -1, 1]);
        assert_eq!(r, (1, 4));
        let mut row = vec![-1i8; 4];
        assert_eq!(set_profile(&mut row, &feats, &[iv(2, 3)], iv(2, 3), |a, b| a == b), (4, 0));
    }
}
