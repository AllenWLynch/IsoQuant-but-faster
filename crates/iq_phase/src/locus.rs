//! Loci and read equivalence classes.
//!
//! Model reads (phased_model/code/analyze.py): final assignment type unique,
//! unique_minor_difference or ambiguous, with a non-empty compatible set C_r (final model
//! ids). Inconsistent, intergenic, noninformative and unassigned reads are excluded. A read
//! name is used once (first occurrence).
//!
//! Unit = gene; genes that share a read are merged into one locus (connected component of the
//! gene graph). The locus transcripts are the transcripts seen in its reads, sorted by id.
//!
//! Classes (gibbs.py `build_unit`): key = (sorted C_r indices, log-lambda ratio d clipped to
//! +-30 and rounded half-to-even to a 0.25-nat bin). Per class: read count, and over reads
//! with d != 0 the sums of l_r(A) = 1/(1 + e^-d) and l_r(B) = 1 - l_r(A) with the exact d.

use rustc_hash::{FxHashMap, FxHashSet};

pub const BIN: f64 = 0.25;
pub const DCLIP: f64 = 30.0;

/// One locus, classes in flat arrays (sorted by key).
#[derive(Clone, Debug, Default)]
pub struct Locus {
    pub chr: String,
    /// Gene ids, sorted.
    pub genes: Vec<String>,
    /// Transcript ids, sorted; class transcript indices refer to this.
    pub tx: Vec<String>,
    /// Gene id per transcript.
    pub tx_gene: Vec<String>,
    /// Start of each class's transcripts in `cls_tx`; one extra entry at the end.
    pub cls_off: Vec<u32>,
    pub cls_tx: Vec<u32>,
    /// Binned log lambda(A) - log lambda(B).
    pub cls_d: Vec<f64>,
    pub cls_n: Vec<u32>,
    pub cls_la: Vec<f64>,
    pub cls_lb: Vec<f64>,
}

impl Locus {
    pub fn n_classes(&self) -> usize {
        self.cls_n.len()
    }
    pub fn n_reads(&self) -> u64 {
        self.cls_n.iter().map(|&n| n as u64).sum()
    }
    pub fn class_tx(&self, c: usize) -> &[u32] {
        &self.cls_tx[self.cls_off[c] as usize..self.cls_off[c + 1] as usize]
    }
    /// Locus name for the gene table: gene ids joined by ','.
    pub fn name(&self) -> String {
        self.genes.join(",")
    }

    /// Builds one locus from reads `(C_r as transcript ids, d)` (tests, TSV replay).
    pub fn from_reads<'a>(
        chr: &str,
        reads: impl IntoIterator<Item = (Vec<&'a str>, f64)>,
        gene_of: impl Fn(&str) -> String,
    ) -> Locus {
        let mut acc = ChrAcc::default();
        for (c, d) in reads {
            let tx: Vec<(String, String)> = c.iter().map(|t| (t.to_string(), gene_of(t))).collect();
            acc.add_read(None, tx.iter().map(|(t, g)| (t.as_str(), g.as_str())), d);
        }
        let mut v = acc.build(chr);
        // all reads in one locus regardless of gene links
        if v.len() == 1 {
            return v.pop().unwrap();
        }
        merge_loci(chr, v)
    }
}

fn merge_loci(chr: &str, v: Vec<Locus>) -> Locus {
    let mut reads: Vec<(Vec<String>, f64, f64, f64, u32)> = Vec::new();
    let mut tg: FxHashMap<String, String> = FxHashMap::default();
    for l in &v {
        for (t, g) in l.tx.iter().zip(&l.tx_gene) {
            tg.insert(t.clone(), g.clone());
        }
        for c in 0..l.n_classes() {
            let tx = l.class_tx(c).iter().map(|&i| l.tx[i as usize].clone()).collect();
            reads.push((tx, l.cls_d[c], l.cls_la[c], l.cls_lb[c], l.cls_n[c]));
        }
    }
    let mut txs: Vec<String> = tg.keys().cloned().collect();
    txs.sort();
    let tix: FxHashMap<&str, u32> = txs.iter().enumerate().map(|(i, t)| (t.as_str(), i as u32)).collect();
    let mut genes: Vec<String> = tg.values().cloned().collect();
    genes.sort();
    genes.dedup();
    let mut cls: Vec<(Vec<u32>, f64, u32, f64, f64)> = reads
        .into_iter()
        .map(|(tx, d, la, lb, n)| {
            let mut k: Vec<u32> = tx.iter().map(|t| tix[t.as_str()]).collect();
            k.sort_unstable();
            (k, d, n, la, lb)
        })
        .collect();
    cls.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.total_cmp(&b.1)));
    let mut l = Locus {
        chr: chr.to_string(),
        genes,
        tx_gene: txs.iter().map(|t| tg[t].clone()).collect(),
        tx: txs,
        cls_off: vec![0],
        ..Default::default()
    };
    for (k, d, n, la, lb) in cls {
        l.cls_tx.extend(k);
        l.cls_off.push(l.cls_tx.len() as u32);
        l.cls_d.push(d);
        l.cls_n.push(n);
        l.cls_la.push(la);
        l.cls_lb.push(lb);
    }
    l
}

/// Model reads of one chromosome, accumulated block by block.
#[derive(Default, Debug)]
pub struct ChrAcc {
    seen: FxHashSet<u64>,
    tx_ix: FxHashMap<String, u32>,
    tx_names: Vec<String>,
    tx_gene: Vec<u32>,
    gene_ix: FxHashMap<String, u32>,
    gene_names: Vec<String>,
    read_off: Vec<u32>,
    read_tx: Vec<u32>,
    read_d: Vec<f64>,
    scratch: Vec<u32>,
}

impl ChrAcc {
    pub fn new() -> Self {
        ChrAcc { read_off: vec![0], ..Default::default() }
    }

    /// Adds a read: its compatible transcripts `(transcript id, gene id)` and d. `hash`: read
    /// name hash for de-duplication (`None` = always add). Returns false if skipped.
    pub fn add_read<'a>(
        &mut self,
        hash: Option<u64>,
        tx: impl IntoIterator<Item = (&'a str, &'a str)>,
        d: f64,
    ) -> bool {
        if self.read_off.is_empty() {
            self.read_off.push(0);
        }
        if let Some(h) = hash {
            if self.seen.contains(&h) {
                return false;
            }
        }
        self.scratch.clear();
        for (t, g) in tx {
            let ti = match self.tx_ix.get(t) {
                Some(&i) => i,
                None => {
                    let gi = match self.gene_ix.get(g) {
                        Some(&i) => i,
                        None => {
                            let i = self.gene_names.len() as u32;
                            self.gene_ix.insert(g.to_string(), i);
                            self.gene_names.push(g.to_string());
                            i
                        }
                    };
                    let i = self.tx_names.len() as u32;
                    self.tx_ix.insert(t.to_string(), i);
                    self.tx_names.push(t.to_string());
                    self.tx_gene.push(gi);
                    i
                }
            };
            if !self.scratch.contains(&ti) {
                self.scratch.push(ti);
            }
        }
        if self.scratch.is_empty() {
            return false;
        }
        if let Some(h) = hash {
            self.seen.insert(h);
        }
        self.read_tx.extend_from_slice(&self.scratch);
        self.read_off.push(self.read_tx.len() as u32);
        self.read_d.push(d);
        true
    }

    pub fn n_reads(&self) -> usize {
        self.read_d.len()
    }

    /// Splits the reads into loci (connected components of genes linked by reads) and builds
    /// their classes. Loci are ordered by their first gene id.
    pub fn build(self, chr: &str) -> Vec<Locus> {
        let ng = self.gene_names.len();
        let mut parent: Vec<u32> = (0..ng as u32).collect();
        fn find(p: &mut [u32], mut x: u32) -> u32 {
            while p[x as usize] != x {
                p[x as usize] = p[p[x as usize] as usize];
                x = p[x as usize];
            }
            x
        }
        let nr = self.read_d.len();
        let rtx = |r: usize| &self.read_tx[self.read_off[r] as usize..self.read_off[r + 1] as usize];
        for r in 0..nr {
            let t = rtx(r);
            let g0 = self.tx_gene[t[0] as usize];
            for &ti in &t[1..] {
                let (a, b) = (find(&mut parent, g0), find(&mut parent, self.tx_gene[ti as usize]));
                if a != b {
                    parent[a.max(b) as usize] = a.min(b);
                }
            }
        }
        // locus id per root
        let mut lid = vec![u32::MAX; ng];
        let mut n_loci = 0u32;
        let mut gene_locus = vec![0u32; ng];
        for g in 0..ng as u32 {
            let r = find(&mut parent, g) as usize;
            if lid[r] == u32::MAX {
                lid[r] = n_loci;
                n_loci += 1;
            }
            gene_locus[g as usize] = lid[r];
        }
        // reads per locus
        let mut by_locus: Vec<Vec<u32>> = vec![Vec::new(); n_loci as usize];
        for r in 0..nr {
            let g = self.tx_gene[rtx(r)[0] as usize];
            by_locus[gene_locus[g as usize] as usize].push(r as u32);
        }
        let mut out: Vec<Locus> = Vec::with_capacity(n_loci as usize);
        let mut local: FxHashMap<u32, u32> = FxHashMap::default();
        let mut classes: FxHashMap<(Vec<u32>, i32), (u32, f64, f64)> = FxHashMap::default();
        let mut key = Vec::new();
        for reads in by_locus {
            if reads.is_empty() {
                continue;
            }
            // transcripts of the locus, sorted by id
            let mut txs: Vec<u32> = reads.iter().flat_map(|&r| rtx(r as usize).iter().copied()).collect();
            txs.sort_unstable();
            txs.dedup();
            txs.sort_by(|a, b| self.tx_names[*a as usize].cmp(&self.tx_names[*b as usize]));
            local.clear();
            for (i, &t) in txs.iter().enumerate() {
                local.insert(t, i as u32);
            }
            classes.clear();
            for &r in &reads {
                key.clear();
                key.extend(rtx(r as usize).iter().map(|t| local[t]));
                key.sort_unstable();
                let d = self.read_d[r as usize];
                let (bin, la, ev) = bin_of(d);
                let e = classes.entry((key.clone(), bin)).or_insert((0, 0.0, 0.0));
                e.0 += 1;
                if ev {
                    e.1 += la;
                    e.2 += 1.0 - la;
                }
            }
            let mut keys: Vec<_> = classes.drain().collect();
            keys.sort_by(|a, b| a.0.cmp(&b.0));
            let mut genes: Vec<String> =
                txs.iter().map(|&t| self.gene_names[self.tx_gene[t as usize] as usize].clone()).collect();
            genes.sort();
            genes.dedup();
            let mut l = Locus {
                chr: chr.to_string(),
                genes,
                tx: txs.iter().map(|&t| self.tx_names[t as usize].clone()).collect(),
                tx_gene: txs.iter().map(|&t| self.gene_names[self.tx_gene[t as usize] as usize].clone()).collect(),
                cls_off: Vec::with_capacity(keys.len() + 1),
                ..Default::default()
            };
            l.cls_off.push(0);
            for ((k, bin), (n, la, lb)) in keys {
                l.cls_tx.extend_from_slice(&k);
                l.cls_off.push(l.cls_tx.len() as u32);
                l.cls_d.push(bin as f64 * BIN);
                l.cls_n.push(n);
                l.cls_la.push(la);
                l.cls_lb.push(lb);
            }
            out.push(l);
        }
        out.sort_by(|a, b| a.genes[0].cmp(&b.genes[0]));
        out
    }
}

/// (bin index, l_r(A), has evidence) of a read's d.
#[inline]
pub fn bin_of(d: f64) -> (i32, f64, bool) {
    let la = 1.0 / (1.0 + (-d.clamp(-700.0, 700.0)).exp());
    let dc = d.clamp(-DCLIP, DCLIP);
    ((dc / BIN).round_ties_even() as i32, la, d != 0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loci_are_gene_components() {
        let mut acc = ChrAcc::new();
        acc.add_read(Some(1), [("t1", "g1")], 0.0);
        acc.add_read(Some(2), [("t2", "g1")], 1.1);
        acc.add_read(Some(2), [("t2", "g1")], 1.1); // duplicate name
        acc.add_read(Some(3), [("t3", "g2"), ("t4", "g3")], -0.2);
        acc.add_read(Some(4), [("t4", "g3")], 0.0);
        acc.add_read(Some(5), [("t9", "g9")], 0.125); // bin rounds half to even -> 0
        let l = acc.build("chr1");
        assert_eq!(l.len(), 3);
        assert_eq!(l[0].genes, vec!["g1"]);
        assert_eq!(l[0].n_reads(), 2);
        assert_eq!(l[1].genes, vec!["g2", "g3"]);
        assert_eq!(l[1].tx, vec!["t3", "t4"]);
        assert_eq!(l[1].n_classes(), 2);
        assert_eq!(l[1].class_tx(0), &[0, 1]);
        assert_eq!(l[1].cls_d[0], -0.25);
        assert_eq!(l[2].cls_d[0], 0.0);
        assert!(l[2].cls_la[0] > 0.5); // exact d used for l(A)
        assert_eq!(bin_of(0.375).0, 2); // 1.5 -> 2
        assert_eq!(bin_of(99.0).0, 120);
    }
}
