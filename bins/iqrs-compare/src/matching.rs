//! Structural matching of models across the two sides (IDs are never compared).
//!
//! Models of both sides live in one global index space (`ref` = `0..nr`, `test` = `nr..`).
//! - Multi-exon models are linked when `(chr, strand, intron chain)` is identical.
//! - Mono-exon models are linked across sides when on the same chr and strand with reciprocal
//!   span overlap >= `MONO_RECIP` (00_scope decision).
//!
//! Links are closed transitively with a union-find; each connected component is one "matched
//! unit". Several models with one chain (alt-end variants) collapse into one unit, and counts
//! are summed per unit and side. A model is *matched* when its unit contains a model of the
//! other side.

use crate::io::Model;
use rustc_hash::FxHashMap;

pub const MONO_RECIP: f64 = 0.8;

pub struct Uf(Vec<u32>);
impl Uf {
    pub fn new(n: usize) -> Uf {
        Uf((0..n as u32).collect())
    }
    pub fn find(&mut self, mut x: u32) -> u32 {
        while self.0[x as usize] != x {
            let p = self.0[x as usize];
            self.0[x as usize] = self.0[p as usize];
            x = p;
        }
        x
    }
    pub fn union(&mut self, a: u32, b: u32) {
        let (a, b) = (self.find(a), self.find(b));
        if a != b {
            let (lo, hi) = if a < b { (a, b) } else { (b, a) };
            self.0[hi as usize] = lo;
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct Comp {
    pub ref_count: f64,
    pub test_count: f64,
    pub n_ref: u32,
    pub n_test: u32,
    /// Global model indices, ref models first.
    pub members: Vec<u32>,
}
impl Comp {
    pub fn matched(&self) -> bool {
        self.n_ref > 0 && self.n_test > 0
    }
}

pub struct Matching {
    pub nr: usize,
    pub comp_of: Vec<u32>,
    pub comps: Vec<Comp>,
}

/// Reciprocal overlap of two closed intervals: `min(ov/len_a, ov/len_b)`.
pub fn recip(a: (i32, i32), b: (i32, i32)) -> f64 {
    let ov = (a.1.min(b.1) - a.0.max(b.0) + 1) as f64;
    if ov <= 0.0 {
        return 0.0;
    }
    let la = (a.1 - a.0 + 1) as f64;
    let lb = (b.1 - b.0 + 1) as f64;
    (ov / la).min(ov / lb)
}

/// Hash key of a multi-exon model: `[chr, strand, i1.l, i1.r, ...]`.
pub fn chain_key(m: &Model) -> Vec<i32> {
    let mut k = Vec::with_capacity(2 + 2 * m.exons.len());
    k.push(m.chr as i32);
    k.push(m.strand as i32);
    for (l, r) in m.introns() {
        k.push(l);
        k.push(r);
    }
    k
}

/// Sorted mono-exon intervals `(chr, strand, s, e, idx)` for window queries.
pub fn mono_list<'a>(models: impl Iterator<Item = (u32, &'a Model)>) -> Vec<(u32, u8, i32, i32, u32)> {
    let mut v: Vec<_> =
        models.filter(|(_, m)| m.is_mono()).map(|(i, m)| (m.chr, m.strand, m.exons[0].0, m.exons[0].1, i)).collect();
    v.sort_unstable();
    v
}

/// Entries of a [`mono_list`] that reciprocally overlap `(s, e)` by >= `MONO_RECIP`.
pub fn mono_hits(list: &[(u32, u8, i32, i32, u32)], chr: u32, strand: u8, s: i32, e: i32, mut f: impl FnMut(u32, f64)) {
    // ov >= 0.8 * len_other and len_other <= len / 0.8 bound the partner's start to
    // [s - len/4, s + len/5]; the window below is a safe superset.
    let len = e - s + 1;
    let lo = (chr, strand, s - len / 2 - 2, i32::MIN, 0u32);
    let i0 = list.partition_point(|x| *x < lo);
    for x in &list[i0..] {
        if x.0 != chr || x.1 != strand || x.2 > e {
            break;
        }
        let r = recip((s, e), (x.2, x.3));
        if r >= MONO_RECIP {
            f(x.4, r);
        }
    }
}

pub fn match_models(refm: &[Model], test: &[Model]) -> Matching {
    let nr = refm.len();
    let n = nr + test.len();
    let mut uf = Uf::new(n);
    let all = |g: usize| if g < nr { &refm[g] } else { &test[g - nr] };

    // multi-exon: identical chains
    let mut chains: FxHashMap<Vec<i32>, u32> = FxHashMap::default();
    for g in 0..n {
        let m = all(g);
        if m.exons.len() < 2 {
            continue;
        }
        match chains.entry(chain_key(m)) {
            std::collections::hash_map::Entry::Occupied(o) => uf.union(*o.get(), g as u32),
            std::collections::hash_map::Entry::Vacant(v) => {
                v.insert(g as u32);
            }
        }
    }
    drop(chains);

    // mono-exon: cross-side reciprocal overlap
    let tl = mono_list(test.iter().enumerate().map(|(i, m)| ((i + nr) as u32, m)));
    for (i, m) in refm.iter().enumerate() {
        if !m.is_mono() {
            continue;
        }
        let (s, e) = m.exons[0];
        let mut hits = Vec::new();
        mono_hits(&tl, m.chr, m.strand, s, e, |j, _| hits.push(j));
        for j in hits {
            uf.union(i as u32, j);
        }
    }

    // components
    let mut comp_of = vec![u32::MAX; n];
    let mut root_comp: FxHashMap<u32, u32> = FxHashMap::default();
    let mut comps: Vec<Comp> = Vec::new();
    for g in 0..n {
        let r = uf.find(g as u32);
        let c = *root_comp.entry(r).or_insert_with(|| {
            comps.push(Comp::default());
            (comps.len() - 1) as u32
        });
        comp_of[g] = c;
        let cp = &mut comps[c as usize];
        cp.members.push(g as u32);
        let m = all(g);
        if g < nr {
            cp.n_ref += 1;
            cp.ref_count += m.count;
        } else {
            cp.n_test += 1;
            cp.test_count += m.count;
        }
    }
    Matching { nr, comp_of, comps }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m(id: &str, strand: u8, exons: &[(i32, i32)], count: f64) -> Model {
        Model {
            id: id.into(),
            gene: String::new(),
            chr: 0,
            strand,
            exons: exons.to_vec(),
            source: "IsoQuant".into(),
            count,
        }
    }

    #[test]
    fn recip_overlap() {
        assert_eq!(recip((1, 100), (1, 100)), 1.0);
        assert!((recip((1, 100), (21, 100)) - 0.8).abs() < 1e-12);
        assert!(recip((1, 100), (22, 100)) < 0.8);
        assert_eq!(recip((1, 10), (11, 20)), 0.0);
    }

    #[test]
    fn chains_and_monos() {
        let r = vec![
            m("a", b'+', &[(1, 10), (20, 30), (40, 50)], 5.0),
            m("a2", b'+', &[(3, 10), (20, 30), (40, 60)], 1.0), // same chain, other ends
            m("b", b'+', &[(1, 10), (21, 30)], 2.0),
            m("mono", b'+', &[(1000, 1099)], 3.0),
            m("mono_minus", b'-', &[(2000, 2099)], 3.0),
        ];
        let t = vec![
            m("x", b'+', &[(5, 10), (20, 30), (40, 45)], 7.0),
            m("y", b'+', &[(1, 11), (21, 30)], 2.0), // intron shifted by 1
            m("z", b'+', &[(1019, 1099)], 4.0),      // recip 0.81
            m("w", b'+', &[(2000, 2099)], 1.0),      // wrong strand
        ];
        let mt = match_models(&r, &t);
        let matched = |g: usize| mt.comps[mt.comp_of[g] as usize].matched();
        assert!(matched(0) && matched(1) && matched(5));
        assert!(!matched(2) && !matched(6));
        assert!(matched(3) && matched(7));
        assert!(!matched(4) && !matched(8));
        let u = &mt.comps[mt.comp_of[0] as usize];
        assert_eq!((u.ref_count, u.test_count, u.n_ref, u.n_test), (6.0, 7.0, 2, 1));
    }
}
