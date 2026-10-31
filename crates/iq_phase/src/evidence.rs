//! Per-read het-site evidence, collected inline during the read pass.
//!
//! Conventions (phased_model/code/het_evidence.py): at every phased het SNV a read aligns a
//! base to (CIGAR M/=/X; deletions and introns do not cover), the read base is compared with
//! the haplotype alleles: equal to A -> A call, equal to B -> B call, anything else skipped.
//! Each call keeps its base quality; eps = 10^(-Q/10) clamped to [1e-3, 0.25] (Q = 20 when the
//! record has no qualities). log lambda(A) - log lambda(B) = sum over calls of
//! +-[log(1 - eps) - log(eps)].
//!
//! Storage: per read (name hash) a run of calls, each 5 bytes (site index u32 + code u8 with
//! bit 7 = B allele, bits 0-6 = min(Q, 127)). Only reads with at least one A/B call are kept.

use crate::vcf::ChrSites;
use iq_types::RecordLite;

pub const EPS_MIN: f64 = 1e-3;
pub const EPS_MAX: f64 = 0.25;
/// Quality used when a record has no quality string (pysam `query_qualities is None`).
pub const DEFAULT_Q: u8 = 20;

/// `log(1 - eps) - log(eps)` per quality code 0..=127.
fn weight_table() -> &'static [f64; 128] {
    static T: std::sync::OnceLock<[f64; 128]> = std::sync::OnceLock::new();
    T.get_or_init(|| {
        let mut t = [0.0; 128];
        for (q, w) in t.iter_mut().enumerate() {
            let eps = 10f64.powf(-(q as f64) / 10.0).clamp(EPS_MIN, EPS_MAX);
            *w = (-eps).ln_1p() - eps.ln();
        }
        t
    })
}

/// log lambda(A) - log lambda(B) of a run of call codes.
#[inline]
pub fn delta_of(codes: &[u8]) -> f64 {
    let t = weight_table();
    let mut d = 0.0;
    for &c in codes {
        let w = t[(c & 0x7f) as usize];
        if c & 0x80 == 0 {
            d += w;
        } else {
            d -= w;
        }
    }
    d
}

/// Growing evidence of one block (fetch order).
#[derive(Default, Debug)]
pub struct EvidenceBuf {
    hashes: Vec<u64>,
    /// Start of each read's calls in `site`/`code`; one extra entry at the end.
    off: Vec<u32>,
    site: Vec<u32>,
    code: Vec<u8>,
}

impl EvidenceBuf {
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds the calls of one record. Secondary, supplementary and unmapped records are skipped
    /// (the read pass keeps primary alignments only). `qual`: raw Phred qualities (empty or
    /// starting with 0xff = absent). Returns the number of A/B calls.
    pub fn push(&mut self, sites: &ChrSites, rec: &RecordLite<'_>, qual: &[u8]) -> usize {
        use iq_types::read::flags;
        if rec.flag & (flags::UNMAPPED | flags::SECONDARY | flags::SUPPLEMENTARY) != 0 || rec.seq.is_empty() {
            return 0;
        }
        if rec.pos0 < 0 || rec.pos0 > u32::MAX as i64 {
            return 0;
        }
        let n = sites.pos.len();
        let mut si = sites.lower_bound(rec.pos0 as u32);
        if si >= n {
            return 0;
        }
        let no_qual = qual.is_empty() || qual[0] == 0xff;
        let start = self.site.len();
        let mut rp = rec.pos0 as u64;
        let mut qp = 0usize;
        for &w in rec.cigar {
            let len = (w >> 4) as u64;
            match w & 0xf {
                0 | 7 | 8 => {
                    let end = rp + len;
                    while si < n && (sites.pos[si] as u64) < end {
                        let p = sites.pos[si] as u64;
                        if p >= rp {
                            let q = qp + (p - rp) as usize;
                            if let Some(&base) = rec.seq.get(q) {
                                let base = base.to_ascii_uppercase();
                                let allele_b = if base == sites.a[si] {
                                    Some(false)
                                } else if base == sites.b[si] {
                                    Some(true)
                                } else {
                                    None
                                };
                                if let Some(isb) = allele_b {
                                    let qv =
                                        if no_qual { DEFAULT_Q } else { qual.get(q).copied().unwrap_or(DEFAULT_Q) };
                                    self.site.push(si as u32);
                                    self.code.push(((isb as u8) << 7) | qv.min(127));
                                }
                            }
                        }
                        si += 1;
                    }
                    rp = end;
                    qp += len as usize;
                }
                1 | 4 => qp += len as usize,
                2 | 3 => {
                    rp += len;
                    while si < n && (sites.pos[si] as u64) < rp {
                        si += 1;
                    }
                }
                _ => {}
            }
            if si >= n {
                break;
            }
        }
        let k = self.site.len() - start;
        if k > 0 {
            if self.off.is_empty() {
                self.off.push(start as u32);
            }
            self.hashes.push(iq_types::name_hash(rec.qname));
            self.off.push(self.site.len() as u32);
        }
        k
    }

    /// Number of reads with calls so far.
    pub fn n_reads(&self) -> usize {
        self.hashes.len()
    }

    /// Sorts by name hash and keeps the reads `keep` accepts (first record of a repeated hash).
    pub fn finish(self, keep: impl Fn(u64) -> bool) -> BlockEvidence {
        let mut idx: Vec<u32> = (0..self.hashes.len() as u32).filter(|&i| keep(self.hashes[i as usize])).collect();
        idx.sort_by_key(|&i| (self.hashes[i as usize], i));
        idx.dedup_by_key(|i| self.hashes[*i as usize]);
        let mut out = BlockEvidence { hashes: Vec::with_capacity(idx.len()), off: vec![0], ..Default::default() };
        for &i in &idx {
            let (s, e) = (self.off[i as usize] as usize, self.off[i as usize + 1] as usize);
            out.hashes.push(self.hashes[i as usize]);
            out.site.extend_from_slice(&self.site[s..e]);
            out.code.extend_from_slice(&self.code[s..e]);
            out.off.push(out.site.len() as u32);
        }
        out.hashes.shrink_to_fit();
        out.site.shrink_to_fit();
        out.code.shrink_to_fit();
        out
    }
}

/// Finished evidence of one block, looked up by read-name hash.
#[derive(Default, Debug, Clone)]
pub struct BlockEvidence {
    hashes: Vec<u64>,
    off: Vec<u32>,
    /// Site index (into the chromosome's [`ChrSites`]) per call.
    pub site: Vec<u32>,
    pub code: Vec<u8>,
}

impl BlockEvidence {
    pub fn n_reads(&self) -> usize {
        self.hashes.len()
    }
    pub fn bytes(&self) -> usize {
        self.hashes.len() * 8 + self.off.len() * 4 + self.site.len() * 5
    }
    /// Call codes of a read (`None` = no A/B call).
    pub fn calls(&self, hash: u64) -> Option<&[u8]> {
        let i = self.hashes.binary_search(&hash).ok()?;
        Some(&self.code[self.off[i] as usize..self.off[i + 1] as usize])
    }
    /// log lambda(A) - log lambda(B) of a read (0 without calls).
    pub fn delta(&self, hash: u64) -> f64 {
        self.calls(hash).map_or(0.0, delta_of)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cig(ops: &[(u32, u32)]) -> Vec<u32> {
        ops.iter().map(|&(l, o)| (l << 4) | o).collect()
    }

    #[test]
    fn calls_follow_cigar() {
        // sites at 0-based 102 (A=G,B=T), 110 (A=C,B=A; in an intron), 205 (A=A,B=G)
        let sites = ChrSites { pos: vec![102, 110, 205], a: b"GCA".to_vec(), b: b"TAG".to_vec() };
        // read at 100: 5M (100-104) 2I 3S? no: 5M 2I 95N 10M (200-209)
        let cigar = cig(&[(5, 0), (2, 1), (95, 3), (10, 0)]);
        // query: 0-4 aligned 100-104, 5-6 insertion, 7-16 aligned 200-209
        let mut seq = b"NNNNNNNNNNNNNNNNN".to_vec();
        seq[2] = b'T'; // pos 102 -> B
        seq[7 + 5] = b'A'; // pos 205 -> A
        let qual = vec![30u8; 17];
        let rec = RecordLite { qname: b"r1", flag: 0, pos0: 100, mapq: 60, cigar: &cigar, seq: &seq, tag: None };
        let mut b = EvidenceBuf::new();
        assert_eq!(b.push(&sites, &rec, &qual), 2);
        let ev = b.finish(|_| true);
        let c = ev.calls(iq_types::name_hash(b"r1")).unwrap();
        assert_eq!(c, &[0x80 | 30, 30]);
        assert_eq!(ev.site, vec![0, 2]);
        assert!(ev.delta(iq_types::name_hash(b"r1")).abs() < 1e-12);
        // a third base is skipped; missing qualities -> Q20
        let mut seq2 = seq.clone();
        seq2[2] = b'C';
        let rec = RecordLite { qname: b"r2", seq: &seq2, ..rec };
        let mut b = EvidenceBuf::new();
        assert_eq!(b.push(&sites, &rec, &[0xff]), 1);
        let ev = b.finish(|_| true);
        let d = ev.delta(iq_types::name_hash(b"r2"));
        assert!((d - ((1.0f64 - 0.01).ln() - 0.01f64.ln())).abs() < 1e-12);
    }

    #[test]
    fn eps_clamped() {
        let t = weight_table();
        assert!((t[60] - ((1.0f64 - 1e-3).ln() - 1e-3f64.ln())).abs() < 1e-12);
        assert!((t[0] - ((0.75f64).ln() - 0.25f64.ln())).abs() < 1e-12);
    }
}
