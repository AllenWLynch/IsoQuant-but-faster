//! Phased heterozygous SNVs from a (bgzipped) VCF.
//!
//! Conventions (phased_model/code/het_evidence.py `het_sites`): one sample; keep records whose
//! GT is phased (`|`), diploid, without missing alleles, heterozygous, with both GT alleles of
//! length 1 (SNVs; indels skipped). Haplotype A = the first GT allele, B = the second. Phase
//! blocks (PS/PB/PI) are ignored: GT is taken as one genome-wide phase (GTEx v8). A repeated
//! position keeps the last record (the prototype's dict assignment).
//!
//! The file is streamed once (bgzip or plain text); the index is not needed for that.

use rustc_hash::FxHashMap;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// Phased het SNVs of one chromosome, sorted by position.
#[derive(Clone, Debug, Default)]
pub struct ChrSites {
    /// 0-based positions, strictly increasing.
    pub pos: Vec<u32>,
    /// Haplotype-A allele (upper case ASCII).
    pub a: Vec<u8>,
    /// Haplotype-B allele.
    pub b: Vec<u8>,
}

impl ChrSites {
    pub fn len(&self) -> usize {
        self.pos.len()
    }
    pub fn is_empty(&self) -> bool {
        self.pos.is_empty()
    }
    /// Index of the first site at or after `pos0`.
    #[inline]
    pub fn lower_bound(&self, pos0: u32) -> usize {
        self.pos.partition_point(|&p| p < pos0)
    }
}

/// All chromosomes' sites.
#[derive(Clone, Debug, Default)]
pub struct HetSites {
    pub chrs: FxHashMap<String, ChrSites>,
    /// Sample column used.
    pub sample: String,
}

impl HetSites {
    pub fn get(&self, chr: &str) -> Option<&ChrSites> {
        self.chrs.get(chr).filter(|c| !c.is_empty())
    }
    pub fn n_sites(&self) -> usize {
        self.chrs.values().map(|c| c.len()).sum()
    }
}

fn open(path: &Path) -> Result<Box<dyn BufRead + Send>, String> {
    let f = std::fs::File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut magic = [0u8; 2];
    let mut f2 = std::fs::File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let n = f2.read(&mut magic).map_err(|e| format!("{}: {e}", path.display()))?;
    if n == 2 && magic == [0x1f, 0x8b] {
        let workers = std::num::NonZero::new(4).unwrap();
        Ok(Box::new(BufReader::with_capacity(
            1 << 20,
            noodles_bgzf::io::MultithreadedReader::with_worker_count(workers, f),
        )))
    } else {
        Ok(Box::new(BufReader::with_capacity(1 << 20, f)))
    }
}

/// `GT` of one sample field: `Some((i, j))` for a phased diploid call without missing alleles.
fn phased_gt(gt: &[u8]) -> Option<(usize, usize)> {
    let bar = gt.iter().position(|&c| c == b'|')?;
    if gt.contains(&b'/') {
        return None;
    }
    let parse = |s: &[u8]| -> Option<usize> {
        if s.is_empty() || !s.iter().all(u8::is_ascii_digit) {
            return None;
        }
        std::str::from_utf8(s).ok()?.parse().ok()
    };
    let rest = &gt[bar + 1..];
    if rest.contains(&b'|') {
        return None; // not diploid
    }
    Some((parse(&gt[..bar])?, parse(rest)?))
}

/// Parses the VCF. `sample`: the column to use (`None` = the single sample of the file).
pub fn load(path: &Path, sample: Option<&str>) -> Result<HetSites, String> {
    let mut r = open(path)?;
    let ctx = |e: std::io::Error| format!("{}: {e}", path.display());
    let mut line = Vec::with_capacity(1 << 12);
    let mut col: Option<usize> = None;
    let mut out = HetSites::default();
    let mut cur_chr: Vec<u8> = Vec::new();
    let mut cur: Option<(String, ChrSites)> = None;
    loop {
        line.clear();
        if r.read_until(b'\n', &mut line).map_err(ctx)? == 0 {
            break;
        }
        while line.last().is_some_and(|&c| c == b'\n' || c == b'\r') {
            line.pop();
        }
        if line.is_empty() || line.starts_with(b"##") {
            continue;
        }
        if line.starts_with(b"#CHROM") {
            let names: Vec<&[u8]> = line.split(|&c| c == b'\t').skip(9).collect();
            let i = match sample {
                Some(s) => names
                    .iter()
                    .position(|n| *n == s.as_bytes())
                    .ok_or_else(|| format!("{}: sample {s:?} not in the VCF header", path.display()))?,
                None if names.len() == 1 => 0,
                None => {
                    return Err(format!(
                        "{}: {} samples in the VCF; choose one with --phase-sample",
                        path.display(),
                        names.len()
                    ))
                }
            };
            out.sample = String::from_utf8_lossy(names[i]).into_owned();
            col = Some(9 + i);
            continue;
        }
        let Some(col) = col else {
            return Err(format!("{}: no #CHROM header line before records", path.display()));
        };
        let mut f = line.split(|&c| c == b'\t');
        let chrom = f.next().unwrap_or_default();
        let pos = f.next().unwrap_or_default();
        let _id = f.next();
        let refa = f.next().unwrap_or_default();
        let alt = f.next().unwrap_or_default();
        let Some(format) = f.nth(3) else { continue };
        let Some(sample_f) = f.nth(col - 9) else { continue };
        let gt_i = if format.starts_with(b"GT") {
            0
        } else {
            match format.split(|&c| c == b':').position(|k| k == b"GT") {
                Some(i) => i,
                None => continue,
            }
        };
        let gt = sample_f.split(|&c| c == b':').nth(gt_i).unwrap_or_default();
        let Some((g0, g1)) = phased_gt(gt) else { continue };
        if g0 == g1 {
            continue;
        }
        let allele = |i: usize| if i == 0 { Some(refa) } else { alt.split(|&c| c == b',').nth(i - 1) };
        let (Some(a), Some(b)) = (allele(g0), allele(g1)) else { continue };
        if a.len() != 1 || b.len() != 1 {
            continue;
        }
        let Some(p1) = std::str::from_utf8(pos).ok().and_then(|s| s.parse::<u64>().ok()) else {
            return Err(format!("{}: bad POS in line {:?}", path.display(), String::from_utf8_lossy(&line)));
        };
        if p1 == 0 || p1 > u32::MAX as u64 {
            continue;
        }
        let p0 = (p1 - 1) as u32;
        if chrom != cur_chr.as_slice() {
            if let Some((n, s)) = cur.take() {
                merge_chr(&mut out, n, s);
            }
            cur_chr = chrom.to_vec();
            cur = Some((String::from_utf8_lossy(chrom).into_owned(), ChrSites::default()));
        }
        let s = &mut cur.as_mut().unwrap().1;
        s.pos.push(p0);
        s.a.push(a[0].to_ascii_uppercase());
        s.b.push(b[0].to_ascii_uppercase());
    }
    if let Some((n, s)) = cur.take() {
        merge_chr(&mut out, n, s);
    }
    for s in out.chrs.values_mut() {
        normalize(s);
    }
    Ok(out)
}

fn merge_chr(out: &mut HetSites, name: String, s: ChrSites) {
    let e = out.chrs.entry(name).or_default();
    e.pos.extend(s.pos);
    e.a.extend(s.a);
    e.b.extend(s.b);
}

/// Sorts by position (stable) and keeps the last record of a repeated position.
fn normalize(s: &mut ChrSites) {
    let mut idx: Vec<usize> = (0..s.pos.len()).collect();
    idx.sort_by_key(|&i| s.pos[i]);
    let mut out = ChrSites::default();
    for (k, &i) in idx.iter().enumerate() {
        if idx.get(k + 1).is_some_and(|&j| s.pos[j] == s.pos[i]) {
            continue;
        }
        out.pos.push(s.pos[i]);
        out.a.push(s.a[i]);
        out.b.push(s.b[i]);
    }
    *s = out;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_small_vcf() {
        let dir = std::env::temp_dir().join(format!("iq_phase_vcf_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("t.vcf");
        std::fs::write(
            &p,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
             chr1\t10\t.\tG\tA\t.\tPASS\t.\tGT:PW\t0|1:0|1\n\
             chr1\t20\t.\tC\tT\t.\tPASS\t.\tGT\t1|0\n\
             chr1\t30\t.\tC\tT\t.\tPASS\t.\tGT\t1|1\n\
             chr1\t40\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\n\
             chr1\t50\t.\tCA\tC\t.\tPASS\t.\tGT\t0|1\n\
             chr1\t60\t.\tC\tA,g\t.\tPASS\t.\tGT\t2|1\n\
             chr1\t70\t.\tC\tT\t.\tPASS\t.\tGT\t.|1\n\
             chr2\t5\t.\ta\tc\t.\tPASS\t.\tGT\t0|1\n",
        )
        .unwrap();
        let s = load(&p, None).unwrap();
        let c1 = &s.chrs["chr1"];
        assert_eq!(c1.pos, vec![9, 19, 59]);
        assert_eq!(c1.a, b"GTG".to_vec());
        assert_eq!(c1.b, b"ACA".to_vec());
        assert_eq!(s.chrs["chr2"].a, b"A".to_vec());
        assert_eq!(s.sample, "S1");
        assert!(load(&p, Some("nope")).is_err());
    }
}
