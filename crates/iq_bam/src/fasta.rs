//! mmap'd FASTA + `.fai`: line-length arithmetic, no case folding.
//!
//! Equals pyfaidx `Fasta(fa)[chr][s-1:e]` (default options: no case change, python slice
//! clipping) for `1 <= s`. A slice inside one FASTA line is borrowed from the map; a slice
//! across lines is copied without the line terminators.

use iq_types::{Iv, RefSlice};
use memmap2::Mmap;
use std::borrow::Cow;
use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};

#[derive(Clone, Debug)]
struct FaiEntry {
    name: String,
    len: u64,
    offset: u64,
    line_bases: u64,
    line_bytes: u64,
}

#[derive(Debug)]
pub struct FastaIndex {
    map: Mmap,
    entries: Vec<FaiEntry>,
}

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

impl FastaIndex {
    /// Opens `fa` and its `fa.fai`.
    pub fn open(fa: &Path) -> io::Result<FastaIndex> {
        let fai = PathBuf::from(format!("{}.fai", fa.display()));
        if fa.extension().is_some_and(|x| x == "gz" || x == "bgz") {
            return Err(invalid(format!(
                "{}: compressed references are not supported; decompress it first",
                fa.display()
            )));
        }
        let txt = std::fs::read_to_string(&fai).map_err(|e| {
            if e.kind() == io::ErrorKind::NotFound && fa.exists() {
                io::Error::new(
                    e.kind(),
                    format!(
                        "reference index {} not found; create it with `samtools faidx {}`",
                        fai.display(),
                        fa.display()
                    ),
                )
            } else {
                io::Error::new(e.kind(), format!("{}: {e}", fai.display()))
            }
        })?;
        let f = File::open(fa).map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", fa.display())))?;
        // SAFETY: the reference is opened read-only and is not modified while mapped.
        let map = unsafe { Mmap::map(&f)? };
        let mut entries = Vec::new();
        for (k, line) in txt.lines().enumerate() {
            if line.is_empty() {
                continue;
            }
            let c: Vec<&str> = line.split('\t').collect();
            let num = |i: usize| -> io::Result<u64> {
                c.get(i)
                    .and_then(|v| v.parse().ok())
                    .ok_or_else(|| invalid(format!("{}:{}: bad fai line", fai.display(), k + 1)))
            };
            let e = FaiEntry {
                name: c[0].to_owned(),
                len: num(1)?,
                offset: num(2)?,
                line_bases: num(3)?,
                line_bytes: num(4)?,
            };
            if e.len > 0 && (e.line_bases == 0 || e.line_bytes < e.line_bases) {
                return Err(invalid(format!("{}:{}: bad line lengths", fai.display(), k + 1)));
            }
            let last = e.offset
                + (e.len.saturating_sub(1) / e.line_bases.max(1)) * e.line_bytes
                + (e.len.saturating_sub(1) % e.line_bases.max(1));
            if e.len > 0 && last >= map.len() as u64 {
                return Err(invalid(format!("{}: {} extends past the end of {}", fai.display(), e.name, fa.display())));
            }
            entries.push(e);
        }
        Ok(FastaIndex { map, entries })
    }

    fn entry(&self, chr: &str) -> Option<&FaiEntry> {
        self.entries.iter().find(|e| e.name == chr)
    }

    /// Sequence names in `.fai` order.
    pub fn names(&self) -> impl Iterator<Item = &str> + '_ {
        self.entries.iter().map(|e| e.name.as_str())
    }

    pub fn chrom_len(&self, chr: &str) -> Option<u64> {
        self.entry(chr).map(|e| e.len)
    }

    /// Bases of the 1-based closed interval `iv`, clipped to `[1, chrom_len]` (empty when
    /// nothing is left). `None` for an unknown `chr`.
    pub fn slice(&self, chr: &str, iv: Iv) -> Option<Cow<'_, [u8]>> {
        let e = self.entry(chr)?;
        let s0 = (iv.s.max(1) - 1) as u64; // 0-based inclusive
        let e0 = (iv.e.max(0) as u64).min(e.len); // 0-based exclusive
        if s0 >= e0 {
            return Some(Cow::Borrowed(&[]));
        }
        let at = |i: u64| (e.offset + (i / e.line_bases) * e.line_bytes + i % e.line_bases) as usize;
        let (first_line, last_line) = (s0 / e.line_bases, (e0 - 1) / e.line_bases);
        if first_line == last_line {
            return Some(Cow::Borrowed(&self.map[at(s0)..at(e0 - 1) + 1]));
        }
        let mut out = Vec::with_capacity((e0 - s0) as usize);
        let mut i = s0;
        while i < e0 {
            let line_end = ((i / e.line_bases) + 1) * e.line_bases;
            let j = line_end.min(e0);
            out.extend_from_slice(&self.map[at(i)..at(j - 1) + 1]);
            i = j;
        }
        Some(Cow::Owned(out))
    }

    /// Convenience: a [`RefSlice`] over an owned buffer.
    pub fn ref_slice<'a>(&self, chr: &str, iv: Iv, buf: &'a mut Vec<u8>) -> Option<RefSlice<'a>> {
        let s = self.slice(chr, iv)?;
        buf.clear();
        buf.extend_from_slice(&s);
        Some(RefSlice { start: iv.s.max(1), seq: buf })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp_fasta(name: &str, body: &str, fai: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("iq_bam_fa_{}_{name}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("t.fa");
        std::fs::write(&fa, body).unwrap();
        std::fs::write(dir.join("t.fa.fai"), fai).unwrap();
        fa
    }

    #[test]
    fn multiline_and_clipping() {
        // chr a: ACGTacgtNN (line width 4), chr b: TTTT with CRLF
        let body = ">a\nACGT\nacgt\nNN\n>b\r\nTT\r\nTT\r\n";
        let fai = "a\t10\t3\t4\t5\nb\t4\t20\t2\t4\n";
        let fa = tmp_fasta("ml", body, fai);
        let f = FastaIndex::open(&fa).unwrap();
        let sl = |c: &str, s, e| f.slice(c, Iv { s, e }).unwrap().into_owned();
        assert_eq!(sl("a", 1, 10), b"ACGTacgtNN");
        assert_eq!(sl("a", 2, 3), b"CG");
        assert!(matches!(f.slice("a", Iv { s: 2, e: 3 }).unwrap(), Cow::Borrowed(_)));
        assert_eq!(sl("a", 4, 5), b"Ta");
        assert_eq!(sl("a", 0, 2), b"AC");
        assert_eq!(sl("a", 9, 50), b"NN");
        assert_eq!(sl("a", 11, 50), b"");
        assert_eq!(sl("a", 5, 4), b"");
        assert_eq!(sl("b", 1, 4), b"TTTT");
        assert_eq!(sl("b", 2, 3), b"TT");
        assert!(f.slice("c", Iv { s: 1, e: 2 }).is_none());
        assert_eq!(f.chrom_len("b"), Some(4));
    }
}
