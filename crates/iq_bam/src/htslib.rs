//! rust-htslib fallback backend (design §6.1), enabled with `--features htslib`.
//!
//! htslib does the region iteration (`hts_itr_next`) and BGZF inflate (`hts_set_threads`);
//! records are copied into the same flat [`RecordBatch`] as the noodles backend.
//!
//! Building hts-sys needs libclang for bindgen and, on some systems, the GCC builtin headers:
//! `LIBCLANG_PATH=<repo>/env/lib BINDGEN_EXTRA_CLANG_ARGS=-I/usr/lib/gcc/x86_64-redhat-linux/11/include`.

use crate::bam::{half_open, Decoder, FetchOpts, Parts, RecordBatch, RecordSource};
use iq_types::Iv0;
use rust_htslib::bam::{self as hbam, Read as _};
use std::io;
use std::path::{Path, PathBuf};

/// An indexed BAM read through htslib. Every fetch opens its own `IndexedReader`.
#[derive(Debug, Clone)]
pub struct HtsBam {
    path: PathBuf,
    refs: Vec<String>,
    pub tag: Option<[u8; 2]>,
    /// htslib decompression threads (0/1 = none).
    pub threads: usize,
}

fn io_err(e: rust_htslib::errors::Error) -> io::Error {
    io::Error::other(e.to_string())
}

impl HtsBam {
    pub fn open(path: &Path) -> io::Result<HtsBam> {
        let r = hbam::IndexedReader::from_path(path).map_err(io_err)?;
        let refs = r.header().target_names().iter().map(|n| String::from_utf8_lossy(n).into_owned()).collect();
        Ok(HtsBam { path: path.to_owned(), refs, tag: None, threads: 0 })
    }

    pub fn references(&self) -> &[String] {
        &self.refs
    }

    /// Region query (same semantics and batches as [`crate::try_fetch`]).
    pub fn try_fetch(&self, chr: &str, region: Iv0, opts: FetchOpts) -> io::Result<HtsFetch> {
        let tid = self
            .refs
            .iter()
            .position(|r| r == chr)
            .ok_or_else(|| io::Error::new(io::ErrorKind::NotFound, format!("unknown reference {chr:?}")))?;
        let mut r = hbam::IndexedReader::from_path(&self.path).map_err(io_err)?;
        if self.threads > 1 {
            r.set_threads(self.threads).map_err(io_err)?;
        }
        let (beg, end) = half_open(region);
        r.fetch((tid as i32, beg, end)).map_err(io_err)?;
        Ok(HtsFetch {
            reader: r,
            rec: hbam::Record::new(),
            dec: Decoder { tag: self.tag, seq: opts.seq },
            batch_records: opts.batch_records.max(1),
            seq_no: 0,
            done: false,
        })
    }
}

/// Streaming htslib region query.
pub struct HtsFetch {
    reader: hbam::IndexedReader,
    rec: hbam::Record,
    dec: Decoder,
    batch_records: usize,
    seq_no: u64,
    done: bool,
}

impl HtsFetch {
    fn push_current(&mut self, b: &mut RecordBatch) {
        let rec = &self.rec;
        let c = &rec.inner().core;
        let l_seq = c.l_qseq as usize;
        let data = {
            let inner = rec.inner();
            // SAFETY: `data` holds `l_data` initialized bytes owned by the record.
            unsafe { std::slice::from_raw_parts(inner.data as *const u8, inner.l_data as usize) }
        };
        let seq_off = c.l_qname as usize + 4 * c.n_cigar as usize;
        let qual_off = seq_off + l_seq.div_ceil(2);
        let aux_off = qual_off + l_seq;
        self.dec.push_parts(
            Parts {
                ref_id: c.tid,
                pos0: c.pos as i32,
                flag: c.flag,
                mapq: c.qual,
                qname: rec.qname(),
                cigar: rec.raw_cigar().iter().copied(),
                packed_seq: &data[seq_off..qual_off],
                l_seq,
                qual0: data.get(qual_off).copied().filter(|_| l_seq > 0),
                aux: &data[aux_off..],
            },
            b,
        );
    }
}

impl Iterator for HtsFetch {
    type Item = io::Result<RecordBatch>;
    fn next(&mut self) -> Option<Self::Item> {
        let mut b = RecordBatch::default();
        while !self.done && b.len() < self.batch_records {
            match self.reader.read(&mut self.rec) {
                None => self.done = true,
                Some(Err(e)) => {
                    self.done = true;
                    return Some(Err(io_err(e)));
                }
                Some(Ok(())) => self.push_current(&mut b),
            }
        }
        if b.is_empty() {
            return None;
        }
        b.seq_no = self.seq_no;
        self.seq_no += 1;
        Some(Ok(b))
    }
}

impl RecordSource for HtsBam {
    fn fetch(&self, chr: &str, region: Iv0) -> Box<dyn Iterator<Item = RecordBatch> + '_> {
        match self.try_fetch(chr, region, FetchOpts::default()) {
            Ok(it) => {
                let path = self.path.clone();
                Box::new(it.map(move |r| r.unwrap_or_else(|e| panic!("{}: {e}", path.display()))))
            }
            Err(e) if e.kind() == io::ErrorKind::NotFound => Box::new(std::iter::empty()),
            Err(e) => panic!("{}: {e}", self.path.display()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Both backends give identical batches on the pysam edge-case BAM.
    #[test]
    fn same_as_noodles() {
        let p = Path::new("tests/data/edge.bam");
        let mut h = HtsBam::open(p).unwrap();
        h.tag = Some(*b"HP");
        let mut n = crate::BamHandle::open(p).unwrap();
        n.set_tag(Some(*b"HP"));
        let key = |b: &RecordBatch| -> Vec<String> {
            (0..b.len())
                .map(|i| {
                    let r = b.get(i);
                    format!("{r:?} {:?} {}", b.tag_type(i), b.qual_len(i))
                })
                .collect()
        };
        for (chr, s, e) in
            [("c1", 0, 2_999_999), ("c1", 100, 110), ("c1", 2000, 2000), ("c2", 0, 99_999), ("c1", 50_000, 50_000)]
        {
            let iv = Iv0 { s, e };
            let a: Vec<String> =
                h.try_fetch(chr, iv, FetchOpts::default()).unwrap().flat_map(|b| key(&b.unwrap())).collect();
            let b: Vec<String> =
                crate::try_fetch(&n, chr, iv, FetchOpts::default()).unwrap().flat_map(|b| key(&b.unwrap())).collect();
            assert!(!a.is_empty());
            assert_eq!(a, b, "{chr}:{s}-{e}");
        }
    }
}
