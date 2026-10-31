//! BAM access (noodles BGZF + BAI, own record decoder).
//!
//! Records are decoded straight from the raw BAM bytes into a flat [`RecordBatch`] (one
//! allocation per batch, none per record): qname bytes, CIGAR words, IUPAC-decoded sequence
//! and the one typed read-group tag. Region semantics follow htslib's `hts_itr_next`, which is
//! what pysam `fetch(chr, start0, end)` uses:
//! - chunks from the BAI bins overlapping `[beg, end)`, pruned by the linear index;
//! - records are read in file order; reading stops at the first record on another reference
//!   or with `pos >= end`;
//! - a record is yielded iff `rec_end > beg && pos < end`, where `rec_end = pos + rlen` and
//!   `rlen` = reference length of the CIGAR, or 1 if that is 0 or the read is unmapped
//!   (`bam_endpos`).
//!
//! CIGARs with more than 65535 ops are stored as `<l_seq>S<rlen>N` plus a `CG:B:I` tag;
//! like htslib (`bam_tag2cigar`) the real CIGAR is taken from the tag.

use iq_types::{Iv0, RecordLite, TagValue};
use noodles_bgzf as bgzf;
use noodles_bgzf::io::Seek as _;
use noodles_bgzf::VirtualPosition;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use std::fs::File;
use std::io::{self, Read};
use std::num::NonZero;
use std::path::{Path, PathBuf};

const BAI_MIN_SHIFT: u32 = 14;
const BAI_DEPTH: u32 = 5;
/// Largest coordinate a BAI can address (2^29).
const BAI_MAX_POS: i64 = 1 << 29;

/// Per-handle options.
#[derive(Clone, Debug, Default)]
pub struct BamOptions {
    /// Read-group tag decoded into [`RecordLite::tag`] (e.g. `*b"HP"`), or none.
    pub tag: Option<[u8; 2]>,
    /// BGZF inflate threads per fetch; 0 or 1 = inline (single-threaded) inflate.
    pub threads: usize,
}

/// Per-fetch options.
#[derive(Clone, Copy, Debug)]
pub struct FetchOpts {
    /// Decode sequences (`RecordLite::seq`); pass A does not need them.
    pub seq: bool,
    /// Target number of records per batch.
    pub batch_records: usize,
}

impl Default for FetchOpts {
    fn default() -> Self {
        FetchOpts { seq: true, batch_records: 4096 }
    }
}

/// An opened, indexed BAM (header, index, reference names, unmapped count). `Sync`: every
/// fetch opens its own file descriptor, so worker threads can fetch concurrently.
#[derive(Debug)]
pub struct BamHandle {
    path: PathBuf,
    refs: Vec<String>,
    ref_lens: Vec<u64>,
    index: noodles_bam::bai::Index,
    opts: BamOptions,
}

impl BamHandle {
    /// Opens `path` and `path.bai` (or `path` with `.bam` replaced by `.bai`).
    pub fn open(path: &Path) -> io::Result<BamHandle> {
        Self::open_with(path, BamOptions::default())
    }

    pub fn open_with(path: &Path, opts: BamOptions) -> io::Result<BamHandle> {
        let (refs, ref_lens) = read_header(path)?;
        let mut bai = PathBuf::from(format!("{}.bai", path.display()));
        if !bai.exists() {
            let alt = path.with_extension("bai");
            if alt.exists() {
                bai = alt;
            }
        }
        let index = noodles_bam::bai::fs::read(&bai)
            .map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", bai.display())))?;
        if index.reference_sequences().len() > refs.len() {
            return Err(invalid(format!("{}: index has more references than the header", bai.display())));
        }
        Ok(BamHandle { path: path.to_owned(), refs, ref_lens, index, opts })
    }

    pub fn path(&self) -> &Path {
        &self.path
    }
    pub fn options(&self) -> &BamOptions {
        &self.opts
    }
    /// Sets the read-group tag decoded into [`RecordLite::tag`].
    pub fn set_tag(&mut self, tag: Option<[u8; 2]>) {
        self.opts.tag = tag;
    }
    /// Sets the number of BGZF inflate threads per fetch.
    pub fn set_threads(&mut self, threads: usize) {
        self.opts.threads = threads;
    }
    /// Reference names in header order.
    pub fn references(&self) -> &[String] {
        &self.refs
    }
    /// Reference lengths in header order.
    pub fn reference_lengths(&self) -> &[u64] {
        &self.ref_lens
    }
    /// Header index of reference `chr`.
    pub fn ref_id(&self, chr: &str) -> Option<usize> {
        self.refs.iter().position(|r| r == chr)
    }
    /// pysam `bam.unmapped`: per-reference unmapped counts from the index metadata plus the
    /// unplaced unmapped count (`nocoordinate`). Used for §8.2 `__not_aligned`.
    pub fn unmapped(&self) -> u64 {
        use noodles_csi::binning_index::ReferenceSequence as _;
        use noodles_csi::BinningIndex as _;
        let placed: u64 = self
            .index
            .reference_sequences()
            .iter()
            .filter_map(|rs| rs.metadata())
            .map(|m| m.unmapped_record_count())
            .sum();
        placed + self.index.unplaced_unmapped_record_count().unwrap_or(0)
    }

    /// Index chunks for `tid`, `[beg, end)` (htslib `hts_itr_query`).
    fn chunks(&self, tid: usize, beg: i64, end: i64) -> Vec<Chunk> {
        let Some(rs) = self.index.reference_sequences().get(tid) else {
            return Vec::new();
        };
        let mut chunks = Vec::new();
        reg2bins(beg, end, |bin| {
            if let Some(b) = rs.bins().get(&(bin as usize)) {
                chunks.extend_from_slice(b.chunks());
            }
        });
        let lin = rs.index();
        let w = (beg >> BAI_MIN_SHIFT) as usize;
        let min_off = lin.get(w).copied().unwrap_or(VirtualPosition::MIN);
        noodles_csi::binning_index::optimize_chunks(&chunks, min_off)
    }
}

fn invalid(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

/// Reads the BAM header: reference names and lengths.
fn read_header(path: &Path) -> io::Result<(Vec<String>, Vec<u64>)> {
    let f = File::open(path).map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", path.display())))?;
    let mut r = bgzf::io::Reader::new(f);
    let mut magic = [0u8; 4];
    r.read_exact(&mut magic)?;
    if &magic != b"BAM\x01" {
        return Err(invalid(format!("{}: not a BAM file", path.display())));
    }
    let l_text = read_u32(&mut r)? as u64;
    io::copy(&mut (&mut r).take(l_text), &mut io::sink())?;
    let n_ref = read_u32(&mut r)? as usize;
    let mut refs = Vec::with_capacity(n_ref);
    let mut lens = Vec::with_capacity(n_ref);
    let mut buf = Vec::new();
    for _ in 0..n_ref {
        let l_name = read_u32(&mut r)? as usize;
        buf.resize(l_name, 0);
        r.read_exact(&mut buf)?;
        if buf.last() == Some(&0) {
            buf.pop();
        }
        refs.push(String::from_utf8(buf.clone()).map_err(|_| invalid("non-UTF-8 reference name"))?);
        lens.push(read_u32(&mut r)? as u64);
    }
    Ok((refs, lens))
}

fn read_u32<R: Read>(r: &mut R) -> io::Result<u32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}

/// htslib `reg2bins` for BAI (min_shift 14, depth 5), `[beg, end)`.
fn reg2bins(beg: i64, end: i64, mut f: impl FnMut(u32)) {
    let end = end - 1;
    if end < beg {
        return;
    }
    f(0);
    let mut t: u32 = 0;
    let mut s = BAI_MIN_SHIFT + BAI_DEPTH * 3;
    for l in 1..=BAI_DEPTH {
        t += 1 << ((l - 1) * 3);
        s -= 3;
        let b = t + (beg >> s) as u32;
        let e = t + (end >> s) as u32;
        for k in b..=e {
            f(k);
        }
    }
}

// ------------------------------------------------------------------------------------------
// RecordBatch
// ------------------------------------------------------------------------------------------

const NO_TAG: u32 = u32::MAX;

#[derive(Clone, Copy, Debug)]
struct Meta {
    pos0: i32,
    ref_id: i32,
    flag: u16,
    mapq: u8,
    /// SAM type letter of the tag (`b'Z'`, `b'C'`, ...), 0 if absent; array subtype for `B`.
    tag_ty: [u8; 2],
    tag: u32,
    name: u32,
    cigar: u32,
    seq: u32,
    l_seq: u32,
    qual_missing: bool,
}

/// A batch of decoded records in file order, with a sequence number for re-ordering after
/// parallel processing (§5.5). Flat storage: no per-record allocation.
#[derive(Debug, Default)]
#[non_exhaustive]
pub struct RecordBatch {
    pub seq_no: u64,
    meta: Vec<Meta>,
    names: Vec<u8>,
    cigars: Vec<u32>,
    seqs: Vec<u8>,
    tags: Vec<TagValue>,
}

impl RecordBatch {
    pub fn len(&self) -> usize {
        self.meta.len()
    }
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
    pub fn clear(&mut self) {
        self.meta.clear();
        self.names.clear();
        self.cigars.clear();
        self.seqs.clear();
        self.tags.clear();
    }
    /// `[start(i), start(i + 1))` of one flat array (`total` = its length).
    #[inline]
    fn range(&self, i: usize, start: fn(&Meta) -> u32, total: usize) -> std::ops::Range<usize> {
        start(&self.meta[i]) as usize..self.meta.get(i + 1).map_or(total, |m| start(m) as usize)
    }
    /// Record `i` as a borrowed view (sequence decoded to IUPAC letters, or empty if the
    /// batch was fetched without sequences or the record has none).
    pub fn get(&self, i: usize) -> RecordLite<'_> {
        let m = &self.meta[i];
        let name = self.range(i, |m| m.name, self.names.len());
        let cig = self.range(i, |m| m.cigar, self.cigars.len());
        let seq = self.range(i, |m| m.seq, self.seqs.len());
        RecordLite {
            qname: &self.names[name],
            flag: m.flag,
            pos0: m.pos0 as i64,
            mapq: m.mapq,
            cigar: &self.cigars[cig],
            seq: &self.seqs[seq],
            tag: (m.tag != NO_TAG).then(|| &self.tags[m.tag as usize]),
        }
    }
    /// Header reference index of record `i` (`-1` = none).
    pub fn ref_id(&self, i: usize) -> i32 {
        self.meta[i].ref_id
    }
    /// `l_seq` of record `i` (also when sequences were not decoded).
    pub fn seq_len(&self, i: usize) -> usize {
        self.meta[i].l_seq as usize
    }
    /// pysam `len(query_qualities)` (0 when the quality string is absent).
    pub fn qual_len(&self, i: usize) -> usize {
        if self.meta[i].qual_missing {
            0
        } else {
            self.meta[i].l_seq as usize
        }
    }
    /// Type of record `i`'s tag as pysam `get_tag(.., with_value_type=True)` reports it:
    /// `"Z"`, `"i"`, `"C"`, ..., or `"B"` + subtype (`"BC"`).
    pub fn tag_type(&self, i: usize) -> Option<&str> {
        let t = &self.meta[i].tag_ty;
        let n = match t {
            [0, _] => return None,
            [_, 0] => 1,
            _ => 2,
        };
        std::str::from_utf8(&t[..n]).ok()
    }
    pub fn iter(&self) -> impl ExactSizeIterator<Item = RecordLite<'_>> + '_ {
        (0..self.len()).map(|i| self.get(i))
    }

    fn intern_tag(&mut self, v: TagValue) -> u32 {
        let n = self.tags.len();
        let lo = n.saturating_sub(16);
        if let Some(k) = self.tags[lo..].iter().position(|t| *t == v) {
            return (lo + k) as u32;
        }
        self.tags.push(v);
        n as u32
    }
}

/// `seq_nt16_str` decode table for byte pairs.
static SEQ_PAIRS: [[u8; 2]; 256] = {
    const NT16: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
    let mut t = [[0u8; 2]; 256];
    let mut i = 0;
    while i < 256 {
        t[i] = [NT16[i >> 4], NT16[i & 15]];
        i += 1;
    }
    t
};

/// Appends the IUPAC letters of a 4-bit packed sequence of `l_seq` bases.
pub fn decode_seq_into(packed: &[u8], l_seq: usize, out: &mut Vec<u8>) {
    let full = l_seq / 2;
    out.reserve(l_seq);
    for &b in &packed[..full] {
        out.extend_from_slice(&SEQ_PAIRS[b as usize]);
    }
    if l_seq % 2 == 1 {
        out.push(SEQ_PAIRS[packed[full] as usize][0]);
    }
}

#[inline]
fn le_u16(b: &[u8], o: usize) -> u16 {
    u16::from_le_bytes([b[o], b[o + 1]])
}
#[inline]
fn le_u32(b: &[u8], o: usize) -> u32 {
    u32::from_le_bytes([b[o], b[o + 1], b[o + 2], b[o + 3]])
}

/// Size in bytes of one aux value of type `ty` starting at `v`; `None` if malformed.
fn aux_value_len(ty: u8, v: &[u8]) -> Option<usize> {
    Some(match ty {
        b'A' | b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        b'Z' | b'H' => v.iter().position(|&c| c == 0)? + 1,
        b'B' => {
            if v.len() < 5 {
                return None;
            }
            let w = match v[0] {
                b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                _ => return None,
            };
            5 + w * le_u32(v, 1) as usize
        }
        _ => return None,
    })
}

/// Finds tag `tag` in the aux block: `(type, value bytes)`.
fn find_aux(aux: &[u8], tag: [u8; 2]) -> Option<(u8, &[u8])> {
    let mut o = 0;
    while o + 3 <= aux.len() {
        let ty = aux[o + 2];
        let v = &aux[o + 3..];
        let n = aux_value_len(ty, v)?;
        if n > v.len() {
            return None;
        }
        if aux[o..o + 2] == tag {
            return Some((ty, &v[..n]));
        }
        o += 3 + n;
    }
    None
}

/// Typed tag value as pysam `get_tag` returns it.
fn tag_value(ty: u8, v: &[u8]) -> Option<TagValue> {
    Some(match ty {
        b'A' => TagValue::Char(v[0] as char),
        b'c' => TagValue::Int(v[0] as i8 as i64),
        b'C' => TagValue::Int(v[0] as i64),
        b's' => TagValue::Int(le_u16(v, 0) as i16 as i64),
        b'S' => TagValue::Int(le_u16(v, 0) as i64),
        b'i' => TagValue::Int(le_u32(v, 0) as i32 as i64),
        b'I' => TagValue::Int(le_u32(v, 0) as i64),
        b'f' => TagValue::Float(f32::from_bits(le_u32(v, 0)) as f64),
        b'Z' | b'H' => TagValue::Str(String::from_utf8_lossy(&v[..v.len() - 1]).into_owned()),
        // pysam returns `array.array`; `str()` of it is what a read grouper would see.
        b'B' => {
            let (code, w) = match v[0] {
                b'c' => ('b', 1),
                b'C' => ('B', 1),
                b's' => ('h', 2),
                b'S' => ('H', 2),
                b'i' => ('i', 4),
                b'I' => ('I', 4),
                b'f' => ('f', 4),
                _ => return None,
            };
            let n = le_u32(v, 1) as usize;
            let el: Vec<String> = (0..n)
                .map(|k| {
                    let e = &v[5 + k * w..];
                    match v[0] {
                        b'c' => (e[0] as i8).to_string(),
                        b'C' => e[0].to_string(),
                        b's' => (le_u16(e, 0) as i16).to_string(),
                        b'S' => le_u16(e, 0).to_string(),
                        b'i' => (le_u32(e, 0) as i32).to_string(),
                        b'I' => le_u32(e, 0).to_string(),
                        _ => format!("{:?}", f32::from_bits(le_u32(e, 0)) as f64),
                    }
                })
                .collect();
            if n == 0 {
                TagValue::Str(format!("array('{code}')"))
            } else {
                TagValue::Str(format!("array('{code}', [{}])", el.join(", ")))
            }
        }
        _ => return None,
    })
}

/// Reference length consumed by CIGAR words (M, D, N, =, X).
#[inline]
fn cigar_rlen(words: impl Iterator<Item = u32>) -> i64 {
    words.filter(|w| matches!(w & 0xf, 0 | 2 | 3 | 7 | 8)).map(|w| (w >> 4) as i64).sum()
}

/// What the region filter needs from a raw record.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Core {
    pub ref_id: i32,
    pub pos: i64,
    /// htslib `bam_endpos`.
    pub end: i64,
}

/// Raw-record decoder state shared by both backends.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Decoder {
    pub tag: Option<[u8; 2]>,
    pub seq: bool,
}

/// Field offsets of one raw BAM record (without the `block_size` prefix).
struct Layout {
    l_read_name: usize,
    n_cigar: usize,
    l_seq: usize,
    cigar_off: usize,
    seq_off: usize,
    qual_off: usize,
    aux_off: usize,
}

fn layout(raw: &[u8]) -> io::Result<Layout> {
    if raw.len() < 32 {
        return Err(invalid("BAM record shorter than 32 bytes"));
    }
    let l_read_name = raw[8] as usize;
    let n_cigar = le_u16(raw, 12) as usize;
    let l_seq = le_u32(raw, 16) as usize;
    let cigar_off = 32 + l_read_name;
    let seq_off = cigar_off + 4 * n_cigar;
    let qual_off = seq_off + l_seq.div_ceil(2);
    let aux_off = qual_off + l_seq;
    if aux_off > raw.len() || l_read_name == 0 {
        return Err(invalid("truncated BAM record"));
    }
    Ok(Layout { l_read_name, n_cigar, l_seq, cigar_off, seq_off, qual_off, aux_off })
}

/// CIGAR from a `CG:B:I` tag when the stored CIGAR is the `<l_seq>S<rlen>N` placeholder
/// (htslib `bam_tag2cigar`).
fn long_cigar<'a>(raw: &'a [u8], l: &Layout) -> Option<&'a [u8]> {
    if l.n_cigar == 0 || (le_u32(raw, 0) as i32) < 0 || (le_u32(raw, 4) as i32) < 0 {
        return None;
    }
    let c0 = le_u32(raw, l.cigar_off);
    if c0 & 0xf != 4 || (c0 >> 4) as usize != l.l_seq {
        return None;
    }
    let (ty, v) = find_aux(&raw[l.aux_off..], *b"CG")?;
    if ty != b'B' || v[0] != b'I' {
        return None;
    }
    let n = le_u32(v, 1) as usize;
    Some(&v[5..5 + 4 * n])
}

impl Decoder {
    /// Region-filter fields of a raw record.
    pub(crate) fn core(&self, raw: &[u8]) -> io::Result<Core> {
        let l = layout(raw)?;
        let ref_id = le_u32(raw, 0) as i32;
        let pos = le_u32(raw, 4) as i32 as i64;
        let flag = le_u16(raw, 14);
        let rlen = if flag & 4 != 0 {
            0
        } else if let Some(cg) = long_cigar(raw, &l) {
            cigar_rlen(cg.chunks_exact(4).map(|c| le_u32(c, 0)))
        } else {
            cigar_rlen(raw[l.cigar_off..l.seq_off].chunks_exact(4).map(|c| le_u32(c, 0)))
        };
        Ok(Core { ref_id, pos, end: pos + rlen.max(1) })
    }

    /// Appends a raw record (BAM layout without `block_size`) to `b`.
    pub(crate) fn push(&self, raw: &[u8], b: &mut RecordBatch) -> io::Result<()> {
        let l = layout(raw)?;
        let words = long_cigar(raw, &l).unwrap_or(&raw[l.cigar_off..l.seq_off]);
        self.push_parts(
            Parts {
                ref_id: le_u32(raw, 0) as i32,
                pos0: le_u32(raw, 4) as i32,
                flag: le_u16(raw, 14),
                mapq: raw[9],
                qname: &raw[32..32 + l.l_read_name - 1],
                cigar: words.chunks_exact(4).map(|c| le_u32(c, 0)),
                packed_seq: &raw[l.seq_off..l.qual_off],
                l_seq: l.l_seq,
                qual0: raw.get(l.qual_off).copied(),
                aux: &raw[l.aux_off..],
            },
            b,
        );
        Ok(())
    }

    /// Appends one record given as its decoded parts (shared with the htslib backend).
    pub(crate) fn push_parts(&self, p: Parts<'_, impl Iterator<Item = u32>>, b: &mut RecordBatch) {
        let name = b.names.len() as u32;
        let cigar = b.cigars.len() as u32;
        let seq = b.seqs.len() as u32;
        b.names.extend_from_slice(p.qname);
        b.cigars.extend(p.cigar);
        if self.seq {
            decode_seq_into(p.packed_seq, p.l_seq, &mut b.seqs);
        }
        let (mut tag_ty, mut tag) = ([0, 0], NO_TAG);
        if let Some(t) = self.tag {
            if let Some((ty, v)) = find_aux(p.aux, t) {
                if let Some(val) = tag_value(ty, v) {
                    tag_ty = [ty, if ty == b'B' { v[0] } else { 0 }];
                    tag = b.intern_tag(val);
                }
            }
        }
        b.meta.push(Meta {
            pos0: p.pos0,
            ref_id: p.ref_id,
            flag: p.flag,
            mapq: p.mapq,
            tag_ty,
            tag,
            name,
            cigar,
            seq,
            l_seq: p.l_seq as u32,
            qual_missing: p.l_seq == 0 || p.qual0 == Some(0xff),
        });
    }
}

/// One record's fields for [`Decoder::push_parts`].
pub(crate) struct Parts<'a, C> {
    pub ref_id: i32,
    pub pos0: i32,
    pub flag: u16,
    pub mapq: u8,
    /// Without the NUL terminator.
    pub qname: &'a [u8],
    pub cigar: C,
    pub packed_seq: &'a [u8],
    pub l_seq: usize,
    /// First quality byte (`0xff` = absent).
    pub qual0: Option<u8>,
    pub aux: &'a [u8],
}

// ------------------------------------------------------------------------------------------
// Fetch
// ------------------------------------------------------------------------------------------

/// Backend abstraction (noodles by default, rust-htslib with `--features htslib`).
pub trait RecordSource {
    fn fetch(&self, chr: &str, region: Iv0) -> Box<dyn Iterator<Item = RecordBatch> + '_>;
}

impl RecordSource for BamHandle {
    fn fetch(&self, chr: &str, region: Iv0) -> Box<dyn Iterator<Item = RecordBatch> + '_> {
        Box::new(fetch(self, chr, region))
    }
}

enum Reader {
    St(bgzf::io::Reader<File>),
    Mt(bgzf::io::MultithreadedReader<File>),
}

impl Reader {
    fn seek(&mut self, vp: VirtualPosition) -> io::Result<()> {
        match self {
            Reader::St(r) => r.seek(vp).map(drop),
            Reader::Mt(r) => r.seek_to_virtual_position(vp).map(drop),
        }
    }
    fn vpos(&self) -> VirtualPosition {
        match self {
            Reader::St(r) => r.virtual_position(),
            Reader::Mt(r) => r.virtual_position(),
        }
    }
    /// Reads one record into `buf`; `false` at EOF.
    fn read_record(&mut self, buf: &mut Vec<u8>) -> io::Result<bool> {
        let r: &mut dyn Read = match self {
            Reader::St(r) => r,
            Reader::Mt(r) => r,
        };
        let mut len = [0u8; 4];
        match r.read_exact(&mut len) {
            Ok(()) => {}
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(false),
            Err(e) => return Err(e),
        }
        let n = u32::from_le_bytes(len) as usize;
        buf.resize(n, 0);
        r.read_exact(buf)?;
        Ok(true)
    }
}

/// Streaming region query; yields `io::Result<RecordBatch>`.
pub struct FetchIter {
    reader: Option<Reader>,
    chunks: Vec<Chunk>,
    next_chunk: usize,
    cur_end: Option<VirtualPosition>,
    tid: i32,
    beg: i64,
    end: i64,
    dec: Decoder,
    batch_records: usize,
    seq_no: u64,
    raw: Vec<u8>,
    done: bool,
}

impl FetchIter {
    fn fill(&mut self, b: &mut RecordBatch) -> io::Result<()> {
        let Some(reader) = self.reader.as_mut() else {
            self.done = true;
            return Ok(());
        };
        while b.len() < self.batch_records {
            // advance to the next chunk once the current one is exhausted
            if self.cur_end.is_none_or(|e| reader.vpos() >= e) {
                let Some(c) = self.chunks.get(self.next_chunk) else {
                    self.done = true;
                    return Ok(());
                };
                self.next_chunk += 1;
                if self.cur_end.is_none() || reader.vpos() != c.start() {
                    reader.seek(c.start())?;
                }
                self.cur_end = Some(c.end());
                continue;
            }
            if !reader.read_record(&mut self.raw)? {
                self.done = true;
                return Ok(());
            }
            let core = self.dec.core(&self.raw)?;
            if core.ref_id != self.tid || core.pos >= self.end {
                self.done = true;
                return Ok(());
            }
            if core.end > self.beg {
                self.dec.push(&self.raw, b)?;
            }
        }
        Ok(())
    }
}

impl Iterator for FetchIter {
    type Item = io::Result<RecordBatch>;
    fn next(&mut self) -> Option<Self::Item> {
        let mut b = RecordBatch::default();
        while !self.done && b.is_empty() {
            if let Err(e) = self.fill(&mut b) {
                self.done = true;
                return Some(Err(e));
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

/// htslib-style region bounds: 0-based closed `region` → `[beg, end)` clipped to the BAI range.
pub(crate) fn half_open(region: Iv0) -> (i64, i64) {
    (region.s.max(0), (region.e + 1).min(BAI_MAX_POS))
}

/// Region query with explicit options. Errors: unknown reference, I/O.
pub fn try_fetch(bam: &BamHandle, chr: &str, region: Iv0, opts: FetchOpts) -> io::Result<FetchIter> {
    let tid = bam.ref_id(chr).ok_or_else(|| io::Error::new(io::ErrorKind::NotFound, format!("unknown reference {chr:?}")))?;
    let (beg, end) = half_open(region);
    let chunks = if beg < end { bam.chunks(tid, beg, end) } else { Vec::new() };
    let reader = if chunks.is_empty() {
        None
    } else {
        let f = File::open(&bam.path)?;
        Some(match NonZero::new(bam.opts.threads).filter(|n| n.get() > 1) {
            Some(n) => Reader::Mt(bgzf::io::MultithreadedReader::with_worker_count(n, f)),
            None => Reader::St(bgzf::io::Reader::new(f)),
        })
    };
    Ok(FetchIter {
        reader,
        chunks,
        next_chunk: 0,
        cur_end: None,
        tid: tid as i32,
        beg,
        end,
        dec: Decoder { tag: bam.opts.tag, seq: opts.seq },
        batch_records: opts.batch_records.max(1),
        seq_no: 0,
        raw: Vec::with_capacity(1 << 16),
        done: false,
    })
}

/// `fn fetch(&BamHandle, chr, Iv0) -> impl Iterator<Item=RecordBatch>` (U2): all records
/// overlapping the 0-based closed `region`, in file order (pysam
/// `fetch(chr, region.s, region.e + 1)`). An unknown `chr` yields nothing.
///
/// # Panics
/// On a corrupt or unreadable BAM (use [`try_fetch`] to handle errors).
pub fn fetch<'a>(bam: &'a BamHandle, chr: &str, region: Iv0) -> impl Iterator<Item = RecordBatch> + 'a {
    fetch_with(bam, chr, region, FetchOpts::default())
}

/// [`fetch`] with options (e.g. no sequences for pass A).
pub fn fetch_with<'a>(bam: &'a BamHandle, chr: &str, region: Iv0, opts: FetchOpts) -> impl Iterator<Item = RecordBatch> + 'a {
    let it = match try_fetch(bam, chr, region, opts) {
        Ok(it) => Some(it),
        Err(e) if e.kind() == io::ErrorKind::NotFound => None,
        Err(e) => panic!("{}: {e}", bam.path.display()),
    };
    let path = bam.path.clone();
    it.into_iter().flatten().map(move |r| r.unwrap_or_else(|e| panic!("{}: {e}", path.display())))
}

pub(crate) mod fixture {
    use super::*;
    use iq_types::fixture::schema::{TagFixture, TagValueFixture, U2FetchIn, U2FetchOut};
    use iq_types::fixture::{FixError, FixResult, FixtureCtx, UnitFixture};
    use std::collections::HashMap;
    use std::sync::{Arc, Mutex, OnceLock};

    pub struct U2Fetch;

    /// `"chr:s-e"` (1-based closed) → `(chr, 0-based closed)`.
    pub(crate) fn parse_region(r: &str) -> FixResult<(&str, Iv0)> {
        let bad = || FixError(format!("bad region {r:?}"));
        let (chr, se) = r.rsplit_once(':').ok_or_else(bad)?;
        let (s, e) = se.split_once('-').ok_or_else(bad)?;
        let num = |x: &str| x.replace(',', "").parse::<i64>().map_err(|_| bad());
        Ok((chr, Iv0 { s: num(s)? - 1, e: num(e)? - 1 }))
    }

    fn to_fixture(bam: &BamHandle, b: &RecordBatch, i: usize, tag: &str) -> U2FetchOut {
        let r = b.get(i);
        let ref_id = b.ref_id(i);
        let mut out = U2FetchOut {
            qname: String::from_utf8_lossy(r.qname).into_owned(),
            flag: r.flag,
            ref_id,
            ref_name: usize::try_from(ref_id).ok().and_then(|t| bam.references().get(t).cloned()),
            pos0: r.pos0,
            // pysam reference_end = bam_endpos: None for unmapped reads or an empty CIGAR, and
            // pos + 1 for a CIGAR without reference-consuming ops (RecordLite::end0 gives pos)
            ref_end0: (!r.is_unmapped() && !r.cigar.is_empty()).then(|| r.end0().max(r.pos0 + 1)),
            mapq: r.mapq,
            cigar: r.cigar.iter().map(|&w| ((w & 0xf) as u8, w >> 4)).collect(),
            seq: (!r.seq.is_empty()).then(|| String::from_utf8_lossy(r.seq).into_owned()),
            qual_len: b.qual_len(i) as i64,
            tags: Default::default(),
        };
        if let (Some(v), Some(ty)) = (r.tag, b.tag_type(i)) {
            let value = match v {
                TagValue::Int(x) => TagValueFixture::Int(*x),
                TagValue::Float(x) => TagValueFixture::Float(*x),
                TagValue::Str(x) => TagValueFixture::Str(x.clone()),
                TagValue::Char(c) => TagValueFixture::Str(c.to_string()),
            };
            out.tags.insert(tag.to_owned(), TagFixture { ty: ty.to_owned(), value });
        }
        out
    }

    type Cache = Mutex<HashMap<(String, String), Arc<Vec<U2FetchOut>>>>;

    /// All records of `(bam, region)`, fetched once per process.
    fn records(i: &U2FetchIn) -> FixResult<Arc<Vec<U2FetchOut>>> {
        static CACHE: OnceLock<Cache> = OnceLock::new();
        let cache = CACHE.get_or_init(Default::default);
        let key = (i.bam.clone(), i.region.clone());
        if let Some(v) = cache.lock().unwrap().get(&key) {
            return Ok(v.clone());
        }
        let (chr, iv) = parse_region(&i.region)?;
        let bam = BamHandle::open_with(Path::new(&i.bam), BamOptions { tag: Some(*b"HP"), threads: 1 })?;
        let mut out = Vec::new();
        for b in try_fetch(&bam, chr, iv, FetchOpts::default())? {
            let b = b?;
            out.extend((0..b.len()).map(|k| to_fixture(&bam, &b, k, "HP")));
        }
        let v = Arc::new(out);
        cache.lock().unwrap().insert(key, v.clone());
        Ok(v)
    }

    impl UnitFixture for U2Fetch {
        const UNIT: &'static str = "U2";
        const FN: &'static str = "fetch";
        type In = U2FetchIn;
        type Out = U2FetchOut;

        fn run(i: &U2FetchIn, _ctx: &FixtureCtx) -> FixResult<U2FetchOut> {
            let recs = records(i)?;
            recs.get(i.n as usize)
                .cloned()
                .ok_or_else(|| FixError(format!("record {} not found: fetch gave {} records", i.n, recs.len())))
        }
    }
}

/// Registers the U2 runner (`fn = "fetch"`).
pub fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    reg.add::<fixture::U2Fetch>();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seq_decode_iupac() {
        // codes 0..15 packed pairwise
        let packed: Vec<u8> = (0..8u8).map(|k| (2 * k) << 4 | (2 * k + 1)).collect();
        let mut out = Vec::new();
        decode_seq_into(&packed, 16, &mut out);
        assert_eq!(out, b"=ACMGRSVTWYHKDBN");
        out.clear();
        decode_seq_into(&packed, 3, &mut out);
        assert_eq!(out, b"=AC");
    }

    #[test]
    fn reg2bins_matches_htslib() {
        let mut v = Vec::new();
        reg2bins(0, 1, |b| v.push(b));
        assert_eq!(v, vec![0, 1, 9, 73, 585, 4681]);
        v.clear();
        reg2bins(16384, 16385, |b| v.push(b));
        assert_eq!(*v.last().unwrap(), 4682);
    }

    #[test]
    fn tag_values() {
        assert_eq!(tag_value(b'c', &[0xff]), Some(TagValue::Int(-1)));
        assert_eq!(tag_value(b'S', &[0x01, 0x01]), Some(TagValue::Int(257)));
        assert_eq!(tag_value(b'Z', b"ab\0"), Some(TagValue::Str("ab".into())));
        assert_eq!(tag_value(b'A', b"x"), Some(TagValue::Char('x')));
        let b = [b'C', 2, 0, 0, 0, 1, 2];
        assert_eq!(tag_value(b'B', &b), Some(TagValue::Str("array('B', [1, 2])".into())));
    }
}
