//! Per-read types: raw record view (U2 → U4/U7), aligned read (U4 out), the logical
//! [`ReadRecord`] (collect output = construct input, design §2.3), reference slices and
//! pass-A blocks (U3).

use crate::assign::Assignment;
use crate::prim::{GeneIdx, Iv, Iv0, PolyAInfo, Pos, Strand};
use serde::{Deserialize, Serialize};

// ------------------------------------------------------------------------------------------
// ReadRecord
// ------------------------------------------------------------------------------------------

/// Logical read record = one 4.1 `ReadAssignment` as dumped after collect (design §2.3).
/// Physical storage is the columnar [`crate::BlockStore`]; convert with
/// [`crate::BlockStore::from_records`] / [`crate::BlockStore::to_records`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ReadRecord {
    /// Position in block, BAM order (= 4.1 dump order, §5.5).
    pub idx: u32,
    /// QNAME. Kept in the BlockStore's cold-name arena.
    pub name: String,
    /// `ReadAssignment.exons` (after polyA/T trimming).
    pub exons: Vec<Iv>,
    /// `ReadAssignment.corrected_exons`.
    pub corrected_exons: Vec<Iv>,
    /// `boundary_deletions`: `(ref_pos, length)` pairs near exon boundaries (U4).
    pub boundary_deletions: Vec<(Pos, u32)>,
    pub polya: PolyAInfo,
    /// `polyA_found`.
    pub polya_found: bool,
    pub strand: Strand,
    pub mapped_strand: Strand,
    /// `mapping_quality`.
    pub mapq: u8,
    /// Set only by U9.
    pub multimapper: bool,
    /// Read-group pool id ([`crate::GroupPool`]); 0 when no `--read-group-tag`.
    pub group: u32,
    /// Reference assignment (`assignment_type`, `gene_assignment_type`, `isoform_matches`).
    pub assignment: Assignment,
    /// `genomic_region`: the block region (kept for the U9 noninformative tie rule).
    pub genomic_region: Iv,
}

impl ReadRecord {
    /// `corrected_introns` = `junctions_from_blocks(corrected_exons)`.
    pub fn corrected_introns(&self) -> Vec<Iv> {
        crate::prim::junctions_from_blocks(&self.corrected_exons)
    }
}

// ------------------------------------------------------------------------------------------
// Raw records (U2) and aligned reads (U4)
// ------------------------------------------------------------------------------------------

/// BAM CIGAR operation codes (SAM spec; pysam `cigartuples` use the same integers).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[repr(u8)]
pub enum CigarOp {
    M = 0,
    I = 1,
    D = 2,
    N = 3,
    S = 4,
    H = 5,
    P = 6,
    Eq = 7,
    X = 8,
}

impl CigarOp {
    pub const fn from_code(c: u8) -> Option<CigarOp> {
        Some(match c {
            0 => CigarOp::M,
            1 => CigarOp::I,
            2 => CigarOp::D,
            3 => CigarOp::N,
            4 => CigarOp::S,
            5 => CigarOp::H,
            6 => CigarOp::P,
            7 => CigarOp::Eq,
            8 => CigarOp::X,
            _ => return None,
        })
    }
    /// Consumes the reference (M, D, N, =, X).
    pub const fn consumes_ref(self) -> bool {
        matches!(self, CigarOp::M | CigarOp::D | CigarOp::N | CigarOp::Eq | CigarOp::X)
    }
    /// Consumes the query (M, I, S, =, X).
    pub const fn consumes_query(self) -> bool {
        matches!(self, CigarOp::M | CigarOp::I | CigarOp::S | CigarOp::Eq | CigarOp::X)
    }
}

/// Decode one raw BAM CIGAR word (`len << 4 | op`). Unknown op codes map to `None`.
#[inline]
pub const fn cigar_decode(raw: u32) -> (Option<CigarOp>, u32) {
    (CigarOp::from_code((raw & 0xf) as u8), raw >> 4)
}
/// Encode `(op, len)` as a raw BAM CIGAR word.
#[inline]
pub const fn cigar_encode(op: CigarOp, len: u32) -> u32 {
    (len << 4) | op as u32
}

/// SAM flag bits used by IsoQuant.
pub mod flags {
    pub const PAIRED: u16 = 0x1;
    pub const UNMAPPED: u16 = 0x4;
    pub const REVERSE: u16 = 0x10;
    pub const SECONDARY: u16 = 0x100;
    pub const SUPPLEMENTARY: u16 = 0x800;
}

/// Value of the read-group BAM tag (any type is accepted, D-O2).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum TagValue {
    Str(String),
    Int(i64),
    Float(f64),
    Char(char),
}

impl TagValue {
    /// Python `str(value)` formatting for integer/string/char values (design §8.3). Floats
    /// use Rust's shortest round-trip form, which equals Python `repr` for finite values
    /// except for exponent formatting; HP tags are never floats in practice.
    pub fn to_group_string(&self) -> String {
        match self {
            TagValue::Str(s) => s.clone(),
            TagValue::Int(i) => i.to_string(),
            TagValue::Float(f) => format!("{f:?}"),
            TagValue::Char(c) => c.to_string(),
        }
    }
}

/// Group string for a read without the tag (`read_groups.py:49-59`).
pub const MISSING_GROUP: &str = "NA";

/// Borrowed view of one BAM record: what U4/U7/U8 need (design U2).
#[derive(Clone, Copy, Debug)]
pub struct RecordLite<'a> {
    pub qname: &'a [u8],
    pub flag: u16,
    /// 0-based leftmost reference position (`reference_start`).
    pub pos0: i64,
    pub mapq: u8,
    /// Raw BAM CIGAR words (`len << 4 | op`), see [`cigar_decode`].
    pub cigar: &'a [u32],
    /// Sequence decoded to ASCII letters `=ACMGRSVTWYHKDBN` (IUPAC kept, design §6.1).
    pub seq: &'a [u8],
    pub tag: Option<&'a TagValue>,
}

impl RecordLite<'_> {
    pub const fn is_unmapped(&self) -> bool {
        self.flag & flags::UNMAPPED != 0
    }
    pub const fn is_reverse(&self) -> bool {
        self.flag & flags::REVERSE != 0
    }
    pub const fn is_secondary(&self) -> bool {
        self.flag & flags::SECONDARY != 0
    }
    pub const fn is_supplementary(&self) -> bool {
        self.flag & flags::SUPPLEMENTARY != 0
    }
    /// 0-based exclusive reference end (`reference_end`), computed from the CIGAR.
    pub fn end0(&self) -> i64 {
        let mut end = self.pos0;
        for &w in self.cigar {
            if let (Some(op), len) = cigar_decode(w) {
                if op.consumes_ref() {
                    end += len as i64;
                }
            }
        }
        end
    }
}

/// Owned record, for fixtures and tests. `as_lite()` gives the borrowed view.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RecordOwned {
    pub qname: String,
    pub flag: u16,
    pub pos0: i64,
    pub mapq: u8,
    pub cigar: Vec<u32>,
    pub seq: Vec<u8>,
    pub tag: Option<TagValue>,
}

impl RecordOwned {
    pub fn as_lite(&self) -> RecordLite<'_> {
        RecordLite {
            qname: self.qname.as_bytes(),
            flag: self.flag,
            pos0: self.pos0,
            mapq: self.mapq,
            cigar: &self.cigar,
            seq: &self.seq,
            tag: self.tag.as_ref(),
        }
    }
}

/// Output of U4 `process_alignment` (`AlignmentInfo` after `add_polya_info`, plus
/// `extract_relevant_deletions`).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct AlignedRead {
    /// `read_exons` (1-based closed, after polyA/T trimming).
    pub exons: Vec<Iv>,
    /// `read_blocks`: 0-based closed query-coordinate ranges, one per exon.
    pub read_blocks: Vec<(i32, i32)>,
    /// `cigar_blocks`: closed CIGAR-index ranges, one per exon.
    pub cigar_blocks: Vec<(i32, i32)>,
    pub polya: PolyAInfo,
    pub exons_changed: bool,
    /// `(ref_pos, length)`.
    pub deletions: Vec<(Pos, u32)>,
}

/// Borrowed reference sequence window. `seq[0]` is the base at 1-based position `start`.
/// No case folding (pyfaidx default, design §6.1).
#[derive(Clone, Copy, Debug)]
pub struct RefSlice<'a> {
    pub start: Pos,
    pub seq: &'a [u8],
}

impl RefSlice<'_> {
    /// Base at 1-based position `pos`, or `None` outside the window.
    #[inline]
    pub fn base(&self, pos: Pos) -> Option<u8> {
        let off = pos.checked_sub(self.start)?;
        if off < 0 {
            return None;
        }
        self.seq.get(off as usize).copied()
    }
    /// Bases of the closed 1-based interval `iv`, clipped to the window.
    pub fn slice(&self, iv: Iv) -> &[u8] {
        let lo = (iv.s - self.start).max(0) as usize;
        let hi = ((iv.e - self.start + 1).max(0) as usize).min(self.seq.len());
        if lo >= hi {
            &[]
        } else {
            &self.seq[lo..hi]
        }
    }
    /// Last 1-based position covered.
    pub fn end(&self) -> Pos {
        self.start + self.seq.len() as Pos - 1
    }
}

// ------------------------------------------------------------------------------------------
// Pass-A regions (U3, consumed by U8/U21)
// ------------------------------------------------------------------------------------------

/// 4.1 `AbstractAlignmentStorage.COVERAGE_BIN` (`alignment_processor.py:116`).
pub const COVERAGE_BIN: i64 = 256;
/// `AlignmentCollector.ABS_COV_VALLEY` (`alignment_processor.py:258`).
pub const ABS_COV_VALLEY: u32 = 1;
/// `AlignmentCollector.REL_COV_VALLEY` (`alignment_processor.py:259`).
pub const REL_COV_VALLEY: f64 = 0.01;
/// `AlignmentCollector.SMALL_CHR_IDS` (`alignment_processor.py:261`).
pub const SMALL_CHR_IDS: [&str; 3] = ["MT", "chrM", "chrMT"];
/// `AlignmentCollector.SMALL_CHR_LEN` (`alignment_processor.py:262`).
pub const SMALL_CHR_LEN: u64 = 500_000;

/// One block (4.1 "split region") produced by pass A.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Block {
    /// 0-based closed region (4.1 tuple).
    pub region: Iv0,
    /// Max coverage from `split_coverage_regions`.
    pub max_cov: u32,
    /// `skip_read_fraction` of the coverage cap (1 = keep every read).
    pub skip: u32,
    /// Genes overlapping the block, sorted by (start, file order).
    pub gene_ids: Vec<GeneIdx>,
    /// Cost estimate: primary records owned by this block (max overlap) plus an upper bound
    /// on the boundary-straddling reads it also processes (4.1 double processing; D-R2 retired).
    pub owned_reads: u64,
}

/// A pass-A cluster (4.1 `alignment_storage` region between gene-separated gaps).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Cluster {
    pub region: Iv0,
    pub blocks: Vec<Block>,
    pub n_records: u64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cigar_roundtrip() {
        for c in 0..9u8 {
            let op = CigarOp::from_code(c).unwrap();
            assert_eq!(cigar_decode(cigar_encode(op, 1234)), (Some(op), 1234));
        }
        assert_eq!(CigarOp::from_code(9), None);
    }

    #[test]
    fn record_end() {
        let r = RecordOwned {
            qname: "r".into(),
            flag: 0,
            pos0: 100,
            mapq: 60,
            cigar: vec![
                cigar_encode(CigarOp::S, 5),
                cigar_encode(CigarOp::M, 10),
                cigar_encode(CigarOp::N, 100),
                cigar_encode(CigarOp::I, 3),
                cigar_encode(CigarOp::D, 2),
                cigar_encode(CigarOp::M, 10),
            ],
            seq: b"A".repeat(28),
            tag: Some(TagValue::Int(1)),
        };
        assert_eq!(r.as_lite().end0(), 222);
        assert_eq!(r.as_lite().tag.unwrap().to_group_string(), "1");
    }

    #[test]
    fn ref_slice() {
        let s = RefSlice { start: 10, seq: b"ACGTACGT" };
        assert_eq!(s.base(10), Some(b'A'));
        assert_eq!(s.base(17), Some(b'T'));
        assert_eq!(s.base(18), None);
        assert_eq!(s.base(9), None);
        assert_eq!(s.slice(Iv::new(8, 11)), b"AC");
        assert_eq!(s.end(), 17);
    }
}
