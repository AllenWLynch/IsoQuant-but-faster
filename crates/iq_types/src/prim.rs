//! Primitive and shared scalar types.
//!
//! Coordinates follow IsoQuant: genomic intervals are **1-based and closed**. BAM positions
//! are converted on decode (0-based `pos` becomes `pos + 1`). The only 0-based intervals in
//! the system are the pass-A cluster/block regions of U3 ([`Iv0`]), which mirror 4.1's
//! 0-based `(pos, end-1)` tuples.

use serde::de::{self, Deserializer, SeqAccess, Visitor};
use serde::ser::{SerializeTuple, Serializer};
use serde::{Deserialize, Serialize};
use std::fmt;

/// 1-based closed genomic coordinate. `-1` means "absent" wherever IsoQuant uses `-1`
/// (e.g. [`PolyAInfo`] fields).
pub type Pos = i32;
/// Index into `AnnotationIndex::chroms` (U1).
pub type ChrIdx = u16;
/// Index into `AnnotationIndex::genes` in production; into the fixture gene [`crate::Interner`]
/// when running fixtures (see `fixture::IdMaps`).
pub type GeneIdx = u32;
/// Index into `AnnotationIndex::transcripts` in production; into the fixture transcript
/// [`crate::Interner`] when running fixtures.
pub type TxIdx = u32;
/// Index into a block's [`crate::ModelStore`] = model creation order within the block.
pub type ModelIdx = u32;
/// Read index within a block (BAM order = 4.1 dump order).
pub type ReadIdx = u32;

/// Closed interval `[s, e]` with 1-based coordinates (IsoQuant's `(start, end)` tuples).
///
/// `Ord` is lexicographic `(s, e)`, identical to Python tuple ordering.
/// Serialized as a two-element JSON array `[s, e]`.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Default)]
pub struct Iv {
    pub s: Pos,
    pub e: Pos,
}

impl Iv {
    #[inline]
    pub const fn new(s: Pos, e: Pos) -> Self {
        Iv { s, e }
    }
    /// Number of bases, `e - s + 1` (IsoQuant `interval_len`).
    #[inline]
    pub const fn len(&self) -> i64 {
        self.e as i64 - self.s as i64 + 1
    }
    /// True if `e < s` (an empty or inverted interval).
    #[inline]
    pub const fn is_empty(&self) -> bool {
        self.e < self.s
    }
    /// IsoQuant `overlaps(a, b)`: `not (a.e < b.s or b.e < a.s)` for closed intervals.
    #[inline]
    pub const fn overlaps(&self, o: &Iv) -> bool {
        !(self.e < o.s || o.e < self.s)
    }
    /// IsoQuant `contains(bigger, smaller)`: `bigger.s <= smaller.s and bigger.e >= smaller.e`.
    #[inline]
    pub const fn contains(&self, o: &Iv) -> bool {
        self.s <= o.s && self.e >= o.e
    }
}

impl fmt::Debug for Iv {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {})", self.s, self.e)
    }
}

impl From<(Pos, Pos)> for Iv {
    fn from(t: (Pos, Pos)) -> Self {
        Iv { s: t.0, e: t.1 }
    }
}
impl From<Iv> for (Pos, Pos) {
    fn from(v: Iv) -> Self {
        (v.s, v.e)
    }
}

impl Serialize for Iv {
    fn serialize<S: Serializer>(&self, ser: S) -> Result<S::Ok, S::Error> {
        let mut t = ser.serialize_tuple(2)?;
        t.serialize_element(&self.s)?;
        t.serialize_element(&self.e)?;
        t.end()
    }
}

impl<'de> Deserialize<'de> for Iv {
    fn deserialize<D: Deserializer<'de>>(de: D) -> Result<Self, D::Error> {
        struct V;
        impl<'de> Visitor<'de> for V {
            type Value = Iv;
            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                f.write_str("a two-element array [start, end]")
            }
            fn visit_seq<A: SeqAccess<'de>>(self, mut a: A) -> Result<Iv, A::Error> {
                let s = a.next_element()?.ok_or_else(|| de::Error::invalid_length(0, &self))?;
                let e = a.next_element()?.ok_or_else(|| de::Error::invalid_length(1, &self))?;
                if a.next_element::<de::IgnoredAny>()?.is_some() {
                    return Err(de::Error::invalid_length(3, &self));
                }
                Ok(Iv { s, e })
            }
        }
        de.deserialize_tuple(2, V)
    }
}

/// 0-based closed interval, used only for pass-A regions (U3) where 4.1 keeps 0-based tuples
/// `(pos0, end0 - 1)`. Same serde shape as [`Iv`].
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Default, Serialize, Deserialize)]
#[serde(from = "Iv", into = "Iv")]
pub struct Iv0 {
    pub s: i64,
    pub e: i64,
}
impl From<Iv> for Iv0 {
    fn from(v: Iv) -> Self {
        Iv0 { s: v.s as i64, e: v.e as i64 }
    }
}
impl From<Iv0> for Iv {
    fn from(v: Iv0) -> Self {
        Iv { s: v.s as Pos, e: v.e as Pos }
    }
}
impl fmt::Debug for Iv0 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "0b({}, {})", self.s, self.e)
    }
}

/// Strand: IsoQuant uses the characters `'+'`, `'-'`, `'.'`. Serialized as those strings.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Default, Debug)]
#[repr(u8)]
pub enum Strand {
    Plus = 0,
    Minus = 1,
    #[default]
    Unknown = 2,
}

impl Strand {
    pub const ALL: [Strand; 3] = [Strand::Plus, Strand::Minus, Strand::Unknown];
    #[inline]
    pub const fn as_char(self) -> char {
        match self {
            Strand::Plus => '+',
            Strand::Minus => '-',
            Strand::Unknown => '.',
        }
    }
    #[inline]
    pub const fn as_str(self) -> &'static str {
        match self {
            Strand::Plus => "+",
            Strand::Minus => "-",
            Strand::Unknown => ".",
        }
    }
    pub fn from_char(c: char) -> Option<Strand> {
        match c {
            '+' => Some(Strand::Plus),
            '-' => Some(Strand::Minus),
            '.' => Some(Strand::Unknown),
            _ => None,
        }
    }
    pub fn parse(s: &str) -> Option<Strand> {
        let mut it = s.chars();
        match (it.next(), it.next()) {
            (Some(c), None) => Strand::from_char(c),
            _ => None,
        }
    }
    /// 2-bit code used in `BlockStore` flags.
    #[inline]
    pub const fn bits(self) -> u16 {
        self as u16
    }
    #[inline]
    pub const fn from_bits(b: u16) -> Strand {
        match b & 3 {
            0 => Strand::Plus,
            1 => Strand::Minus,
            _ => Strand::Unknown,
        }
    }
    /// The opposite strand; `.` stays `.`.
    pub const fn flip(self) -> Strand {
        match self {
            Strand::Plus => Strand::Minus,
            Strand::Minus => Strand::Plus,
            Strand::Unknown => Strand::Unknown,
        }
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl Serialize for Strand {
    fn serialize<S: Serializer>(&self, ser: S) -> Result<S::Ok, S::Error> {
        ser.serialize_str(self.as_str())
    }
}
impl<'de> Deserialize<'de> for Strand {
    fn deserialize<D: Deserializer<'de>>(de: D) -> Result<Self, D::Error> {
        let s = <std::borrow::Cow<'de, str>>::deserialize(de)?;
        Strand::parse(&s).ok_or_else(|| de::Error::custom(format!("invalid strand {s:?}")))
    }
}

/// IsoQuant `PolyAInfo` (`terminal_prediction/polya_finder.py:81-86`). `-1` = absent.
///
/// Field order = constructor order: external polyA, external polyT, internal polyA, internal
/// polyT. Serialized as a 4-element array `[ext_a, ext_t, int_a, int_t]`.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Serialize, Deserialize)]
#[serde(from = "[Pos; 4]", into = "[Pos; 4]")]
pub struct PolyAInfo {
    pub ext_a: Pos,
    pub ext_t: Pos,
    pub int_a: Pos,
    pub int_t: Pos,
}

/// All positions absent: `PolyAInfo(-1, -1, -1, -1)`.
pub const NO_POLYA: PolyAInfo = PolyAInfo { ext_a: -1, ext_t: -1, int_a: -1, int_t: -1 };

impl Default for PolyAInfo {
    fn default() -> Self {
        NO_POLYA
    }
}
impl From<[Pos; 4]> for PolyAInfo {
    fn from(a: [Pos; 4]) -> Self {
        PolyAInfo { ext_a: a[0], ext_t: a[1], int_a: a[2], int_t: a[3] }
    }
}
impl From<PolyAInfo> for [Pos; 4] {
    fn from(p: PolyAInfo) -> Self {
        [p.ext_a, p.ext_t, p.int_a, p.int_t]
    }
}

/// Introns of a sorted exon list, IsoQuant `junctions_from_blocks` (`common.py:462-468`):
/// `(e_i + 1, s_{i+1} - 1)` only when `e_i + 1 < s_{i+1}`.
pub fn junctions_from_blocks(exons: &[Iv]) -> Vec<Iv> {
    let mut out = Vec::with_capacity(exons.len().saturating_sub(1));
    for w in exons.windows(2) {
        if w[0].e + 1 < w[1].s {
            out.push(Iv::new(w[0].e + 1, w[1].s - 1));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iv_serde_is_pair() {
        let v = Iv::new(10, 20);
        assert_eq!(serde_json::to_string(&v).unwrap(), "[10,20]");
        assert_eq!(serde_json::from_str::<Iv>("[10,20]").unwrap(), v);
        assert!(serde_json::from_str::<Iv>("[10]").is_err());
        assert!(serde_json::from_str::<Iv>("[10,2,3]").is_err());
    }

    #[test]
    fn strand_codes() {
        for s in Strand::ALL {
            assert_eq!(Strand::from_bits(s.bits()), s);
            assert_eq!(Strand::from_char(s.as_char()), Some(s));
            let j = serde_json::to_string(&s).unwrap();
            assert_eq!(serde_json::from_str::<Strand>(&j).unwrap(), s);
        }
        assert_eq!(Strand::parse("x"), None);
    }

    #[test]
    fn polya_serde() {
        let p = PolyAInfo { ext_a: 5, ext_t: -1, int_a: 7, int_t: -1 };
        assert_eq!(serde_json::to_string(&p).unwrap(), "[5,-1,7,-1]");
        assert_eq!(serde_json::from_str::<PolyAInfo>("[5,-1,7,-1]").unwrap(), p);
    }

    #[test]
    fn junctions() {
        let ex = [Iv::new(1, 10), Iv::new(11, 20), Iv::new(30, 40)];
        assert_eq!(junctions_from_blocks(&ex), vec![Iv::new(21, 29)]);
    }
}
