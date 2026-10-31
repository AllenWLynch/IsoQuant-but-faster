//! Columnar per-block read storage (design §2.3): [`BlockStore`], [`ChainTable`] interning
//! and the [`ReadView`] accessor that construct-side algorithms (U10-U19) are written
//! against.
//!
//! ## Exon-list encoding
//! An exon list `[(s0,e0), (s1,e1), …, (sk,ek)]` is stored as `(chain, start=s0, end=ek)`,
//! where `chain` interns the **gaps** `[(e0+1, s1-1), …, (e_{k-1}+1, sk-1)]`. For ordinary
//! reads the gaps are exactly the introns (`junctions_from_blocks`). If two consecutive
//! exons touch or overlap, the gap is empty/inverted and is *not* an intron; the encoding
//! is still lossless, and [`ReadView::corrected_introns`] filters such gaps out. A chain
//! with no such gap is *pure* ([`ChainTable::is_pure`]) and then
//! [`ReadView::corrected_introns_slice`] returns it without copying.
//! Mono-exon reads use chain [`ChainTable::MONO`] (= 0, empty). An empty exon list uses
//! [`ChainTable::NONE`].

use crate::assign::{Assignment, GeneKey, IsoRef, IsoformMatch, MatchEvent};
use crate::enums::{AssignmentType, MatchClass};
use crate::prim::{Iv, PolyAInfo, Pos, Strand};
use crate::read::ReadRecord;
use rustc_hash::FxHashMap;
use smallvec::SmallVec;

/// xxh3-64 of a QNAME (duplicate-name detection, §2.3, §5.5).
#[inline]
pub fn name_hash(qname: &[u8]) -> u64 {
    xxhash_rust::xxh3::xxh3_64(qname)
}

fn chain_hash(chain: &[Iv]) -> u64 {
    let mut h = xxhash_rust::xxh3::Xxh3::new();
    for iv in chain {
        h.update(&iv.s.to_le_bytes());
        h.update(&iv.e.to_le_bytes());
    }
    h.digest()
}

// ------------------------------------------------------------------------------------------
// ChainTable
// ------------------------------------------------------------------------------------------

/// Per-block hash-consed table of gap chains (flat `Vec<Iv>` + offsets). Reads sharing an
/// intron chain store it once.
#[derive(Clone, Debug)]
pub struct ChainTable {
    flat: Vec<Iv>,
    /// `offs[id]..offs[id+1]` is chain `id`; `offs.len() == n_chains + 1`.
    offs: Vec<u32>,
    pure: Vec<bool>,
    index: FxHashMap<u64, SmallVec<[u32; 1]>>,
}

impl Default for ChainTable {
    fn default() -> Self {
        Self::new()
    }
}

impl ChainTable {
    /// Chain id of the empty chain (mono-exon lists).
    pub const MONO: u32 = 0;
    /// Sentinel chain id for an empty exon list (never stored in the table).
    pub const NONE: u32 = u32::MAX;

    pub fn new() -> Self {
        let mut index = FxHashMap::default();
        index.insert(chain_hash(&[]), SmallVec::from_slice(&[0u32]));
        ChainTable { flat: Vec::new(), offs: vec![0, 0], pure: vec![true], index }
    }

    /// Number of distinct chains (including the empty chain 0).
    pub fn len(&self) -> usize {
        self.offs.len() - 1
    }
    /// Always false: the empty chain is always present.
    pub fn is_empty(&self) -> bool {
        false
    }
    /// Total stored gap intervals.
    pub fn n_intervals(&self) -> usize {
        self.flat.len()
    }

    /// Interns `chain`, returning its id; identical chains get identical ids.
    pub fn intern(&mut self, chain: &[Iv]) -> u32 {
        let h = chain_hash(chain);
        if let Some(cands) = self.index.get(&h) {
            for &id in cands {
                if self.get(id) == chain {
                    return id;
                }
            }
        }
        let id = self.len() as u32;
        self.flat.extend_from_slice(chain);
        self.offs.push(self.flat.len() as u32);
        self.pure.push(chain.iter().all(|g| g.s <= g.e));
        self.index.entry(h).or_default().push(id);
        id
    }

    /// The stored gaps of chain `id`.
    #[inline]
    pub fn get(&self, id: u32) -> &[Iv] {
        let i = id as usize;
        &self.flat[self.offs[i] as usize..self.offs[i + 1] as usize]
    }

    /// True if every gap of chain `id` is a real intron (`s <= e`).
    #[inline]
    pub fn is_pure(&self, id: u32) -> bool {
        self.pure[id as usize]
    }

    /// Encodes a sorted exon list as `(chain, start, end)`.
    pub fn encode(&mut self, exons: &[Iv]) -> (u32, Pos, Pos) {
        match exons {
            [] => (Self::NONE, 0, 0),
            [only] => (Self::MONO, only.s, only.e),
            _ => {
                let gaps: SmallVec<[Iv; 16]> = exons.windows(2).map(|w| Iv::new(w[0].e + 1, w[1].s - 1)).collect();
                (self.intern(&gaps), exons[0].s, exons[exons.len() - 1].e)
            }
        }
    }

    /// Decodes `(chain, start, end)` back into the exon list, appending to `out`.
    pub fn decode_into(&self, chain: u32, start: Pos, end: Pos, out: &mut impl Extend<Iv>) {
        if chain == Self::NONE {
            return;
        }
        let gaps = self.get(chain);
        let mut s = start;
        let mut v: SmallVec<[Iv; 16]> = SmallVec::with_capacity(gaps.len() + 1);
        for g in gaps {
            v.push(Iv::new(s, g.s - 1));
            s = g.e + 1;
        }
        v.push(Iv::new(s, end));
        out.extend(v);
    }

    /// Number of exons of an encoded list.
    #[inline]
    pub fn n_exons(&self, chain: u32) -> usize {
        if chain == Self::NONE {
            0
        } else {
            self.get(chain).len() + 1
        }
    }
}

// ------------------------------------------------------------------------------------------
// BlockStore
// ------------------------------------------------------------------------------------------

/// Bit layout of [`BlockStore`] flags.
pub mod flag_bits {
    pub const STRAND_SHIFT: u16 = 0; // 2 bits
    pub const MAPPED_STRAND_SHIFT: u16 = 2; // 2 bits
    pub const POLYA_FOUND: u16 = 1 << 4;
    pub const MULTIMAPPER: u16 = 1 << 5;
    pub const HAS_DELS: u16 = 1 << 6;
}
use flag_bits::*;

/// One stored isoform match (row of the match table).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MatchRow {
    pub gene: Option<GeneKey>,
    pub tx: Option<IsoRef>,
    pub tx_strand: Strand,
    pub class: MatchClass,
    /// Stored as f64 (exact; the §2.3 f32 compaction is a later optimisation behind this API).
    pub penalty: f64,
    pub ev_off: u32,
    pub n_ev: u32,
}

/// Columnar storage of one block's reads (design §2.3). Row `i` is the read with block
/// index `i` unless the records given to [`BlockStore::from_records`] had non-dense `idx`
/// values (sampled fixtures), in which case the original `idx` is kept in a side column.
#[derive(Clone, Debug, Default)]
pub struct BlockStore {
    /// Block region; every record's `genomic_region` unless listed in `region_exc`.
    pub region: Iv,
    region_exc: Vec<(u32, Iv)>,
    idx: Option<Vec<u32>>,

    name_hash: Vec<u64>,
    names: String,
    name_off: Vec<u32>,

    pub chains: ChainTable,
    chain: Vec<u32>,
    start: Vec<Pos>,
    end: Vec<Pos>,
    /// 0 = raw exons equal corrected exons; otherwise `1 + index` into `raw_side`.
    raw: Vec<u32>,
    raw_side: Vec<(u32, Pos, Pos)>,

    polya: Vec<PolyAInfo>,
    flags: Vec<u16>,
    mapq: Vec<u8>,
    group: Vec<u32>,
    atype: Vec<AssignmentType>,
    gene_atype: Vec<AssignmentType>,

    match_off: Vec<u32>,
    matches: Vec<MatchRow>,
    events: Vec<MatchEvent>,

    dels_off: Vec<u32>,
    dels: Vec<(Pos, u32)>,
}

impl BlockStore {
    /// Empty store for a block with region `region`.
    pub fn new(region: Iv) -> Self {
        BlockStore { region, name_off: vec![0], match_off: vec![0], dels_off: vec![0], ..Default::default() }
    }

    /// Number of reads.
    #[inline]
    pub fn len(&self) -> usize {
        self.flags.len()
    }
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.flags.is_empty()
    }

    /// Builds a store from records in block order. The block region is taken from the
    /// first record's `genomic_region` (or `(0,0)` if empty).
    pub fn from_records(records: &[ReadRecord]) -> Self {
        let region = records.first().map(|r| r.genomic_region).unwrap_or_default();
        let mut st = BlockStore::new(region);
        for r in records {
            st.push(r);
        }
        st
    }

    /// Materializes every record (inverse of [`Self::from_records`]).
    pub fn to_records(&self) -> Vec<ReadRecord> {
        (0..self.len()).map(|i| self.record(i)).collect()
    }

    /// Appends one record. Records must be pushed in block order.
    pub fn push(&mut self, r: &ReadRecord) {
        let i = self.len() as u32;
        // idx: keep a side column only once a non-dense idx is seen
        match &mut self.idx {
            Some(col) => col.push(r.idx),
            None if r.idx != i => {
                let mut col: Vec<u32> = (0..i).collect();
                col.push(r.idx);
                self.idx = Some(col);
            }
            None => {}
        }
        if r.genomic_region != self.region {
            self.region_exc.push((i, r.genomic_region));
        }
        self.name_hash.push(name_hash(r.name.as_bytes()));
        self.names.push_str(&r.name);
        self.name_off.push(self.names.len() as u32);

        let (c, s, e) = self.chains.encode(&r.corrected_exons);
        self.chain.push(c);
        self.start.push(s);
        self.end.push(e);
        if r.exons == r.corrected_exons {
            self.raw.push(0);
        } else {
            let enc = self.chains.encode(&r.exons);
            self.raw_side.push(enc);
            self.raw.push(self.raw_side.len() as u32);
        }

        self.polya.push(r.polya);
        let mut f = (r.strand.bits() << STRAND_SHIFT) | (r.mapped_strand.bits() << MAPPED_STRAND_SHIFT);
        if r.polya_found {
            f |= POLYA_FOUND;
        }
        if r.multimapper {
            f |= MULTIMAPPER;
        }
        if !r.boundary_deletions.is_empty() {
            f |= HAS_DELS;
        }
        self.flags.push(f);
        self.mapq.push(r.mapq);
        self.group.push(r.group);
        self.atype.push(r.assignment.atype);
        self.gene_atype.push(r.assignment.gene_atype);

        for m in &r.assignment.matches {
            let ev_off = self.events.len() as u32;
            self.events.extend_from_slice(&m.events);
            self.matches.push(MatchRow {
                gene: m.gene,
                tx: m.tx,
                tx_strand: m.tx_strand,
                class: m.class,
                penalty: m.penalty,
                ev_off,
                n_ev: m.events.len() as u32,
            });
        }
        self.match_off.push(self.matches.len() as u32);

        self.dels.extend_from_slice(&r.boundary_deletions);
        self.dels_off.push(self.dels.len() as u32);
    }

    /// Accessor for read `i`.
    #[inline]
    pub fn view(&self, i: usize) -> ReadView<'_> {
        debug_assert!(i < self.len());
        ReadView { st: self, i }
    }

    /// All reads in block order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = ReadView<'_>> + '_ {
        (0..self.len()).map(move |i| self.view(i))
    }

    /// Materializes read `i` as a [`ReadRecord`].
    pub fn record(&self, i: usize) -> ReadRecord {
        let v = self.view(i);
        ReadRecord {
            idx: v.idx(),
            name: v.name().to_owned(),
            exons: v.exons().into_vec(),
            corrected_exons: v.corrected_exons().into_vec(),
            boundary_deletions: v.boundary_deletions().to_vec(),
            polya: v.polya(),
            polya_found: v.polya_found(),
            strand: v.strand(),
            mapped_strand: v.mapped_strand(),
            mapq: v.mapq(),
            multimapper: v.multimapper(),
            group: v.group(),
            assignment: v.assignment(),
            genomic_region: v.genomic_region(),
        }
    }

    // ---- barrier patches (U9 / §5.4) ----

    /// Overwrites `assignment_type` / `gene_assignment_type` of read `i`.
    pub fn set_types(&mut self, i: usize, atype: AssignmentType, gene_atype: AssignmentType) {
        self.atype[i] = atype;
        self.gene_atype[i] = gene_atype;
    }
    /// Sets the multimapper flag of read `i`.
    pub fn set_multimapper(&mut self, i: usize, v: bool) {
        if v {
            self.flags[i] |= MULTIMAPPER;
        } else {
            self.flags[i] &= !MULTIMAPPER;
        }
    }

    /// Approximate heap bytes held (for admission-control calibration, §5.4).
    pub fn heap_bytes(&self) -> usize {
        use std::mem::size_of;
        self.name_hash.capacity() * 8
            + self.names.capacity()
            + self.name_off.capacity() * 4
            + self.chains.flat.capacity() * size_of::<Iv>()
            + self.chains.offs.capacity() * 4
            + (self.chain.capacity() + self.raw.capacity()) * 4
            + (self.start.capacity() + self.end.capacity()) * 4
            + self.raw_side.capacity() * 12
            + self.polya.capacity() * size_of::<PolyAInfo>()
            + self.flags.capacity() * 2
            + self.mapq.capacity()
            + self.group.capacity() * 4
            + self.atype.capacity() * 2
            + self.match_off.capacity() * 4
            + self.matches.capacity() * size_of::<MatchRow>()
            + self.events.capacity() * size_of::<MatchEvent>()
            + self.dels_off.capacity() * 4
            + self.dels.capacity() * 8
    }

    /// Resident bytes by element count (not capacity). Use for admission control: `heap_bytes` counts
    /// allocated capacity and overstates resident memory by 1.35-2x (docs/vtg_proto.md).
    pub fn len_bytes(&self) -> usize {
        use std::mem::size_of;
        self.name_hash.len() * 8
            + self.names.len()
            + self.name_off.len() * 4
            + self.chains.flat.len() * size_of::<Iv>()
            + self.chains.offs.len() * 4
            + (self.chain.len() + self.raw.len()) * 4
            + (self.start.len() + self.end.len()) * 4
            + self.raw_side.len() * 12
            + self.polya.len() * size_of::<PolyAInfo>()
            + self.flags.len() * 2
            + self.mapq.len()
            + self.group.len() * 4
            + self.atype.len() * 2
            + self.match_off.len() * 4
            + self.matches.len() * size_of::<MatchRow>()
            + self.events.len() * size_of::<MatchEvent>()
            + self.dels_off.len() * 4
            + self.dels.len() * 8
    }
}

/// Borrowed accessor for one read of a [`BlockStore`].
#[derive(Clone, Copy)]
pub struct ReadView<'a> {
    st: &'a BlockStore,
    i: usize,
}

/// Small materialized exon/intron list.
pub type IvVec = SmallVec<[Iv; 8]>;

impl<'a> ReadView<'a> {
    /// Row in the store (= block read index for dense stores).
    #[inline]
    pub fn row(&self) -> usize {
        self.i
    }
    /// `ReadRecord.idx`.
    #[inline]
    pub fn idx(&self) -> u32 {
        match &self.st.idx {
            Some(c) => c[self.i],
            None => self.i as u32,
        }
    }
    pub fn name(&self) -> &'a str {
        let o = &self.st.name_off;
        &self.st.names[o[self.i] as usize..o[self.i + 1] as usize]
    }
    #[inline]
    pub fn name_hash(&self) -> u64 {
        self.st.name_hash[self.i]
    }
    pub fn genomic_region(&self) -> Iv {
        match self.st.region_exc.binary_search_by_key(&(self.i as u32), |x| x.0) {
            Ok(k) => self.st.region_exc[k].1,
            Err(_) => self.st.region,
        }
    }

    // ---- exons ----

    /// Interned chain id of `corrected_exons` (equal ids ⇔ equal gap chains).
    #[inline]
    pub fn chain_id(&self) -> u32 {
        self.st.chain[self.i]
    }
    /// First corrected exon start (`corrected_exons[0][0]`); 0 if there are no exons.
    #[inline]
    pub fn start(&self) -> Pos {
        self.st.start[self.i]
    }
    /// Last corrected exon end (`corrected_exons[-1][1]`); 0 if there are no exons.
    #[inline]
    pub fn end(&self) -> Pos {
        self.st.end[self.i]
    }
    pub fn n_corrected_exons(&self) -> usize {
        self.st.chains.n_exons(self.chain_id())
    }
    pub fn corrected_exons(&self) -> IvVec {
        let mut v = IvVec::new();
        self.st.chains.decode_into(self.chain_id(), self.start(), self.end(), &mut v);
        v
    }
    /// `corrected_introns` = `junctions_from_blocks(corrected_exons)`.
    pub fn corrected_introns(&self) -> IvVec {
        let c = self.chain_id();
        if c == ChainTable::NONE {
            return IvVec::new();
        }
        self.st.chains.get(c).iter().copied().filter(|g| g.s <= g.e).collect()
    }
    /// Zero-copy intron slice when the chain is pure (the usual case).
    pub fn corrected_introns_slice(&self) -> Option<&'a [Iv]> {
        let c = self.chain_id();
        if c == ChainTable::NONE {
            Some(&[])
        } else if self.st.chains.is_pure(c) {
            Some(self.st.chains.get(c))
        } else {
            None
        }
    }
    /// True if raw `exons` differ from `corrected_exons`.
    #[inline]
    pub fn exons_differ(&self) -> bool {
        self.st.raw[self.i] != 0
    }
    /// `ReadAssignment.exons` (raw, post polyA trimming).
    pub fn exons(&self) -> IvVec {
        match self.st.raw[self.i] {
            0 => self.corrected_exons(),
            k => {
                let (c, s, e) = self.st.raw_side[k as usize - 1];
                let mut v = IvVec::new();
                self.st.chains.decode_into(c, s, e, &mut v);
                v
            }
        }
    }
    /// `(exons[0][0], exons[-1][1])` without materializing; `(0, 0)` if empty.
    pub fn exons_span(&self) -> (Pos, Pos) {
        match self.st.raw[self.i] {
            0 => (self.start(), self.end()),
            k => {
                let (_, s, e) = self.st.raw_side[k as usize - 1];
                (s, e)
            }
        }
    }

    // ---- scalars ----

    #[inline]
    pub fn polya(&self) -> PolyAInfo {
        self.st.polya[self.i]
    }
    #[inline]
    pub fn polya_found(&self) -> bool {
        self.st.flags[self.i] & POLYA_FOUND != 0
    }
    #[inline]
    pub fn multimapper(&self) -> bool {
        self.st.flags[self.i] & MULTIMAPPER != 0
    }
    #[inline]
    pub fn strand(&self) -> Strand {
        Strand::from_bits(self.st.flags[self.i] >> STRAND_SHIFT)
    }
    #[inline]
    pub fn mapped_strand(&self) -> Strand {
        Strand::from_bits(self.st.flags[self.i] >> MAPPED_STRAND_SHIFT)
    }
    #[inline]
    pub fn mapq(&self) -> u8 {
        self.st.mapq[self.i]
    }
    #[inline]
    pub fn group(&self) -> u32 {
        self.st.group[self.i]
    }
    pub fn boundary_deletions(&self) -> &'a [(Pos, u32)] {
        let o = &self.st.dels_off;
        &self.st.dels[o[self.i] as usize..o[self.i + 1] as usize]
    }

    // ---- assignment ----

    #[inline]
    pub fn atype(&self) -> AssignmentType {
        self.st.atype[self.i]
    }
    #[inline]
    pub fn gene_atype(&self) -> AssignmentType {
        self.st.gene_atype[self.i]
    }
    pub fn n_matches(&self) -> usize {
        (self.st.match_off[self.i + 1] - self.st.match_off[self.i]) as usize
    }
    /// The isoform matches, in order.
    pub fn matches(&self) -> impl ExactSizeIterator<Item = MatchRef<'a>> + 'a {
        let st = self.st;
        let lo = st.match_off[self.i] as usize;
        let hi = st.match_off[self.i + 1] as usize;
        st.matches[lo..hi].iter().map(move |row| MatchRef { row, st })
    }
    /// `isoform_matches[0]` if any.
    pub fn first_match(&self) -> Option<MatchRef<'a>> {
        self.matches().next()
    }
    /// Materialized [`Assignment`].
    pub fn assignment(&self) -> Assignment {
        Assignment {
            atype: self.atype(),
            gene_atype: self.gene_atype(),
            matches: self.matches().map(|m| m.to_match()).collect(),
        }
    }
}

impl std::fmt::Debug for ReadView<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ReadView").field("idx", &self.idx()).field("name", &self.name()).finish()
    }
}

/// Borrowed isoform match of a stored read.
#[derive(Clone, Copy)]
pub struct MatchRef<'a> {
    row: &'a MatchRow,
    st: &'a BlockStore,
}

impl<'a> MatchRef<'a> {
    #[inline]
    pub fn gene(&self) -> Option<GeneKey> {
        self.row.gene
    }
    #[inline]
    pub fn tx(&self) -> Option<IsoRef> {
        self.row.tx
    }
    #[inline]
    pub fn tx_strand(&self) -> Strand {
        self.row.tx_strand
    }
    #[inline]
    pub fn class(&self) -> MatchClass {
        self.row.class
    }
    #[inline]
    pub fn penalty(&self) -> f64 {
        self.row.penalty
    }
    pub fn events(&self) -> &'a [MatchEvent] {
        let lo = self.row.ev_off as usize;
        &self.st.events[lo..lo + self.row.n_ev as usize]
    }
    pub fn to_match(&self) -> IsoformMatch {
        IsoformMatch {
            gene: self.gene(),
            tx: self.tx(),
            tx_strand: self.tx_strand(),
            class: self.class(),
            penalty: self.penalty(),
            events: self.events().iter().copied().collect(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enums::EventType;

    #[test]
    fn chain_encoding_lossless_with_touching_exons() {
        let mut t = ChainTable::new();
        let ex = [Iv::new(1, 10), Iv::new(11, 20), Iv::new(15, 30), Iv::new(40, 50)];
        let (c, s, e) = t.encode(&ex);
        assert!(!t.is_pure(c));
        let mut out = Vec::new();
        t.decode_into(c, s, e, &mut out);
        assert_eq!(out, ex);
        assert_eq!(t.encode(&[]).0, ChainTable::NONE);
        assert_eq!(t.encode(&[Iv::new(3, 4)]).0, ChainTable::MONO);
    }

    #[test]
    fn view_basics() {
        let mut a = Assignment::empty(AssignmentType::Unique);
        a.matches.push(IsoformMatch::new(
            MatchClass::FullSpliceMatch,
            Some(GeneKey::Ref(3)),
            Some(IsoRef::Ref(7)),
            [MatchEvent::new(EventType::Fsm)],
            Strand::Plus,
            0.0,
        ));
        let r = ReadRecord {
            idx: 0,
            name: "read1".into(),
            exons: vec![Iv::new(100, 200), Iv::new(300, 400)],
            corrected_exons: vec![Iv::new(101, 200), Iv::new(300, 400)],
            boundary_deletions: vec![(150, 3)],
            polya: PolyAInfo { ext_a: 400, ext_t: -1, int_a: -1, int_t: -1 },
            polya_found: true,
            strand: Strand::Plus,
            mapped_strand: Strand::Minus,
            mapq: 60,
            multimapper: false,
            group: 2,
            assignment: a,
            genomic_region: Iv::new(1, 1000),
        };
        let st = BlockStore::from_records(std::slice::from_ref(&r));
        let v = st.view(0);
        assert_eq!(v.corrected_introns_slice().unwrap(), &[Iv::new(201, 299)]);
        assert_eq!(v.exons_span(), (100, 400));
        assert!(v.exons_differ());
        assert_eq!(v.first_match().unwrap().events()[0].ty, EventType::Fsm);
        assert_eq!(st.to_records(), vec![r]);
    }
}
