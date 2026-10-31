//! Work-unit formation (design §5.3).
//!
//! A *legal* work unit is a maximal run of consecutive blocks of one chromosome in which
//! adjacent blocks share at least one gene. Blocks sharing a gene must be constructed in one
//! unit, in genomic order, because `detected_known` (a known isoform is emitted by the first
//! block in file order) is threaded through them. Legal units are independent.
//!
//! On top of the legal split the scheduler may *pack* consecutive small legal units into one
//! scheduling unit (fewer tasks, less bookkeeping). Packing never changes outputs: a packed
//! unit constructs its blocks in genomic order exactly like separate units would, and
//! `detected_known` cannot couple blocks that share no gene.

use iq_types::Block;
use std::ops::Range;

/// How blocks are grouped into scheduling units (`IQRS_DEBUG_UNIT_SPLIT`).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum SplitMode {
    /// Legal units, consecutive small ones packed up to [`PACK_TARGET_READS`].
    #[default]
    Packed,
    /// Legal units, no packing (the finest legal split): `IQRS_DEBUG_UNIT_SPLIT=max`.
    Max,
    /// One unit per chromosome (4.1's granularity): `IQRS_DEBUG_UNIT_SPLIT=whole-chr`.
    WholeChr,
    /// One unit per block, ignoring shared genes. **Illegal** (changes outputs when blocks
    /// share a gene); exists only so tests can show the determinism check has teeth.
    UnsafePerBlock,
}

impl SplitMode {
    /// Parses `IQRS_DEBUG_UNIT_SPLIT` values.
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "" | "packed" | "default" => Ok(SplitMode::Packed),
            "max" => Ok(SplitMode::Max),
            "whole-chr" => Ok(SplitMode::WholeChr),
            "unsafe-per-block" => Ok(SplitMode::UnsafePerBlock),
            v => Err(format!("IQRS_DEBUG_UNIT_SPLIT: unknown mode {v:?} (whole-chr|max|packed)")),
        }
    }

    /// Reads `IQRS_DEBUG_UNIT_SPLIT` (unset = [`SplitMode::Packed`]).
    pub fn from_env() -> Result<Self, String> {
        match std::env::var("IQRS_DEBUG_UNIT_SPLIT") {
            Ok(v) => Self::parse(&v),
            Err(_) => Ok(SplitMode::Packed),
        }
    }
}

/// Packing stops once a packed unit owns this many reads (upper bound; see [`pack_target`]).
pub const PACK_TARGET_READS: u64 = 5_000;
/// Lower bound of [`pack_target`].
pub const PACK_MIN_READS: u64 = 200;

/// Packing target for a chromosome with `chr_reads` owned reads on `threads` workers: about
/// 32 units per worker, clamped to [`PACK_MIN_READS`, `PACK_TARGET_READS`]. (Outputs never
/// depend on it.)
pub fn pack_target(chr_reads: u64, threads: usize) -> u64 {
    (chr_reads / (32 * threads.max(1) as u64)).clamp(PACK_MIN_READS, PACK_TARGET_READS)
}
/// ... or has this many blocks.
pub const PACK_MAX_BLOCKS: usize = 512;

/// `true` when two blocks share a gene (gene lists are sorted by (start, file order), not by
/// index, so this is a small set intersection rather than a merge).
fn share_gene(a: &Block, b: &Block) -> bool {
    if a.gene_ids.is_empty() || b.gene_ids.is_empty() {
        return false;
    }
    // Gene lists are short (typically 0-5); quadratic is fastest here.
    a.gene_ids.iter().any(|g| b.gene_ids.contains(g))
}

/// Legal work units: block index ranges `[lo, hi)` in genomic order (design §5.3).
pub fn form_units(blocks: &[Block]) -> Vec<Range<usize>> {
    let mut out = Vec::new();
    if blocks.is_empty() {
        return out;
    }
    let mut lo = 0;
    for i in 1..blocks.len() {
        if !share_gene(&blocks[i - 1], &blocks[i]) {
            out.push(lo..i);
            lo = i;
        }
    }
    out.push(lo..blocks.len());
    out
}

/// Scheduling units of one chromosome under `mode`.
pub fn schedule_units(blocks: &[Block], mode: SplitMode, target: u64) -> Vec<Range<usize>> {
    match mode {
        SplitMode::WholeChr => {
            if blocks.is_empty() {
                vec![]
            } else {
                std::iter::once(0..blocks.len()).collect()
            }
        }
        SplitMode::UnsafePerBlock => (0..blocks.len()).map(|i| i..i + 1).collect(),
        SplitMode::Max => form_units(blocks),
        SplitMode::Packed => pack(blocks, form_units(blocks), target),
    }
}

fn pack(blocks: &[Block], legal: Vec<Range<usize>>, target: u64) -> Vec<Range<usize>> {
    let reads = |r: &Range<usize>| blocks[r.clone()].iter().map(|b| b.owned_reads).sum::<u64>();
    let mut out: Vec<Range<usize>> = Vec::with_capacity(legal.len());
    let mut cur: Option<(Range<usize>, u64)> = None;
    for r in legal {
        let n = reads(&r);
        cur = match cur {
            Some((c, cn)) if cn + n <= target && (r.end - c.start) <= PACK_MAX_BLOCKS => {
                Some((c.start..r.end, cn + n))
            }
            Some((c, _)) => {
                out.push(c);
                Some((r, n))
            }
            None => Some((r, n)),
        };
    }
    if let Some((c, _)) = cur {
        out.push(c);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use iq_types::Iv0;

    fn blk(genes: &[u32], reads: u64) -> Block {
        Block { region: Iv0 { s: 0, e: 1 }, max_cov: 1, skip: 1, gene_ids: genes.to_vec(), owned_reads: reads }
    }

    #[test]
    fn legal_units() {
        let b = vec![blk(&[], 5), blk(&[1], 5), blk(&[1, 2], 5), blk(&[2], 5), blk(&[3], 5), blk(&[], 5), blk(&[], 5)];
        assert_eq!(form_units(&b), vec![0..1, 1..4, 4..5, 5..6, 6..7]);
        assert_eq!(form_units(&[]), Vec::<Range<usize>>::new());
        // gene lists are not index-sorted
        let b = vec![blk(&[7, 3], 1), blk(&[3], 1)];
        assert_eq!(form_units(&b), vec![0..2]);
    }

    #[test]
    fn packing_respects_legal_boundaries() {
        let b = vec![blk(&[], 10), blk(&[1], 30_000), blk(&[1], 1), blk(&[], 10), blk(&[], 10), blk(&[2], 25_000)];
        let p = schedule_units(&b, SplitMode::Packed, 20_000);
        assert_eq!(p, vec![0..1, 1..3, 3..5, 5..6]);
        // every packed unit is a union of consecutive legal units
        let legal = form_units(&b);
        for u in &p {
            assert!(legal.iter().any(|l| l.start == u.start));
            assert!(legal.iter().any(|l| l.end == u.end));
        }
        assert_eq!(schedule_units(&b, SplitMode::WholeChr, 20_000), vec![0..6]);
        assert_eq!(schedule_units(&b, SplitMode::Max, 20_000), legal);
        assert_eq!(schedule_units(&b, SplitMode::UnsafePerBlock, 20_000).len(), 6);
    }

    #[test]
    fn parse_modes() {
        assert_eq!(SplitMode::parse("whole-chr").unwrap(), SplitMode::WholeChr);
        assert_eq!(SplitMode::parse("max").unwrap(), SplitMode::Max);
        assert!(SplitMode::parse("bogus").is_err());
    }
}
