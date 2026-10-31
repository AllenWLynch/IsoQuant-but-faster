//! `common.get_read_blocks` (`common.py:553-621`), ported as a state machine over raw ops.
//!
//! Semantics kept from 4.1:
//! - a block opens at the first M/I/D/=/X after a gap (so a leading D/I is part of the block);
//! - N and S close the open block; it is emitted only if it contained a match op;
//! - H and P are ignored (but still advance the CIGAR index);
//! - trailing D/I stay inside the last block (`ref_pos - 1` includes a trailing deletion).

use iq_types::read::cigar_decode;
use iq_types::{CigarOp, Iv, Pos};

/// Reference blocks (1-based closed), read blocks (0-based closed query ranges) and
/// CIGAR-index blocks (closed index ranges), one entry per exon.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ReadBlocks {
    pub ref_blocks: Vec<Iv>,
    pub read_blocks: Vec<(i32, i32)>,
    pub cigar_blocks: Vec<(i32, i32)>,
}

impl ReadBlocks {
    pub fn clear(&mut self) {
        self.ref_blocks.clear();
        self.read_blocks.clear();
        self.cigar_blocks.clear();
    }
    pub fn len(&self) -> usize {
        self.ref_blocks.len()
    }
    pub fn is_empty(&self) -> bool {
        self.ref_blocks.is_empty()
    }
}

/// `get_read_blocks(ref_start, cigar_tuples)` with `pos0 = reference_start`.
pub fn read_blocks(pos0: i64, cigar: &[u32]) -> ReadBlocks {
    let mut b = ReadBlocks::default();
    read_blocks_into(pos0, cigar, &mut b.ref_blocks, &mut b.read_blocks, &mut b.cigar_blocks);
    b
}

/// Same as [`read_blocks`], but appends into caller-owned buffers (cleared first).
pub fn read_blocks_into(
    pos0: i64,
    cigar: &[u32],
    ref_blocks: &mut Vec<Iv>,
    read_blocks: &mut Vec<(i32, i32)>,
    cigar_blocks: &mut Vec<(i32, i32)>,
) {
    ref_blocks.clear();
    read_blocks.clear();
    cigar_blocks.clear();
    let mut read_pos: i64 = 0;
    let mut ref_pos: i64 = pos0 + 1;
    // (ref start, read start, cigar index start) of the open block
    let mut open: Option<(i64, i64, usize)> = None;
    let mut has_match = false;

    let mut close =
        |open: &mut Option<(i64, i64, usize)>, has_match: &mut bool, ref_pos: i64, read_pos: i64, idx: usize| {
            if let Some((rs, qs, cs)) = open.take() {
                if *has_match {
                    ref_blocks.push(Iv::new(rs as Pos, (ref_pos - 1) as Pos));
                    read_blocks.push((qs as i32, (read_pos - 1) as i32));
                    cigar_blocks.push((cs as i32, idx as i32 - 1));
                }
            }
            *has_match = false;
        };

    for (idx, &w) in cigar.iter().enumerate() {
        let (op, len) = cigar_decode(w);
        let len = len as i64;
        let Some(op) = op else { continue };
        match op {
            CigarOp::M | CigarOp::Eq | CigarOp::X => {
                if open.is_none() {
                    open = Some((ref_pos, read_pos, idx));
                }
                read_pos += len;
                ref_pos += len;
                has_match = true;
            }
            CigarOp::I => {
                if open.is_none() {
                    open = Some((ref_pos, read_pos, idx));
                }
                read_pos += len;
            }
            CigarOp::D => {
                if open.is_none() {
                    open = Some((ref_pos, read_pos, idx));
                }
                ref_pos += len;
            }
            CigarOp::N => {
                close(&mut open, &mut has_match, ref_pos, read_pos, idx);
                ref_pos += len;
            }
            CigarOp::S => {
                close(&mut open, &mut has_match, ref_pos, read_pos, idx);
                read_pos += len;
            }
            CigarOp::H | CigarOp::P => {}
        }
    }
    if open.is_some() && has_match {
        close(&mut open, &mut has_match, ref_pos, read_pos, cigar.len());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use iq_types::read::cigar_encode;
    use CigarOp::*;

    fn cig(ops: &[(CigarOp, u32)]) -> Vec<u32> {
        ops.iter().map(|&(o, l)| cigar_encode(o, l)).collect()
    }

    #[test]
    fn spliced_with_clips() {
        let b = read_blocks(99, &cig(&[(S, 5), (M, 10), (N, 100), (M, 3), (I, 2), (M, 4), (S, 1)]));
        assert_eq!(b.ref_blocks, vec![Iv::new(100, 109), Iv::new(210, 216)]);
        assert_eq!(b.read_blocks, vec![(5, 14), (15, 23)]);
        assert_eq!(b.cigar_blocks, vec![(1, 1), (3, 5)]);
    }

    #[test]
    fn leading_and_trailing_indels_and_dropped_blocks() {
        // leading D belongs to the block; trailing D extends it; a D-only block is dropped
        let b = read_blocks(0, &cig(&[(H, 3), (D, 2), (M, 5), (D, 3), (N, 10), (D, 4), (N, 10), (M, 2)]));
        assert_eq!(b.ref_blocks, vec![Iv::new(1, 10), Iv::new(35, 36)]);
        assert_eq!(b.cigar_blocks, vec![(1, 3), (7, 7)]);
        assert!(read_blocks(0, &cig(&[(S, 4), (I, 3)])).is_empty());
    }
}
