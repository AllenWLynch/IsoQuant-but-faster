//! `extract_relevant_deletions` (`transcript_splice_site_corrector.py:75-100`).
//!
//! 4.1 checks every deletion against every exon boundary (O(D·E)). Here it is one sweep:
//! boundaries `[s0, e0, s1, e1, ...]` are non-decreasing (exons are sorted and disjoint),
//! and so are the deletion windows' lower bounds `ref_pos - margin`, so one boundary pointer
//! only moves forward.

use iq_types::read::cigar_decode;
use iq_types::{CigarOp, Iv, Pos};

/// `BOUNDARY_MARGIN = WINDOW_SIZE + max(ACCEPTED_DEL_CASES)` (`:59`) = 14.
pub const BOUNDARY_MARGIN: i32 = 14;

/// `extract_relevant_deletions(cigartuples, ref_start, exons, margin)`, appending
/// `(ref_pos, length)` pairs to `out` (not cleared).
pub fn extract_relevant_deletions_into(
    cigar: &[u32],
    ref_start: Pos,
    exons: &[Iv],
    margin: i32,
    out: &mut Vec<(Pos, u32)>,
) {
    if cigar.is_empty() || exons.is_empty() {
        return;
    }
    debug_assert!(exons.windows(2).all(|w| w[0].e <= w[1].s) && exons.iter().all(|e| e.s <= e.e));
    let nb = exons.len() * 2;
    let boundary = |k: usize| if k & 1 == 0 { exons[k >> 1].s } else { exons[k >> 1].e };
    let mut k = 0usize;
    let mut ref_pos = ref_start as i64;
    for &w in cigar {
        let (op, len) = cigar_decode(w);
        match op {
            Some(CigarOp::D) => {
                let lo = ref_pos - margin as i64;
                let hi = ref_pos + len as i64 - 1 + margin as i64;
                while k < nb && (boundary(k) as i64) < lo {
                    k += 1;
                }
                if k < nb && boundary(k) as i64 <= hi {
                    out.push((ref_pos as Pos, len));
                }
                ref_pos += len as i64;
            }
            Some(o) if o.consumes_ref() => ref_pos += len as i64,
            _ => {}
        }
    }
}

/// Allocating convenience wrapper of [`extract_relevant_deletions_into`].
pub fn extract_relevant_deletions(cigar: &[u32], ref_start: Pos, exons: &[Iv], margin: i32) -> Vec<(Pos, u32)> {
    let mut v = Vec::new();
    extract_relevant_deletions_into(cigar, ref_start, exons, margin, &mut v);
    v
}

// ------------------------------------------------------------------------------------------
// Fixture runner
// ------------------------------------------------------------------------------------------

pub(crate) mod fixture {
    use super::*;
    use iq_types::fixture::schema::{U4ExtractDeletionsIn, U4ExtractDeletionsOut};
    use iq_types::fixture::{FixResult, FixtureCtx, UnitFixture};

    pub struct U4Deletions;

    impl UnitFixture for U4Deletions {
        const UNIT: &'static str = "U4";
        const FN: &'static str = "extract_relevant_deletions";
        type In = U4ExtractDeletionsIn;
        type Out = U4ExtractDeletionsOut;

        fn run(i: &U4ExtractDeletionsIn, _ctx: &FixtureCtx) -> FixResult<U4ExtractDeletionsOut> {
            let cigar: Vec<u32> = i.cigar.iter().map(|&(op, len)| (len << 4) | op as u32).collect();
            Ok(U4ExtractDeletionsOut {
                deletions: extract_relevant_deletions(&cigar, i.ref_start as Pos, &i.exons, i.margin),
            })
        }
    }
}

pub(crate) fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    reg.add::<fixture::U4Deletions>();
}

#[cfg(test)]
mod tests {
    use super::*;
    use iq_types::read::cigar_encode;
    use CigarOp::*;

    fn cig(ops: &[(CigarOp, u32)]) -> Vec<u32> {
        ops.iter().map(|&(o, l)| cigar_encode(o, l)).collect()
    }

    /// 4.1's O(D·E) definition.
    fn naive(cigar: &[u32], ref_start: Pos, exons: &[Iv], margin: i32) -> Vec<(Pos, u32)> {
        let b: Vec<Pos> = exons.iter().flat_map(|e| [e.s, e.e]).collect();
        let mut out = vec![];
        let mut r = ref_start;
        for &w in cigar {
            let (op, len) = cigar_decode(w);
            let len_i = len as Pos;
            match op {
                Some(D) => {
                    if b.iter().any(|&x| r - margin <= x && x <= r + len_i - 1 + margin) {
                        out.push((r, len));
                    }
                    r += len_i;
                }
                Some(o) if o.consumes_ref() => r += len_i,
                _ => {}
            }
        }
        out
    }

    #[test]
    fn sweep_equals_naive() {
        let c = cig(&[
            (S, 3),
            (M, 20),
            (D, 1),
            (M, 30),
            (D, 2),
            (M, 5),
            (N, 100),
            (M, 3),
            (D, 5),
            (M, 40),
            (D, 1),
            (M, 2),
            (I, 2),
        ]);
        let rb = crate::cigar::read_blocks(999, &c).ref_blocks;
        for margin in [0, 1, 5, 14, 30] {
            assert_eq!(extract_relevant_deletions(&c, 1000, &rb, margin), naive(&c, 1000, &rb, margin));
        }
        assert!(extract_relevant_deletions(&c, 1000, &[], 14).is_empty());
    }
}
