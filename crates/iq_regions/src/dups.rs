//! Duplicate-name candidates (design §5.1, §5.5): 64-bit name hashes of primary records
//! that occur more than once on a chromosome.
//!
//! A scalable blocked Bloom filter replaces a `Vec<u64>` of every hash (8 B/read): ~1.5 B/read.
//! A Bloom filter has no false negatives, so every repeated hash is reported; false
//! positives (≈0.5-1%) only add candidates, which U9 checks against real names anyway.

/// Bits per expected element of each filter stage.
const BITS_PER_ELEM: usize = 12;
/// Bits set per element inside one 512-bit block.
const K: u32 = 6;
const FIRST_CAPACITY: usize = 1 << 16;

struct Stage {
    /// 512-bit blocks.
    blocks: Vec<[u64; 8]>,
    capacity: usize,
    len: usize,
}

impl Stage {
    fn new(capacity: usize) -> Stage {
        let n_blocks = (capacity * BITS_PER_ELEM).div_ceil(512).next_power_of_two();
        Stage { blocks: vec![[0u64; 8]; n_blocks], capacity, len: 0 }
    }
    #[inline]
    fn probe(&self, h: u64) -> (usize, [(usize, u64); K as usize]) {
        let b = (h as usize) & (self.blocks.len() - 1);
        // bit positions from a second, independent mix of the hash
        let mut x = h.rotate_left(29).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        let mut bits = [(0usize, 0u64); K as usize];
        for slot in bits.iter_mut() {
            let bit = (x >> 55) as usize; // 9 bits: 0..512
            *slot = (bit >> 6, 1u64 << (bit & 63));
            x = x.wrapping_mul(0x2545_F491_4F6C_DD1D).rotate_left(17);
        }
        (b, bits)
    }
    #[inline]
    fn contains(&self, h: u64) -> bool {
        let (b, bits) = self.probe(h);
        let blk = &self.blocks[b];
        bits.iter().all(|&(w, m)| blk[w] & m != 0)
    }
    #[inline]
    fn insert(&mut self, h: u64) {
        let (b, bits) = self.probe(h);
        let blk = &mut self.blocks[b];
        for (w, m) in bits {
            blk[w] |= m;
        }
        self.len += 1;
    }
}

/// Streaming detector of repeated 64-bit hashes.
pub struct DupFilter {
    stages: Vec<Stage>,
    candidates: Vec<u64>,
}

impl Default for DupFilter {
    fn default() -> Self {
        Self::with_capacity(FIRST_CAPACITY)
    }
}

impl DupFilter {
    /// `expected`: expected number of inserted hashes (e.g. mapped count from the BAM index).
    pub fn with_capacity(expected: usize) -> Self {
        DupFilter { stages: vec![Stage::new(expected.max(1024))], candidates: Vec::new() }
    }

    /// Inserts `h`; records it as a candidate if it may have been seen before.
    #[inline]
    pub fn insert(&mut self, h: u64) {
        if self.stages.iter().any(|s| s.contains(h)) {
            self.candidates.push(h);
            return;
        }
        let last = self.stages.last_mut().expect("stage");
        if last.len >= last.capacity {
            let cap = last.capacity * 2;
            self.stages.push(Stage::new(cap));
        }
        self.stages.last_mut().expect("stage").insert(h);
    }

    /// Sorted, deduplicated candidate hashes.
    pub fn finish(mut self) -> Vec<u64> {
        self.candidates.sort_unstable();
        self.candidates.dedup();
        self.candidates
    }

    /// Heap bytes in use.
    pub fn heap_bytes(&self) -> usize {
        self.stages.iter().map(|s| s.blocks.len() * 64).sum::<usize>() + self.candidates.capacity() * 8
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_all_repeats_few_false_positives() {
        let mut f = DupFilter::with_capacity(1000);
        let h = |i: u64| iq_types::name_hash(format!("read{i}").as_bytes());
        let n = 200_000u64;
        for i in 0..n {
            f.insert(h(i));
        }
        let dups: Vec<u64> = (0..n).step_by(97).map(h).collect();
        for &d in &dups {
            f.insert(d);
        }
        let c = f.finish();
        for d in dups {
            assert!(c.binary_search(&d).is_ok());
        }
        let fp = c.len() as f64 / n as f64;
        assert!(fp < 0.03, "candidate rate {fp}");
    }
}
