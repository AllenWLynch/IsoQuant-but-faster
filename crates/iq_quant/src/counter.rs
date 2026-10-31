//! `AssignedFeatureCounter` for transcript models, strategy `unique_only`
//! (`quantification/long_read_counter.py:327-394`, weights `:151-183`), ungrouped and
//! grouped in one pass.
//!
//! [`FeatureCounter`] is the per-block partial. Features are dense `u32`s, which are
//! `ModelIdx` (block-local creation order) in production, and groups are `GroupPool` ids.
//! Each read costs O(#features) with no allocation in the common case and no string work.
//! The per-chromosome merge keyed by final feature id lives in [`crate::output::ChrCounts`].

use iq_types::smallvec::SmallVec;
use iq_types::{AssignmentType, FinalAssign, ModelIdx};

/// Group list of one feature: `(group, count)` in first-touch order (4.1 `IncrementalDict`,
/// including weight-0 touches by ambiguous reads).
pub type Groups = SmallVec<[(u32, f64); 4]>;

/// One feature's state.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FeatState {
    /// Ungrouped count (4.1 group 0 of the ungrouped counter), summed in read order.
    pub total: f64,
    /// Per-group counts (the grouped counter), first-touch order.
    pub groups: Groups,
    /// Ever touched (has an entry in 4.1 `feature_counter`).
    pub touched: bool,
    /// In 4.1 `all_features` (some read added a positive weight).
    pub in_all: bool,
}

/// `unique_only` read weight (`ReadWeightCounter`, `long_read_counter.py:151-183`) and the
/// per-type branch of `add_read_info` (`:356-394`).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum ReadFate {
    /// `__no_feature += 1` (unassigned type, no matches, or no feature).
    NoFeature,
    /// `ambiguous`: `__ambiguous += 1`; every feature is touched with `weight`.
    Ambiguous { weight_one: bool },
    /// Inconsistent types: weight 0 under `unique_only` (nothing touched); under
    /// `unique_splicing_consistent` a single-feature `inconsistent_non_intronic` read
    /// counts 1.0 (`weight_one`).
    Inconsistent { weight_one: bool },
    /// `unique` / `umd`: +1.0 to the first feature.
    Unique,
    /// Any other type with features (`discarded`/`suspended` with features): counted for
    /// TPM only, nothing touched.
    Other,
}

/// Counting strategy (`CountingStrategy`). Transcript models use `UniqueOnly` (the 4.1
/// default `--transcript_quantification`); `UniqueSplicingConsistent` is the gene default,
/// kept only so the gene-model counter fixtures can be checked.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum Strategy {
    #[default]
    UniqueOnly,
    UniqueSplicingConsistent,
}

/// Branch of `add_read_info` for assignment type `atype` and `n_features` distinct features.
/// `has_matches` = `isoform_matches` non-empty.
pub fn read_fate(strategy: Strategy, atype: AssignmentType, has_matches: bool, n_features: usize) -> ReadFate {
    if atype.is_unassigned() || !has_matches || n_features == 0 {
        ReadFate::NoFeature
    } else if atype == AssignmentType::Ambiguous {
        // process_ambiguous: 1 feature -> 1.0; several -> 0.0 (no use_ambiguous)
        ReadFate::Ambiguous { weight_one: n_features == 1 }
    } else if atype.is_inconsistent() {
        // process_inconsistent: no strategy here has use_inconsistent; only
        // unique_splicing_consistent has use_inconsistent_minor
        let minor = strategy == Strategy::UniqueSplicingConsistent
            && atype == AssignmentType::InconsistentNonIntronic
            && n_features == 1;
        ReadFate::Inconsistent { weight_one: minor }
    } else if atype.is_unique() {
        ReadFate::Unique
    } else {
        ReadFate::Other
    }
}

/// Per-block counter partial over dense feature ids.
#[derive(Clone, Debug, Default)]
pub struct FeatureCounter {
    /// Indexed by feature id; grown on demand.
    pub feats: Vec<FeatState>,
    /// `__ambiguous`.
    pub ambiguous: u64,
    /// `__no_feature` (`not_assigned_reads`).
    pub no_feature: u64,
}

impl FeatureCounter {
    pub fn new() -> Self {
        Self::default()
    }

    #[inline]
    fn slot(&mut self, f: u32) -> &mut FeatState {
        let i = f as usize;
        if i >= self.feats.len() {
            self.feats.resize_with(i + 1, FeatState::default);
        }
        &mut self.feats[i]
    }

    #[inline]
    fn inc(&mut self, f: u32, group: u32, w: f64, positive: bool) {
        let s = self.slot(f);
        s.touched = true;
        s.total += w;
        match s.groups.iter_mut().find(|(g, _)| *g == group) {
            Some(e) => e.1 += w,
            None => s.groups.push((group, w)),
        }
        if positive {
            s.in_all = true;
        }
    }

    /// `add_read_info` for one read: `features` are the read's **distinct** feature ids
    /// (match order). For `unique`, 4.1 takes `list(set)[0]`; with one feature (the only
    /// case seen) that is `features[0]`. Returns the fate.
    pub fn add_features(
        &mut self,
        strategy: Strategy,
        atype: AssignmentType,
        has_matches: bool,
        features: &[u32],
        group: u32,
    ) -> ReadFate {
        let fate = read_fate(strategy, atype, has_matches, features.len());
        match fate {
            ReadFate::NoFeature => self.no_feature += 1,
            ReadFate::Ambiguous { weight_one } => {
                let w = if weight_one { 1.0 } else { 0.0 };
                for &f in features {
                    self.inc(f, group, w, weight_one);
                }
                self.ambiguous += 1;
            }
            ReadFate::Unique => self.inc(features[0], group, 1.0, true),
            ReadFate::Inconsistent { weight_one: true } => {
                for &f in features {
                    self.inc(f, group, 1.0, true);
                }
            }
            ReadFate::Inconsistent { weight_one: false } | ReadFate::Other => {}
        }
        fate
    }

    /// Number of feature slots.
    pub fn len(&self) -> usize {
        self.feats.len()
    }
    pub fn is_empty(&self) -> bool {
        self.feats.is_empty()
    }
}

/// The transcript-model counter of one block (ungrouped + grouped), keyed by `ModelIdx`.
#[derive(Clone, Debug, Default)]
pub struct TxCounter {
    pub inner: FeatureCounter,
}

impl TxCounter {
    pub fn new() -> Self {
        Self::default()
    }

    /// `add_read_info` for one final assignment of a read in read group `group`
    /// (`GroupPool` id; `GroupPool::NONE` without grouping). Must be called in read order
    /// (f64 sums in 4.1 order).
    pub fn add(&mut self, fa: &FinalAssign, group: u32) -> ReadFate {
        // distinct models in match order (Python set; order only matters for `unique`
        // with >1 model, never observed, where 4.1 order is hash order)
        let mut uniq: SmallVec<[u32; 4]> = SmallVec::new();
        for &m in &fa.tx {
            if !uniq.contains(&m) {
                uniq.push(m);
            }
        }
        // `has_matches`: U19 always produces matches for reads that name a model, and a
        // read without models is `__no_feature` either way.
        self.inner.add_features(Strategy::UniqueOnly, fa.atype, !uniq.is_empty(), &uniq, group)
    }

    /// Ungrouped count of a model.
    pub fn count(&self, m: ModelIdx) -> f64 {
        self.inner.feats.get(m as usize).map_or(0.0, |s| s.total)
    }

    /// State of a model, if it was ever touched.
    pub fn state(&self, m: ModelIdx) -> Option<&FeatState> {
        self.inner.feats.get(m as usize).filter(|s| s.touched)
    }

    pub fn ambiguous(&self) -> u64 {
        self.inner.ambiguous
    }
    pub fn no_feature(&self) -> u64 {
        self.inner.no_feature
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use AssignmentType as A;

    fn fa(atype: A, tx: &[u32]) -> FinalAssign {
        FinalAssign { atype, gene_atype: atype, tx: tx.iter().copied().collect() }
    }

    #[test]
    fn unique_only_weights() {
        let mut c = TxCounter::new();
        assert_eq!(c.add(&fa(A::Unique, &[3]), 1), ReadFate::Unique);
        c.add(&fa(A::UniqueMinorDifference, &[3]), 2);
        c.add(&fa(A::Ambiguous, &[3, 4]), 3); // weight 0, touches both, group 3 appended
        c.add(&fa(A::Ambiguous, &[4, 4]), 1); // dedup -> one feature, weight 1
        c.add(&fa(A::Inconsistent, &[5]), 1);
        c.add(&fa(A::Noninformative, &[]), 1);
        c.add(&fa(A::Discarded, &[]), 1);
        c.add(&fa(A::Intergenic, &[3]), 1);
        assert_eq!(c.count(3), 2.0);
        assert_eq!(c.count(4), 1.0);
        assert_eq!(c.state(3).unwrap().groups, Groups::from_vec(vec![(1, 1.0), (2, 1.0), (3, 0.0)]));
        let s4 = c.state(4).unwrap();
        assert_eq!(s4.groups, Groups::from_vec(vec![(3u32, 0.0), (1, 1.0)]));
        assert!(s4.in_all);
        assert!(c.state(5).is_none());
        assert_eq!(c.ambiguous(), 2);
        assert_eq!(c.no_feature(), 3);
    }

    #[test]
    fn ambiguous_zero_not_in_all() {
        let mut c = TxCounter::new();
        c.add(&fa(A::Ambiguous, &[0, 1]), 0);
        assert!(c.state(0).unwrap().touched && !c.state(0).unwrap().in_all);
    }
}
