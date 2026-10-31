//! Id-string ↔ numeric-key maps used when converting fixture DTOs to core types.
//!
//! Conventions (fixture mode only; production uses the `AnnotationIndex`):
//! - annotated genes → `GeneKey::Ref(genes.intern(id))`;
//! - ids starting with `novel_gene_` → `GeneKey::Novel(novel_genes.intern(id))`;
//! - reference transcripts → `IsoRef::Ref(txs.intern(id))`;
//! - models → `IsoRef::Model(models.intern(id))`: `ModelIdx` = first-seen order of the 4.1
//!   transcript id within the block, which equals creation order when snapshots are
//!   converted in call order. Call [`IdMaps::reset_block`] between blocks.
//!   A known model's id is also its reference transcript id.
//! - read groups: the `read_group` list of a record → [`crate::GroupPool`] id (empty list
//!   → 0; several names are joined with `\x1f`).

use crate::assign::{GeneKey, IsoRef};
use crate::intern::{GroupPool, Interner};

/// Prefix of 4.1 novel gene ids (`novel_gene_<chr>_<N>`).
pub const NOVEL_GENE_PREFIX: &str = "novel_gene_";

/// Which id space transcript strings belong to in a given payload.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum TxSpace {
    /// Reference transcripts (collect side, reference `GeneInfo`).
    Ref,
    /// Block models (`from_models` indices, store snapshots).
    Model,
}

#[derive(Clone, Debug, Default)]
pub struct IdMaps {
    pub chroms: Interner,
    pub genes: Interner,
    pub novel_genes: Interner,
    pub txs: Interner,
    pub models: Interner,
    pub groups: GroupPool,
}

const GROUP_SEP: char = '\x1f';

impl IdMaps {
    pub fn new() -> Self {
        Self::default()
    }
    /// Forget block-local model ids (call when moving to another block).
    pub fn reset_block(&mut self) {
        self.models = Interner::new();
    }

    pub fn gene_key(&mut self, id: &str) -> GeneKey {
        if id.starts_with(NOVEL_GENE_PREFIX) {
            GeneKey::Novel(self.novel_genes.intern(id))
        } else {
            GeneKey::Ref(self.genes.intern(id))
        }
    }
    pub fn gene_name(&self, k: GeneKey) -> &str {
        match k {
            GeneKey::Ref(i) => self.genes.resolve(i),
            GeneKey::Novel(i) => self.novel_genes.resolve(i),
        }
    }

    pub fn iso_ref(&mut self, id: &str, space: TxSpace) -> IsoRef {
        match space {
            TxSpace::Ref => IsoRef::Ref(self.txs.intern(id)),
            TxSpace::Model => IsoRef::Model(self.models.intern(id)),
        }
    }
    pub fn iso_name(&self, r: IsoRef) -> &str {
        match r {
            IsoRef::Ref(i) => self.txs.resolve(i),
            IsoRef::Model(i) => self.models.resolve(i),
        }
    }

    pub fn group_id(&mut self, names: &[String]) -> u32 {
        match names {
            [] => GroupPool::NONE,
            [one] => self.groups.intern(one),
            many => self.groups.intern(&many.join(&GROUP_SEP.to_string())),
        }
    }
    pub fn group_names(&self, id: u32) -> Vec<String> {
        match self.groups.resolve(id) {
            None => Vec::new(),
            Some(s) => s.split(GROUP_SEP).map(str::to_owned).collect(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spaces() {
        let mut m = IdMaps::new();
        let g = m.gene_key("ENSDARG1");
        let n = m.gene_key("novel_gene_8_3");
        assert!(matches!(g, GeneKey::Ref(0)) && matches!(n, GeneKey::Novel(0)));
        assert_eq!(m.gene_name(n), "novel_gene_8_3");
        let r = m.iso_ref("ENSDART1", TxSpace::Ref);
        let k = m.iso_ref("ENSDART1", TxSpace::Model);
        assert_eq!(m.iso_name(r), m.iso_name(k));
        let gid = m.group_id(&["A".into(), "x".into()]);
        assert_eq!(m.group_names(gid), vec!["A".to_string(), "x".to_string()]);
        assert_eq!(m.group_id(&[]), 0);
        assert!(m.group_names(0).is_empty());
    }
}
