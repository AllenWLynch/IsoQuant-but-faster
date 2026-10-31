//! String interning (replaces 4.1 `utils/string_pools.py`).
//!
//! - [`Interner`]: append-only `str ↔ u32`, ids dense from 0 in first-intern order.
//! - [`GroupPool`]: read-group strings; id 0 is reserved for "no grouping" (no
//!   `--read-group-tag`), real groups start at 1.

use rustc_hash::FxHashMap;
use std::sync::Arc;

/// Append-only string interner. Ids are dense and assigned in first-intern order, so they
/// are deterministic given a deterministic intern order.
#[derive(Clone, Debug, Default)]
pub struct Interner {
    map: FxHashMap<Arc<str>, u32>,
    strs: Vec<Arc<str>>,
}

impl Interner {
    pub fn new() -> Self {
        Self::default()
    }
    /// Id of `s`, interning it if new.
    pub fn intern(&mut self, s: &str) -> u32 {
        if let Some(&id) = self.map.get(s) {
            return id;
        }
        let id = self.strs.len() as u32;
        let a: Arc<str> = Arc::from(s);
        self.strs.push(a.clone());
        self.map.insert(a, id);
        id
    }
    /// Id of `s` if already interned.
    pub fn get(&self, s: &str) -> Option<u32> {
        self.map.get(s).copied()
    }
    /// String of `id`. Panics on an unknown id.
    pub fn resolve(&self, id: u32) -> &str {
        &self.strs[id as usize]
    }
    pub fn try_resolve(&self, id: u32) -> Option<&str> {
        self.strs.get(id as usize).map(|s| &**s)
    }
    pub fn len(&self) -> usize {
        self.strs.len()
    }
    pub fn is_empty(&self) -> bool {
        self.strs.is_empty()
    }
    /// All strings in id order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = &str> {
        self.strs.iter().map(|s| &**s)
    }
}

/// Read-group pool. Id 0 = no grouping; groups are interned from id 1.
#[derive(Clone, Debug, Default)]
pub struct GroupPool {
    inner: Interner,
}

impl GroupPool {
    /// Id meaning "no read grouping configured".
    pub const NONE: u32 = 0;

    pub fn new() -> Self {
        Self::default()
    }
    /// Id of group string `g` (≥ 1).
    pub fn intern(&mut self, g: &str) -> u32 {
        self.inner.intern(g) + 1
    }
    pub fn get(&self, g: &str) -> Option<u32> {
        self.inner.get(g).map(|i| i + 1)
    }
    /// Group string of `id`; `None` for [`GroupPool::NONE`].
    pub fn resolve(&self, id: u32) -> Option<&str> {
        if id == Self::NONE {
            None
        } else {
            Some(self.inner.resolve(id - 1))
        }
    }
    /// Number of real groups.
    pub fn len(&self) -> usize {
        self.inner.len()
    }
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interner() {
        let mut i = Interner::new();
        assert_eq!(i.intern("b"), 0);
        assert_eq!(i.intern("a"), 1);
        assert_eq!(i.intern("b"), 0);
        assert_eq!(i.resolve(1), "a");
        assert_eq!(i.get("c"), None);
        assert_eq!(i.iter().collect::<Vec<_>>(), vec!["b", "a"]);
    }

    #[test]
    fn groups() {
        let mut g = GroupPool::new();
        assert_eq!(g.intern("A"), 1);
        assert_eq!(g.intern("NA"), 2);
        assert_eq!(g.resolve(0), None);
        assert_eq!(g.resolve(2), Some("NA"));
    }
}
