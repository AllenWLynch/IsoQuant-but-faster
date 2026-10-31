//! Python container semantics the intron graph depends on.
//!
//! - [`PyDict`]: insertion-ordered map with O(1) `del` (tombstones). Deleting and later
//!   re-inserting a key moves it to the end, and assigning to an existing key keeps its
//!   position, exactly like a Python `dict`. [`PyDict::entry_default`] is a `defaultdict`
//!   read (`d[k]`), which inserts a default value when the key is missing; 4.1 output order
//!   depends on those implicit insertions, so every such read is mirrored explicitly.
//! - [`VSet`]: a small sorted vertex set (Python `set`).

use iq_types::smallvec::SmallVec;
use iq_types::Vertex;
use std::collections::HashMap;
use std::hash::Hash;

/// Insertion-ordered map with Python `dict` semantics (see module docs).
#[derive(Clone, Debug)]
pub(crate) struct PyDict<K, V> {
    slots: Vec<Option<(K, V)>>,
    pos: HashMap<K, usize>,
}

impl<K, V> Default for PyDict<K, V> {
    fn default() -> Self {
        PyDict { slots: Vec::new(), pos: HashMap::new() }
    }
}

impl<K: Copy + Eq + Hash, V> PyDict<K, V> {
    #[inline]
    pub(crate) fn len(&self) -> usize {
        self.pos.len()
    }
    #[inline]
    pub(crate) fn is_empty(&self) -> bool {
        self.pos.is_empty()
    }
    #[inline]
    pub(crate) fn contains(&self, k: &K) -> bool {
        self.pos.contains_key(k)
    }
    #[inline]
    pub(crate) fn get(&self, k: &K) -> Option<&V> {
        self.pos.get(k).map(|&i| &self.slots[i].as_ref().expect("live slot").1)
    }
    #[inline]
    pub(crate) fn get_mut(&mut self, k: &K) -> Option<&mut V> {
        match self.pos.get(k) {
            Some(&i) => Some(&mut self.slots[i].as_mut().expect("live slot").1),
            None => None,
        }
    }
    /// `d[k] = v`: overwrites in place, or appends.
    pub(crate) fn insert(&mut self, k: K, v: V) {
        match self.pos.get(&k) {
            Some(&i) => self.slots[i].as_mut().expect("live slot").1 = v,
            None => {
                self.pos.insert(k, self.slots.len());
                self.slots.push(Some((k, v)));
            }
        }
    }
    /// `del d[k]` (no-op when absent).
    pub(crate) fn remove(&mut self, k: &K) -> Option<V> {
        let i = self.pos.remove(k)?;
        let v = self.slots[i].take().map(|(_, v)| v);
        if self.slots.len() > 32 && self.slots.len() > 2 * self.pos.len() {
            self.compact();
        }
        v
    }
    fn compact(&mut self) {
        self.slots.retain(Option::is_some);
        for (i, s) in self.slots.iter().enumerate() {
            let (k, _) = s.as_ref().expect("live slot");
            *self.pos.get_mut(k).expect("indexed key") = i;
        }
    }
    /// Entries in dict order.
    pub(crate) fn iter(&self) -> impl Iterator<Item = (&K, &V)> + '_ {
        self.slots.iter().filter_map(|s| s.as_ref().map(|(k, v)| (k, v)))
    }
    /// Keys in dict order.
    pub(crate) fn keys(&self) -> impl Iterator<Item = K> + '_ {
        self.iter().map(|(k, _)| *k)
    }
    /// `sorted(d.keys())`.
    pub(crate) fn keys_sorted(&self) -> Vec<K>
    where
        K: Ord,
    {
        let mut v: Vec<K> = self.keys().collect();
        v.sort_unstable();
        v
    }
}

impl<K: Copy + Eq + Hash, V: Default> PyDict<K, V> {
    /// `defaultdict.__getitem__`: inserts `V::default()` at the end when missing.
    pub(crate) fn entry_default(&mut self, k: K) -> &mut V {
        let i = match self.pos.get(&k) {
            Some(&i) => i,
            None => {
                let i = self.slots.len();
                self.pos.insert(k, i);
                self.slots.push(Some((k, V::default())));
                i
            }
        };
        &mut self.slots[i].as_mut().expect("live slot").1
    }
}

/// A small sorted vertex set (Python `set`). Iteration is ascending; 4.1 iterates in hash
/// order, which only changes the order in which `defaultdict` reads insert adjacency keys,
/// never the graph content (see `graph.rs` module docs).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct VSet(SmallVec<[Vertex; 4]>);

impl VSet {
    #[inline]
    pub(crate) fn len(&self) -> usize {
        self.0.len()
    }
    #[inline]
    pub(crate) fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
    pub(crate) fn insert(&mut self, v: Vertex) -> bool {
        match self.0.binary_search(&v) {
            Ok(_) => false,
            Err(i) => {
                self.0.insert(i, v);
                true
            }
        }
    }
    /// `set.discard` (4.1 uses `set.remove`, which can only fail on an asymmetric graph).
    pub(crate) fn remove(&mut self, v: &Vertex) -> bool {
        match self.0.binary_search(v) {
            Ok(i) => {
                self.0.remove(i);
                true
            }
            Err(_) => false,
        }
    }
    /// `set.update`.
    pub(crate) fn extend(&mut self, other: &VSet) {
        for v in other.iter() {
            self.insert(*v);
        }
    }
    /// Ascending (= Python `sorted`).
    #[inline]
    pub(crate) fn iter(&self) -> std::slice::Iter<'_, Vertex> {
        self.0.iter()
    }
    #[inline]
    pub(crate) fn to_vec(&self) -> Vec<Vertex> {
        self.0.to_vec()
    }
}

/// `common.find_closest` (`common.py:119`): first element (in list order) with the strictly
/// smallest distance; `None` for an empty list.
#[inline]
pub(crate) fn find_closest(value: i64, list: &[i64]) -> Option<(i64, i64)> {
    let mut best: Option<(i64, i64)> = None;
    for &v in list {
        let d = (v - value).abs();
        if best.is_none_or(|(_, bd)| d < bd) {
            best = Some((v, d));
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pydict_order() {
        let mut d: PyDict<i32, u32> = PyDict::default();
        d.insert(3, 1);
        d.insert(1, 1);
        *d.entry_default(2) += 5;
        d.insert(3, 9); // keeps position
        assert_eq!(d.keys().collect::<Vec<_>>(), vec![3, 1, 2]);
        d.remove(&3);
        *d.entry_default(3) += 1; // re-inserted at the end
        assert_eq!(d.iter().map(|(k, v)| (*k, *v)).collect::<Vec<_>>(), vec![(1, 1), (2, 5), (3, 1)]);
        for k in 10..100 {
            d.insert(k, 0);
        }
        for k in 10..95 {
            d.remove(&k);
        }
        assert_eq!(d.keys().collect::<Vec<_>>(), vec![1, 2, 3, 95, 96, 97, 98, 99]);
        assert_eq!(d.get(&96), Some(&0));
        assert_eq!(d.len(), 8);
    }

    #[test]
    fn closest_ties() {
        assert_eq!(find_closest(10, &[12, 8]), Some((12, 2)));
        assert_eq!(find_closest(10, &[8, 12]), Some((8, 2)));
        assert_eq!(find_closest(10, &[]), None);
    }
}
