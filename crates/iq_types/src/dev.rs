//! Deviation IDs (design §9) and the runtime quirk switches of the `quirks41` feature.
//!
//! ## Quirk mechanism
//! Each dropped 4.1 behaviour lives behind exactly one call site:
//! ```ignore
//! if iq_types::quirk(DevId::P1) { /* 4.1 behaviour */ } else { /* iqrs behaviour */ }
//! ```
//! - Without the `quirks41` feature, [`quirk`] is a `const fn` returning `false`, so the
//!   4.1 branch compiles away. Release builds never enable the feature.
//! - With `quirks41`, [`quirk`] reads a thread-local [`QuirkSet`]. The fixture runner uses
//!   [`with_quirks`] to run a case with all quirks on (parity: must be exact), with none (=
//!   default-build behaviour), or with exactly one on (to attribute a mismatch to an ID).
//!
//! Default inside a `quirks41` build when nothing was set: **all quirks on**, so plain
//! `cargo test --features quirks41` checks parity.
//!
//! **`IQRS_QUIRKS` (test-only, `quirks41` builds):** overrides that default for every thread
//! (including rayon workers), for deviation-effect attribution of whole runs. Syntax:
//! `all` or `none`, then comma-separated `+D-XX` / `-D-XX`, e.g. `IQRS_QUIRKS=all,-D-M6`
//! (parity minus one quirk) or `IQRS_QUIRKS=none,+D-M6` (default plus one). Read once;
//! an invalid spec panics. [`with_quirks`] still overrides it per scope.

use serde::{Deserialize, Serialize};

macro_rules! dev_ids {
    ($( $(#[$m:meta])* $v:ident = $s:literal, $verdict:ident; )+) => {
        /// A deviation ID of design §9. Serialized as its label, e.g. `"D-P1"`.
        #[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
        #[repr(u8)]
        pub enum DevId { $( $(#[$m])* $v ),+ }

        impl DevId {
            /// All IDs in §9 table order.
            pub const ALL: &'static [DevId] = &[ $( DevId::$v ),+ ];
            /// Label as written in the design doc (`"D-P1"`).
            pub const fn label(self) -> &'static str {
                match self { $( DevId::$v => $s ),+ }
            }
            pub fn from_label(s: &str) -> Option<DevId> {
                match s { $( $s => Some(DevId::$v), )+ _ => None }
            }
            /// §9 verdict.
            pub const fn verdict(self) -> Verdict9 {
                match self { $( DevId::$v => Verdict9::$verdict ),+ }
            }
        }
    };
}

/// Verdict column of the §9 deviation log.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Verdict9 {
    /// Behaviour changed; `quirks41` restores it.
    Drop,
    /// 4.1 semantics kept (scope narrowed only).
    Keep,
    /// Deliberate user decision.
    User,
    /// No behavioural content (dead code not ported, or no output effect).
    None,
}

dev_ids! {
    /// Genes = all gene types by default (`--gene-types isoquant` = 4.1).
    A1 = "D-A1", User;
    /// gffutils `region()` queried with 0-based coordinates (1-bp shift).
    R1 = "D-R1", Drop;
    /// Valley-straddling reads processed in every block they overlap (duplicate resolver keeps the best copy).
    /// Retired 2026-10-29 (user): 4.1 behaviour in both builds, no quirk site; the single-owner rule cost reads
    /// their gene (`validation/deviation_effects.md`: ONT 1.2 pp assigned→no_feature, Revio read-through loci).
    R2 = "D-R2", Keep;
    /// `find_polya` never tests the last window.
    P1 = "D-P1", Drop;
    /// `count_polya_exons` counts non-contiguously.
    P2 = "D-P2", Drop;
    /// `extract_relevant_deletions` anchored at the post-trim first exon.
    P3 = "D-P3", Drop;
    /// polyA detection gated differently for intergenic reads (no effect in scope).
    P4 = "D-P4", Drop;
    /// No common exon → Python `[-1]` index in `categorize_exon_elongation_subtype`.
    A3 = "D-A3", Drop;
    /// `score * 1.5 >= best` with negative best drops every candidate.
    A4 = "D-A4", Drop;
    /// Dead `check_read_ends` assignment-type writes (not ported).
    A5 = "D-A5", None;
    /// `BasicReadAssignment.penalty_score` always 0.
    C1 = "D-C1", Drop;
    /// `cluster_introns` substitute = largest coordinate tuple.
    G1 = "D-G1", Drop;
    /// min-count test skipped for introns with unclustered similar neighbours.
    G2 = "D-G2", Drop;
    /// An ID consumed per rejected FL path.
    M1 = "D-M1", Drop;
    /// `detected_known_isoforms` process-wide across chromosomes.
    M2 = "D-M2", Keep;
    /// FL-path `PolyAInfo` branch dead (not ported).
    M3 = "D-M3", None;
    /// `detect_similar_isoforms` polyA branch dead (not ported).
    M4 = "D-M4", None;
    /// ModelStore keyed by read name.
    M5 = "D-M5", Drop;
    /// Gene joiner stops when one pair remains (`len(scores) > 1`).
    M6 = "D-M6", Drop;
    /// `correct_novel_transcript_ends` captures the last exon's start before editing the first exon, so on a
    /// mono-exon model the end write undoes a start correction (regression in fe84ccdc/0d6e47f2; 0bd9b1e2 read it live).
    E1 = "D-E1", Drop;
    /// Peak tie order = numpy argsort (iqrs emulates the non-SIMD introsort; SIMD CPUs may differ).
    T1 = "D-T1", User;
    /// Zero-peak fallback (ported, unreachable).
    T2 = "D-T2", None;
    /// `set_junction_profiles` redundant per-gene loop.
    X1 = "D-X1", None;
    /// Novel `exon_id` bare `N` first, then `chr.N`.
    O1 = "D-O1", Drop;
    /// Integer read-group tag crashes 4.1.
    O2 = "D-O2", Drop;
    /// Grouped output name override.
    O3 = "D-O3", None;
    /// GTF header version line.
    O4 = "D-O4", None;
}

impl DevId {
    #[inline]
    const fn bit(self) -> u64 {
        1u64 << (self as u8)
    }
}

impl std::fmt::Display for DevId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.label())
    }
}
impl Serialize for DevId {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.serialize_str(self.label())
    }
}
impl<'de> Deserialize<'de> for DevId {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let s = <std::borrow::Cow<'de, str>>::deserialize(d)?;
        DevId::from_label(&s).ok_or_else(|| serde::de::Error::custom(format!("unknown deviation id {s:?}")))
    }
}

/// A set of enabled quirks.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub struct QuirkSet(u64);

impl QuirkSet {
    pub const NONE: QuirkSet = QuirkSet(0);
    /// Every ID (IDs without a code site are harmless no-ops).
    pub const fn all() -> QuirkSet {
        let mut b = 0u64;
        let mut i = 0;
        while i < DevId::ALL.len() {
            b |= DevId::ALL[i].bit();
            i += 1;
        }
        QuirkSet(b)
    }
    pub const fn only(id: DevId) -> QuirkSet {
        QuirkSet(id.bit())
    }
    pub const fn with(self, id: DevId) -> QuirkSet {
        QuirkSet(self.0 | id.bit())
    }
    pub const fn contains(self, id: DevId) -> bool {
        self.0 & id.bit() != 0
    }
    pub const fn without(self, id: DevId) -> QuirkSet {
        QuirkSet(self.0 & !id.bit())
    }
    /// Parses an `IQRS_QUIRKS` spec: `all|none[,+D-XX|,-D-XX]...`.
    pub fn parse(spec: &str) -> Result<QuirkSet, String> {
        let mut toks = spec.split(',').map(str::trim).filter(|t| !t.is_empty());
        let mut set = match toks.next() {
            Some("all") => QuirkSet::all(),
            Some("none") => QuirkSet::NONE,
            other => return Err(format!("quirk spec must start with all|none, got {other:?}")),
        };
        for t in toks {
            let (add, label) = match t.split_at(1) {
                ("+", l) => (true, l),
                ("-", l) => (false, l),
                _ => return Err(format!("quirk spec token {t:?} must start with + or -")),
            };
            let id = DevId::from_label(label).ok_or_else(|| format!("unknown deviation id {label:?}"))?;
            set = if add { set.with(id) } else { set.without(id) };
        }
        Ok(set)
    }
}

/// True if `quirks41` is compiled in (runtime switching is available).
pub const QUIRKS41: bool = cfg!(feature = "quirks41");

/// Process-wide default quirk set: `IQRS_QUIRKS` if set, else all.
#[cfg(feature = "quirks41")]
fn env_default() -> QuirkSet {
    static S: std::sync::OnceLock<QuirkSet> = std::sync::OnceLock::new();
    *S.get_or_init(|| match std::env::var("IQRS_QUIRKS") {
        Ok(v) => QuirkSet::parse(&v).unwrap_or_else(|e| panic!("IQRS_QUIRKS: {e}")),
        Err(_) => QuirkSet::all(),
    })
}

#[cfg(feature = "quirks41")]
thread_local! {
    static QUIRKS: std::cell::Cell<QuirkSet> = std::cell::Cell::new(env_default());
}

/// Is the 4.1 behaviour for `id` enabled on this thread? Always `false` without `quirks41`.
#[cfg(not(feature = "quirks41"))]
#[inline(always)]
pub const fn quirk(_id: DevId) -> bool {
    false
}

/// Is the 4.1 behaviour for `id` enabled on this thread?
#[cfg(feature = "quirks41")]
#[inline]
pub fn quirk(id: DevId) -> bool {
    QUIRKS.with(|q| q.get().contains(id))
}

/// The quirk set active on this thread (`NONE` without `quirks41`).
pub fn current_quirks() -> QuirkSet {
    #[cfg(feature = "quirks41")]
    {
        QUIRKS.with(|q| q.get())
    }
    #[cfg(not(feature = "quirks41"))]
    {
        QuirkSet::NONE
    }
}

/// Runs `f` with `set` active on this thread, restoring the previous set afterwards (also
/// on panic). Without `quirks41` the set is ignored (always none).
pub fn with_quirks<R>(set: QuirkSet, f: impl FnOnce() -> R) -> R {
    #[cfg(feature = "quirks41")]
    {
        struct Restore(QuirkSet);
        impl Drop for Restore {
            fn drop(&mut self) {
                QUIRKS.with(|q| q.set(self.0));
            }
        }
        let _g = Restore(QUIRKS.with(|q| q.replace(set)));
        f()
    }
    #[cfg(not(feature = "quirks41"))]
    {
        let _ = set;
        f()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn labels_roundtrip() {
        assert_eq!(DevId::ALL.len(), 27);
        for &d in DevId::ALL {
            assert_eq!(DevId::from_label(d.label()), Some(d));
            let j = serde_json::to_string(&d).unwrap();
            assert_eq!(serde_json::from_str::<DevId>(&j).unwrap(), d);
        }
    }

    #[test]
    fn parse_spec() {
        assert_eq!(QuirkSet::parse("all").unwrap(), QuirkSet::all());
        assert_eq!(QuirkSet::parse("none, +D-M6").unwrap(), QuirkSet::only(DevId::M6));
        let s = QuirkSet::parse("all,-D-M6,-D-R2").unwrap();
        assert!(!s.contains(DevId::M6) && !s.contains(DevId::R2) && s.contains(DevId::P1));
        assert!(QuirkSet::parse("D-M6").is_err());
        assert!(QuirkSet::parse("all,D-M6").is_err());
        assert!(QuirkSet::parse("all,-D-ZZ").is_err());
    }

    #[test]
    fn quirk_switching() {
        let inner = with_quirks(QuirkSet::only(DevId::P1), || (quirk(DevId::P1), quirk(DevId::P2)));
        if QUIRKS41 {
            assert_eq!(inner, (true, false));
            assert!(quirk(DevId::G1)); // default: all on
        } else {
            assert_eq!(inner, (false, false));
            assert!(!quirk(DevId::G1));
        }
    }
}
