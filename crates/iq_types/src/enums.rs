//! IsoQuant 4.1 enums with their **exact numeric codes**, so fixtures can carry raw integers
//! (design §2.1). Every enum here:
//! - serializes/deserializes as its integer code (unknown codes are a hard error);
//! - has `code()`, `from_code()`, `py_name()` (the Python member name) and `ALL`.
//!
//! Abbreviation in citations: `IA` = `assignment/isoform_assignment.py`.

use serde::{Deserialize, Deserializer, Serialize, Serializer};

/// Declares an integer-coded enum mirroring a Python `Enum`.
macro_rules! py_enum {
    (
        $(#[$m:meta])*
        pub enum $name:ident : $repr:ty {
            $( $(#[$vm:meta])* $var:ident = $code:literal => $py:literal ),+ $(,)?
        }
    ) => {
        $(#[$m])*
        #[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
        #[repr($repr)]
        pub enum $name { $( $(#[$vm])* $var = $code ),+ }

        impl $name {
            /// Every member, in Python declaration order.
            pub const ALL: &'static [$name] = &[ $( $name::$var ),+ ];
            /// The IsoQuant numeric code.
            #[inline]
            pub const fn code(self) -> $repr { self as $repr }
            /// Inverse of [`Self::code`]; `None` for codes IsoQuant does not define.
            pub const fn from_code(c: $repr) -> Option<$name> {
                match c { $( $code => Some($name::$var), )+ _ => None }
            }
            /// The Python member name (for readable diffs and logs).
            pub const fn py_name(self) -> &'static str {
                match self { $( $name::$var => $py ),+ }
            }
            /// Look up by Python member name.
            pub fn from_py_name(s: &str) -> Option<$name> {
                match s { $( $py => Some($name::$var), )+ _ => None }
            }
        }

        impl Serialize for $name {
            fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
                (self.code() as i64).serialize(s)
            }
        }
        impl<'de> Deserialize<'de> for $name {
            fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
                let v = i64::deserialize(d)?;
                <$repr>::try_from(v)
                    .ok()
                    .and_then($name::from_code)
                    .ok_or_else(|| serde::de::Error::custom(
                        format!(concat!("invalid ", stringify!($name), " code {}"), v)))
            }
        }
        impl std::fmt::Display for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                f.write_str(self.py_name())
            }
        }
    };
}

py_enum! {
    /// `ReadAssignmentType`, IA:39-53.
    pub enum AssignmentType: u8 {
        Noninformative = 0 => "noninformative",
        Unique = 1 => "unique",
        UniqueMinorDifference = 2 => "unique_minor_difference",
        Ambiguous = 10 => "ambiguous",
        Intergenic = 20 => "intergenic",
        Inconsistent = 30 => "inconsistent",
        InconsistentNonIntronic = 31 => "inconsistent_non_intronic",
        InconsistentAmbiguous = 32 => "inconsistent_ambiguous",
        /// A single overlapping gene, no isoform resemblance (mirrors `inconsistent`).
        InconsistentGenic = 33 => "inconsistent_genic",
        /// Several overlapping genes (mirrors `inconsistent_ambiguous`).
        InconsistentMultigenic = 34 => "inconsistent_multigenic",
        Discarded = 40 => "discarded",
        Suspended = 255 => "suspended",
    }
}

impl AssignmentType {
    /// IA:55-60.
    pub const fn is_inconsistent(self) -> bool {
        matches!(
            self,
            Self::Inconsistent
                | Self::InconsistentAmbiguous
                | Self::InconsistentNonIntronic
                | Self::InconsistentGenic
                | Self::InconsistentMultigenic
        )
    }
    /// IA:62-65.
    pub const fn is_consistent(self) -> bool {
        matches!(self, Self::Unique | Self::UniqueMinorDifference | Self::Ambiguous)
    }
    /// IA:67-69.
    pub const fn is_unassigned(self) -> bool {
        matches!(self, Self::Noninformative | Self::Intergenic)
    }
    /// IA:71-73.
    pub const fn is_unique(self) -> bool {
        matches!(self, Self::UniqueMinorDifference | Self::Unique)
    }
    /// IA:75-78.
    pub const fn is_ambiguous(self) -> bool {
        matches!(self, Self::Ambiguous | Self::InconsistentAmbiguous | Self::InconsistentMultigenic)
    }
}

py_enum! {
    /// `MatchClassification` (SQANTI-like), IA:83-95.
    pub enum MatchClass: u8 {
        Undefined = 0 => "undefined",
        FullSpliceMatch = 10 => "full_splice_match",
        IncompleteSpliceMatch = 11 => "incomplete_splice_match",
        MonoExonMatch = 12 => "mono_exon_match",
        NovelInCatalog = 20 => "novel_in_catalog",
        NovelNotInCatalog = 21 => "novel_not_in_catalog",
        Genic = 30 => "genic",
        GenicIntron = 31 => "genic_intron",
        Intergenic = 32 => "intergenic",
        Antisense = 40 => "antisense",
        Fusion = 50 => "fusion",
    }
}

py_enum! {
    /// `MatchEventSubtype`, IA:128-210. Python orders members by value (`__lt__`, IA:212);
    /// [`EventType::code`] gives that order.
    pub enum EventType: u16 {
        None = 0 => "none",
        Undefined = 1 => "undefined",
        // non-contradictory
        MonoExonic = 10 => "mono_exonic",
        Antisense = 9 => "antisense",
        Fsm = 11 => "fsm",
        IsmLeft = 15 => "ism_left",
        IsmRight = 13 => "ism_right",
        IsmInternal = 14 => "ism_internal",
        MonoExonMatch = 19 => "mono_exon_match",
        // alignment artifacts
        IntronShift = 21 => "intron_shift",
        ExonMisalignment = 22 => "exon_misalignment",
        FakeTerminalExonLeft = 24 => "fake_terminal_exon_left",
        FakeTerminalExonRight = 26 => "fake_terminal_exon_right",
        TerminalExonMisalignmentRight = 28 => "terminal_exon_misalignment_right",
        TerminalExonMisalignmentLeft = 29 => "terminal_exon_misalignment_left",
        // minor alternations
        ExonElongationLeft = 25 => "exon_elongation_left",
        ExonElongationRight = 23 => "exon_elongation_right",
        // intron retentions
        IntronRetention = 31 => "intron_retention",
        UnsplicedIntronRetention = 32 => "unspliced_intron_retention",
        IncompleteIntronRetentionLeft = 38 => "incomplete_intron_retention_left",
        IncompleteIntronRetentionRight = 39 => "incomplete_intron_retention_right",
        FakeMicroIntronRetention = 34 => "fake_micro_intron_retention",
        // alternative donor/acceptor sites
        AltLeftSiteKnown = 101 => "alt_left_site_known",
        AltRightSiteKnown = 102 => "alt_right_site_known",
        AltLeftSiteNovel = 103 => "alt_left_site_novel",
        AltRightSiteNovel = 104 => "alt_right_site_novel",
        // additional introns in the middle
        ExtraIntronNovel = 1012 => "extra_intron_novel",
        ExtraIntronKnown = 1011 => "extra_intron_known",
        // extra introns on the sides
        ExtraIntronFlankingLeft = 1013 => "extra_intron_flanking_left",
        ExtraIntronFlankingRight = 1015 => "extra_intron_flanking_right",
        // significant exon elongation
        MajorExonElongationLeft = 1005 => "major_exon_elongation_left",
        MajorExonElongationRight = 1003 => "major_exon_elongation_right",
        // other intron modifications
        IntronMigration = 114 => "intron_migration",
        IntronAlternationNovel = 115 => "intron_alternation_novel",
        IntronAlternationKnown = 116 => "intron_alternation_known",
        // mutually exclusive
        MutuallyExclusiveExonsNovel = 121 => "mutually_exclusive_exons_novel",
        MutuallyExclusiveExonsKnown = 122 => "mutually_exclusive_exons_known",
        // exon skipping
        ExonSkippingKnown = 133 => "exon_skipping_known",
        ExonSkippingNovel = 134 => "exon_skipping_novel",
        ExonMergeKnown = 135 => "exon_merge_known",
        ExonMergeNovel = 136 => "exon_merge_novel",
        // exon gain
        ExonGainKnown = 145 => "exon_gain_known",
        ExonGainNovel = 146 => "exon_gain_novel",
        ExonDetachKnown = 147 => "exon_detach_known",
        ExonDetachNovel = 148 => "exon_detach_novel",
        // terminal exon shift
        TerminalExonShiftKnown = 151 => "terminal_exon_shift_known",
        TerminalExonShiftNovel = 152 => "terminal_exon_shift_novel",
        // other
        AlternativeStructureNovel = 181 => "alternative_structure_novel",
        AlternativeStructureKnown = 182 => "alternative_structure_known",
        // TTS and TSS
        AlternativePolyaSiteRight = 200 => "alternative_polya_site_right",
        AlternativePolyaSiteLeft = 201 => "alternative_polya_site_left",
        InternalPolyaLeft = 202 => "internal_polya_left",
        InternalPolyaRight = 203 => "internal_polya_right",
        AlternativeTssLeft = 204 => "alternative_tss_left",
        AlternativeTssRight = 205 => "alternative_tss_right",
        CorrectPolyaSiteLeft = 222 => "correct_polya_site_left",
        CorrectPolyaSiteRight = 223 => "correct_polya_site_right",
        AlignedPolyaTail = 300 => "aligned_polya_tail",
        TerminalSiteMatchLeft = 250 => "terminal_site_match_left",
        TerminalSiteMatchRight = 251 => "terminal_site_match_right",
        TerminalSiteMatchLeftPrecise = 252 => "terminal_site_match_left_precise",
        TerminalSiteMatchRightPrecise = 253 => "terminal_site_match_right_precise",
    }
}

impl EventType {
    /// `event_subtype_cost`, IA:271-346. `None` for `antisense`, which has no entry in the
    /// Python dict (a lookup would raise `KeyError`).
    pub const fn cost(self) -> Option<f64> {
        use EventType::*;
        Some(match self {
            None | Undefined | MonoExonic | Fsm | IsmLeft | IsmRight | IsmInternal | MonoExonMatch => 0.0,
            IntronShift | ExonMisalignment | TerminalExonMisalignmentRight | TerminalExonMisalignmentLeft => 0.1,
            FakeTerminalExonLeft | FakeTerminalExonRight => 0.2,
            ExonElongationLeft | ExonElongationRight => 0.1,
            IntronRetention | UnsplicedIntronRetention => 0.6,
            IncompleteIntronRetentionLeft | IncompleteIntronRetentionRight => 0.7,
            FakeMicroIntronRetention => 0.1,
            AltLeftSiteKnown | AltRightSiteKnown | AltLeftSiteNovel | AltRightSiteNovel => 1.0,
            ExtraIntronNovel | ExtraIntronKnown => 1.0,
            ExtraIntronFlankingLeft | ExtraIntronFlankingRight => 1.0,
            MajorExonElongationLeft | MajorExonElongationRight => 0.6,
            IntronMigration | IntronAlternationNovel | IntronAlternationKnown => 1.0,
            MutuallyExclusiveExonsNovel | MutuallyExclusiveExonsKnown => 0.8,
            ExonSkippingKnown | ExonSkippingNovel => 1.0,
            ExonMergeKnown | ExonMergeNovel => 0.75,
            ExonGainKnown | ExonGainNovel => 1.0,
            ExonDetachKnown | ExonDetachNovel => 0.75,
            TerminalExonShiftKnown | TerminalExonShiftNovel => 0.5,
            AlternativeStructureNovel | AlternativeStructureKnown => 1.0,
            AlternativePolyaSiteLeft | AlternativePolyaSiteRight => 0.6,
            InternalPolyaLeft | InternalPolyaRight => 0.5,
            AlternativeTssLeft | AlternativeTssRight => 0.6,
            CorrectPolyaSiteLeft | CorrectPolyaSiteRight | AlignedPolyaTail => 0.0,
            TerminalSiteMatchLeft
            | TerminalSiteMatchRight
            | TerminalSiteMatchLeftPrecise
            | TerminalSiteMatchRightPrecise => 0.0,
            Antisense => return Option::None,
        })
    }

    /// `MatchEventSubtype.is_alignment_artifact`, IA:215-222.
    pub const fn is_alignment_artifact(self) -> bool {
        use EventType::*;
        matches!(
            self,
            IntronShift
                | ExonMisalignment
                | TerminalExonMisalignmentRight
                | TerminalExonMisalignmentLeft
                | FakeTerminalExonLeft
                | FakeTerminalExonRight
                | FakeMicroIntronRetention
        )
    }
    /// `MatchEventSubtype.is_minor_error`, IA:224-234.
    pub const fn is_minor_error(self) -> bool {
        use EventType::*;
        matches!(
            self,
            ExonElongationLeft
                | ExonElongationRight
                | IntronShift
                | ExonMisalignment
                | FakeMicroIntronRetention
                | FakeTerminalExonLeft
                | FakeTerminalExonRight
                | TerminalExonMisalignmentLeft
                | TerminalExonMisalignmentRight
        )
    }
    /// `MatchEventSubtype.is_consistent`, IA:236-250.
    pub const fn is_consistent(self) -> bool {
        use EventType::*;
        matches!(
            self,
            None | MonoExonic
                | MonoExonMatch
                | Fsm
                | IsmLeft
                | IsmRight
                | IsmInternal
                | CorrectPolyaSiteLeft
                | CorrectPolyaSiteRight
                | TerminalSiteMatchLeft
                | TerminalSiteMatchRight
                | TerminalSiteMatchLeftPrecise
                | TerminalSiteMatchRightPrecise
        )
    }
    /// IA:252-255.
    pub const fn is_major_elongation(self) -> bool {
        matches!(self, EventType::MajorExonElongationLeft | EventType::MajorExonElongationRight)
    }
    /// IA:257-260.
    pub const fn is_minor_elongation(self) -> bool {
        matches!(self, EventType::ExonElongationLeft | EventType::ExonElongationRight)
    }
    /// Member of `nnic_event_types`, IA:362-370.
    pub const fn is_nnic(self) -> bool {
        use EventType::*;
        matches!(
            self,
            AltLeftSiteNovel
                | AltRightSiteNovel
                | ExtraIntronNovel
                | ExtraIntronFlankingLeft
                | ExtraIntronFlankingRight
                | MutuallyExclusiveExonsNovel
                | ExonGainNovel
                | ExonSkippingNovel
                | ExonDetachNovel
                | ExonMergeNovel
                | TerminalExonShiftNovel
                | AlternativeStructureNovel
                | IntronAlternationNovel
        )
    }
    /// Member of `nic_event_types`, IA:376-390.
    pub const fn is_nic(self) -> bool {
        use EventType::*;
        matches!(
            self,
            UnsplicedIntronRetention
                | IntronRetention
                | AltLeftSiteKnown
                | AltRightSiteKnown
                | ExtraIntronKnown
                | IntronMigration
                | MutuallyExclusiveExonsKnown
                | ExonSkippingKnown
                | ExonDetachKnown
                | ExonMergeKnown
                | TerminalExonShiftKnown
                | ExonGainKnown
                | AlternativeStructureKnown
                | IntronAlternationKnown
                | MajorExonElongationLeft
                | MajorExonElongationRight
                | IncompleteIntronRetentionLeft
                | IncompleteIntronRetentionRight
                | InternalPolyaRight
                | InternalPolyaLeft
                | AlternativePolyaSiteLeft
                | AlternativePolyaSiteRight
                | AlternativeTssRight
                | AlternativeTssLeft
        )
    }
    /// Member of `nonintronic_events`, IA:392-398.
    pub const fn is_nonintronic(self) -> bool {
        use EventType::*;
        matches!(
            self,
            AlternativePolyaSiteLeft
                | AlternativePolyaSiteRight
                | AlternativeTssLeft
                | AlternativeTssRight
                | InternalPolyaRight
                | InternalPolyaLeft
                | MajorExonElongationLeft
                | MajorExonElongationRight
                | ExonElongationLeft
                | ExonElongationRight
        )
    }
    /// `is_major_inconsistency`: member of `all_major_events = nic ∪ nnic`, IA:262-264, 400.
    pub const fn is_major_inconsistency(self) -> bool {
        self.is_nic() || self.is_nnic()
    }
    /// `is_intronic_inconsistency`: `all_major_events − nonintronic_events`, IA:266-268, 402.
    pub const fn is_intronic_inconsistency(self) -> bool {
        self.is_major_inconsistency() && !self.is_nonintronic()
    }
}

/// `elongation_cost(params, elongation_len)`, IA:349-359, with `minor_exon_extension` /
/// `major_exon_extension` passed explicitly. Operation order matches Python for f64 parity.
pub fn elongation_cost(minor_exon_extension: i64, major_exon_extension: i64, elongation_len: i64) -> f64 {
    let min_cost = 0.1_f64; // exon_elongation_left
    let max_cost = 0.6_f64; // major_exon_elongation_left
    if elongation_len <= minor_exon_extension {
        min_cost
    } else if elongation_len >= major_exon_extension {
        max_cost
    } else {
        min_cost
            + (max_cost - min_cost) * ((elongation_len - minor_exon_extension) as f64)
                / ((major_exon_extension - minor_exon_extension) as f64)
    }
}

py_enum! {
    /// `TranscriptModelType`, `gene_info.py:28-31`.
    pub enum TranscriptModelType: u8 {
        Known = 1 => "known",
        NovelInCatalog = 2 => "novel_in_catalog",
        NovelNotInCatalog = 10 => "novel_not_in_catalog",
    }
}

py_enum! {
    /// `AmbiguityResolvingMethod`, `assignment/long_read_assigner.py:48-52`. (The Python
    /// class also holds the constants `minimal_score = -0.5` and `top_scored_factor = 1.5`,
    /// `:54-55`; see [`crate::params::MINIMAL_SCORE`], [`crate::params::TOP_SCORED_FACTOR`].)
    pub enum AmbiguityResolving: i8 {
        None = -10 => "none",
        MonoexonOnly = 10 => "monoexon_only",
        MonoexonAndFsm = 20 => "monoexon_and_fsm",
        All = 30 => "all",
    }
}

py_enum! {
    /// `StrandnessReportingLevel`, `model_construction/context.py:29-33`.
    pub enum StrandnessReporting: u8 {
        OnlyCanonical = 1 => "only_canonical",
        OnlyStranded = 2 => "only_stranded",
        All = 3 => "all",
        Auto = 10 => "auto",
    }
}

py_enum! {
    /// `PolyAUsageStrategies`, `dataset_processor.py:86-89`.
    pub enum PolyAUsage: u8 {
        Auto = 1 => "auto",
        Never = 2 => "never",
        Always = 3 => "always",
    }
}

py_enum! {
    /// `MultimapResolvingStrategy`, `assignment/multimap_resolver.py:19-22`.
    pub enum MultimapStrategy: u8 {
        IgnoreMultimapper = 0 => "ignore_multimapper",
        Merge = 1 => "merge",
        TakeBest = 2 => "take_best",
    }
}

py_enum! {
    /// `PolyATrimmed`, `alignment/alignment_processor.py:48-55`. Only `none` is in scope.
    pub enum PolyATrimmed: u8 {
        None = 1 => "none",
        All = 2 => "all",
        Stranded = 3 => "stranded",
        Tag = 4 => "tag",
        List = 5 => "list",
        Flnc = 6 => "flnc",
    }
}

/// `SupplementaryMatchConstants`, IA:412-419. Kept as `i64` because `undefined_position`
/// (`1 << 31`) does not fit in `i32`.
pub mod match_consts {
    /// `extra_left_mod_position = (1 << 30) - 1`.
    pub const EXTRA_LEFT_MOD_POSITION: i64 = (1 << 30) - 1;
    /// `extra_right_mod_position = (1 << 30) + 1`.
    pub const EXTRA_RIGHT_MOD_POSITION: i64 = (1 << 30) + 1;
    /// `undefined_position = 1 << 31`.
    pub const UNDEFINED_POSITION: i64 = 1 << 31;
    /// `absent_position = (1 << 31) - 1`.
    pub const ABSENT_POSITION: i64 = (1 << 31) - 1;
    pub const UNDEFINED_REGION: (i64, i64) = (UNDEFINED_POSITION, UNDEFINED_POSITION);
    pub const EXTRA_LEFT_REGION: (i64, i64) = (EXTRA_LEFT_MOD_POSITION, EXTRA_LEFT_MOD_POSITION);
    pub const EXTRA_RIGHT_REGION: (i64, i64) = (EXTRA_RIGHT_MOD_POSITION, EXTRA_RIGHT_MOD_POSITION);
}
pub use match_consts::ABSENT_POSITION;

#[cfg(test)]
mod tests {
    use super::*;

    /// Every Python member, by (name, code), copied from the 4.1 source. Guards against
    /// typos in the macro tables above: each must round-trip and the counts must match.
    #[test]
    fn assignment_type_codes() {
        let py: &[(&str, u8)] = &[
            ("unique", 1),
            ("noninformative", 0),
            ("intergenic", 20),
            ("ambiguous", 10),
            ("unique_minor_difference", 2),
            ("inconsistent", 30),
            ("inconsistent_non_intronic", 31),
            ("inconsistent_ambiguous", 32),
            ("inconsistent_genic", 33),
            ("inconsistent_multigenic", 34),
            ("discarded", 40),
            ("suspended", 255),
        ];
        check(py, AssignmentType::ALL, AssignmentType::from_py_name, |x| x.code(), AssignmentType::from_code);
    }

    #[test]
    fn match_class_codes() {
        let py: &[(&str, u8)] = &[
            ("undefined", 0),
            ("full_splice_match", 10),
            ("incomplete_splice_match", 11),
            ("mono_exon_match", 12),
            ("novel_in_catalog", 20),
            ("novel_not_in_catalog", 21),
            ("genic", 30),
            ("antisense", 40),
            ("fusion", 50),
            ("intergenic", 32),
            ("genic_intron", 31),
        ];
        check(py, MatchClass::ALL, MatchClass::from_py_name, |x| x.code(), MatchClass::from_code);
    }

    #[test]
    fn event_type_codes() {
        let py: &[(&str, u16)] = &[
            ("none", 0),
            ("undefined", 1),
            ("mono_exonic", 10),
            ("antisense", 9),
            ("fsm", 11),
            ("ism_left", 15),
            ("ism_right", 13),
            ("ism_internal", 14),
            ("mono_exon_match", 19),
            ("intron_shift", 21),
            ("exon_misalignment", 22),
            ("fake_terminal_exon_left", 24),
            ("fake_terminal_exon_right", 26),
            ("terminal_exon_misalignment_right", 28),
            ("terminal_exon_misalignment_left", 29),
            ("exon_elongation_left", 25),
            ("exon_elongation_right", 23),
            ("intron_retention", 31),
            ("unspliced_intron_retention", 32),
            ("incomplete_intron_retention_left", 38),
            ("incomplete_intron_retention_right", 39),
            ("fake_micro_intron_retention", 34),
            ("alt_left_site_known", 101),
            ("alt_right_site_known", 102),
            ("alt_left_site_novel", 103),
            ("alt_right_site_novel", 104),
            ("extra_intron_novel", 1012),
            ("extra_intron_known", 1011),
            ("extra_intron_flanking_left", 1013),
            ("extra_intron_flanking_right", 1015),
            ("major_exon_elongation_left", 1005),
            ("major_exon_elongation_right", 1003),
            ("intron_migration", 114),
            ("intron_alternation_novel", 115),
            ("intron_alternation_known", 116),
            ("mutually_exclusive_exons_novel", 121),
            ("mutually_exclusive_exons_known", 122),
            ("exon_skipping_known", 133),
            ("exon_skipping_novel", 134),
            ("exon_merge_known", 135),
            ("exon_merge_novel", 136),
            ("exon_gain_known", 145),
            ("exon_gain_novel", 146),
            ("exon_detach_known", 147),
            ("exon_detach_novel", 148),
            ("terminal_exon_shift_known", 151),
            ("terminal_exon_shift_novel", 152),
            ("alternative_structure_novel", 181),
            ("alternative_structure_known", 182),
            ("alternative_polya_site_right", 200),
            ("alternative_polya_site_left", 201),
            ("internal_polya_left", 202),
            ("internal_polya_right", 203),
            ("alternative_tss_left", 204),
            ("alternative_tss_right", 205),
            ("correct_polya_site_left", 222),
            ("correct_polya_site_right", 223),
            ("aligned_polya_tail", 300),
            ("terminal_site_match_left", 250),
            ("terminal_site_match_right", 251),
            ("terminal_site_match_left_precise", 252),
            ("terminal_site_match_right_precise", 253),
        ];
        check(py, EventType::ALL, EventType::from_py_name, |x| x.code(), EventType::from_code);
    }

    #[test]
    fn small_enum_codes() {
        check(
            &[("known", 1u8), ("novel_in_catalog", 2), ("novel_not_in_catalog", 10)],
            TranscriptModelType::ALL,
            TranscriptModelType::from_py_name,
            |x| x.code(),
            TranscriptModelType::from_code,
        );
        check(
            &[("none", -10i8), ("monoexon_only", 10), ("monoexon_and_fsm", 20), ("all", 30)],
            AmbiguityResolving::ALL,
            AmbiguityResolving::from_py_name,
            |x| x.code(),
            AmbiguityResolving::from_code,
        );
        check(
            &[("only_canonical", 1u8), ("only_stranded", 2), ("all", 3), ("auto", 10)],
            StrandnessReporting::ALL,
            StrandnessReporting::from_py_name,
            |x| x.code(),
            StrandnessReporting::from_code,
        );
        check(
            &[("auto", 1u8), ("never", 2), ("always", 3)],
            PolyAUsage::ALL,
            PolyAUsage::from_py_name,
            |x| x.code(),
            PolyAUsage::from_code,
        );
        check(
            &[("ignore_multimapper", 0u8), ("merge", 1), ("take_best", 2)],
            MultimapStrategy::ALL,
            MultimapStrategy::from_py_name,
            |x| x.code(),
            MultimapStrategy::from_code,
        );
        check(
            &[("none", 1u8), ("all", 2), ("stranded", 3), ("tag", 4), ("list", 5), ("flnc", 6)],
            PolyATrimmed::ALL,
            PolyATrimmed::from_py_name,
            |x| x.code(),
            PolyATrimmed::from_code,
        );
    }

    fn check<T, C>(
        py: &[(&str, C)],
        all: &[T],
        by_name: fn(&str) -> Option<T>,
        code: fn(T) -> C,
        by_code: fn(C) -> Option<T>,
    ) where
        T: Copy + PartialEq + std::fmt::Debug + Serialize + for<'de> Deserialize<'de>,
        C: Copy + PartialEq + std::fmt::Debug + Into<i64>,
    {
        assert_eq!(py.len(), all.len(), "member count differs from 4.1");
        for &(name, c) in py {
            let v = by_name(name).unwrap_or_else(|| panic!("missing member {name}"));
            assert_eq!(code(v), c, "{name}");
            assert_eq!(by_code(c), Some(v));
            let j = serde_json::to_string(&v).unwrap();
            assert_eq!(j, c.into().to_string());
            assert_eq!(serde_json::from_str::<T>(&j).unwrap(), v);
        }
        assert!(serde_json::from_str::<T>("12345").is_err());
    }

    #[test]
    fn event_cost_table_complete() {
        // every member except antisense has a cost (IA:271-346 has 62 entries)
        let with_cost = EventType::ALL.iter().filter(|e| e.cost().is_some()).count();
        assert_eq!(with_cost, EventType::ALL.len() - 1);
        assert_eq!(EventType::Antisense.cost(), None);
        assert_eq!(EventType::ExonMergeKnown.cost(), Some(0.75));
        assert_eq!(EventType::IncompleteIntronRetentionRight.cost(), Some(0.7));
    }

    #[test]
    fn event_sets() {
        // sizes of the Python sets
        let n = |f: fn(EventType) -> bool| EventType::ALL.iter().filter(|e| f(**e)).count();
        assert_eq!(n(EventType::is_nnic), 13);
        assert_eq!(n(EventType::is_nic), 24);
        assert_eq!(n(EventType::is_nonintronic), 10);
        assert_eq!(n(EventType::is_major_inconsistency), 37);
        assert_eq!(n(EventType::is_intronic_inconsistency), 37 - 8);
        assert_eq!(n(EventType::is_alignment_artifact), 7);
        assert_eq!(n(EventType::is_minor_error), 9);
        assert_eq!(n(EventType::is_consistent), 13);
    }

    #[test]
    fn assignment_type_predicates() {
        use AssignmentType::*;
        assert!(InconsistentGenic.is_inconsistent() && !InconsistentGenic.is_ambiguous());
        assert!(InconsistentMultigenic.is_ambiguous());
        assert!(Intergenic.is_unassigned() && !Discarded.is_unassigned());
        assert!(UniqueMinorDifference.is_unique() && UniqueMinorDifference.is_consistent());
    }

    #[test]
    fn elongation() {
        assert_eq!(elongation_cost(50, 300, 10), 0.1);
        assert_eq!(elongation_cost(50, 300, 300), 0.6);
        // Python: 0.1 + (0.6 - 0.1) * (175 - 50) / (300 - 50)
        assert_eq!(elongation_cost(50, 300, 175), 0.1 + (0.6 - 0.1) * 125.0 / 250.0);
    }
}
