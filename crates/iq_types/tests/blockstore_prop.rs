//! BlockStore property tests: `from_records → to_records` is the identity, and chain
//! interning stores each distinct gap chain once.

use iq_types::*;
use proptest::prelude::*;
use std::collections::HashSet;

fn exon_list(max_exons: usize) -> impl Strategy<Value = Vec<Iv>> {
    // sorted exons; gaps may be 0 or negative (touching/overlapping exons) to exercise the
    // lossless encoding of non-intron gaps
    (1i32..10_000, prop::collection::vec((1i32..300, -3i32..400), 0..max_exons)).prop_map(|(start, parts)| {
        let mut out = Vec::new();
        let mut s = start;
        for (len, gap) in parts {
            let e = s + len;
            out.push(Iv::new(s, e));
            s = e + 1 + gap;
        }
        out
    })
}

fn event() -> impl Strategy<Value = MatchEvent> {
    (
        prop::sample::select(EventType::ALL.to_vec()),
        prop_oneof![Just(match_consts::UNDEFINED_REGION), (-1i64..50, -1i64..50)],
        (-1i64..50, -1i64..(1i64 << 31)),
        -100i64..100,
    )
        .prop_map(|(ty, iso_region, read_region, info)| MatchEvent { ty, iso_region, read_region, info })
}

fn strand() -> impl Strategy<Value = Strand> {
    prop::sample::select(Strand::ALL.to_vec())
}

fn gene_key() -> impl Strategy<Value = Option<GeneKey>> {
    prop_oneof![
        Just(None),
        (0u32..5).prop_map(|g| Some(GeneKey::Ref(g))),
        (0u32..3).prop_map(|g| Some(GeneKey::Novel(g)))
    ]
}
fn iso_ref() -> impl Strategy<Value = Option<IsoRef>> {
    prop_oneof![
        Just(None),
        (0u32..9).prop_map(|t| Some(IsoRef::Ref(t))),
        (0u32..9).prop_map(|t| Some(IsoRef::Model(t)))
    ]
}

fn isoform_match() -> impl Strategy<Value = IsoformMatch> {
    (
        gene_key(),
        iso_ref(),
        strand(),
        prop::sample::select(MatchClass::ALL.to_vec()),
        -5.0f64..5.0,
        prop::collection::vec(event(), 0..4),
    )
        .prop_map(|(gene, tx, tx_strand, class, penalty, ev)| IsoformMatch {
            gene,
            tx,
            tx_strand,
            class,
            penalty,
            events: ev.into_iter().collect(),
        })
}

fn record(chains: Vec<Vec<Iv>>) -> impl Strategy<Value = ReadRecord> {
    (
        (prop::sample::select(chains), prop::option::of(exon_list(6))),
        prop::collection::vec((1i32..10_000, 1u32..20), 0..3),
        prop::array::uniform4(-1i32..10_000),
        any::<bool>(),
        strand(),
        strand(),
        any::<u8>(),
        any::<bool>(),
        0u32..4,
        (
            prop::sample::select(AssignmentType::ALL.to_vec()),
            prop::sample::select(AssignmentType::ALL.to_vec()),
            prop::collection::vec(isoform_match(), 0..3),
        ),
        (prop::bool::weighted(0.9), "[a-zA-Z0-9:_/]{1,20}"),
    )
        .prop_map(
            |(
                (corrected, raw),
                dels,
                pa,
                polya_found,
                st,
                mst,
                mapq,
                mm,
                group,
                (at, gat, matches),
                (same_region, name),
            )| {
                ReadRecord {
                    idx: 0,
                    name,
                    exons: raw.unwrap_or_else(|| corrected.clone()),
                    corrected_exons: corrected,
                    boundary_deletions: dels,
                    polya: pa.into(),
                    polya_found,
                    strand: st,
                    mapped_strand: mst,
                    mapq,
                    multimapper: mm,
                    group,
                    assignment: Assignment { atype: at, gene_atype: gat, matches: matches.into_iter().collect() },
                    genomic_region: if same_region { Iv::new(1, 100_000) } else { Iv::new(5, 50) },
                }
            },
        )
}

fn block(dense_idx: bool) -> impl Strategy<Value = Vec<ReadRecord>> {
    // a small pool of chains so that reads share them
    prop::collection::vec(exon_list(8), 1..6).prop_flat_map(|pool| prop::collection::vec(record(pool), 0..40)).prop_map(
        move |mut v| {
            for (i, r) in v.iter_mut().enumerate() {
                r.idx = if dense_idx { i as u32 } else { (i as u32) * 3 + 1 };
            }
            v
        },
    )
}

fn gaps(ex: &[Iv]) -> Vec<Iv> {
    ex.windows(2).map(|w| Iv::new(w[0].e + 1, w[1].s - 1)).collect()
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 256, ..ProptestConfig::default() })]

    #[test]
    fn roundtrip_identity(recs in block(true)) {
        let st = BlockStore::from_records(&recs);
        prop_assert_eq!(st.len(), recs.len());
        prop_assert_eq!(st.to_records(), recs);
    }

    #[test]
    fn roundtrip_identity_sparse_idx(recs in block(false)) {
        let st = BlockStore::from_records(&recs);
        prop_assert_eq!(st.to_records(), recs);
    }

    #[test]
    fn views_agree_with_records(recs in block(true)) {
        let st = BlockStore::from_records(&recs);
        for (v, r) in st.iter().zip(&recs) {
            prop_assert_eq!(v.corrected_introns().into_vec(), junctions_from_blocks(&r.corrected_exons));
            if let Some(s) = v.corrected_introns_slice() {
                prop_assert_eq!(s.to_vec(), r.corrected_introns());
            }
            prop_assert_eq!(v.n_corrected_exons(), r.corrected_exons.len());
            prop_assert_eq!(v.exons_differ(), r.exons != r.corrected_exons);
            prop_assert_eq!(v.name_hash(), name_hash(r.name.as_bytes()));
            let span = r.exons.first().map(|a| (a.s, r.exons.last().unwrap().e)).unwrap_or((0, 0));
            prop_assert_eq!(v.exons_span(), span);
            prop_assert_eq!(v.n_matches(), r.assignment.matches.len());
        }
    }

    #[test]
    fn chain_interning_dedups(recs in block(true)) {
        let st = BlockStore::from_records(&recs);
        let mut distinct: HashSet<Vec<Iv>> = HashSet::new();
        distinct.insert(Vec::new()); // the empty chain is always present
        for r in &recs {
            for ex in [&r.corrected_exons, &r.exons] {
                if ex.len() > 1 { distinct.insert(gaps(ex)); }
            }
        }
        prop_assert_eq!(st.chains.len(), distinct.len());
        // equal corrected chains ⇔ equal chain ids
        for (a, ra) in st.iter().zip(&recs) {
            for (b, rb) in st.iter().zip(&recs) {
                if ra.corrected_exons.len() > 1 && rb.corrected_exons.len() > 1 {
                    prop_assert_eq!(a.chain_id() == b.chain_id(), gaps(&ra.corrected_exons) == gaps(&rb.corrected_exons));
                }
            }
        }
    }
}

#[test]
fn shared_chains_stored_once() {
    let mk = |i: u32, ex: Vec<Iv>| ReadRecord {
        idx: i,
        name: format!("r{i}"),
        exons: ex.clone(),
        corrected_exons: ex,
        boundary_deletions: vec![],
        polya: NO_POLYA,
        polya_found: false,
        strand: Strand::Plus,
        mapped_strand: Strand::Plus,
        mapq: 60,
        multimapper: false,
        group: 0,
        assignment: Assignment::empty(AssignmentType::Noninformative),
        genomic_region: Iv::new(1, 1000),
    };
    // same introns, different terminal exon ends → one chain
    let recs: Vec<_> = (0..1000)
        .map(|i| mk(i, vec![Iv::new(100 + (i % 7) as i32, 200), Iv::new(300, 400), Iv::new(500, 600 + (i % 5) as i32)]))
        .collect();
    let st = BlockStore::from_records(&recs);
    assert_eq!(st.chains.len(), 2); // empty + one
    assert_eq!(st.chains.n_intervals(), 2);
    assert_eq!(st.to_records(), recs);
}
