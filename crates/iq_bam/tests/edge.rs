//! Edge cases against pysam: `tests/data/edge.bam` and its U2 fixture file are produced by
//! `tests/gen_edge_bam.py` (IUPAC bases, every HP tag type, `*` sequence, missing quality,
//! placed-unmapped and zero-reference-length reads, a 70k-op CIGAR stored as `CG:B:I`,
//! 150 random regions). cargo runs tests with the crate dir as cwd.

use iq_bam::{try_fetch, BamHandle, FetchOpts, RecordBatch};
use iq_types::fixture::{read_lines, run_path, Case, Registry, Summary};
use iq_types::Iv0;
use std::collections::BTreeMap;
use std::path::Path;

const BAM: &str = "tests/data/edge.bam";

#[test]
fn edge_fixture_exact() {
    let mut reg = Registry::new();
    iq_bam::register_fixtures(&mut reg);
    let mut s = Summary::new(5);
    run_path(&reg, "U2", Path::new("tests/data/edge"), &mut s).unwrap();
    eprintln!("{s}");
    let t = &s.by_key[&("U2".to_owned(), "fetch".to_owned())];
    assert!(t.exact > 10_000 && t.mismatch == 0 && t.decode_error == 0, "{s}");
}

/// The runner only checks record `n` exists; also check no extra records are returned.
#[test]
fn edge_record_counts() {
    let mut want: BTreeMap<String, u64> = BTreeMap::new();
    for c in read_lines::<Case>(Path::new("tests/data/edge/units/U2.jsonl.gz")).unwrap() {
        let c = c.unwrap();
        *want.entry(c.input.get("region").unwrap().as_str().unwrap().to_owned()).or_default() += 1;
    }
    let bam = BamHandle::open(Path::new(BAM)).unwrap();
    for (region, n) in &want {
        let (chr, se) = region.rsplit_once(':').unwrap();
        let (s, e) = se.split_once('-').unwrap();
        let iv = Iv0 { s: s.parse::<i64>().unwrap() - 1, e: e.parse::<i64>().unwrap() - 1 };
        let got: usize = try_fetch(&bam, chr, iv, FetchOpts::default()).unwrap().map(|b| b.unwrap().len()).sum();
        assert_eq!(got as u64, *n, "{region}");
    }
    // regions without any record are absent from the fixture; spot-check one
    let got: usize = try_fetch(&bam, "c2", Iv0 { s: 99_999, e: 99_999 }, FetchOpts::default()).unwrap().count();
    assert_eq!(got, 0);
}

#[test]
fn unmapped_matches_pysam() {
    let bam = BamHandle::open(Path::new(BAM)).unwrap();
    let want: u64 = std::fs::read_to_string("tests/data/edge.unmapped.txt").unwrap().trim().parse().unwrap();
    assert_eq!(bam.unmapped(), want);
    assert_eq!(bam.references(), ["c1", "c2"]);
    assert_eq!(bam.reference_lengths(), [3_000_000, 100_000]);
}

#[test]
fn no_seq_option_and_batching() {
    let mut bam = BamHandle::open(Path::new(BAM)).unwrap();
    bam.set_tag(Some(*b"HP"));
    let all = |opts: FetchOpts, bam: &BamHandle| -> Vec<RecordBatch> {
        try_fetch(bam, "c1", Iv0 { s: 0, e: 2_999_999 }, opts).unwrap().map(|b| b.unwrap()).collect()
    };
    let a = all(FetchOpts { seq: false, batch_records: 7 }, &bam);
    let b = all(FetchOpts::default(), &bam);
    assert!(a.len() > 1);
    assert!(a.iter().enumerate().all(|(k, x)| x.seq_no == k as u64 && x.len() <= 7));
    let ra: Vec<_> = a
        .iter()
        .flat_map(|x| x.iter().map(|r| (r.qname.to_vec(), r.pos0, r.cigar.to_vec(), r.tag.cloned(), r.seq.len())))
        .collect();
    let rb: Vec<_> = b
        .iter()
        .flat_map(|x| x.iter().map(|r| (r.qname.to_vec(), r.pos0, r.cigar.to_vec(), r.tag.cloned(), 0)))
        .collect();
    assert_eq!(ra, rb);
    // multithreaded inflate gives the same records
    bam.set_threads(3);
    let c = all(FetchOpts::default(), &bam);
    let rc: Vec<_> = c.iter().flat_map(|x| x.iter().map(|r| (r.qname.to_vec(), r.seq.to_vec()))).collect();
    let rb2: Vec<_> = b.iter().flat_map(|x| x.iter().map(|r| (r.qname.to_vec(), r.seq.to_vec()))).collect();
    assert_eq!(rc, rb2);
}
