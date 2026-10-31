//! # iq_bam — U2: BAM and FASTA IO
//!
//! noodles-based `fetch` with parallel BGZF inflate yielding
//! [`RecordBatch`]es of lazily parsed records ([`iq_types::RecordLite`] views), and an
//! mmap'd FASTA ([`FastaIndex`]). `--features htslib` swaps in a rust-htslib
//! [`RecordSource`] implementation.
//!
//! Fixture runners: unit `U2`, fn `fetch` ([`register_fixtures`]).

pub mod bam;
pub mod fasta;
#[cfg(feature = "htslib")]
pub mod htslib;
#[cfg(feature = "htslib")]
pub use htslib::HtsBam;

pub use bam::{
    decode_seq_into, fetch, fetch_with, try_fetch, BamHandle, BamOptions, FetchIter, FetchOpts, RecordBatch,
    RecordSource,
};
pub use fasta::FastaIndex;

/// Registers this crate's fixture runners.
pub fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    bam::register_fixtures(reg);
}

#[cfg(test)]
mod fixture_tests {
    use iq_types::fixture::{run_path, Registry, Summary};
    use std::path::PathBuf;

    /// Runs the U2 fixtures (`$IQRS_FIXTURES` or `../fixtures`) when they exist.
    #[test]
    fn u2_fixtures_exact() {
        let root = iq_types::fixture::fixtures_root()
            .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../../fixtures"));
        if !root.is_dir() {
            eprintln!("no fixtures at {}; skipped", root.display());
            return;
        }
        let mut reg = Registry::new();
        super::register_fixtures(&mut reg);
        let mut s = Summary::new(5);
        for e in std::fs::read_dir(&root).unwrap() {
            let u2 = e.unwrap().path().join("u2");
            if u2.is_dir() {
                run_path(&reg, "U2", &u2, &mut s).unwrap();
            }
        }
        eprintln!("{s}");
        assert!(s.all_ok(), "{s}");
        let t = &s.by_key[&("U2".to_owned(), "fetch".to_owned())];
        assert!(t.exact > 0 && t.mismatch == 0);
    }
}
