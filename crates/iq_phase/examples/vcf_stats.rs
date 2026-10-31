//! `vcf_stats VCF [SAMPLE]`: phased het SNVs per contig.
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let s = iq_phase::vcf::load(std::path::Path::new(&a[1]), a.get(2).map(|s| s.as_str())).unwrap();
    let mut v: Vec<_> = s.chrs.iter().map(|(k, c)| (k.clone(), c.len())).collect();
    v.sort();
    for (k, n) in v {
        println!("{k}\t{n}");
    }
}
