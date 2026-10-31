//! Phase W and the final writers.
//!
//! Per chromosome, blocks are consumed in genomic order by [`ChrPhaseW`]: models are named
//! ([`crate::ids::NovelNamer`]), the GTF fragment is printed ([`crate::gtf::GtfChrWriter`]),
//! and block counters are merged by final feature id into [`ChrCounts`] (4.1's
//! per-chromosome `AssignedFeatureCounter` is keyed by feature string, so a known transcript
//! that appears in several blocks is summed). Chromosomes are independent, so U21 may run
//! them in parallel. [`write_outputs`] then concatenates the chromosome outputs in 4.1's
//! fragment order (`merge_files`, `utils/file_utils.py:157-182`) and writes the counts
//! stats rows (`merge_counts`, `:185-213`).

use crate::counter::{Groups, TxCounter};
use crate::gtf::{print_models, AnnotBlock, GtfChrWriter};
use crate::ids::{BlockNames, NovelNamer};
use iq_annot::AnnotationIndex;
use iq_types::{GroupPool, IsoformIndex, TranscriptModel};
use rustc_hash::FxHashMap;
use std::fmt::Write as _;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

/// `'%.2f' % x`. Python formats the exact binary value with round-half-even, as Rust's
/// `{:.2}` does (exact decimal expansion; checked in tests for exact binary ties).
pub fn fmt_fixed2(x: f64) -> String {
    format!("{x:.2}")
}

/// Appends `'%.2f' % x`.
fn push_fixed2(out: &mut String, x: f64) {
    let _ = write!(out, "{x:.2}");
}

// ------------------------------------------------------------------------------------------
// Per-chromosome counts
// ------------------------------------------------------------------------------------------

/// One feature of a chromosome counter.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChrFeat {
    pub id: String,
    pub total: f64,
    pub groups: Groups,
    /// In `all_features` (printed).
    pub in_all: bool,
    /// In `confirmed_features`; unconfirmed features print 0 (`dump`, `:398-405`). Every
    /// model of a block is confirmed (`model_read_counter.py:169`), so production entries
    /// are always confirmed.
    pub confirmed: bool,
}

/// 4.1's per-chromosome transcript-model counter state (ungrouped + grouped).
#[derive(Clone, Debug, Default)]
pub struct ChrCounts {
    pub feats: Vec<ChrFeat>,
    index: FxHashMap<String, u32>,
    pub ambiguous: u64,
    pub no_feature: u64,
}

impl ChrCounts {
    pub fn new() -> Self {
        Self::default()
    }

    /// Entry of feature `id`, created if new.
    pub fn entry(&mut self, id: &str) -> &mut ChrFeat {
        let k = match self.index.get(id) {
            Some(&k) => k,
            None => {
                let k = self.feats.len() as u32;
                self.index.insert(id.to_owned(), k);
                self.feats.push(ChrFeat { id: id.to_owned(), confirmed: true, ..Default::default() });
                k
            }
        };
        &mut self.feats[k as usize]
    }

    /// Merges one block's counter (blocks must come in genomic order). `models` is the
    /// block's final storage, `names` its final names (storage order).
    ///
    /// Sums across blocks are `chr_total + block_total`; 4.1 adds read by read. Under
    /// `unique_only` every weight is 0 or 1, so both are exact integers.
    pub fn add_block(&mut self, counter: &TxCounter, models: &[TranscriptModel], names: &BlockNames) {
        for (i, m) in models.iter().enumerate() {
            let Some(s) = counter.state(m.key) else {
                continue;
            };
            let e = self.entry(&names.tx[i]);
            e.total += s.total;
            e.in_all |= s.in_all;
            for &(g, v) in &s.groups {
                match e.groups.iter_mut().find(|(h, _)| *h == g) {
                    Some(x) => x.1 += v,
                    None => e.groups.push((g, v)),
                }
            }
        }
        self.ambiguous += counter.ambiguous();
        self.no_feature += counter.no_feature();
    }

    /// Printed features (`all_features`) sorted by Python string order.
    fn printed(&self) -> Vec<&ChrFeat> {
        let mut v: Vec<&ChrFeat> = self.feats.iter().filter(|f| f.in_all).collect();
        v.sort_unstable_by(|a, b| a.id.as_bytes().cmp(b.id.as_bytes()));
        v
    }

    /// `dump_ungrouped` rows (no header): `"%s\t%.2f\n"`.
    pub fn write_ungrouped_rows(&self, out: &mut String) {
        for f in self.printed() {
            out.push_str(&f.id);
            out.push('\t');
            push_fixed2(out, if f.confirmed { f.total } else { 0.0 });
            out.push('\n');
        }
    }

    /// `dump_grouped` rows (no header): `"%s\t%s\t%.2f\n"` for groups with count != 0, in
    /// first-touch order. `group_name` resolves a group id.
    pub fn write_grouped_rows<'a>(&self, group_name: impl Fn(u32) -> &'a str, out: &mut String) {
        for f in self.printed() {
            for &(g, v) in &f.groups {
                let v = if f.confirmed { v } else { 0.0 };
                if v != 0.0 {
                    out.push_str(&f.id);
                    out.push('\t');
                    out.push_str(group_name(g));
                    out.push('\t');
                    push_fixed2(out, v);
                    out.push('\n');
                }
            }
        }
    }

    /// The per-chromosome `.stats` rows.
    pub fn stats_rows(&self, not_aligned: u64) -> String {
        format!("__ambiguous\t{}\n__no_feature\t{}\n__not_aligned\t{}\n", self.ambiguous, self.no_feature, not_aligned)
    }
}

/// Group display name: the `GroupPool` string, `NA` for `GroupPool::NONE`
/// (`AbstractReadGrouper.default_group_id`).
pub fn group_name(pool: &GroupPool, g: u32) -> &str {
    pool.resolve(g).unwrap_or("NA")
}

// ------------------------------------------------------------------------------------------
// Phase W per chromosome
// ------------------------------------------------------------------------------------------

/// One block's construct output as Phase W consumes it.
pub struct BlockOutput<'a> {
    /// The block's `gene_info` (its `genes` are the block genes).
    pub index: &'a IsoformIndex,
    /// Final model storage (storage order).
    pub models: &'a [TranscriptModel],
    pub counter: &'a TxCounter,
}

/// Finished output of one chromosome.
#[derive(Clone, Debug, Default)]
pub struct ChrOutput {
    pub chr: String,
    /// The chromosome's `transcript_models.gtf` fragment (no header).
    pub gtf: String,
    pub counts: ChrCounts,
}

/// Phase W state of one chromosome.
pub struct ChrPhaseW<'a> {
    ann: &'a AnnotationIndex,
    chr_idx: Option<iq_types::ChrIdx>,
    namer: NovelNamer,
    gtf: GtfChrWriter,
    out: ChrOutput,
}

impl<'a> ChrPhaseW<'a> {
    /// `chr` = chromosome name as in the BAM/FASTA (used in ids and GTF column 1).
    pub fn new(ann: &'a AnnotationIndex, chr: &str) -> Self {
        let chr_idx = ann.chrom(chr);
        let forbidden = chr_idx.map(|c| ann.forbidden_ids(c)).unwrap_or(&[]);
        ChrPhaseW {
            ann,
            chr_idx,
            namer: NovelNamer::new(chr, forbidden),
            gtf: GtfChrWriter::new(chr),
            out: ChrOutput { chr: chr.to_owned(), ..Default::default() },
        }
    }

    /// Consumes the next block (genomic order). Returns the block's final names.
    pub fn add_block(&mut self, b: &BlockOutput<'_>) -> BlockNames {
        let ann = self.ann;
        let names = self.namer.name_block(b.models, |t| ann.tx_id(t), |g| ann.gene_id(g));
        let pms = print_models(ann, b.models, &names);
        let annot = AnnotBlock { ann, chr: self.chr_idx, index: b.index, models: b.models };
        self.gtf.dump_block(&annot, &pms, &mut self.out.gtf);
        self.out.counts.add_block(b.counter, b.models, &names);
        names
    }

    pub fn finish(self) -> ChrOutput {
        self.out
    }
}

// ------------------------------------------------------------------------------------------
// Final files
// ------------------------------------------------------------------------------------------

/// 4.1 `merge_files` order key: `[int(t) if t.isdigit() else t.lower() for t in
/// re.split(r"(\d+)", s)]`. Digit runs compare numerically, the rest lowercased.
pub fn natural_key_cmp(a: &str, b: &str) -> std::cmp::Ordering {
    use std::cmp::Ordering;
    fn split(s: &str) -> Vec<(bool, &str)> {
        // alternating: text, digits, text, ... (text pieces may be empty, as in re.split)
        let mut v = Vec::new();
        let bytes = s.as_bytes();
        let mut i = 0;
        let mut start = 0;
        let mut in_digits = false;
        while i < bytes.len() {
            let d = bytes[i].is_ascii_digit();
            if d != in_digits {
                v.push((in_digits, &s[start..i]));
                start = i;
                in_digits = d;
            }
            i += 1;
        }
        v.push((in_digits, &s[start..]));
        if v.first().is_some_and(|x| x.0) {
            v.insert(0, (false, ""));
        }
        v
    }
    fn num_cmp(a: &str, b: &str) -> Ordering {
        let a = a.trim_start_matches('0');
        let b = b.trim_start_matches('0');
        a.len().cmp(&b.len()).then_with(|| a.cmp(b))
    }
    let (sa, sb) = (split(a), split(b));
    for (x, y) in sa.iter().zip(sb.iter()) {
        let o = if x.0 && y.0 { num_cmp(x.1, y.1) } else { x.1.to_lowercase().cmp(&y.1.to_lowercase()) };
        if o != Ordering::Equal {
            return o;
        }
    }
    sa.len().cmp(&sb.len())
}

/// Chromosome order of the fragments of the file whose basename is `<prefix><suffix>`
/// (fragments are `<prefix>_<chr with '/'→'_'><suffix>`). Returns indices into `chrs`.
pub fn fragment_order(prefix: &str, suffix: &str, chrs: &[&str]) -> Vec<usize> {
    let names: Vec<String> = chrs.iter().map(|c| format!("{prefix}_{}{suffix}", c.replace('/', "_"))).collect();
    let mut order: Vec<usize> = (0..chrs.len()).collect();
    order.sort_by(|&a, &b| natural_key_cmp(&names[a], &names[b]));
    order
}

/// Output configuration.
#[derive(Clone, Debug)]
pub struct OutputSpec {
    pub dir: PathBuf,
    /// File prefix (4.1 default `OUT`).
    pub prefix: String,
    /// `args._cmd_line` (argv joined by ' ').
    pub cmd_line: String,
    /// Version line value (D-O4: `4.1.0-iqrs-<ver>`).
    pub version: String,
    /// `--read-group-tag` (e.g. `HP`); `None` = no grouped output.
    pub group_tag: Option<String>,
    /// `--grouped-out PATH` (D-O3); default is the 4.1 name.
    pub grouped_out: Option<PathBuf>,
    /// gzip the three files (4.1 writes them plain; the default is `false`).
    pub gzip: bool,
}

impl OutputSpec {
    pub fn gtf_suffix() -> &'static str {
        ".transcript_models.gtf"
    }
    pub fn counts_suffix() -> &'static str {
        ".discovered_transcript_counts.tsv"
    }
    pub fn grouped_suffix(tag: &str) -> String {
        format!(".discovered_transcript_grouped_tag_{tag}_counts.linear.tsv")
    }
    fn path(&self, suffix: &str) -> PathBuf {
        let mut p = self.dir.join(format!("{}{suffix}", self.prefix));
        if self.gzip {
            p.as_mut_os_string().push(".gz");
        }
        p
    }
    pub fn gtf_path(&self) -> PathBuf {
        self.path(Self::gtf_suffix())
    }
    pub fn counts_path(&self) -> PathBuf {
        self.path(Self::counts_suffix())
    }
    pub fn grouped_path(&self) -> Option<PathBuf> {
        let tag = self.group_tag.as_ref()?;
        Some(self.grouped_out.clone().unwrap_or_else(|| self.path(&Self::grouped_suffix(tag))))
    }
}

fn open_out(path: &Path, gzip: bool) -> io::Result<Box<dyn Write>> {
    let f = BufWriter::with_capacity(1 << 20, File::create(path)?);
    Ok(if gzip { Box::new(flate2::write::GzEncoder::new(f, flate2::Compression::default())) } else { Box::new(f) })
}

/// The GTF header (`GFFPrinter.__init__` + `common_header`).
pub fn gtf_header(spec: &OutputSpec) -> String {
    format!(
        "# {} IsoQuant generated GTF\n# Command line: {}\n# IsoQuant version: {}\n",
        spec.prefix, spec.cmd_line, spec.version
    )
}

/// Text of the three files, in memory (`write_outputs` writes these).
pub struct FinalTexts {
    pub gtf: Vec<u8>,
    pub counts: String,
    pub grouped: Option<String>,
}

/// Assembles the final file contents from per-chromosome outputs (any order; 4.1 fragment
/// order is applied here). `not_aligned` = BAM unmapped count (`__not_aligned` when > 0).
pub fn assemble(spec: &OutputSpec, chrs: &[ChrOutput], pool: &GroupPool, not_aligned: u64) -> FinalTexts {
    let names: Vec<&str> = chrs.iter().map(|c| c.chr.as_str()).collect();
    let p = spec.prefix.as_str();

    let mut gtf = gtf_header(spec).into_bytes();
    for i in fragment_order(p, OutputSpec::gtf_suffix(), &names) {
        gtf.extend_from_slice(chrs[i].gtf.as_bytes());
    }

    // merge_counts: header of the first fragment, rows, then stats (dict insertion order of
    // the first stats file; __not_aligned replaced by the unaligned count when > 0)
    let mut counts = String::new();
    if !chrs.is_empty() {
        counts.push_str("feature_id\tcount\n");
    }
    for i in fragment_order(p, OutputSpec::counts_suffix(), &names) {
        chrs[i].counts.write_ungrouped_rows(&mut counts);
    }
    let amb: u64 = chrs.iter().map(|c| c.counts.ambiguous).sum();
    let nf: u64 = chrs.iter().map(|c| c.counts.no_feature).sum();
    if !chrs.is_empty() {
        let _ = write!(counts, "__ambiguous\t{amb}\n__no_feature\t{nf}\n__not_aligned\t{not_aligned}\n");
    } else if not_aligned > 0 {
        let _ = writeln!(counts, "__not_aligned\t{not_aligned}");
    }

    let grouped = spec.group_tag.as_ref().map(|tag| {
        let mut s = String::new();
        if !chrs.is_empty() {
            s.push_str("feature_id\tgroup_id\tcount\n");
        }
        for i in fragment_order(p, &OutputSpec::grouped_suffix(tag), &names) {
            chrs[i].counts.write_grouped_rows(|g| group_name(pool, g), &mut s);
        }
        s
    });
    FinalTexts { gtf, counts, grouped }
}

/// Writes the three output files.
pub fn write_outputs(spec: &OutputSpec, chrs: &[ChrOutput], pool: &GroupPool, not_aligned: u64) -> io::Result<()> {
    std::fs::create_dir_all(&spec.dir)?;
    let t = assemble(spec, chrs, pool, not_aligned);
    let mut w = open_out(&spec.gtf_path(), spec.gzip)?;
    w.write_all(&t.gtf)?;
    w.flush()?;
    drop(w);
    let mut w = open_out(&spec.counts_path(), spec.gzip)?;
    w.write_all(t.counts.as_bytes())?;
    w.flush()?;
    drop(w);
    if let (Some(path), Some(g)) = (spec.grouped_path(), &t.grouped) {
        let gz = spec.gzip || path.extension().is_some_and(|e| e == "gz");
        let mut w = open_out(&path, gz)?;
        w.write_all(g.as_bytes())?;
        w.flush()?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cmp::Ordering::*;

    #[test]
    fn fixed2_half_even_on_exact_ties() {
        // exact binary ties: Python '%.2f' rounds half to even
        assert_eq!(fmt_fixed2(0.125), "0.12");
        assert_eq!(fmt_fixed2(0.375), "0.38");
        assert_eq!(fmt_fixed2(2.5), "2.50");
        assert_eq!(fmt_fixed2(1.0 / 3.0), "0.33");
        assert_eq!(fmt_fixed2(2.675), "2.67"); // binary value below the tie
        assert_eq!(fmt_fixed2(0.0), "0.00");
        assert_eq!(fmt_fixed2(1e17), "100000000000000000.00");
        for n in 1..=16u32 {
            for k in 0..=(4 * n) {
                let x = k as f64 / n as f64;
                let s = fmt_fixed2(x);
                // reference: exact rational rounding half-even of k/n*100
                let num = k as u64 * 100;
                let (q, r) = (num / n as u64, num % n as u64);
                let q = if 2 * r > n as u64 || (2 * r == n as u64 && q % 2 == 1) { q + 1 } else { q };
                // only valid where k/n is exactly representable or not a tie; ties need exactness
                if 2 * r == n as u64 && !(n.is_power_of_two()) {
                    continue;
                }
                assert_eq!(s, format!("{}.{:02}", q / 100, q % 100), "{k}/{n}");
            }
        }
    }

    #[test]
    fn natural_order() {
        assert_eq!(natural_key_cmp("OUT_2.x", "OUT_10.x"), Less);
        assert_eq!(natural_key_cmp("OUT_X.x", "OUT_10.x"), Greater);
        assert_eq!(natural_key_cmp("OUT_chrM.x", "OUT_chrm.x"), Equal);
        let chrs = ["10", "2", "MT", "1", "chrUn/1", "KI270.1"];
        let o = fragment_order("OUT", ".transcript_models.gtf", &chrs);
        let got: Vec<&str> = o.iter().map(|&i| chrs[i]).collect();
        assert_eq!(got, ["1", "2", "10", "chrUn/1", "KI270.1", "MT"]);
        // the suffix takes part in the comparison
        let chrs = ["a", "a.t"];
        let g: Vec<&str> = fragment_order("OUT", ".transcript_models.gtf", &chrs).iter().map(|&i| chrs[i]).collect();
        assert_eq!(g, ["a.t", "a"]);
        let c: Vec<&str> =
            fragment_order("OUT", ".discovered_transcript_counts.tsv", &chrs).iter().map(|&i| chrs[i]).collect();
        assert_eq!(c, ["a", "a.t"]);
    }

    #[test]
    fn counts_rows() {
        let mut c = ChrCounts::new();
        c.entry("transcript:ENS1").total = 3.0;
        c.entry("transcript:ENS1").in_all = true;
        c.entry("transcript:ENS1").groups = Groups::from_vec(vec![(2, 0.0), (1, 3.0)]);
        let e = c.entry("transcript12.8.nic");
        e.total = 1.0;
        e.in_all = true;
        e.groups = Groups::from_vec(vec![(1, 1.0)]);
        c.entry("transcript9.8.nic").groups.push((1, 0.0)); // touched, not in all_features
        let mut s = String::new();
        c.write_ungrouped_rows(&mut s);
        assert_eq!(s, "transcript12.8.nic\t1.00\ntranscript:ENS1\t3.00\n");
        let mut pool = GroupPool::new();
        pool.intern("A");
        pool.intern("B");
        let mut g = String::new();
        c.write_grouped_rows(|x| group_name(&pool, x), &mut g);
        assert_eq!(g, "transcript12.8.nic\tA\t1.00\ntranscript:ENS1\tA\t3.00\n");
    }
}

#[cfg(test)]
mod phase_w_tests {
    use super::*;
    use iq_annot::{index_from_genes, parse_gff3_reader, GeneTypes};
    use iq_types::smallvec::smallvec;
    use iq_types::{AssignmentType, FinalAssign, GeneKey, Iv, ModelKind, Strand};

    const GFF: &str = "##gff-version 3\n\
8\tens\tgene\t100\t900\t.\t+\t.\tID=gene:G1;Name=g1;biotype=pc;gene_id=G1\n\
8\tens\tmRNA\t100\t900\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;Name=g1-201;transcript_id=T1\n\
8\tens\texon\t100\t200\t.\t+\t.\tParent=transcript:T1;exon_id=E1\n\
8\tens\tCDS\t150\t200\t.\t+\t0\tParent=transcript:T1\n\
8\tens\texon\t700\t900\t.\t+\t.\tParent=transcript:T1;exon_id=E2\n\
8\tens\tgene\t5000\t6000\t.\t-\t.\tID=gene:G2;Name=g2\n\
8\tens\tmRNA\t5000\t6000\t.\t-\t.\tID=transcript3.8.nic;Parent=gene:G2\n\
8\tens\texon\t5000\t6000\t.\t-\t.\tParent=transcript3.8.nic\n";

    fn m(key: u32, kind: ModelKind, gene: GeneKey, strand: Strand, ex: &[(i32, i32)]) -> TranscriptModel {
        TranscriptModel {
            key,
            strand,
            gene,
            kind,
            exons: ex.iter().map(|&(s, e)| Iv::new(s as _, e as _)).collect(),
            intron_path: Box::new([]),
            polya_confirmed: false,
            tss_confirmed: false,
            source_known: None,
        }
    }

    #[test]
    fn production_block() {
        let ann = parse_gff3_reader(GFF.as_bytes(), GeneTypes::Isoquant).unwrap();
        let chr = ann.chrom("8").unwrap();
        let g1 = ann.chrom_genes(chr).start;
        let t1 = ann.gene_txs(g1).start;
        let index = index_from_genes(&ann, &[g1], 6);
        let models = vec![
            m(0, ModelKind::Known(t1), GeneKey::Ref(g1), Strand::Plus, &[(100, 200), (700, 900)]),
            m(3, ModelKind::Nic, GeneKey::Ref(g1), Strand::Plus, &[(50, 200), (700, 950)]),
            m(5, ModelKind::Nnic, GeneKey::Novel(9), Strand::Minus, &[(2000, 2100), (2300, 2400)]),
        ];
        let mut c = TxCounter::new();
        c.add(&FinalAssign { atype: AssignmentType::Unique, gene_atype: AssignmentType::Unique, tx: smallvec![3] }, 1);
        c.add(&FinalAssign { atype: AssignmentType::Unique, gene_atype: AssignmentType::Unique, tx: smallvec![0] }, 2);
        c.add(
            &FinalAssign { atype: AssignmentType::Ambiguous, gene_atype: AssignmentType::Unique, tx: smallvec![0, 3] },
            1,
        );
        let mut w = ChrPhaseW::new(&ann, "8");
        let names = w.add_block(&BlockOutput { index: &index, models: &models, counter: &c });
        // forbidden: transcript3 from the annotation
        assert_eq!(names.tx, ["transcript:T1", "transcript1.8.nic", "transcript2.8.nnic"]);
        assert_eq!(names.gene, ["gene:G1", "gene:G1", "novel_gene_8_4"]);
        let out = w.finish();
        let want = "8\tens\tgene\t50\t950\t.\t+\t.\tgene_id \"gene:G1\"; transcripts \"2\"; Name \"g1\"; biotype \"pc\"; \n\
8\tens\ttranscript\t100\t900\t.\t+\t.\tgene_id \"gene:G1\"; transcript_id \"transcript:T1\"; exons \"2\"; Name \"g1-201\"; \n\
8\tens\texon\t100\t200\t.\t+\t.\tgene_id \"gene:G1\"; transcript_id \"transcript:T1\"; exon_number \"1\"; exon_id \"E1\";  Name \"g1-201\"; \n\
8\tens\tCDS\t150\t200\t.\t+\t.\tgene_id \"gene:G1\"; transcript_id \"transcript:T1\"; exon_number \"2\"; exon_id \"1\"; \n\
8\tens\texon\t700\t900\t.\t+\t.\tgene_id \"gene:G1\"; transcript_id \"transcript:T1\"; exon_number \"3\"; exon_id \"E2\";  Name \"g1-201\"; \n\
8\tIsoQuant\ttranscript\t50\t950\t.\t+\t.\tgene_id \"gene:G1\"; transcript_id \"transcript1.8.nic\"; exons \"2\";\n\
8\tIsoQuant\texon\t50\t200\t.\t+\t.\tgene_id \"gene:G1\"; transcript_id \"transcript1.8.nic\"; exon_number \"1\"; exon_id \"2\"; \n\
8\tIsoQuant\texon\t700\t950\t.\t+\t.\tgene_id \"gene:G1\"; transcript_id \"transcript1.8.nic\"; exon_number \"2\"; exon_id \"3\"; \n\
8\tIsoQuant\tgene\t2000\t2400\t.\t-\t.\tgene_id \"novel_gene_8_4\"; transcripts \"1\"; \n\
8\tIsoQuant\ttranscript\t2000\t2400\t.\t-\t.\tgene_id \"novel_gene_8_4\"; transcript_id \"transcript2.8.nnic\"; exons \"2\";\n\
8\tIsoQuant\texon\t2300\t2400\t.\t-\t.\tgene_id \"novel_gene_8_4\"; transcript_id \"transcript2.8.nnic\"; exon_number \"1\"; exon_id \"4\"; \n\
8\tIsoQuant\texon\t2000\t2100\t.\t-\t.\tgene_id \"novel_gene_8_4\"; transcript_id \"transcript2.8.nnic\"; exon_number \"2\"; exon_id \"5\"; \n";
        assert_eq!(out.gtf, want);
        let mut pool = GroupPool::new();
        pool.intern("A");
        pool.intern("B");
        let spec = OutputSpec {
            dir: PathBuf::new(),
            prefix: "OUT".into(),
            cmd_line: "iqrs x".into(),
            version: "4.1.0-iqrs-0.1.0".into(),
            group_tag: Some("HP".into()),
            grouped_out: None,
            gzip: false,
        };
        let t = assemble(&spec, &[out], &pool, 7);
        assert!(String::from_utf8(t.gtf)
            .unwrap()
            .starts_with("# OUT IsoQuant generated GTF\n# Command line: iqrs x\n"));
        assert_eq!(
            t.counts,
            "feature_id\tcount\ntranscript1.8.nic\t1.00\ntranscript:T1\t1.00\n__ambiguous\t1\n__no_feature\t0\n__not_aligned\t7\n"
        );
        assert_eq!(
            t.grouped.unwrap(),
            "feature_id\tgroup_id\tcount\ntranscript1.8.nic\tA\t1.00\ntranscript:T1\tB\t1.00\n"
        );
    }

    #[test]
    fn write_files_plain_and_gzip() {
        use std::io::Read;
        let dir = std::env::temp_dir().join(format!("iq_quant_test_{}", std::process::id()));
        let mut c = ChrCounts::new();
        let e = c.entry("transcript1.8.nic");
        e.total = 2.0;
        e.in_all = true;
        e.groups.push((1, 2.0));
        let chrs = [ChrOutput { chr: "8".into(), gtf: "x\n".into(), counts: c }];
        let mut pool = GroupPool::new();
        pool.intern("A");
        for gz in [false, true] {
            let spec = OutputSpec {
                dir: dir.clone(),
                prefix: "OUT".into(),
                cmd_line: "c".into(),
                version: "v".into(),
                group_tag: Some("HP".into()),
                grouped_out: None,
                gzip: gz,
            };
            write_outputs(&spec, &chrs, &pool, 0).unwrap();
            let read = |p: PathBuf| -> String {
                let mut s = String::new();
                if gz {
                    flate2::read::GzDecoder::new(File::open(p).unwrap()).read_to_string(&mut s).unwrap();
                } else {
                    File::open(p).unwrap().read_to_string(&mut s).unwrap();
                }
                s
            };
            assert!(read(spec.gtf_path()).ends_with("# IsoQuant version: v\nx\n"));
            assert_eq!(
                read(spec.counts_path()),
                "feature_id\tcount\ntranscript1.8.nic\t2.00\n__ambiguous\t0\n__no_feature\t0\n__not_aligned\t0\n"
            );
            let gp = spec.grouped_path().unwrap();
            assert!(gp.to_str().unwrap().contains("OUT.discovered_transcript_grouped_tag_HP_counts.linear.tsv"));
            assert_eq!(read(gp), "feature_id\tgroup_id\tcount\ntranscript1.8.nic\tA\t2.00\n");
        }
        let _ = std::fs::remove_dir_all(&dir);
    }
}
