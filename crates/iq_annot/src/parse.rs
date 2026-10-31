//! GFF3 / GTF streaming parsers.
//!
//! One pass over the file keeps two kinds of rows:
//! - *leaf* rows (exon, CDS, codons, UTRs): parent row, interval, type; exons with an
//!   `exon_id` also feed the exon-id table;
//! - *feature* rows (everything else that has an `ID` or a `Parent`): coordinates, ids,
//!   the raw attribute column.
//!
//! [`finish`] then resolves the hierarchy, selects genes and transcripts ([`GeneTypes`]),
//! applies the CDS-as-exon rule and lays out the flat [`AnnotationIndex`].

use crate::attrs;
use crate::index::{AnnotationIndex, GeneRow, GeneTypes, StrRef, TxRow};
use iq_types::{ChrIdx, GeneIdx, Iv, Pos, Strand};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Annotation parse error (malformed line, duplicate gene/transcript ID, ...).
#[derive(Debug)]
pub struct AnnotError(pub String);

impl std::fmt::Display for AnnotError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for AnnotError {}

/// 4.1 `GeneInfo.DEFAULT_RNA_TYPES` (`gene_info.py:230-231`).
pub const DEFAULT_RNA_TYPES: &[&str] = &[
    "transcript",
    "RNA",
    "mRNA",
    "miRNA",
    "ncRNA",
    "tRNA",
    "rRNA",
    "snRNA",
    "lnc_RNA",
    "snoRNA",
    "telomerase_RNA",
    "Y_RNA",
    "scRNA",
];
/// 4.1 `GeneInfo.OTHER_FEATURES` (`gene_info.py:229`; exact type match, so Ensembl's
/// `five_prime_UTR`/`three_prime_UTR` are not included).
pub const OTHER_FEATURES: &[&str] = &["CDS", "start_codon", "stop_codon", "UTR"];
/// Types kept as leaf rows (never parents).
const LEAF_TYPES: &[&str] = &["exon", "CDS", "start_codon", "stop_codon", "UTR", "five_prime_UTR", "three_prime_UTR"];

const NONE: u32 = u32::MAX;

/// `fn parse_gff3(path, GeneTypes) -> AnnotationIndex` (U1).
pub fn parse_gff3(path: &Path, gene_types: GeneTypes) -> Result<AnnotationIndex, AnnotError> {
    let f = File::open(path).map_err(|e| AnnotError(format!("{}: {e}", path.display())))?;
    parse_gff3_reader(BufReader::with_capacity(1 << 20, f), gene_types)
}

/// GFF3 from any reader.
pub fn parse_gff3_reader<R: BufRead>(r: R, gene_types: GeneTypes) -> Result<AnnotationIndex, AnnotError> {
    let mut raw = Raw::default();
    raw.read(r, Format::Gff3)?;
    raw.finish(gene_types)
}

/// GTF input (parity is only claimed for GFF3).
pub fn parse_gtf(path: &Path, gene_types: GeneTypes) -> Result<AnnotationIndex, AnnotError> {
    let f = File::open(path).map_err(|e| AnnotError(format!("{}: {e}", path.display())))?;
    parse_gtf_reader(BufReader::with_capacity(1 << 20, f), gene_types)
}

/// GTF from any reader.
pub fn parse_gtf_reader<R: BufRead>(r: R, gene_types: GeneTypes) -> Result<AnnotationIndex, AnnotError> {
    let mut raw = Raw::default();
    raw.read(r, Format::Gtf)?;
    raw.finish(gene_types)
}

#[derive(Copy, Clone, PartialEq, Eq)]
enum Format {
    Gff3,
    Gtf,
}

#[derive(Clone, Debug)]
struct Feat {
    chr: ChrIdx,
    iv: Iv,
    strand: Strand,
    ty: u16,
    source: u16,
    /// Empty (`len == 0`) for rows without `ID`.
    id: StrRef,
    attrs: StrRef,
    has_parent: bool,
    /// gffutils primary key (ordering tie-breaks).
    db_id: DbId,
}

#[derive(Copy, Clone, Debug)]
struct Leaf {
    /// Feature row, or `NONE` while unresolved / unknown.
    parent: u32,
    iv: Iv,
    ty: u16,
    db_id: DbId,
}

/// The primary key gffutils gives a row (`create.py:164-237, 375-377`): the `ID`
/// attribute; `<ID>_<k>` for the k-th duplicate (`merge_strategy='create_unique'`);
/// `<featuretype>_<n>` (n-th ID-less row of that type) otherwise. SQLite returns a
/// parent's children in `(parent, child)` primary-key order, i.e. sorted by this string,
/// which makes it the tie-break of every `children(...)` query 4.1 issues.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum DbId {
    Attr(StrRef),
    Dup(StrRef, u32),
    Auto(u16, u32),
}

struct IdEntry {
    /// Feature row holding the ID (`NONE` if a leaf row holds it).
    feat: u32,
    base: StrRef,
    n_dup: u32,
}

#[derive(Default)]
struct Interner {
    names: Vec<String>,
    map: HashMap<String, u16>,
}

impl Interner {
    fn get(&mut self, s: &str) -> Result<u16, AnnotError> {
        if let Some(&i) = self.map.get(s) {
            return Ok(i);
        }
        let i = u16::try_from(self.names.len())
            .map_err(|_| AnnotError("too many distinct types/sources/chromosomes".into()))?;
        self.names.push(s.to_owned());
        self.map.insert(s.to_owned(), i);
        Ok(i)
    }
}

#[derive(Default)]
struct Raw {
    arena: String,
    chroms: Interner,
    seq_region: HashMap<String, u64>,
    types: Interner,
    sources: Interner,
    feats: Vec<Feat>,
    id_map: HashMap<String, IdEntry>,
    /// `_autoincrements` per feature type.
    auto: HashMap<u16, u32>,
    /// Feature rows whose ID was already taken.
    dup_ids: Vec<u32>,
    /// `(child feature row, parent feature row)`.
    feat_links: Vec<(u32, u32)>,
    leaves: Vec<Leaf>,
    /// Unresolved parents: `(is_leaf, row, parent id)`.
    pending: Vec<(bool, u32, String)>,
    /// `(chr, s, e, strand)` → exon_id. **Last exon row in file order wins**: 4.1
    /// `FeatureIdStorage.__init__` (`utils/id_policy.py:56-66`) assigns into a dict while
    /// iterating `genedb.region(seqid, start=1, featuretype="exon")`, which SQLite returns in
    /// rowid (= file) order. Verified on GRCz11: 12,177 keys have differing exon_ids; 4.1
    /// keeps the last in all of them (320,082/320,082 entries equal).
    exon_ids: HashMap<(ChrIdx, Pos, Pos, Strand), StrRef>,
    /// GTF: synthetic gene/transcript rows by id.
    gtf_genes: HashMap<String, u32>,
    gtf_txs: HashMap<String, u32>,
}

fn strand_of(s: &str) -> Strand {
    match s {
        "+" => Strand::Plus,
        "-" => Strand::Minus,
        _ => Strand::Unknown,
    }
}

impl Raw {
    fn push_str(&mut self, s: &str) -> StrRef {
        let off = self.arena.len() as u32;
        self.arena.push_str(s);
        StrRef { off, len: s.len() as u32 }
    }
    fn s(&self, r: StrRef) -> &str {
        &self.arena[r.off as usize..(r.off + r.len) as usize]
    }

    /// Feature row a `Parent` value refers to.
    fn parent_row(&self, p: &str) -> Option<u32> {
        self.id_map.get(p).map(|e| e.feat).filter(|&f| f != NONE)
    }

    /// Assigns the gffutils primary key of a new row (`feat` = its feature row, if any).
    fn db_id(&mut self, id: Option<&str>, ty: u16, feat: u32) -> DbId {
        match id {
            None => {
                let n = self.auto.entry(ty).or_insert(0);
                *n += 1;
                DbId::Auto(ty, *n)
            }
            Some(id) => {
                if let Some(e) = self.id_map.get_mut(id) {
                    e.n_dup += 1;
                    if feat != NONE {
                        self.dup_ids.push(feat);
                    }
                    return DbId::Dup(e.base, e.n_dup);
                }
                let base = self.push_str(id);
                self.id_map.insert(id.to_owned(), IdEntry { feat, base, n_dup: 0 });
                DbId::Attr(base)
            }
        }
    }

    fn render(&self, d: DbId, types: &[String]) -> String {
        match d {
            DbId::Attr(r) => self.s(r).to_owned(),
            DbId::Dup(r, k) => format!("{}_{k}", self.s(r)),
            DbId::Auto(t, n) => format!("{}_{n}", types[t as usize]),
        }
    }

    fn read<R: BufRead>(&mut self, mut r: R, fmt: Format) -> Result<(), AnnotError> {
        let mut buf = Vec::with_capacity(1 << 12);
        let mut line_no = 0usize;
        loop {
            buf.clear();
            let n = r.read_until(b'\n', &mut buf).map_err(|e| AnnotError(e.to_string()))?;
            if n == 0 {
                break;
            }
            line_no += 1;
            while matches!(buf.last(), Some(b'\n' | b'\r')) {
                buf.pop();
            }
            if buf.is_empty() {
                continue;
            }
            let line = std::str::from_utf8(&buf).map_err(|e| AnnotError(format!("line {line_no}: {e}")))?;
            if line.starts_with('#') {
                if line.starts_with("##FASTA") {
                    break;
                }
                if let Some(rest) = line.strip_prefix("##sequence-region") {
                    let f: Vec<&str> = rest.split_whitespace().collect();
                    if let [name, _, end] = f[..] {
                        if let Ok(e) = end.parse() {
                            self.seq_region.insert(name.to_owned(), e);
                        }
                    }
                }
                continue;
            }
            self.row(line, fmt).map_err(|e| AnnotError(format!("line {line_no}: {}", e.0)))?;
        }
        Ok(())
    }

    fn row(&mut self, line: &str, fmt: Format) -> Result<(), AnnotError> {
        let mut it = line.split('\t');
        let mut f: [&str; 9] = [""; 9];
        for (k, slot) in f.iter_mut().enumerate() {
            *slot = it.next().ok_or_else(|| AnnotError(format!("expected 9 tab-separated columns, got {k}")))?;
        }
        let [seqid, source, ty, start, end, _score, strand, _phase, attr] = f;
        let start: Pos = start.trim().parse().map_err(|_| AnnotError(format!("bad start {start:?}")))?;
        let end: Pos = end.trim().parse().map_err(|_| AnnotError(format!("bad end {end:?}")))?;
        let iv = Iv::new(start, end);
        let strand = strand_of(strand);
        match fmt {
            Format::Gff3 => self.gff3_row(seqid, source, ty, iv, strand, attr),
            Format::Gtf => self.gtf_row(seqid, source, ty, iv, strand, attr),
        }
    }

    fn gff3_row(
        &mut self,
        seqid: &str,
        source: &str,
        ty: &str,
        iv: Iv,
        strand: Strand,
        attr: &str,
    ) -> Result<(), AnnotError> {
        let parents = attrs::attr_values(attr, "Parent");
        let id = attrs::attr_first(attr, "ID");
        let tyi = self.types.get(ty)?;
        if LEAF_TYPES.contains(&ty) {
            let db_id = self.db_id(id.as_deref(), tyi, NONE);
            if parents.is_empty() || !(ty == "exon" || OTHER_FEATURES.contains(&ty)) {
                // orphans cannot belong to a transcript; UTR rows are never read (their
                // primary-key counter still advances above)
                return Ok(());
            }
            let chr = self.chroms.get(seqid)?;
            if ty == "exon" {
                if let Some(eid) = attrs::attr_first(attr, "exon_id") {
                    let r = self.push_str(&eid);
                    self.exon_ids.insert((chr, iv.s, iv.e, strand), r);
                }
            }
            for p in parents {
                let row = self.leaves.len() as u32;
                let parent = match self.parent_row(&p) {
                    Some(x) => x,
                    None => {
                        self.pending.push((true, row, p));
                        NONE
                    }
                };
                self.leaves.push(Leaf { parent, iv, ty: tyi, db_id });
            }
            return Ok(());
        }
        if id.is_none() && parents.is_empty() {
            self.db_id(None, tyi, NONE);
            return Ok(()); // cannot take part in the hierarchy (e.g. biological_region)
        }
        let chr = self.chroms.get(seqid)?;
        let row = self.feats.len() as u32;
        let db_id = self.db_id(id.as_deref(), tyi, row);
        let idr = match db_id {
            DbId::Attr(r) | DbId::Dup(r, _) => r,
            DbId::Auto(..) => StrRef::default(),
        };
        let attrs = self.push_str(attr);
        let feat = Feat {
            chr,
            iv,
            strand,
            ty: tyi,
            source: self.sources.get(source)?,
            id: idr,
            attrs,
            has_parent: !parents.is_empty(),
            db_id,
        };
        self.feats.push(feat);
        for p in parents {
            match self.parent_row(&p) {
                Some(x) => self.feat_links.push((row, x)),
                None => self.pending.push((false, row, p)),
            }
        }
        Ok(())
    }

    /// GTF row: `gene_id`/`transcript_id` attributes; gene and transcript rows are
    /// used when present, otherwise synthesized from their children's spans.
    fn gtf_row(
        &mut self,
        seqid: &str,
        source: &str,
        ty: &str,
        iv: Iv,
        strand: Strand,
        attr: &str,
    ) -> Result<(), AnnotError> {
        let kv = gtf_attributes(attr);
        let get = |k: &str| kv.iter().find(|(kk, _)| kk == k).map(|(_, v)| v.clone());
        let gid = get("gene_id");
        let tid = get("transcript_id");
        let chr = self.chroms.get(seqid)?;
        let src = self.sources.get(source)?;
        // GTF attributes are re-expressed as a GFF3 column so the shared printing path works.
        let as_gff3 = |extra: &[(&str, &str)]| -> String {
            let mut s = String::new();
            for (k, v) in extra.iter().copied().chain(kv.iter().map(|(k, v)| (k.as_str(), v.as_str()))) {
                if !s.is_empty() {
                    s.push(';');
                }
                s.push_str(k);
                s.push('=');
                s.push_str(&gff3_escape(v));
            }
            s
        };
        let Some(gid) = gid else { return Ok(()) };
        let gene_row = |this: &mut Raw, explicit: bool| -> Result<u32, AnnotError> {
            if let Some(&g) = this.gtf_genes.get(&gid) {
                if explicit {
                    let a = as_gff3(&[("ID", &gid)]);
                    let ar = this.push_str(&a);
                    let tyi = this.types.get(ty)?;
                    let f = &mut this.feats[g as usize];
                    f.iv = iv;
                    f.strand = strand;
                    f.ty = tyi;
                    f.source = src;
                    f.attrs = ar;
                } else {
                    let f = &mut this.feats[g as usize];
                    f.iv = Iv::new(f.iv.s.min(iv.s), f.iv.e.max(iv.e));
                }
                return Ok(g);
            }
            let row = this.feats.len() as u32;
            let idr = this.push_str(&gid);
            let a = if explicit { as_gff3(&[("ID", &gid)]) } else { format!("ID={}", gff3_escape(&gid)) };
            let ar = this.push_str(&a);
            let tyi = this.types.get("gene")?;
            this.feats.push(Feat {
                chr,
                iv,
                strand,
                ty: tyi,
                source: src,
                id: idr,
                attrs: ar,
                has_parent: false,
                db_id: DbId::Attr(idr),
            });
            this.gtf_genes.insert(gid.clone(), row);
            Ok(row)
        };
        if ty == "gene" {
            gene_row(self, true)?;
            return Ok(());
        }
        let Some(tid) = tid else { return Ok(()) };
        let g = gene_row(self, false)?;
        let t = match self.gtf_txs.get(&tid) {
            Some(&t) => {
                if ty == "transcript" {
                    let a = as_gff3(&[("ID", &tid), ("Parent", &gid)]);
                    let ar = self.push_str(&a);
                    let f = &mut self.feats[t as usize];
                    f.iv = iv;
                    f.strand = strand;
                    f.attrs = ar;
                } else {
                    let f = &mut self.feats[t as usize];
                    f.iv = Iv::new(f.iv.s.min(iv.s), f.iv.e.max(iv.e));
                }
                t
            }
            None => {
                let row = self.feats.len() as u32;
                let idr = self.push_str(&tid);
                let a = if ty == "transcript" {
                    as_gff3(&[("ID", &tid), ("Parent", &gid)])
                } else {
                    format!("ID={};Parent={}", gff3_escape(&tid), gff3_escape(&gid))
                };
                let ar = self.push_str(&a);
                let tyi = self.types.get("transcript")?;
                self.feats.push(Feat {
                    chr,
                    iv,
                    strand,
                    ty: tyi,
                    source: src,
                    id: idr,
                    attrs: ar,
                    has_parent: true,
                    db_id: DbId::Attr(idr),
                });
                self.gtf_txs.insert(tid.clone(), row);
                self.feat_links.push((row, g));
                row
            }
        };
        if ty != "transcript" && LEAF_TYPES.contains(&ty) {
            let tyi = self.types.get(ty)?;
            if ty == "exon" {
                if let Some(eid) = get("exon_id") {
                    let r = self.push_str(&eid);
                    self.exon_ids.insert((chr, iv.s, iv.e, strand), r);
                }
            }
            let db_id = self.db_id(None, tyi, NONE);
            self.leaves.push(Leaf { parent: t, iv, ty: tyi, db_id });
        }
        Ok(())
    }

    fn finish(mut self, gene_types: GeneTypes) -> Result<AnnotationIndex, AnnotError> {
        // 1. resolve forward parent references (unknown parents are ignored, as gffutils
        //    keeps a relation to a feature that does not exist)
        for (is_leaf, row, p) in std::mem::take(&mut self.pending) {
            if let Some(x) = self.parent_row(&p) {
                if is_leaf {
                    self.leaves[row as usize].parent = x;
                } else {
                    self.feat_links.push((row, x));
                }
            }
        }
        let nf = self.feats.len();
        // 2. CSR children lists (file order)
        let feat_kids = csr(nf, self.feat_links.iter().map(|&(c, p)| (p, c)));
        let leaf_kids = csr(
            nf,
            self.leaves.iter().enumerate().filter(|(_, l)| l.parent != NONE).map(|(i, l)| (l.parent, i as u32)),
        );
        let ty_exon = self.types.map.get("exon").copied();
        let ty_cds = self.types.map.get("CDS").copied();
        let is_exon_cds = |t: u16| Some(t) == ty_exon || Some(t) == ty_cds;
        let has_exons = |f: u32| leaf_kids.get(f).iter().any(|&l| is_exon_cds(self.leaves[l as usize].ty));
        let other_tys: Vec<u16> = OTHER_FEATURES.iter().filter_map(|t| self.types.map.get(*t).copied()).collect();
        let ty_gene = self.types.map.get("gene").copied();
        let rna_tys: Vec<u16> = DEFAULT_RNA_TYPES.iter().filter_map(|t| self.types.map.get(*t).copied()).collect();

        // 3. genes and their transcripts (feature rows, unordered = file order)
        let mut sel: Vec<(u32, Vec<u32>)> = Vec::new();
        for g in 0..nf as u32 {
            let f = &self.feats[g as usize];
            let txs: Vec<u32> = match gene_types {
                GeneTypes::Isoquant => {
                    if Some(f.ty) != ty_gene {
                        continue;
                    }
                    // gffutils children(level=None): descendants at any depth, DISTINCT
                    let mut seen = vec![g];
                    let mut stack: Vec<u32> = feat_kids.get(g).to_vec();
                    let mut out = Vec::new();
                    while let Some(c) = stack.pop() {
                        if seen.contains(&c) {
                            continue;
                        }
                        seen.push(c);
                        if rna_tys.contains(&self.feats[c as usize].ty) {
                            out.push(c);
                        }
                        stack.extend_from_slice(feat_kids.get(c));
                    }
                    out.sort_unstable();
                    out
                }
                GeneTypes::All => {
                    if f.has_parent {
                        continue;
                    }
                    let mut out: Vec<u32> = feat_kids.get(g).iter().copied().filter(|&c| has_exons(c)).collect();
                    if out.is_empty() {
                        continue;
                    }
                    out.dedup();
                    out
                }
            };
            sel.push((g, txs));
        }
        // duplicate IDs on selected gene/transcript rows are fatal (gtf2db.py:285-329)
        if !self.dup_ids.is_empty() {
            let mut used = vec![false; nf];
            for (g, txs) in &sel {
                used[*g as usize] = true;
                for &t in txs {
                    used[t as usize] = true;
                }
            }
            for &d in &self.dup_ids {
                let id = self.s(self.feats[d as usize].id);
                let first = self.id_map[id].feat;
                let is_gt = |r: u32| {
                    used[r as usize]
                        || Some(self.feats[r as usize].ty) == ty_gene
                        || rna_tys.contains(&self.feats[r as usize].ty)
                };
                if is_gt(d) || is_gt(first) {
                    return Err(AnnotError(format!("duplicate gene/transcript ID {id:?} in annotation")));
                }
            }
        }
        sel.sort_by_key(|(g, _)| {
            let f = &self.feats[*g as usize];
            (f.chr, f.iv.s, *g)
        });

        // 4. lay out
        let mut idx = AnnotationIndex { gene_types, ..Default::default() };
        let mut genes = Vec::with_capacity(sel.len());
        let mut txs_out: Vec<TxRow> = Vec::new();
        let mut exons: Vec<Iv> = Vec::new();
        let mut other: Vec<(Iv, u16)> = Vec::new();
        let mut feats_buf: Vec<(Iv, bool, u32)> = Vec::new();
        let mut kids_buf: Vec<(String, u32)> = Vec::new();
        for (gi, (g, txs)) in sel.iter().enumerate() {
            let gf = &self.feats[*g as usize];
            let tx_off = txs_out.len() as u32;
            // unordered children(): primary-key order; order_by='start': (start, key)
            let mut keyed: Vec<(String, u32)> =
                txs.iter().map(|&t| (self.render(self.feats[t as usize].db_id, &self.types.names), t)).collect();
            keyed.sort();
            let mut ordered: Vec<(u32, &str, u32)> =
                keyed.iter().enumerate().map(|(k, (key, t))| (k as u32, key.as_str(), *t)).collect();
            ordered.sort_by_key(|&(_, key, t)| (self.feats[t as usize].iv.s, key));
            for (rank, _, t) in ordered {
                let tf = &self.feats[t as usize];
                feats_buf.clear();
                kids_buf.clear();
                let exon_off = exons.len() as u32;
                let other_off = other.len() as u32;
                for &l in leaf_kids.get(t) {
                    let lf = self.leaves[l as usize];
                    if is_exon_cds(lf.ty) || other_tys.contains(&lf.ty) {
                        kids_buf.push((self.render(lf.db_id, &self.types.names), l));
                    }
                }
                kids_buf.sort();
                for (k, (_, l)) in kids_buf.iter().enumerate() {
                    let lf = self.leaves[*l as usize];
                    if is_exon_cds(lf.ty) {
                        feats_buf.push((lf.iv, Some(lf.ty) == ty_exon, k as u32));
                    }
                    if other_tys.contains(&lf.ty) {
                        other.push((lf.iv, lf.ty));
                    }
                }
                feats_buf.sort_by_key(|&(iv, _, k)| (iv.s, k));
                cds_as_exon(&feats_buf, &mut exons);
                txs_out.push(TxRow {
                    gene: gi as GeneIdx,
                    iv: tf.iv,
                    strand: tf.strand,
                    ty: tf.ty,
                    source: tf.source,
                    id: tf.id,
                    attrs: tf.attrs,
                    exon_off,
                    n_exon: exons.len() as u32 - exon_off,
                    other_off,
                    n_other: other.len() as u32 - other_off,
                    file_rank: rank,
                });
            }
            genes.push(GeneRow {
                chr: gf.chr,
                iv: gf.iv,
                strand: gf.strand,
                ty: gf.ty,
                source: gf.source,
                id: gf.id,
                attrs: gf.attrs,
                tx_off,
                n_tx: txs_out.len() as u32 - tx_off,
            });
        }

        // release parse-time tables before the arena is compacted (peak memory)
        drop(leaf_kids);
        drop(feat_kids);
        self.leaves = Vec::new();
        self.id_map = HashMap::new();
        self.feat_links = Vec::new();
        drop(sel);

        // forbidden novel ids (id_policy.py:20-45)
        let nchr = self.chroms.names.len();
        let mut forbidden = vec![Vec::new(); nchr];
        let ty_tx = self.types.map.get("transcript").copied();
        let ty_mrna = self.types.map.get("mRNA").copied();
        for f in &self.feats {
            let id = self.s(f.id);
            if Some(f.ty) == ty_gene {
                if let Some(n) =
                    id.strip_prefix("novel_gene_").and_then(|_| py_int(id.rsplit('_').next().unwrap_or("")))
                {
                    forbidden[f.chr as usize].push(n);
                }
            }
            if (Some(f.ty) == ty_tx || Some(f.ty) == ty_mrna) && id.starts_with("transcript") {
                if let Some(n) = py_int(id.split('.').next().unwrap_or("").get(10..).unwrap_or("")) {
                    forbidden[f.chr as usize].push(n);
                }
            }
        }
        for v in &mut forbidden {
            v.sort_unstable();
            v.dedup();
        }

        // per-chromosome gene ranges + prefix-max ends
        let mut chr_genes = vec![(0u32, 0u32); nchr];
        let mut prefix = Vec::with_capacity(genes.len());
        let mut k = 0usize;
        for (c, slot) in chr_genes.iter_mut().enumerate() {
            let lo = k;
            let mut m = Pos::MIN;
            while k < genes.len() && genes[k].chr as usize == c {
                m = m.max(genes[k].iv.e);
                prefix.push(m);
                k += 1;
            }
            *slot = (lo as u32, k as u32);
        }

        // compact the arena: keep only strings the index references
        let mut arena = String::new();
        let mut keep = |r: StrRef, src: &str| -> StrRef {
            let off = arena.len() as u32;
            arena.push_str(&src[r.off as usize..(r.off + r.len) as usize]);
            StrRef { off, len: r.len }
        };
        for g in &mut genes {
            g.id = keep(g.id, &self.arena);
            g.attrs = keep(g.attrs, &self.arena);
        }
        for t in &mut txs_out {
            t.id = keep(t.id, &self.arena);
            t.attrs = keep(t.attrs, &self.arena);
        }
        let mut exon_ids = HashMap::with_capacity(self.exon_ids.len());
        for (k, r) in self.exon_ids.drain() {
            exon_ids.insert(k, keep(r, &self.arena));
        }

        idx.chrom_len = self.chroms.names.iter().map(|n| self.seq_region.get(n).copied()).collect();
        idx.chrom_map = self.chroms.map.iter().map(|(k, &v)| (k.clone(), v as ChrIdx)).collect();
        idx.chroms = self.chroms.names;
        idx.types = self.types.names;
        idx.sources = self.sources.names;
        arena.shrink_to_fit();
        idx.arena = arena;
        genes.shrink_to_fit();
        txs_out.shrink_to_fit();
        exons.shrink_to_fit();
        other.shrink_to_fit();
        idx.genes = genes;
        idx.txs = txs_out;
        idx.exons = exons;
        idx.other = other;
        idx.chr_genes = chr_genes;
        idx.prefix_max_end = prefix;
        idx.exon_ids = exon_ids;
        idx.forbidden = forbidden;
        Ok(idx)
    }
}

/// The CDS-as-exon rule (`set_introns_and_exons`, `gene_info.py:600-626`): `feats` are the
/// transcript's exon/CDS children ordered by `(start, gffutils key)`; an exon is always kept,
/// a CDS only if it overlaps neither neighbour in that list. KEEP (`bb4e729e`, #309).
pub fn cds_as_exon(feats: &[(Iv, bool, u32)], out: &mut Vec<Iv>) {
    for (i, &(iv, is_exon, _)) in feats.iter().enumerate() {
        if is_exon {
            out.push(iv);
            continue;
        }
        let prev = i.checked_sub(1).map(|j| feats[j].0);
        let next = feats.get(i + 1).map(|f| f.0);
        if prev.is_some_and(|p| iv.overlaps(&p)) || next.is_some_and(|n| iv.overlaps(&n)) {
            continue;
        }
        out.push(iv);
    }
}

/// Python `int(s)` for the digit strings `id_policy.py` parses (sign/whitespace allowed;
/// negative values cannot collide with the counter and are dropped).
fn py_int(s: &str) -> Option<u64> {
    let t = s.trim();
    let t = t.strip_prefix('+').unwrap_or(t);
    if t.is_empty()
        || !t.bytes().all(|b| b.is_ascii_digit() || b == b'_')
        || t.starts_with('_')
        || t.ends_with('_')
        || t.contains("__")
    {
        return None;
    }
    t.replace('_', "").parse().ok()
}

/// Compressed children lists: `get(p)` = children of `p` in insertion order.
struct Csr {
    off: Vec<u32>,
    items: Vec<u32>,
}

impl Csr {
    fn get(&self, p: u32) -> &[u32] {
        &self.items[self.off[p as usize] as usize..self.off[p as usize + 1] as usize]
    }
}

fn csr(n: usize, pairs: impl Iterator<Item = (u32, u32)> + Clone) -> Csr {
    let mut off = vec![0u32; n + 1];
    for (p, _) in pairs.clone() {
        off[p as usize + 1] += 1;
    }
    for i in 0..n {
        off[i + 1] += off[i];
    }
    let mut pos = off.clone();
    let mut items = vec![0u32; off[n] as usize];
    let mut pairs: Vec<(u32, u32)> = pairs.collect();
    // children in file order (row index) regardless of link creation order
    pairs.sort_by_key(|&(p, c)| (p, c));
    for (p, c) in pairs {
        items[pos[p as usize] as usize] = c;
        pos[p as usize] += 1;
    }
    Csr { off, items }
}

/// GTF attributes `key "value"; key value;` → pairs (quotes stripped).
fn gtf_attributes(raw: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    for part in raw.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        let (k, v) = match part.find(char::is_whitespace) {
            Some(i) => (&part[..i], part[i..].trim()),
            None => (part, ""),
        };
        let v = v.strip_prefix('"').and_then(|v| v.strip_suffix('"')).unwrap_or(v);
        out.push((k.to_owned(), v.to_owned()));
    }
    out
}

/// Escapes the characters that are structural in a GFF3 attribute column.
fn gff3_escape(v: &str) -> String {
    let mut s = String::with_capacity(v.len());
    for c in v.chars() {
        match c {
            ';' => s.push_str("%3B"),
            '=' => s.push_str("%3D"),
            ',' => s.push_str("%2C"),
            '%' => s.push_str("%25"),
            '&' => s.push_str("%26"),
            '\t' => s.push_str("%09"),
            c => s.push(c),
        }
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(s: Pos, e: Pos) -> Iv {
        Iv::new(s, e)
    }

    #[test]
    fn cds_rule() {
        // exon, CDS inside it (dropped), lone CDS (kept), CDS touching next exon (dropped)
        let f = [
            (iv(1, 100), true, 0),
            (iv(50, 100), false, 1),
            (iv(200, 300), false, 2),
            (iv(400, 450), false, 3),
            (iv(450, 500), true, 4),
        ];
        let mut out = Vec::new();
        cds_as_exon(&f, &mut out);
        assert_eq!(out, vec![iv(1, 100), iv(200, 300), iv(450, 500)]);
    }

    const GFF: &str = "##gff-version 3
##sequence-region   1 1 1000
1\t.\tbiological_region\t1\t5\t.\t.\t.\tlogic_name=cpg
1\tsrc\tgene\t100\t900\t.\t+\t.\tID=gene:G1;Name=g1;biotype=protein_coding;description=a%3Bb
1\tsrc\tmRNA\t300\t900\t.\t+\t.\tID=transcript:T2;Parent=gene:G1;Name=t2;tag=x,y
1\tsrc\texon\t300\t400\t.\t+\t.\tParent=transcript:T2;exon_id=E3
1\tsrc\tmRNA\t100\t800\t.\t+\t.\tID=transcript:T1;Parent=gene:G1
1\tsrc\texon\t100\t200\t.\t+\t.\tParent=transcript:T1;exon_id=E1
1\tsrc\tCDS\t150\t200\t.\t+\t0\tID=CDS:P1;Parent=transcript:T1
1\tsrc\tCDS\t500\t600\t.\t+\t0\tID=CDS:P1;Parent=transcript:T1
1\tsrc\texon\t700\t800\t.\t+\t.\tParent=transcript:T1;exon_id=E2
1\tsrc\tfive_prime_UTR\t100\t149\t.\t+\t.\tParent=transcript:T1
1\tsrc\tlnc_RNA\t850\t900\t.\t+\t.\tID=transcript:T3;Parent=gene:G1
1\tsrc\tncRNA_gene\t2000\t3000\t.\t-\t.\tID=gene:N1
1\tsrc\tprocessed_transcript\t2000\t3000\t.\t-\t.\tID=transcript:N1a;Parent=gene:N1
1\tsrc\texon\t2000\t3000\t.\t-\t.\tParent=transcript:N1a
1\tsrc\tmRNA\t2\t3\t.\t+\t.\tID=transcript:T2;Parent=gene:G1
";

    #[test]
    fn duplicate_transcript_id_is_fatal() {
        let r = parse_gff3_reader(GFF.as_bytes(), GeneTypes::Isoquant);
        assert!(r.is_err());
    }

    #[test]
    fn small_gff() {
        let gff: String = GFF.lines().filter(|l| !l.contains("\t2\t3\t")).map(|l| format!("{l}\n")).collect();
        let a = parse_gff3_reader(gff.as_bytes(), GeneTypes::Isoquant).unwrap();
        assert_eq!(a.n_genes(), 1);
        assert_eq!(a.chrom_len(0), Some(1000));
        let txs: Vec<&str> = a.gene_txs(0).map(|t| a.tx_id(t)).collect();
        assert_eq!(txs, vec!["transcript:T1", "transcript:T2", "transcript:T3"]);
        assert_eq!(a.tx_exons(0), &[iv(100, 200), iv(500, 600), iv(700, 800)]);
        assert!(a.tx_exons(2).is_empty());
        assert_eq!(a.tx_file_rank(0), 0);
        let of: Vec<(Iv, &str)> = a.tx_other_features(0).collect();
        assert_eq!(of, vec![(iv(150, 200), "CDS"), (iv(500, 600), "CDS")]);
        assert_eq!(a.gene_attrs(0), "Name \"g1\"; biotype \"protein_coding\"; description \"a;b\"; ");
        assert_eq!(a.tx_attrs(1), "Name \"t2\"; tag \"x\"; ");
        assert_eq!(a.exon_id(0, iv(700, 800), Strand::Plus), Some("E2"));
        let all = parse_gff3_reader(gff.as_bytes(), GeneTypes::All).unwrap();
        assert_eq!(all.n_genes(), 2);
        assert_eq!(all.gene_type(1), "ncRNA_gene");
        // T3 has no exons → not a transcript in `all` mode
        assert_eq!(all.gene_txs(0).len(), 2);
    }

    /// Tie-breaks follow the gffutils primary key: `CDS:P_10` < `CDS:P_2`, `exon_10` < `exon_9`.
    #[test]
    fn gffutils_key_order() {
        let mut gff = String::from("1\ts\tgene\t1\t5000\t.\t+\t.\tID=G\n1\ts\tmRNA\t1\t5000\t.\t+\t.\tID=T;Parent=G\n");
        // 11 CDS segments sharing an ID → CDS:P, CDS:P_1 .. CDS:P_10; also exon_1..exon_10 ids
        for k in 0..11 {
            let s = 100 * (k + 1);
            gff += &format!("1\ts\tCDS\t{s}\t{}\t.\t+\t0\tID=CDS:P;Parent=T\n", s + 10);
        }
        for k in 0..10 {
            let s = 2000 + 100 * k;
            gff += &format!("1\ts\texon\t{s}\t{}\t.\t+\t.\tParent=T\n", s + 10);
        }
        let a = parse_gff3_reader(gff.as_bytes(), GeneTypes::Isoquant).unwrap();
        let of: Vec<Pos> = a.tx_other_features(0).map(|(iv, _)| iv.s).collect();
        // CDS:P, CDS:P_1, CDS:P_10, CDS:P_2, ..
        assert_eq!(&of[..4], &[100, 200, 1100, 300]);
        assert_eq!(a.tx_exons(0).len(), 21);
    }

    #[test]
    fn gtf_minimal() {
        let gtf = "1\ts\texon\t10\t20\t.\t+\t.\tgene_id \"G\"; transcript_id \"T\"; exon_id \"E\";\n\
                   1\ts\texon\t30\t40\t.\t+\t.\tgene_id \"G\"; transcript_id \"T\";\n";
        let a = parse_gtf_reader(gtf.as_bytes(), GeneTypes::Isoquant).unwrap();
        assert_eq!(a.n_genes(), 1);
        assert_eq!(a.gene_span(0).1, iv(10, 40));
        assert_eq!(a.tx_exons(0), &[iv(10, 20), iv(30, 40)]);
        assert_eq!(a.exon_id(0, iv(10, 20), Strand::Plus), Some("E"));
    }

    #[test]
    fn py_int_rules() {
        assert_eq!(py_int(" 12 "), Some(12));
        assert_eq!(py_int("1_0"), Some(10));
        assert_eq!(py_int("x1"), None);
        assert_eq!(py_int(""), None);
    }
}
