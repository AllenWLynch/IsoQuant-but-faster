//! gffutils-compatible attribute handling (design §4.1).
//!
//! The pipeline's DB was built by gffutils 0.14 with the inferred Ensembl GFF3 dialect
//! (`field separator ';'`, `keyval separator '='`, no quoting, no repeated keys). gffutils
//! `_split_keyvals` (`parser.py:96-322`) then gives, per `key=value` part:
//! - the key is the text before the first `=` (the value keeps any further `=`);
//! - a value containing `", "` is kept whole, otherwise it is split on `,`;
//! - every value is URL-unquoted (`urllib.parse.unquote`); keys are not;
//! - a repeated key appends to the list at its first position.
//!
//! `sort_attribute_values=True` only affects gffutils' *printing* of a feature
//! (`_reconstruct`), not `Feature.attributes`, so value lists keep their file order.

/// Attributes skipped by `set_gene_attributes` for genes (`gene_info.py:579`).
pub const GENE_SKIP: &[&str] = &["gene_id", "ID", "level", "Parent"];
/// ... for transcripts (`gene_info.py:586`).
pub const TX_SKIP: &[&str] = &["transcript_id", "gene_id", "ID", "level", "exons", "Parent"];
/// ... for the per-exon keys (`gene_info.py:593`).
pub const EXON_SKIP: &[&str] = &["transcript_id", "gene_id", "ID", "Parent", "level", "exon_id", "exon", "exon_number"];

/// Iterates the raw `key=value` parts of a GFF3 attribute column (leading/trailing `;`
/// stripped, as the inferred dialect does for files that have them).
fn parts(raw: &str) -> impl Iterator<Item = (&str, &str)> {
    let raw = if raw == "." { "" } else { raw };
    let raw = raw.trim_end_matches(';').trim_start_matches(';');
    let empty = raw.is_empty();
    raw.split(';').filter(move |_| !empty).map(|p| match p.find('=') {
        Some(i) => (&p[..i], &p[i + 1..]),
        None => (p, ""),
    })
}

/// Splits one raw value the gffutils way (before unquoting) and pushes the unquoted pieces.
fn push_values(val: &str, out: &mut Vec<String>) {
    if val.is_empty() {
        return;
    }
    if val.contains(", ") {
        out.push(unquote(val));
    } else {
        out.extend(val.split(',').map(unquote));
    }
}

/// gffutils `Feature.attributes` for a GFF3 attribute column.
pub fn parse_gff3_attributes(raw: &str) -> Vec<(String, Vec<String>)> {
    let mut out: Vec<(String, Vec<String>)> = Vec::new();
    for (k, v) in parts(raw) {
        let pos = match out.iter().position(|(kk, _)| kk == k) {
            Some(p) => p,
            None => {
                out.push((k.to_owned(), Vec::new()));
                out.len() - 1
            }
        };
        push_values(v, &mut out[pos].1);
    }
    out
}

/// All values of `key` (gffutils semantics), without building the full map.
pub fn attr_values(raw: &str, key: &str) -> Vec<String> {
    let mut out = Vec::new();
    for (k, v) in parts(raw) {
        if k == key {
            push_values(v, &mut out);
        }
    }
    out
}

/// First value of `key` (gffutils `attributes[key][0]`), if any.
pub fn attr_first(raw: &str, key: &str) -> Option<String> {
    for (k, v) in parts(raw) {
        if k == key && !v.is_empty() {
            let first = if v.contains(", ") { v } else { v.split(',').next().unwrap_or("") };
            return Some(unquote(first));
        }
    }
    None
}

/// `'%s "%s"; ' % (attr, values[0])` for every attribute not in `skip` with a non-empty
/// value list (`set_gene_attributes`, `gene_info.py:576-596`).
pub fn format_attrs(attrs: &[(String, Vec<String>)], skip: &[&str]) -> String {
    let mut s = String::new();
    for (k, vals) in attrs {
        if skip.contains(&k.as_str()) {
            continue;
        }
        if let Some(v) = vals.first() {
            s.push_str(k);
            s.push_str(" \"");
            s.push_str(v);
            s.push_str("\"; ");
        }
    }
    s
}

/// `urllib.parse.unquote` (UTF-8, `errors='replace'`): valid `%XX` escapes become bytes,
/// invalid ones stay literal.
pub fn unquote(s: &str) -> String {
    if !s.contains('%') {
        return s.to_owned();
    }
    let b = s.as_bytes();
    let mut out = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'%' && i + 2 < b.len() {
            if let (Some(h), Some(l)) = (hex(b[i + 1]), hex(b[i + 2])) {
                out.push(h << 4 | l);
                i += 3;
                continue;
            }
        }
        out.push(b[i]);
        i += 1;
    }
    match String::from_utf8(out) {
        Ok(s) => s,
        Err(e) => String::from_utf8_lossy(e.as_bytes()).into_owned(),
    }
}

fn hex(c: u8) -> Option<u8> {
    match c {
        b'0'..=b'9' => Some(c - b'0'),
        b'a'..=b'f' => Some(c - b'a' + 10),
        b'A'..=b'F' => Some(c - b'A' + 10),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gffutils_split_rules() {
        let raw = "ID=gene:X;Alias=A,B;description=kinase, subunit 1 [Source:x%3BAcc:1];k=a=b;flag;ID=Y";
        let a = parse_gff3_attributes(raw);
        assert_eq!(a[0], ("ID".into(), vec!["gene:X".into(), "Y".into()]));
        assert_eq!(a[1].1, vec!["A".to_string(), "B".into()]);
        assert_eq!(a[2].1, vec!["kinase, subunit 1 [Source:x;Acc:1]".to_string()]);
        assert_eq!(a[3].1, vec!["a=b".to_string()]);
        assert_eq!(a[4], ("flag".into(), vec![]));
        assert_eq!(attr_first(raw, "ID").as_deref(), Some("gene:X"));
        assert_eq!(attr_values(raw, "Alias"), vec!["A".to_string(), "B".into()]);
        assert_eq!(format_attrs(&a, GENE_SKIP), "Alias \"A\"; description \"kinase, subunit 1 [Source:x;Acc:1]\"; k \"a=b\"; ");
    }

    #[test]
    fn unquote_like_python() {
        assert_eq!(unquote("a%3Db%2c%"), "a=b,%");
        assert_eq!(unquote("%zz%4"), "%zz%4");
        assert_eq!(unquote("%C3%A9"), "é");
        assert_eq!(unquote("%FF"), "\u{FFFD}");
        // %2C is decoded after the comma split, so it stays inside one value
        assert_eq!(attr_values("x=a%2Cb,c", "x"), vec!["a,b".to_string(), "c".into()]);
    }
}
