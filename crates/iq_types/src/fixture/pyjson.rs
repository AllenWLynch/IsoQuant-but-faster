//! JSON writer that reproduces the bytes of Python's
//! `json.dumps(obj, separators=(',', ':'), ensure_ascii=True)` (schema §0.2/§0.4), so that
//! fixture lines can be round-tripped byte-for-byte:
//! - floats as `float.__repr__` (`3.0`, `1e-05`, `1.5e+16`, `-0.0`);
//! - non-ASCII characters and DEL escaped as `\uXXXX` (surrogate pairs above the BMP);
//! - compact separators.
//!
//! Key order is the serialization order of the value (struct field order, `IndexMap` order,
//! or `serde_json::Map` order — this crate does not enable `preserve_order`, so a
//! `serde_json::Value` object re-serializes with **sorted** keys; use typed DTOs for
//! byte-equal round trips, or [`to_py_string_sorted`] for `params.json`).
//!
//! Non-finite floats: plain `f64` fields cannot hold them in fixtures; fields documented as
//! `f64*` use [`PyF64`], which (de)serializes `"NaN"`, `"Infinity"`, `"-Infinity"`.

use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::io::{self, Write};

/// Python `repr(float)` for finite values (shortest round-trip digits; positional when
/// `-4 <= exp10 < 16`, else `d[.ddd]e±XX`). NaN/inf are rendered `NaN`/`inf`/`-inf`, as
/// Python `repr` does (JSON writers must not emit those; see [`PyF64`]).
pub fn py_float_repr(x: f64) -> String {
    if x.is_nan() {
        return "nan".into();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf".into() } else { "-inf".into() };
    }
    // Rust's `{:e}` gives the shortest round-trip digit *count* n ("-1.2345e2", "1e-5",
    // "0e0"), but on exact decimal ties it may pick the upper candidate. Python's repr (David
    // Gay dtoa, mode 0) returns the correctly rounded n-digit string (half-even), which Rust's
    // exact `{:.prec$e}` formatting also produces; the closest n-digit string always
    // round-trips when any n-digit string does.
    let short = format!("{x:e}");
    let n_digits = short.split_once('e').map(|(m, _)| m.bytes().filter(u8::is_ascii_digit).count()).unwrap_or(1);
    let s = format!("{:.*e}", n_digits.saturating_sub(1), x);
    debug_assert_eq!(s.parse::<f64>().ok(), Some(x));
    let (mant, exp) = s.split_once('e').expect("LowerExp has an exponent");
    let exp: i32 = exp.parse().expect("integer exponent");
    let (neg, mant) = match mant.strip_prefix('-') {
        Some(m) => (true, m),
        None => (false, mant),
    };
    let digits: String = mant.chars().filter(|c| *c != '.').collect();
    let mut out = String::with_capacity(24);
    if neg {
        out.push('-');
    }
    if (-4..16).contains(&exp) {
        if exp >= 0 {
            let int_len = exp as usize + 1;
            if digits.len() <= int_len {
                out.push_str(&digits);
                out.extend(std::iter::repeat_n('0', int_len - digits.len()));
                out.push_str(".0");
            } else {
                out.push_str(&digits[..int_len]);
                out.push('.');
                out.push_str(&digits[int_len..]);
            }
        } else {
            out.push_str("0.");
            out.extend(std::iter::repeat_n('0', (-exp - 1) as usize));
            out.push_str(&digits);
        }
    } else {
        out.push_str(&digits[..1]);
        if digits.len() > 1 {
            out.push('.');
            out.push_str(&digits[1..]);
        }
        out.push('e');
        out.push(if exp < 0 { '-' } else { '+' });
        let a = exp.unsigned_abs();
        if a < 10 {
            out.push('0');
        }
        out.push_str(&a.to_string());
    }
    out
}

/// serde_json formatter emulating Python's compact, ASCII-only `json.dumps`.
#[derive(Clone, Copy, Debug, Default)]
pub struct PyFormatter;

impl serde_json::ser::Formatter for PyFormatter {
    fn write_f64<W: ?Sized + Write>(&mut self, w: &mut W, v: f64) -> io::Result<()> {
        w.write_all(py_float_repr(v).as_bytes())
    }
    fn write_f32<W: ?Sized + Write>(&mut self, w: &mut W, v: f32) -> io::Result<()> {
        w.write_all(py_float_repr(v as f64).as_bytes())
    }
    fn write_string_fragment<W: ?Sized + Write>(&mut self, w: &mut W, frag: &str) -> io::Result<()> {
        let mut start = 0;
        for (i, c) in frag.char_indices() {
            if c.is_ascii() && c != '\x7f' {
                continue;
            }
            w.write_all(&frag.as_bytes()[start..i])?;
            let mut buf = [0u16; 2];
            for u in c.encode_utf16(&mut buf) {
                write!(w, "\\u{:04x}", u)?;
            }
            start = i + c.len_utf8();
        }
        w.write_all(&frag.as_bytes()[start..])
    }
}

/// Serializes `v` exactly as the Python capture harness writes it.
pub fn to_py_writer<W: Write, T: ?Sized + Serialize>(w: W, v: &T) -> serde_json::Result<()> {
    let mut ser = serde_json::Serializer::with_formatter(w, PyFormatter);
    v.serialize(&mut ser)
}

/// Serializes `v` to a `String` in Python byte format.
pub fn to_py_string<T: ?Sized + Serialize>(v: &T) -> serde_json::Result<String> {
    let mut buf = Vec::with_capacity(256);
    to_py_writer(&mut buf, v)?;
    Ok(String::from_utf8(buf).expect("ASCII output"))
}

/// Like [`to_py_string`] with object keys sorted (Python `sort_keys=True`, used for
/// `params.json`). Goes through `serde_json::Value`, whose maps are sorted.
pub fn to_py_string_sorted<T: ?Sized + Serialize>(v: &T) -> serde_json::Result<String> {
    to_py_string(&serde_json::to_value(v)?)
}

/// `f64` that may be non-finite: written as a number when finite, else as the strings
/// `"NaN"`, `"Infinity"`, `"-Infinity"` (schema §0.4 `f64*`). Equality is bitwise-ish:
/// NaN equals NaN, so fixture comparisons work.
#[derive(Copy, Clone, Debug, Default)]
pub struct PyF64(pub f64);

impl PartialEq for PyF64 {
    fn eq(&self, o: &Self) -> bool {
        (self.0.is_nan() && o.0.is_nan()) || self.0 == o.0
    }
}

impl Serialize for PyF64 {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let x = self.0;
        if x.is_nan() {
            s.serialize_str("NaN")
        } else if x == f64::INFINITY {
            s.serialize_str("Infinity")
        } else if x == f64::NEG_INFINITY {
            s.serialize_str("-Infinity")
        } else {
            s.serialize_f64(x)
        }
    }
}

impl<'de> Deserialize<'de> for PyF64 {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        #[derive(Deserialize)]
        #[serde(untagged)]
        enum R<'a> {
            N(f64),
            #[serde(borrow)]
            S(std::borrow::Cow<'a, str>),
        }
        match R::deserialize(d)? {
            R::N(x) => Ok(PyF64(x)),
            R::S(s) => match &*s {
                "NaN" => Ok(PyF64(f64::NAN)),
                "Infinity" => Ok(PyF64(f64::INFINITY)),
                "-Infinity" => Ok(PyF64(f64::NEG_INFINITY)),
                o => Err(serde::de::Error::custom(format!("invalid f64* string {o:?}"))),
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn float_repr_matches_python() {
        let cases: &[(f64, &str)] = &[
            (0.0, "0.0"),
            (-0.0, "-0.0"),
            (3.0, "3.0"),
            (0.1, "0.1"),
            (1e-05, "1e-05"),
            (0.0001, "0.0001"),
            (0.00012345, "0.00012345"),
            (1e16, "1e+16"),
            (1.5e16, "1.5e+16"),
            (1234567890123456.0, "1234567890123456.0"),
            (-0.5, "-0.5"),
            (123.456, "123.456"),
            (2.5e-300, "2.5e-300"),
            (1e100, "1e+100"),
            (0.30000000000000004, "0.30000000000000004"),
            (100.0, "100.0"),
            (1.7976931348623157e308, "1.7976931348623157e+308"),
            (5e-324, "5e-324"),
            // exact decimal tie at 17 digits: Python rounds half-even
            (f64::from_bits(0x3fc2_b020_0000_0000), "0.14599990844726562"),
        ];
        for &(x, want) in cases {
            assert_eq!(py_float_repr(x), want, "{x:e}");
        }
    }

    #[test]
    fn strings_ascii_escaped() {
        assert_eq!(to_py_string("a\u{e9}\u{7f}\u{1F600}\n\"").unwrap(), r#""a\u00e9\u007f\ud83d\ude00\n\"""#);
        assert_eq!(to_py_string("\u{1}").unwrap(), r#""\u0001""#);
    }

    #[test]
    fn compact_and_nonfinite() {
        #[derive(Serialize, Deserialize, PartialEq, Debug)]
        struct S {
            a: Vec<f64>,
            b: PyF64,
            c: PyF64,
            d: i64,
        }
        let s = S { a: vec![1.0, 2.5], b: PyF64(f64::NAN), c: PyF64(-1e-7), d: 1 << 31 };
        let j = to_py_string(&s).unwrap();
        assert_eq!(j, r#"{"a":[1.0,2.5],"b":"NaN","c":-1e-07,"d":2147483648}"#);
        assert_eq!(serde_json::from_str::<S>(&j).unwrap(), s);
    }
}
