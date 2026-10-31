//! Read grouping from one BAM tag (`read_groups.py:49-59`): `str(value)`, missing → `NA`
//! (D-O2). See `iq_types::TagValue::to_group_string`.

use iq_types::read::MISSING_GROUP;
use iq_types::TagValue;
use std::borrow::Cow;

/// Group string of one record: the tag value as Python `str(value)`, or `NA` when the tag is
/// missing. String tags are borrowed, not copied.
pub fn group_string(tag: Option<&TagValue>) -> Cow<'_, str> {
    match tag {
        None => Cow::Borrowed(MISSING_GROUP),
        Some(TagValue::Str(s)) => Cow::Borrowed(s.as_str()),
        Some(t) => Cow::Owned(t.to_group_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn groups() {
        assert_eq!(group_string(None), "NA");
        assert_eq!(group_string(Some(&TagValue::Int(2))), "2");
        assert_eq!(group_string(Some(&TagValue::Str("HP1".into()))), "HP1");
        assert_eq!(group_string(Some(&TagValue::Char('x'))), "x");
    }
}
