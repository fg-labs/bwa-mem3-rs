//! Test-only helpers shared by the workspace's test suites, so each suite does
//! not carry its own copy of the same decoder. Not published; every crate uses
//! it as a path dev-dependency.

/// One aux tag of a packed BAM record: its key, its type byte, and its raw
/// value bytes (for `Z`/`H`, without the trailing NUL).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AuxTag<'a> {
    pub key: [u8; 2],
    pub ty: u8,
    pub value: &'a [u8],
}

/// The byte offset at which a packed BAM record body's aux tags start: the
/// fixed fields, then the read name, CIGAR, packed SEQ and QUAL. `body` is the
/// record without its leading `block_size`.
///
/// # Panics
///
/// If the body is shorter than its fixed fields or than the lengths those
/// fields declare: a truncated record, not one without tags (`start ==
/// body.len()` is a complete record with none).
#[must_use]
pub fn aux_start(body: &[u8]) -> usize {
    assert!(
        body.len() >= 32,
        "record body shorter than its fixed fields"
    );
    let l_read_name = usize::from(body[8]);
    let n_cigar = usize::from(u16::from_le_bytes([body[12], body[13]]));
    let l_seq = u32::from_le_bytes(body[16..20].try_into().unwrap()) as usize;
    let start = [l_read_name, 4 * n_cigar, l_seq.div_ceil(2), l_seq]
        .into_iter()
        .try_fold(32usize, usize::checked_add)
        .expect("record field lengths overflow usize");
    assert!(
        start <= body.len(),
        "record truncated: its fields need {start} bytes but the body has {}",
        body.len()
    );
    start
}

/// The aux tags of a packed BAM record body (without its leading
/// `block_size`), in order.
///
/// # Panics
///
/// On a truncated body (see [`aux_start`]), an unknown type byte, or a tag
/// whose value runs past the end of the body, rather than reporting a
/// truncated tag as present.
#[must_use]
pub fn aux_tags(body: &[u8]) -> Vec<AuxTag<'_>> {
    let mut o = aux_start(body);
    let mut tags = Vec::new();
    while o < body.len() {
        assert!(o + 3 <= body.len(), "aux tag header truncated at byte {o}");
        let (key, ty) = ([body[o], body[o + 1]], body[o + 2]);
        o += 3;
        let rest = &body[o..];
        let (value_len, skip) = match ty {
            b'A' | b'c' | b'C' => (1, 1),
            b's' | b'S' => (2, 2),
            b'i' | b'I' | b'f' => (4, 4),
            b'Z' | b'H' => {
                let nul = rest
                    .iter()
                    .position(|&b| b == 0)
                    .unwrap_or_else(|| panic!("{key:?}: string value has no NUL"));
                (nul, nul + 1)
            }
            b'B' => {
                assert!(rest.len() >= 5, "{key:?}: B-array header truncated");
                let width = match rest[0] {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    _ => 4,
                };
                let count = u32::from_le_bytes(rest[1..5].try_into().unwrap()) as usize;
                let n = 5 + width * count;
                (n, n)
            }
            other => panic!("unknown aux type {:?}", other as char),
        };
        assert!(
            skip <= rest.len(),
            "{key:?}: value needs {skip} bytes but {} remain",
            rest.len()
        );
        tags.push(AuxTag {
            key,
            ty,
            value: &rest[..value_len],
        });
        o += skip;
    }
    tags
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A body with no read name, CIGAR or SEQ, followed by `aux`.
    fn body_with(aux: &[u8]) -> Vec<u8> {
        let mut b = vec![0u8; 32];
        b.extend_from_slice(aux);
        b
    }

    #[test]
    fn decodes_each_type_and_strips_the_nul() {
        let mut aux = Vec::new();
        aux.extend_from_slice(b"NMi");
        aux.extend_from_slice(&3i32.to_le_bytes());
        aux.extend_from_slice(b"MDZ10A5\0");
        aux.extend_from_slice(b"XSC\x07");
        let body = body_with(&aux);
        let tags = aux_tags(&body);
        assert_eq!(
            tags,
            vec![
                AuxTag {
                    key: *b"NM",
                    ty: b'i',
                    value: &3i32.to_le_bytes()
                },
                AuxTag {
                    key: *b"MD",
                    ty: b'Z',
                    value: b"10A5"
                },
                AuxTag {
                    key: *b"XS",
                    ty: b'C',
                    value: b"\x07"
                },
            ]
        );
    }

    #[test]
    fn a_complete_record_without_tags_has_none() {
        assert!(aux_tags(&body_with(&[])).is_empty());
    }

    #[test]
    #[should_panic(expected = "record truncated")]
    fn rejects_a_body_shorter_than_its_declared_fields() {
        let mut body = vec![0u8; 32];
        body[8] = 5; // a 5-byte read name that is not there
        let _ = aux_tags(&body);
    }

    #[test]
    #[should_panic(expected = "value needs 4 bytes")]
    fn rejects_a_value_running_past_the_body() {
        let _ = aux_tags(&body_with(b"NMi\x01\x00"));
    }

    #[test]
    #[should_panic(expected = "has no NUL")]
    fn rejects_an_unterminated_string() {
        let _ = aux_tags(&body_with(b"MDZ10"));
    }
}
