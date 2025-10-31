use sha3::{
    Shake128, Shake256,
    digest::{ExtendableOutput, Update, XofReader},
};

/// Hash helper that captures the shared mechanics of SHAKE-based XOFs.
fn squeeze_xof<Hash>(out_len: usize, input: &[u8]) -> Vec<u8>
where
    Hash: Default + Update + ExtendableOutput,
{
    let mut hasher = Hash::default();
    hasher.update(input);
    read_xof(hasher.finalize_xof(), out_len)
}

fn read_xof<R>(mut reader: R, out_len: usize) -> Vec<u8>
where
    R: XofReader,
{
    let mut out = vec![0u8; out_len];
    reader.read(&mut out);
    out
}

pub fn shake128(out_len: usize, input: &[u8]) -> Vec<u8> {
    squeeze_xof::<Shake128>(out_len, input)
}

pub fn shake256(out_len: usize, input: &[u8]) -> Vec<u8> {
    squeeze_xof::<Shake256>(out_len, input)
}

#[cfg(test)]
mod tests {

    macro_rules! xof_smoke_tests {
        ($modname:ident, $f:path) => {
            mod $modname {
                #[test]
                fn deterministic_and_length() {
                    let m = b"determinism";
                    let a = $f(64, m);
                    let b = $f(64, m);
                    assert_eq!(a.len(), 64);
                    assert_eq!(b.len(), 64);
                    assert_eq!(a, b);
                }
                #[test]
                fn prefix_property() {
                    let m = b"prefix test";
                    let a = $f(64, m);
                    let b = $f(128, m);
                    assert_eq!(&b[..64], &a[..]);
                }
                #[test]
                fn input_sensitivity() {
                    let a = $f(48, b"abc");
                    let b = $f(48, b"abd");
                    assert_ne!(a, b);
                }
                #[test]
                fn empty_message_is_ok() {
                    let out = $f(32, b"");
                    assert_eq!(out.len(), 32);
                }
            }
        };
    }

    xof_smoke_tests!(shake128_suite, crate::matrix::hash::shake128);
    xof_smoke_tests!(shake256_suite, crate::matrix::hash::shake256);
}
