use sha3::{
    Shake128, Shake256,
    digest::{ExtendableOutput, Update, XofReader},
};

pub fn shake128(out_len: usize, input: &[u8]) -> Vec<u8> {
    let mut hasher = Shake128::default();
    hasher.update(input);
    let mut reader = hasher.finalize_xof();
    let mut out = vec![0u8; out_len];
    reader.read(&mut out);
    out
}

pub fn shake256(out_len: usize, input: &[u8]) -> Vec<u8> {
    let mut hasher = Shake256::default();
    hasher.update(input);
    let mut reader = hasher.finalize_xof();
    let mut out = vec![0u8; out_len];
    reader.read(&mut out);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shake128_deterministic_and_length() {
        let m = b"determinism";
        let a = shake128(64, m);
        let b = shake128(64, m);
        assert_eq!(a.len(), 64);
        assert_eq!(b.len(), 64);
        assert_eq!(a, b);
    }

    #[test]
    fn shake256_deterministic_and_length() {
        let m = b"determinism";
        let a = shake256(96, m);
        let b = shake256(96, m);
        assert_eq!(a.len(), 96);
        assert_eq!(b.len(), 96);
        assert_eq!(a, b);
    }

    #[test]
    fn shake_xof_prefix_property() {
        let m = b"prefix test";
        let a128 = shake128(64, m);
        let b128 = shake128(128, m);
        assert_eq!(&b128[..64], &a128[..]);

        let a256 = shake256(80, m);
        let b256 = shake256(160, m);
        assert_eq!(&b256[..80], &a256[..]);
    }

    #[test]
    fn shake_domain_separation_and_input_sensitivity() {
        let m = b"same input";
        let a = shake128(32, m);
        let b = shake256(32, m);
        assert_ne!(a, b, "SHAKE128 and SHAKE256 should differ");

        let m2 = b"same inpuu";
        let a1 = shake128(48, m);
        let a2 = shake128(48, m2);
        assert_ne!(a1, a2, "Different inputs should give different outputs");
    }
}
