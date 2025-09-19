#[cfg(test)]
mod tests {
    use super::*;
    use crate::keypair::keygen;
    use crate::sign::{sign, verify};

    #[test]
    fn round_trip_sign_verify() {
        let (pk, sk) = keygen();

        let msg = b"The lattice rocks: Dilithium / ML-DSA demo!";
        let sig = sign(&sk.a, &sk.s1, &sk.s2, &pk.t, msg);

        assert!(verify(&pk.a, &pk.t, msg, &sig));
    }

    #[test]
    fn negative_case_wrong_message() {
        let (pk, sk) = keygen();
        let msg = b"hello world";
        let sig = sign(&sk.a, &sk.s1, &sk.s2, &pk.t, msg);

        let wrong = b"hello wurld";
        assert!(!verify(&pk.a, &pk.t, wrong, &sig));
    }
}
