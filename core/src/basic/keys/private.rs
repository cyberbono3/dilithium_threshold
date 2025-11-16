use crate::basic::sign::{
    DilithiumSignature, SigningEngine, SigningEngineConfig, SigningField,
};
use crate::dilithium::error::DilithiumError;
use crate::dilithium::params::{K, L};
use math::{Matrix, poly::Polynomial, traits::FiniteField};

#[derive(Clone, Debug)]
pub(crate) struct PrivateKey<'a, FF: FiniteField> {
    pub(crate) a: Matrix<'a, FF>, // include A here for convenience
    pub(crate) s1: [Polynomial<'a, FF>; L],
    pub(crate) s2: [Polynomial<'a, FF>; K],
}

impl<'a, FF: FiniteField> PrivateKey<'a, FF> {
    /// Construct a secret key from its parts.
    pub(crate) fn new(
        a: Matrix<'a, FF>,
        s1: [Polynomial<'a, FF>; L],
        s2: [Polynomial<'a, FF>; K],
    ) -> Self {
        Self { a, s1, s2 }
    }
}

impl<FF> PrivateKey<'static, FF>
where
    FF: SigningField,
    i64: From<FF>,
{
    /// Produce a Dilithium signature using the embedded secret key.
    pub(crate) fn sign(
        &self,
        msg: &[u8],
    ) -> Result<DilithiumSignature<'static, FF>, DilithiumError> {
        let engine = SigningEngine::new(self, msg);

        SigningEngineConfig::<FF>::rejection_attempts()
            .find_map(|ctr| engine.try_with_counter(ctr))
            .ok_or(DilithiumError::SignatureGenerationFailed)
    }
}

#[cfg(test)]
mod tests {
    use crate::basic::keys::test_helpers;
    use crate::basic::{KeyPair, keygen_with_seeds};
    use crate::dilithium::params::{K, L};
    use crate::matrix::expand_a_from_rho;
    use math::field_element::FieldElement;

    #[test]
    fn private_key_has_expected_dimensions() {
        let KeyPair { private: sk, .. } = test_helpers::fixture_keypair();

        assert_eq!(sk.s1.len(), L);
        assert_eq!(sk.s2.len(), K);
    }

    #[test]
    fn private_matrix_matches_public_and_rho() {
        let KeyPair {
            public: pk,
            private: sk,
        } = test_helpers::fixture_keypair();

        assert_eq!(sk.a.as_slice(), pk.a.as_slice());
        let from_rho = expand_a_from_rho(pk.rho);
        assert_eq!(sk.a.as_slice(), from_rho.as_slice());
    }

    #[test]
    fn private_key_is_deterministic_given_seeds() {
        let seeds = test_helpers::fixture_seeds();
        let KeyPair { private: sk1, .. } =
            keygen_with_seeds::<FieldElement>(seeds)
                .expect("key generation should succeed");
        let KeyPair { private: sk2, .. } =
            keygen_with_seeds::<FieldElement>(seeds)
                .expect("key generation should succeed");

        assert_eq!(sk1.a.as_slice(), sk2.a.as_slice());
        assert_eq!(sk1.s1, sk2.s1);
        assert_eq!(sk1.s2, sk2.s2);
    }
}
