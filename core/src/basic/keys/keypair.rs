use super::{PrivateKey, PublicKey};
use crate::basic::sign::{DilithiumSignature, SigningField};
use crate::basic::utils::expand_secret_vector;
use crate::dilithium::error::DilithiumError;
use crate::dilithium::params::{K, L};
use crate::dilithium::utils::random_bytes;
use crate::matrix::expand_a_from_rho;
use math::{Matrix, error::MatrixError, poly::Polynomial, traits::FiniteField};

pub type KeyPairResult<FF> = Result<KeyPair<'static, FF>, MatrixError>;

/// Convenience container bundling a Dilithium keypair.
#[derive(Clone, Debug)]
pub struct KeyPair<'a, FF: FiniteField> {
    pub public: PublicKey<'a, FF>,
    pub(crate) private: PrivateKey<'a, FF>,
}

/// Read-only view into the secret key components.
#[derive(Clone, Copy)]
pub struct SecretKeyView<'a, FF: FiniteField> {
    pub a: &'a Matrix<'a, FF>,
    pub s1: &'a [Polynomial<'a, FF>; L],
    pub s2: &'a [Polynomial<'a, FF>; K],
}

impl<'a, FF: FiniteField> From<&'a PrivateKey<'a, FF>>
    for SecretKeyView<'a, FF>
{
    fn from(value: &'a PrivateKey<'a, FF>) -> Self {
        Self {
            a: &value.a,
            s1: &value.s1,
            s2: &value.s2,
        }
    }
}

impl<'a, FF: FiniteField> KeyPair<'a, FF> {
    pub(crate) fn new(
        public: PublicKey<'a, FF>,
        private: PrivateKey<'a, FF>,
    ) -> Self {
        Self { public, private }
    }

    /// Expose the secret key fields for verification/testing purposes.
    pub fn secret(&self) -> SecretKeyView<'_, FF> {
        SecretKeyView::from(&self.private)
    }
}

impl<FF> KeyPair<'static, FF>
where
    FF: SigningField,
    i64: From<FF>,
{
    /// Sign `msg` with the embedded private key.
    pub fn sign(
        &self,
        msg: &[u8],
    ) -> Result<DilithiumSignature<'static, FF>, DilithiumError> {
        self.private.sign(msg)
    }

    /// Verify `sig` against `msg` using the embedded public key.
    pub fn verify(&self, msg: &[u8], sig: &DilithiumSignature<'_, FF>) -> bool {
        self.public.verify(msg, sig)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct KeypairSeeds {
    pub rho: [u8; 32],
    pub s1: [u8; 32],
    pub s2: [u8; 32],
}

impl KeypairSeeds {
    #[inline]
    /// Generate fresh random seeds for key generation.
    pub fn random() -> Self {
        Self {
            rho: random_bytes(),
            s1: random_bytes(),
            s2: random_bytes(),
        }
    }

    #[inline]
    /// Assemble seeds from explicit 32-byte inputs.
    pub const fn new(rho: [u8; 32], s1: [u8; 32], s2: [u8; 32]) -> Self {
        Self { rho, s1, s2 }
    }

    /// Consume the seeds to deterministically derive a Dilithium keypair.
    fn into_keypair<FF>(self) -> KeyPairResult<FF>
    where
        FF: math::traits::FiniteField + From<i64>,
    {
        use crate::matrix::traits::MatrixMulExt;

        let KeypairSeeds { rho, s1, s2 } = self;
        let a = expand_a_from_rho(rho);
        let s1 = expand_secret_vector::<L, FF>(&s1);
        let s2 = expand_secret_vector::<K, FF>(&s2);

        let mut t: [Polynomial<'static, FF>; K] = a.mul_polynomials(&s1)?;

        t.iter_mut()
            .zip(s2.iter())
            .for_each(|(dest, addend)| *dest += addend);

        let public_key = PublicKey::new(a.clone(), t, rho);
        let private_key = PrivateKey::new(a, s1, s2);

        Ok(KeyPair::new(public_key, private_key))
    }
}

impl<FF> From<KeypairSeeds> for KeyPairResult<FF>
where
    FF: math::traits::FiniteField + From<i64>,
{
    fn from(seeds: KeypairSeeds) -> Self {
        seeds.into_keypair()
    }
}

/// Generate a random Dilithium keypair.
pub fn keygen<FF: FiniteField + From<i64>>() -> KeyPairResult<FF> {
    KeypairSeeds::random().into()
}

/// Derive a deterministic keypair from caller supplied seeds.
pub fn keygen_with_seeds<FF: FiniteField + From<i64>>(
    seeds: KeypairSeeds,
) -> KeyPairResult<FF> {
    seeds.into()
}
