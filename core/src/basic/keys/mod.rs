pub mod keypair;
mod private;
pub mod public;

pub use keypair::{
    KeyPair, KeyPairResult, KeypairSeeds, keygen, keygen_with_seeds,
};
pub(crate) use private::PrivateKey;
pub use public::PublicKey;

#[cfg(test)]
pub(crate) mod test_helpers {
    use super::keypair::{KeyPair, KeypairSeeds, keygen_with_seeds};
    use math::field_element::FieldElement;

    pub(crate) fn fixture_seeds() -> KeypairSeeds {
        KeypairSeeds::new([0x11; 32], [0x22; 32], [0x33; 32])
    }

    pub(crate) fn fixture_keypair() -> KeyPair<'static, FieldElement> {
        keygen_with_seeds::<FieldElement>(fixture_seeds())
            .expect("deterministic key generation should succeed")
    }
}
