use crate::basic::keypair::PublicKey;
use crate::dilithium::shamir::ShamirShare;
use math::traits::FiniteField;

/// Threshold participant's key material and associated public key.
#[derive(Clone, Debug)]
pub struct ThresholdKeyShare<'a, FF: FiniteField> {
    pub participant_id: usize,
    pub s1_share: ShamirShare<'a, FF>,
    pub s2_share: ShamirShare<'a, FF>,
    pub public_key: PublicKey<'a, FF>,
}

impl<FF: FiniteField> ThresholdKeyShare<'static, FF> {
    /// Package the participant identifier, Shamir shares, and public key into a threshold share.
    pub fn new(
        participant_id: usize,
        s1_share: ShamirShare<'static, FF>,
        s2_share: ShamirShare<'static, FF>,
        public_key: PublicKey<'static, FF>,
    ) -> Self {
        Self {
            participant_id,
            s1_share,
            s2_share,
            public_key,
        }
    }
}

impl<FF: FiniteField> std::fmt::Display for ThresholdKeyShare<'static, FF> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "ThresholdKeyShare(id={})", self.participant_id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::keypair::keygen;
    use math::{field_element::FieldElement, prelude::*};

    #[test]
    fn test_key_share_creation() {
        let poly1: Polynomial<'_, FieldElement> = poly![1, 2, 3];
        let poly2: Polynomial<'_, FieldElement> = poly![4, 5, 6];
        let s1_share = ShamirShare::new(1, vec![poly1]).unwrap();
        let s2_share = ShamirShare::new(1, vec![poly2]).unwrap();

        let pub_key = keygen::<FieldElement>()
            .expect("key generation should succeed")
            .public;

        let key_share = ThresholdKeyShare::new(1, s1_share, s2_share, pub_key);

        assert_eq!(key_share.participant_id, 1);
        assert_eq!(key_share.s1_share.participant_id, 1);
        assert_eq!(key_share.s2_share.participant_id, 1);
    }

    #[test]
    fn test_key_share_display() {
        let poly1: Polynomial<'_, FieldElement> = poly![vec![1]];
        let poly2: Polynomial<'_, FieldElement> = poly![vec![2]];
        let s1_share = ShamirShare::new(5, vec![poly1]).unwrap();
        let s2_share = ShamirShare::new(5, vec![poly2]).unwrap();

        let pub_key = keygen::<FieldElement>()
            .expect("key generation should succeed")
            .public;
        let key_share = ThresholdKeyShare::new(5, s1_share, s2_share, pub_key);

        let display_str = format!("{}", key_share);
        assert_eq!(display_str, "ThresholdKeyShare(id=5)");
    }
}
