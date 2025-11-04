use crate::traits::PolyVectorSource;
use math::{prelude::*, traits::FiniteField};

/// Represents a partial signature from one participant.
///
/// Contains the participant's contribution to the threshold signature
/// along with verification information.
#[derive(Clone, Debug)]
pub struct PartialSignature<'sign, FF: FiniteField> {
    pub participant_id: usize,
    pub z_partial: Vec<Polynomial<'sign, FF>>,
    pub commitment: Vec<Polynomial<'sign, FF>>,
    pub challenge: Polynomial<'sign, FF>,
}

impl<FF: FiniteField> PartialSignature<'static, FF> {
    /// Bundle the per-participant signature artefacts into a single structure.
    pub fn new(
        participant_id: usize,
        z_partial: Vec<Polynomial<'static, FF>>,
        commitment: Vec<Polynomial<'static, FF>>,
        challenge: Polynomial<'static, FF>,
    ) -> Self {
        Self {
            participant_id,
            z_partial,
            commitment,
            challenge,
        }
    }
}

impl<FF: FiniteField> std::fmt::Display for PartialSignature<'static, FF> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "PartialSignature(id={})", self.participant_id)
    }
}

impl<FF: FiniteField> PolyVectorSource<FF> for PartialSignature<'static, FF> {
    fn participant_id(&self) -> usize {
        self.participant_id
    }

    fn polynomials(&self) -> &[Polynomial<'static, FF>] {
        &self.z_partial
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use math::field_element::FieldElement;
    use num_traits::Zero;

    #[test]
    fn test_partial_signature_creation() {
        let z_partial = vec![poly![100, 200, 300]];
        let commitment = vec![poly![10, 20, 30]];
        let challenge: Polynomial<'static, FieldElement> = poly![1, -1, 0, 1];

        let partial_sig =
            PartialSignature::new(3, z_partial, commitment, challenge.clone());

        assert_eq!(partial_sig.participant_id, 3);
        assert_eq!(partial_sig.z_partial.len(), 1);
        assert_eq!(partial_sig.commitment.len(), 1);
        assert_eq!(partial_sig.challenge, challenge);
    }

    #[test]
    fn test_partial_signature_display() {
        let zero_poly = Polynomial::<FieldElement>::zero();
        let z_partial = vec![zero_poly.clone()];
        let commitment = vec![zero_poly.clone()];
        let challenge = zero_poly;

        let partial_sig =
            PartialSignature::new(7, z_partial, commitment, challenge);

        let display_str = format!("{}", partial_sig);
        assert_eq!(display_str, "PartialSignature(id=7)");
    }
}
