use math::{prelude::*, traits::FiniteField};

use crate::dilithium::error::DilithiumResult;
use crate::dilithium::shamir::error::ShamirError;
use crate::traits::PolyVectorSource;

/// Adapted Shamir's Secret Sharing share wrapper.
#[derive(Clone, Debug, PartialEq)]
pub struct ShamirShare<'a, FF: FiniteField> {
    pub participant_id: usize,
    pub share_vector: Vec<Polynomial<'a, FF>>,
}

impl<FF: FiniteField> ShamirShare<'static, FF> {
    pub fn new(
        participant_id: usize,
        share_vector: Vec<Polynomial<'static, FF>>,
    ) -> DilithiumResult<Self> {
        if participant_id == 0 {
            return Err(
                ShamirError::InvalidParticipantId(participant_id).into()
            );
        }

        Ok(ShamirShare {
            participant_id,
            share_vector,
        })
    }

    pub fn vector_length(&self) -> usize {
        self.share_vector.len()
    }
}

impl<FF: FiniteField> PolyVectorSource<FF> for ShamirShare<'static, FF> {
    fn participant_id(&self) -> usize {
        self.participant_id
    }

    fn polynomials(&self) -> &[Polynomial<'static, FF>] {
        &self.share_vector
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_share_creation() {
        let poly1: Polynomial<'_, FieldElement> = poly![1, 2, 3, 4, 5];
        let poly2: Polynomial<'_, FieldElement> = poly![6, 7, 8, 9, 10];
        let share_vector = vec![poly1, poly2];
        let share = ShamirShare::new(1, share_vector.clone()).unwrap();

        assert_eq!(share.participant_id, 1);
        assert_eq!(share.vector_length(), 2);
        assert_eq!(share.share_vector, share_vector);
    }

    #[test]
    fn test_invalid_participant_id() {
        let poly1: Polynomial<'_, FieldElement> = poly![1, 2, 3, 4, 5];
        let share_vector = vec![poly1];

        assert!(ShamirShare::new(0, share_vector.clone()).is_err());
    }

    #[test]
    fn test_share_debug_representation() {
        let poly1: Polynomial<'_, FieldElement> = poly![1, 2, 3];
        let share_vector = vec![poly1];
        let share = ShamirShare::new(1, share_vector).unwrap();

        let debug_str = format!("{:?}", share);
        assert!(debug_str.contains("ShamirShare"));
        assert!(debug_str.contains("participant_id: 1"));
    }
}
