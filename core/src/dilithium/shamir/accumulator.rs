use math::{prelude::*, traits::FiniteField};

use crate::dilithium::error::DilithiumResult;

use super::share::ShamirShare;

pub(super) struct ShareAccumulator<FF: FiniteField> {
    pub(super) participant_id: usize,
    pub(super) buffers: Vec<Vec<FF>>,
}

impl<FF: FiniteField> ShareAccumulator<FF> {
    pub(super) fn new(participant_id: usize, lengths: &[usize]) -> Self {
        let buffers = lengths
            .iter()
            .map(|&len| vec![FF::default(); len])
            .collect();
        Self {
            participant_id,
            buffers,
        }
    }

    pub(super) fn insert(
        &mut self,
        poly_idx: usize,
        coeff_idx: usize,
        value: FF,
    ) {
        if poly_idx >= self.buffers.len() {
            return;
        }
        let poly = &mut self.buffers[poly_idx];
        if coeff_idx >= poly.len() {
            return;
        }
        poly[coeff_idx] = value;
    }

    pub(super) fn finalize(self) -> DilithiumResult<ShamirShare<'static, FF>> {
        let polynomials =
            self.buffers.into_iter().map(Polynomial::from).collect();
        ShamirShare::new(self.participant_id, polynomials)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dilithium::error::DilithiumError;
    use crate::dilithium::shamir::error::ShamirError;
    use math::prelude::*;
    use num_traits::Zero;

    #[test]
    fn accumulators_initialize_with_zero_buffers() {
        let lengths = [2usize, 3, 0];
        let accumulator = ShareAccumulator::<FieldElement>::new(7, &lengths);
        assert_eq!(accumulator.participant_id, 7);
        assert_eq!(accumulator.buffers.len(), lengths.len());
        for (buffer, &expected_len) in accumulator.buffers.iter().zip(&lengths)
        {
            assert_eq!(buffer.len(), expected_len);
            assert!(buffer.iter().all(|value| value.is_zero()));
        }
    }

    #[test]
    fn insert_updates_target_coefficient_only() {
        let lengths = [3usize, 2];
        let mut accumulator =
            ShareAccumulator::<FieldElement>::new(2, &lengths);
        let value = FieldElement::from(123i32);
        accumulator.insert(0, 2, value);

        assert_eq!(accumulator.buffers[0][2], value);
        assert!(accumulator.buffers[0][0].is_zero());
        assert!(accumulator.buffers[0][1].is_zero());
        assert!(accumulator.buffers[1].iter().all(|coeff| coeff.is_zero()));
    }

    #[test]
    fn insert_overwrites_existing_entry() {
        let mut accumulator =
            ShareAccumulator::<FieldElement>::new(9, &[2usize]);
        accumulator.insert(0, 1, FieldElement::from(40i32));
        accumulator.insert(0, 1, FieldElement::from(55i32));

        assert_eq!(accumulator.buffers[0][1], FieldElement::from(55i32));
    }

    #[test]
    fn insert_handles_multiple_polynomials() {
        let mut accumulator =
            ShareAccumulator::<FieldElement>::new(5, &[2usize, 3, 1]);
        accumulator.insert(0, 0, FieldElement::from(1i32));
        accumulator.insert(0, 1, FieldElement::from(2i32));
        accumulator.insert(1, 0, FieldElement::from(10i32));
        accumulator.insert(1, 2, FieldElement::from(30i32));
        accumulator.insert(2, 0, FieldElement::from(-5i32));

        assert_eq!(accumulator.buffers[0], vec![fe!(1), fe!(2)]);
        assert_eq!(
            accumulator.buffers[1],
            vec![fe!(10), FieldElement::default(), fe!(30)]
        );
        assert_eq!(accumulator.buffers[2], vec![fe!(-5)]);
    }

    #[test]
    fn insert_ignores_out_of_bounds_poly_indices() {
        let mut accumulator =
            ShareAccumulator::<FieldElement>::new(4, &[2usize, 1]);
        let snapshot = accumulator.buffers.clone();

        accumulator.insert(5, 0, FieldElement::from(99i32));

        assert_eq!(accumulator.buffers, snapshot);
    }

    #[test]
    fn insert_ignores_out_of_bounds_coeff_indices() {
        let mut accumulator =
            ShareAccumulator::<FieldElement>::new(6, &[1usize]);
        accumulator.insert(0, 0, FieldElement::from(33i32));
        let snapshot = accumulator.buffers.clone();

        accumulator.insert(0, 3, FieldElement::from(77i32));

        assert_eq!(accumulator.buffers, snapshot);
    }

    #[test]
    fn insert_ignores_zero_length_polynomial() {
        let mut accumulator =
            ShareAccumulator::<FieldElement>::new(7, &[0usize, 2]);
        let snapshot = accumulator.buffers.clone();

        accumulator.insert(0, 0, FieldElement::from(99i32));

        assert_eq!(accumulator.buffers, snapshot);
    }

    #[test]
    fn finalize_produces_shamir_share() {
        let mut accumulator = ShareAccumulator::<FieldElement>::new(3, &[1, 2]);
        accumulator.insert(0, 0, FieldElement::from(7));
        accumulator.insert(1, 1, FieldElement::from(11));

        let share = accumulator.finalize().expect("finalization succeeds");
        assert_eq!(share.participant_id, 3);
        assert_eq!(share.vector_length(), 2);
        assert_eq!(
            share.share_vector[0].coefficients()[0],
            FieldElement::from(7)
        );
        assert_eq!(
            share.share_vector[1].coefficients()[1],
            FieldElement::from(11)
        );
    }

    #[test]
    fn finalize_validates_participant_id() {
        let accumulator = ShareAccumulator::<FieldElement>::new(1, &[0]);
        let share = accumulator.finalize();
        assert!(share.is_ok());

        let accumulator_invalid =
            ShareAccumulator::<FieldElement>::new(0, &[0]);
        assert!(matches!(
            accumulator_invalid.finalize(),
            Err(DilithiumError::Shamir(ShamirError::InvalidParticipantId(0)))
        ));
    }

    #[test]
    fn finalize_preserves_zero_length_polynomials() {
        let accumulator =
            ShareAccumulator::<FieldElement>::new(8, &[0usize, 3]);
        let share = accumulator.finalize().expect("share creation succeeds");

        assert_eq!(share.vector_length(), 2);
        assert!(share.share_vector[0].coefficients().is_empty());
        assert!(share.share_vector[1].coefficients().is_empty());
    }
}
