use math::{prelude::*, traits::FiniteField};

use crate::dilithium::error::{DilithiumError, DilithiumResult};
use crate::dilithium::params::validate_threshold_config;
use crate::dilithium::shamir::error::ShamirError;
use crate::dilithium::utils::reconstruct_vector_from_points;

use super::accumulator::ShareAccumulator;
use super::share::ShamirShare;

/// Adapted Shamir secret sharing scheme tailored for polynomial vectors.
#[derive(Debug)]
pub struct AdaptedShamirSSS {
    threshold: usize,
    participant_number: usize,
}

impl AdaptedShamirSSS {
    /// Allocate accumulators for each participant with buffers sized to `lengths`.
    fn share_accumulators<FF: FiniteField>(
        &self,
        lengths: &[usize],
    ) -> Vec<ShareAccumulator<FF>> {
        (1..=self.participant_number)
            .map(|pid| ShareAccumulator::new(pid, lengths))
            .collect()
    }

    /// Create a sharing scheme for the provided threshold/participant configuration.
    pub fn new(
        threshold: usize,
        participant_number: usize,
    ) -> DilithiumResult<Self> {
        if !validate_threshold_config(threshold, participant_number) {
            return Err(DilithiumError::InvalidThreshold(
                threshold,
                participant_number,
            ));
        }

        Ok(AdaptedShamirSSS {
            threshold,
            participant_number,
        })
    }

    /// Split a vector of polynomials into participant shares.
    pub fn split_secret<FF>(
        &self,
        secret_vector: &[Polynomial<'static, FF>],
    ) -> DilithiumResult<Vec<ShamirShare<'static, FF>>>
    where
        FF: FiniteField,
        rand::distributions::Standard: rand::distributions::Distribution<FF>,
    {
        if secret_vector.is_empty() {
            return Ok(Vec::new());
        }

        let lengths: Vec<usize> = secret_vector
            .iter()
            .map(|poly| poly.coefficients().len())
            .collect();

        let mut accumulators: Vec<ShareAccumulator<FF>> =
            self.share_accumulators(&lengths);

        let mut rng = rand::thread_rng();

        for (poly_idx, polynomial) in secret_vector.iter().enumerate() {
            for (coeff_idx, &secret_coeff) in
                polynomial.coefficients().iter().enumerate()
            {
                let shamir_poly =
                    self.create_shamir_polynomial_with(&mut rng, secret_coeff);
                for (pid, accumulator) in accumulators.iter_mut().enumerate() {
                    let point = FF::from((pid + 1) as u64);
                    let value = shamir_poly.evaluate_in_same_field(point);
                    accumulator.insert(poly_idx, coeff_idx, value);
                }
            }
        }

        accumulators
            .into_iter()
            .map(ShareAccumulator::finalize)
            .collect()
    }

    /// Reconstruct the entire secret vector from the first `threshold` shares.
    pub fn reconstruct_secret<FF: FiniteField>(
        &self,
        shares: &[ShamirShare<'static, FF>],
    ) -> DilithiumResult<Vec<Polynomial<'static, FF>>> {
        let (active_shares, vector_length) =
            self.select_active_shares(shares)?;
        let poly_indices: Vec<usize> = (0..vector_length).collect();
        Self::reconstruct_poly_vector(active_shares, &poly_indices)
    }

    #[cfg(test)]
    /// Reconstruct only the requested polynomial indices (test helper).
    pub fn partial_reconstruct<FF: FiniteField>(
        &self,
        shares: &[ShamirShare<'static, FF>],
        poly_indices: &[usize],
    ) -> DilithiumResult<Vec<Polynomial<'static, FF>>> {
        let (active_shares, vector_length) =
            self.select_active_shares(shares)?;
        Self::reconstruct_poly_vector(active_shares, poly_indices)
    }

    /// Helper to invoke the shared polynomial-vector reconstruction utility.
    fn reconstruct_poly_vector<FF: FiniteField>(
        shares: &[ShamirShare<'static, FF>],
        poly_indices: &[usize],
    ) -> DilithiumResult<Vec<Polynomial<'static, FF>>> {
        reconstruct_vector_from_points::<FF, _>(shares, poly_indices)
    }

    /// Take the first `threshold` shares and ensure they are consistent in length.
    fn select_active_shares<'a, FF>(
        &self,
        shares: &'a [ShamirShare<'static, FF>],
    ) -> DilithiumResult<(&'a [ShamirShare<'static, FF>], usize)>
    where
        FF: FiniteField,
    {
        if shares.len() < self.threshold {
            return Err(DilithiumError::InsufficientShares(
                self.threshold,
                shares.len(),
            ));
        }

        let active = &shares[..self.threshold];
        let vector_length = Self::ensure_consistent_lengths(active)?;
        Ok((active, vector_length))
    }

    /// Validate that all shares expose the same polynomial vector length.
    fn ensure_consistent_lengths<FF>(
        shares: &[ShamirShare<'static, FF>],
    ) -> DilithiumResult<usize>
    where
        FF: FiniteField,
    {
        let Some(first_share) = shares.first() else {
            return Ok(0);
        };

        let expected_length = first_share.vector_length();
        if shares
            .iter()
            .any(|share| share.vector_length() != expected_length)
        {
            return Err(ShamirError::InconsistentShareLengths.into());
        }

        Ok(expected_length)
    }

    #[cfg(test)]
    /// Convenience helper used in tests to generate a random sharing polynomial.
    fn create_shamir_polynomial<FF>(
        &self,
        secret: &FF,
    ) -> Polynomial<'static, FF>
    where
        FF: FiniteField,
        rand::distributions::Standard: rand::distributions::Distribution<FF>,
    {
        let mut rng = rand::thread_rng();
        self.create_shamir_polynomial_with(&mut rng, *secret)
    }

    /// Build a random polynomial of degree threshold-1 with the constant term set to `secret`.
    fn create_shamir_polynomial_with<FF, R>(
        &self,
        rng: &mut R,
        secret: FF,
    ) -> Polynomial<'static, FF>
    where
        FF: FiniteField,
        R: rand::Rng + ?Sized,
        rand::distributions::Standard: rand::distributions::Distribution<FF>,
    {
        let coefficients: Vec<FF> = std::iter::once(secret)
            .chain(
                std::iter::repeat_with(|| rand::Rng::r#gen(rng))
                    .take(self.threshold.saturating_sub(1)),
            )
            .collect();
        Polynomial::from(coefficients)
    }

    #[cfg(test)]
    fn evaluate_polynomial(coeffs: &[i32], x: i32) -> i32 {
        let mut result = 0i64;
        let mut x_power = 1i64;

        for &coeff in coeffs {
            result = (result + (coeff as i64 * x_power)) % (Q as i64);
            x_power = (x_power * x as i64) % (Q as i64);
        }

        result as i32
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::Zero;

    impl Default for AdaptedShamirSSS {
        fn default() -> Self {
            Self {
                threshold: 3,
                participant_number: 5,
            }
        }
    }

    mod adapted_shamir_sss_tests {
        use super::*;
        use num_traits::Zero;

        fn setup_adapted_shamir(
            threshold: usize,
            participants: usize,
        ) -> DilithiumResult<AdaptedShamirSSS> {
            AdaptedShamirSSS::new(threshold, participants)
        }

        fn share_with_length(
            id: usize,
            polys: usize,
        ) -> ShamirShare<'static, FieldElement> {
            let vector = (0..polys)
                .map(|offset| {
                    Polynomial::from(vec![FieldElement::from(
                        (id + offset + 1) as i64,
                    )])
                })
                .collect();
            ShamirShare::new(id, vector).expect("valid share")
        }

        #[test]
        fn ensure_consistent_lengths_empty_slice() {
            let len = AdaptedShamirSSS::ensure_consistent_lengths::<
                FieldElement,
            >(&[])
            .expect("empty input should succeed");
            assert_eq!(len, 0);
        }

        #[test]
        fn ensure_consistent_lengths_reports_error() {
            let shares = vec![share_with_length(1, 2), share_with_length(2, 3)];

            let err =
                AdaptedShamirSSS::ensure_consistent_lengths::<FieldElement>(
                    &shares,
                )
                .expect_err("mismatched lengths must fail");
            assert!(matches!(
                err,
                DilithiumError::Shamir(ShamirError::InconsistentShareLengths)
            ));
        }

        #[test]
        fn ensure_consistent_lengths_tracks_length() {
            let shares = vec![share_with_length(1, 4), share_with_length(2, 4)];

            let len =
                AdaptedShamirSSS::ensure_consistent_lengths::<FieldElement>(
                    &shares,
                )
                .expect("consistent lengths should succeed");
            assert_eq!(len, 4);
        }

        #[test]
        fn ensure_consistent_lengths_single_share() {
            let shares = vec![share_with_length(1, 5)];

            let len =
                AdaptedShamirSSS::ensure_consistent_lengths::<FieldElement>(
                    &shares,
                )
                .expect("single share should succeed");
            assert_eq!(len, 5);
        }

        #[test]
        fn test_invalid_threshold_config() {
            assert!(matches!(
                setup_adapted_shamir(1, 5),
                Err(DilithiumError::InvalidThreshold(1, 5))
            ));
            assert!(matches!(
                setup_adapted_shamir(0, 5),
                Err(DilithiumError::InvalidThreshold(0, 5))
            ));

            assert!(matches!(
                setup_adapted_shamir(6, 5),
                Err(DilithiumError::InvalidThreshold(6, 5))
            ));

            assert!(matches!(
                setup_adapted_shamir(3, 300),
                Err(DilithiumError::InvalidThreshold(3, 300))
            ));
        }

        #[test]
        fn test_secret_splitting() {
            let shamir = AdaptedShamirSSS::default();

            let coeffs = fe_vec!(1, 2, 3, 4, 5);
            let secret_poly1: Polynomial<'_, FieldElement> =
                poly!(coeffs.clone());
            assert_eq!(secret_poly1.coefficients().len(), coeffs.len());
            let secret_poly2: Polynomial<'_, FieldElement> =
                poly![10, 20, 30, 40, 50];
            assert_eq!(secret_poly2.coefficients().len(), coeffs.len());
            let secret_vector = poly_vec!(secret_poly1, secret_poly2);

            let shares = shamir.split_secret(&secret_vector).unwrap();

            assert_eq!(shares.len(), shamir.participant_number);

            for (i, share) in shares.iter().enumerate() {
                assert_eq!(share.participant_id, i + 1);
                assert_eq!(share.vector_length(), secret_vector.len());
            }
        }

        #[test]
        fn test_secret_reconstruction() {
            let shamir = AdaptedShamirSSS::default();

            let secret_vector: Vec<Polynomial<'_, FieldElement>> =
                vec![poly![1, 2, 3, 4, 5], poly![10, 20, 30, 40, 50]];

            let shares = shamir.split_secret(&secret_vector).unwrap();

            let reconstructed = shamir
                .reconstruct_secret(&shares[..shamir.threshold])
                .unwrap();

            assert_eq!(reconstructed, secret_vector);
        }

        #[test]
        fn test_reconstruction_with_more_shares() {
            let threshold = 3;
            let shamir = setup_adapted_shamir(threshold, 5).unwrap();

            let secret_poly: Polynomial<'_, FieldElement> = poly![42, 17, 99];
            let secret_vector = vec![secret_poly];

            let shares = shamir.split_secret(&secret_vector).unwrap();

            let reconstructed = shamir.reconstruct_secret(&shares).unwrap();
            assert_eq!(reconstructed, secret_vector);

            let reconstructed2 =
                shamir.reconstruct_secret(&shares[..threshold + 1]).unwrap();
            assert_eq!(reconstructed2, secret_vector);
        }

        #[test]
        fn test_insufficient_shares() {
            let shamir = AdaptedShamirSSS::default();

            let secret_poly: Polynomial<'_, FieldElement> = poly![100, 200];
            let secret_vector = poly_vec!(vec![secret_poly]);

            let shares = shamir.split_secret(&secret_vector).unwrap();

            let result =
                shamir.reconstruct_secret(&shares[..shamir.threshold - 1]);
            assert!(matches!(
                result,
                Err(DilithiumError::InsufficientShares(required, provided))
                    if required == shamir.threshold
                        && provided == shamir.threshold - 1
            ));

            let empty_shares: &[ShamirShare<FieldElement>] = &[];
            assert!(matches!(
                shamir.reconstruct_secret(empty_shares),
                Err(DilithiumError::InsufficientShares(required, 0))
                    if required == shamir.threshold
            ));
        }

        #[test]
        fn test_partial_reconstruction() {
            let threshold = 3;
            let shamir = setup_adapted_shamir(threshold, 5).unwrap();

            let secret_poly1: Polynomial<'_, FieldElement> =
                poly![1, 2, 3, 4, 5];
            let secret_poly2: Polynomial<'_, FieldElement> =
                poly![10, 20, 30, 40, 50];
            let secret_vector =
                poly_vec!(secret_poly1.clone(), secret_poly2.clone());

            let shares = shamir.split_secret(&secret_vector).unwrap();

            let partial = shamir
                .partial_reconstruct(&shares[..threshold], &[0])
                .unwrap();
            assert_eq!(partial.len(), 1);
            assert_eq!(partial.first(), Some(&secret_poly1));

            let partial2 = shamir
                .partial_reconstruct(&shares[..threshold], &[1])
                .unwrap();
            assert_eq!(partial2.len(), 1);
            assert_eq!(partial2.first(), Some(&secret_poly2));
        }

        #[test]
        fn test_different_vector_lengths() {
            let threshold = 3;
            let shamir = setup_adapted_shamir(threshold, 5).unwrap();

            let single_poly: Polynomial<'_, FieldElement> = poly![42, 17];
            let single_vector = vec![single_poly];
            let shares = shamir.split_secret(&single_vector).unwrap();
            let reconstructed =
                shamir.reconstruct_secret(&shares[..threshold]).unwrap();
            assert_eq!(reconstructed, single_vector);

            let poly1: Polynomial<'_, FieldElement> = poly![1, 2, 3];
            let poly2: Polynomial<'_, FieldElement> = poly![10, 20, 30];
            let poly3: Polynomial<'_, FieldElement> = poly![100, 200, 300];
            let long_vector = vec![poly1, poly2, poly3];
            let shares2 = shamir.split_secret(&long_vector).unwrap();
            let reconstructed2 =
                shamir.reconstruct_secret(&shares2[..threshold]).unwrap();
            assert_eq!(reconstructed2, long_vector);
        }

        #[test]
        fn test_zero_polynomial_handling() {
            let threshold = 2;
            let shamir = setup_adapted_shamir(threshold, 3).unwrap();

            let zero_poly: Polynomial<'_, FieldElement> =
                poly![FieldElement::default(); N];
            let zero_vector = vec![zero_poly];

            let shares: Vec<ShamirShare<'static, FieldElement>> =
                shamir.split_secret(&zero_vector).unwrap();
            let reconstructed =
                shamir.reconstruct_secret(&shares[..threshold]).unwrap();

            assert_eq!(reconstructed, zero_vector);

            let reconstructed_poly = reconstructed.first().cloned().unwrap();
            assert!(
                reconstructed_poly
                    .coefficients()
                    .iter()
                    .all(|fe| fe.is_zero())
            );
        }

        #[test]
        fn test_edge_case_thresholds() {
            let shamir_min = AdaptedShamirSSS::new(2, 2).unwrap();
            let secret_poly: Polynomial<'_, FieldElement> = poly![100, 200];
            let secret_vector = vec![secret_poly];

            let shares = shamir_min.split_secret(&secret_vector).unwrap();
            let reconstructed = shamir_min.reconstruct_secret(&shares).unwrap();
            assert_eq!(reconstructed, secret_vector);

            let shamir_large = AdaptedShamirSSS::new(10, 15).unwrap();
            let shares_large =
                shamir_large.split_secret(&secret_vector).unwrap();
            let reconstructed_large = shamir_large
                .reconstruct_secret(&shares_large[..10])
                .unwrap();
            assert_eq!(reconstructed_large, secret_vector);
        }

        #[test]
        fn test_modular_arithmetic() {
            let threshold = 2;
            let shamir = setup_adapted_shamir(threshold, 3).unwrap();

            let poly: Polynomial<'_, FieldElement> = poly![Q - 1, 1, Q - 2];
            let secret_vector = vec![poly];

            let shares = shamir.split_secret(&secret_vector).unwrap();
            let reconstructed =
                shamir.reconstruct_secret(&shares[..threshold]).unwrap();

            assert_eq!(reconstructed, secret_vector);
        }

        #[test]
        fn test_consistency_across_reconstructions() {
            let shamir = AdaptedShamirSSS::default();

            let poly: Polynomial<'_, FieldElement> = poly!(42, 17, 99, 13);
            let secret_vector = vec![poly];

            let shares = shamir.split_secret(&secret_vector).unwrap();

            let combo1 = shamir
                .reconstruct_secret(&shares[0..shamir.threshold])
                .unwrap();

            let combo2 = shamir
                .reconstruct_secret(&shares[1..shamir.threshold + 1])
                .unwrap();
            let combo3 = shamir
                .reconstruct_secret(&shares[2..shamir.threshold + 2])
                .unwrap();

            assert_eq!(combo1, secret_vector);
            assert_eq!(combo2, secret_vector);
            assert_eq!(combo3, secret_vector);
        }
    }

    #[test]
    fn test_create_shamir_polynomial() {
        let shamir = AdaptedShamirSSS::new(2, 3).unwrap();
        let secret = FieldElement::new(42u32);
        let poly = shamir.create_shamir_polynomial(&secret);

        assert_eq!(
            poly.coefficients().iter().filter(|&c| !c.is_zero()).count(),
            2
        );
        assert_eq!(poly.coefficients()[0], secret);

        let shamir = AdaptedShamirSSS::new(5, 7).unwrap();
        let secret = FieldElement::new(123);
        let poly = shamir.create_shamir_polynomial(&secret);

        assert_eq!(
            poly.coefficients().iter().filter(|&c| !c.is_zero()).count(),
            5
        );
        assert_eq!(poly.coefficients()[0], secret);

        let secret = FieldElement::zero();
        let poly = shamir.create_shamir_polynomial(&secret);
        assert_eq!(poly.coefficients()[0], FieldElement::zero());
        assert_eq!(
            poly.coefficients().iter().filter(|&c| !c.is_zero()).count(),
            4
        );

        let secret = FieldElement::new(777);
        let poly1 = shamir.create_shamir_polynomial(&secret);
        let poly2 = shamir.create_shamir_polynomial(&secret);
        let poly3 = shamir.create_shamir_polynomial(&secret);

        assert_eq!(poly1.coefficients()[0], secret);
        assert_eq!(poly2.coefficients()[0], secret);
        assert_eq!(poly3.coefficients()[0], secret);

        let coeffs_differ = (1..5)
            .any(|i| poly1.coefficients()[i] != poly2.coefficients()[i])
            || (1..5)
                .any(|i| poly1.coefficients()[i] != poly3.coefficients()[i]);
        assert!(
            coeffs_differ,
            "Random coefficients should differ between polynomials"
        );

        let shamir = AdaptedShamirSSS::new(10, 15).unwrap();
        let secret = FieldElement::new(999);
        let poly = shamir.create_shamir_polynomial(&secret);

        assert_eq!(
            poly.coefficients().iter().filter(|&c| !c.is_zero()).count(),
            10
        );
        assert_eq!(poly.coefficients()[0], secret);

        let secret = FieldElement::from(-42i32);
        let poly = shamir.create_shamir_polynomial(&secret);
        assert_eq!(poly.coefficients()[0], secret);
    }

    #[test]
    fn test_evaluate_polynomial() {
        let coeffs = vec![42];
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 0), 42);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 1), 42);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 10), 42);

        let coeffs = vec![3, 2];
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 0), 3);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 1), 5);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 2), 7);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 5), 13);

        let coeffs = vec![1, 2, 3];
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 0), 1);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 1), 6);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 2), 17);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 3), 34);

        let coeffs = vec![10, -5, 2];
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 0), 10);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 1), 7);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 2), 8);

        let coeffs = vec![5, 0, 3];
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 0), 5);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 1), 8);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 2), 17);

        let coeffs = vec![];
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 5), 0);

        let large_coeff = Q - 1;
        let coeffs = vec![large_coeff, 2];
        let result = AdaptedShamirSSS::evaluate_polynomial(&coeffs, 1);
        assert_eq!(result, (large_coeff + 2) % Q);

        let coeffs = vec![1, 1, 1];
        let large_x = 1000;
        let expected = (1 + large_x + large_x * large_x) % Q;
        assert_eq!(
            AdaptedShamirSSS::evaluate_polynomial(&coeffs, large_x),
            expected
        );

        let coeffs = vec![5, 3, 2, 4];
        let x = 7;
        let manual_result = (5 + 7 * (3 + 7 * (2 + 7 * 4))) % Q;
        assert_eq!(
            AdaptedShamirSSS::evaluate_polynomial(&coeffs, x),
            manual_result
        );

        let coeffs = vec![0, 0, 0, 0];
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 0), 0);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 100), 0);

        let coeffs = vec![0, 0, 0, 7];
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 0), 0);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 1), 7);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 2), 56);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 3), 189);

        let secret = 42;
        let coeffs = vec![secret, 17, 23, 31];
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 0), secret);
    }
}
