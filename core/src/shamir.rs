use math::{
    polynomial::{Polynomial, N, Q},
    poly_vector::PolynomialVector,
};

use crate::config::validate_threshold_config;
use crate::error::{Result, ThresholdError};
use crate::utils::lagrange_interpolation;

/// Adapted Shamir's Secret Sharing
#[derive(Clone, Debug, PartialEq)]
pub struct ShamirShare {
    pub participant_id: usize,
    pub share_vector: PolynomialVector,
}

impl ShamirShare {
    pub fn new(
        participant_id: usize,
        share_vector: PolynomialVector,
    ) -> Result<Self> {
        if participant_id == 0 {
            return Err(ThresholdError::InvalidParticipantId(participant_id));
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

#[derive(Clone, Copy, Debug)]
struct Share {
    poly_idx: usize,
    coeff_idx: usize,
    value: i32,
}

impl Share {
    pub fn new(poly_idx: usize, coeff_idx: usize, value: i32) -> Self {
        Self {
            poly_idx,
            coeff_idx,
            value,
        }
    }
}

#[derive(Debug)]
pub struct AdaptedShamirSSS {
    threshold: usize,
    participant_number: usize,
}

impl AdaptedShamirSSS {
    /// Initialize the adapted Shamir scheme.
    pub fn new(threshold: usize, participant_number: usize) -> Result<Self> {
        if !validate_threshold_config(threshold, participant_number) {
            return Err(ThresholdError::InvalidThreshold {
                threshold,
                participant_number,
            });
        }

        Ok(AdaptedShamirSSS {
            threshold,
            participant_number,
        })
    }

    /// Split a polynomial vector secret into shares.
    pub fn split_secret(
        &self,
        secret_vector: &PolynomialVector,
    ) -> Result<Vec<ShamirShare>> {
        let vector_length = secret_vector.len();
        let mut participant_shares: Vec<Vec<Share>> = vec![
                Vec::with_capacity(vector_length * N);
                self.participant_number
            ];

        for (poly_idx, poly) in secret_vector.as_slice().iter().enumerate() {
            // Process each coefficient
            for (coeff_idx, &secret_coeff) in poly.coeffs().iter().enumerate() {
                let shamir_poly = self.create_shamir_polynomial(secret_coeff);

                // Evaluate for each participant
                for pid in 1..=self.participant_number {
                    let share_value =
                        Self::evaluate_polynomial(&shamir_poly, pid as i32);
                    let share = Share::new(poly_idx, coeff_idx, share_value);
                    participant_shares[pid - 1].push(share);
                }
            }
        }

        // Organize into shares
        self.organize_shares(participant_shares, vector_length)
    }

    /// Reconstruct the secret from a sufficient number of shares.
    pub fn reconstruct_secret(
        &self,
        shares: &[ShamirShare],
    ) -> Result<PolynomialVector> {
        self.validate_shares(shares)?;

        let active_shares = &shares[..self.threshold];
        let vector_length = active_shares[0].vector_length();

        let poly_indices: Vec<usize> = (0..vector_length).collect();
        Self::reconstruct_polynomials(active_shares, &poly_indices)
    }

    #[cfg(test)]
    /// Partially reconstruct only specified polynomials from the vector.
    pub fn partial_reconstruct(
        &self,
        shares: &[ShamirShare],
        poly_indices: &[usize],
    ) -> Result<PolynomialVector> {
        self.validate_shares(shares)?;

        let active_shares = &shares[..self.threshold];
        Self::reconstruct_polynomials(active_shares, poly_indices)
    }

    /// Common logic for reconstructing polynomials
    fn reconstruct_polynomials(
        active_shares: &[ShamirShare],
        poly_indices: &[usize],
    ) -> Result<PolynomialVector> {
        let mut reconstructed_polys = Vec::with_capacity(poly_indices.len());

        for poly_idx in poly_indices {
            let mut coeffs = vec![0i32; N];

            for (i, c) in coeffs.iter_mut().enumerate().take(N) {
                let points = Self::collect_points(active_shares, *poly_idx, i)?;
                *c = lagrange_interpolation(&points, 0)?;
            }

            reconstructed_polys.push(Polynomial::from(coeffs));
        }

        Ok(PolynomialVector::new(reconstructed_polys))
    }

    /// Validate shares for reconstruction
    fn validate_shares(&self, shares: &[ShamirShare]) -> Result<()> {
        if shares.len() < self.threshold {
            return Err(ThresholdError::InsufficientShares {
                required: self.threshold,
                provided: shares.len(),
            });
        }

        // Check if all shares have the same vector length
        if let Some(first) = shares.first() {
            let expected_length = first.vector_length();
            if shares.iter().any(|s| s.vector_length() != expected_length) {
                return Err(ThresholdError::InconsistentShareLengths);
            }
        }

        Ok(())
    }

    /// Collect points for Lagrange interpolation
    fn collect_points(
        shares: &[ShamirShare],
        poly_idx: usize,
        coeff_idx: usize,
    ) -> Result<Vec<(i32, i32)>> {
        shares
            .iter()
            .map(|share| {
                let poly = share.share_vector.get(poly_idx).ok_or(
                    ThresholdError::InvalidIndex {
                        index: poly_idx,
                        length: share.vector_length(),
                    },
                )?;
                let coeff = poly.coeffs().get(coeff_idx).copied().ok_or(
                    ThresholdError::InvalidIndex {
                        index: coeff_idx,
                        length: share.vector_length(),
                    },
                )?;
                Ok((share.participant_id as i32, coeff))
            })
            .collect()
    }

    /// Organize participant shares into ShamirShare objects
    fn organize_shares(
        &self,
        participant_shares: Vec<Vec<Share>>,
        vector_length: usize,
    ) -> Result<Vec<ShamirShare>> {
        let mut shares = Vec::with_capacity(self.participant_number);

        for pid in 1..=self.participant_number {
            let mut share_polys = vec![vec![0i32; N]; vector_length];

            // Fill in the coefficients for each polynomial
            for share in &participant_shares[pid - 1] {
                share_polys[share.poly_idx][share.coeff_idx] = share.value;
            }

            let share_vector = PolynomialVector::new(
                share_polys.into_iter().map(Polynomial::new).collect(),
            );
            shares.push(ShamirShare::new(pid, share_vector)?);
        }

        Ok(shares)
    }

    /// Create a Shamir polynomial with given secret as constant term.
    fn create_shamir_polynomial(&self, secret: i32) -> Vec<i32> {
        std::iter::once(secret)
            .chain((1..self.threshold).map(|_| rand::random_range(0..Q)))
            .collect()
    }

    /// Evaluate polynomial at given point using Horner's method
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
    use math::{
        polynomial::{Polynomial, N, Q},
        poly_vector::PolynomialVector,
    };

    mod shamir_share_tests {
        use super::*;

        #[test]
        fn test_share_creation() {
            let poly1 = Polynomial::from(vec![1, 2, 3, 4, 5]);
            let poly2 = Polynomial::from(vec![6, 7, 8, 9, 10]);
            let share_vector = PolynomialVector::new(vec![poly1, poly2]);
            let share = ShamirShare::new(1, share_vector.clone()).unwrap();

            assert_eq!(share.participant_id, 1);
            assert_eq!(share.vector_length(), 2);
            assert_eq!(share.share_vector, share_vector);
        }

        #[test]
        fn test_invalid_participant_id() {
            let poly1 = Polynomial::from(vec![1, 2, 3, 4, 5]);
            let share_vector = PolynomialVector::new(vec![poly1]);

            // Test that participant ID 0 is rejected
            assert!(ShamirShare::new(0, share_vector.clone()).is_err());
        }

        #[test]
        fn test_share_debug_representation() {
            let poly1 = Polynomial::from(vec![1, 2, 3]);
            let share_vector = PolynomialVector::new(vec![poly1]);
            let share = ShamirShare::new(1, share_vector).unwrap();

            let debug_str = format!("{:?}", share);
            assert!(debug_str.contains("ShamirShare"));
            assert!(debug_str.contains("participant_id: 1"));
        }
    }

    mod adapted_shamir_sss_tests {
        use super::*;

        #[test]
        fn test_shamir_initialization() {
            let threshold = 3;
            let participants = 5;
            let shamir =
                AdaptedShamirSSS::new(threshold, participants).unwrap();

            // We can't directly access fields in Rust due to privacy,
            // but we can test that creation succeeds
            assert!(AdaptedShamirSSS::new(threshold, participants).is_ok());
        }

        #[test]
        fn test_invalid_threshold_config() {
            // Threshold too small
            assert!(matches!(
                AdaptedShamirSSS::new(1, 5),
                Err(ThresholdError::InvalidThreshold {
                    threshold: 1,
                    participant_number: 5
                })
            ));
            assert!(matches!(
                AdaptedShamirSSS::new(0, 5),
                Err(ThresholdError::InvalidThreshold {
                    threshold: 0,
                    participant_number: 5
                })
            ));

            // Threshold larger than participants
            assert!(matches!(
                AdaptedShamirSSS::new(6, 5),
                Err(ThresholdError::InvalidThreshold {
                    threshold: 6,
                    participant_number: 5
                })
            ));

            // Too many participants (assuming max is 255 based on typical limits)
            assert!(matches!(
                AdaptedShamirSSS::new(3, 300),
                Err(ThresholdError::InvalidThreshold {
                    threshold: 3,
                    participant_number: 300
                })
            ));
        }

        #[test]
        fn test_secret_splitting() {
            let threshold = 3;
            let participants = 5;
            let shamir =
                AdaptedShamirSSS::new(threshold, participants).unwrap();

            // Create test secret vector
            let secret_poly1 = Polynomial::from(vec![1, 2, 3, 4, 5]);
            let secret_poly2 = Polynomial::from(vec![10, 20, 30, 40, 50]);
            let secret_vector =
                PolynomialVector::new(vec![secret_poly1, secret_poly2]);

            let shares = shamir.split_secret(&secret_vector).unwrap();

            // Check we get correct number of shares
            assert_eq!(shares.len(), participants);

            // Check each share has correct participant ID
            for (i, share) in shares.iter().enumerate() {
                assert_eq!(share.participant_id, i + 1);
                assert_eq!(share.vector_length(), 2);
            }
        }

        #[test]
        fn test_secret_reconstruction() {
            let threshold = 3;
            let participants = 5;
            let shamir =
                AdaptedShamirSSS::new(threshold, participants).unwrap();

            // Create test secret vector
            let secret_poly1 = Polynomial::from(vec![1, 2, 3, 4, 5]);
            let secret_poly2 = Polynomial::from(vec![10, 20, 30, 40, 50]);
            let secret_vector =
                PolynomialVector::new(vec![secret_poly1, secret_poly2]);

            // Split secret
            let shares = shamir.split_secret(&secret_vector).unwrap();

            // Reconstruct using exactly threshold shares
            let reconstructed =
                shamir.reconstruct_secret(&shares[..threshold]).unwrap();

            // Check reconstruction is correct
            assert_eq!(reconstructed, secret_vector);
        }

        #[test]
        fn test_reconstruction_with_more_shares() {
            let threshold = 3;
            let participants = 5;
            let shamir =
                AdaptedShamirSSS::new(threshold, participants).unwrap();

            let secret_poly = Polynomial::from(vec![42, 17, 99]);
            let secret_vector = PolynomialVector::new(vec![secret_poly]);

            let shares = shamir.split_secret(&secret_vector).unwrap();

            // Use all shares
            let reconstructed = shamir.reconstruct_secret(&shares).unwrap();
            assert_eq!(reconstructed, secret_vector);

            // Use threshold + 1 shares
            let reconstructed2 =
                shamir.reconstruct_secret(&shares[..threshold + 1]).unwrap();
            assert_eq!(reconstructed2, secret_vector);
        }

        #[test]
        fn test_insufficient_shares() {
            let threshold = 3;
            let participants = 5;
            let shamir =
                AdaptedShamirSSS::new(threshold, participants).unwrap();

            let secret_poly = Polynomial::from(vec![100, 200]);
            let secret_vector = PolynomialVector::new(vec![secret_poly]);

            let shares = shamir.split_secret(&secret_vector).unwrap();

            // Try with threshold - 1 shares
            assert!(shamir
                .reconstruct_secret(&shares[..threshold - 1])
                .is_err());

            // Try with empty slice
            let empty_shares: &[ShamirShare] = &[];
            assert!(shamir.reconstruct_secret(empty_shares).is_err());
        }

        #[test]
        fn test_partial_reconstruction() {
            let threshold = 3;
            let participants = 5;
            let shamir =
                AdaptedShamirSSS::new(threshold, participants).unwrap();

            // Create test secret vector
            let secret_poly1 = Polynomial::from(vec![1, 2, 3, 4, 5]);
            let secret_poly2 = Polynomial::from(vec![10, 20, 30, 40, 50]);
            let secret_vector = PolynomialVector::new(vec![
                secret_poly1.clone(),
                secret_poly2.clone(),
            ]);

            let shares = shamir.split_secret(&secret_vector).unwrap();

            // Reconstruct only first polynomial
            let partial = shamir
                .partial_reconstruct(&shares[..threshold], &[0])
                .unwrap();
            assert_eq!(partial.len(), 1);
            assert_eq!(partial.get(0).unwrap(), &secret_poly1);

            // Reconstruct only second polynomial
            let partial2 = shamir
                .partial_reconstruct(&shares[..threshold], &[1])
                .unwrap();
            assert_eq!(partial2.len(), 1);
            assert_eq!(partial2.get(0).unwrap(), &secret_poly2);
        }

        #[test]
        fn test_different_vector_lengths() {
            let threshold = 3;
            let participants = 5;
            let shamir =
                AdaptedShamirSSS::new(threshold, participants).unwrap();

            // Test with single polynomial
            let single_poly = Polynomial::from(vec![42, 17]);
            let single_vector = PolynomialVector::new(vec![single_poly]);
            let shares = shamir.split_secret(&single_vector).unwrap();
            let reconstructed =
                shamir.reconstruct_secret(&shares[..threshold]).unwrap();
            assert_eq!(reconstructed, single_vector);

            // Test with longer vector
            let poly1 = Polynomial::from(vec![1, 2, 3]);
            let poly2 = Polynomial::from(vec![10, 20, 30]);
            let poly3 = Polynomial::from(vec![100, 200, 300]);
            let long_vector = PolynomialVector::new(vec![poly1, poly2, poly3]);
            let shares2 = shamir.split_secret(&long_vector).unwrap();
            let reconstructed2 =
                shamir.reconstruct_secret(&shares2[..threshold]).unwrap();
            assert_eq!(reconstructed2, long_vector);
        }

        #[test]
        fn test_zero_polynomial_handling() {
            let threshold = 2;
            let participants = 3;
            let shamir =
                AdaptedShamirSSS::new(threshold, participants).unwrap();

            // Create zero polynomial
            let zero_poly = Polynomial::new(vec![0; N]);
            let zero_vector = PolynomialVector::new(vec![zero_poly]);

            let shares = shamir.split_secret(&zero_vector).unwrap();
            let reconstructed =
                shamir.reconstruct_secret(&shares[..threshold]).unwrap();

            assert_eq!(reconstructed, zero_vector);

            // Check all coefficients are zero
            let reconstructed_poly = reconstructed.get(0).unwrap();
            assert!(reconstructed_poly.coeffs().iter().all(|&c| c == 0));
        }

        #[test]
        fn test_random_polynomial_reconstruction() {
            use rand::{prelude::*, rng};

            let threshold = 3;
            let participants = 5;
            let shamir =
                AdaptedShamirSSS::new(threshold, participants).unwrap();

            // Generate random polynomial vector
            let mut random_polys = Vec::new();

            for _ in 0..3 {
                let mut coeffs = vec![0i32; 10]; // Use only first 10 coefficients
                for c in coeffs.iter_mut() {
                    *c = rng().random_range(0..1000);
                }
                random_polys.push(Polynomial::from(coeffs));
            }

            let random_vector = PolynomialVector::new(random_polys);

            // Split and reconstruct
            let shares = shamir.split_secret(&random_vector).unwrap();
            let reconstructed =
                shamir.reconstruct_secret(&shares[..threshold]).unwrap();

            assert_eq!(reconstructed, random_vector);
        }

        #[test]
        fn test_edge_case_thresholds() {
            // Minimum threshold (2 out of 2)
            let shamir_min = AdaptedShamirSSS::new(2, 2).unwrap();
            let secret_poly = Polynomial::from(vec![100, 200]);
            let secret_vector = PolynomialVector::new(vec![secret_poly]);

            let shares = shamir_min.split_secret(&secret_vector).unwrap();
            let reconstructed = shamir_min.reconstruct_secret(&shares).unwrap();
            assert_eq!(reconstructed, secret_vector);

            // Large threshold
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
            // Since mod_inverse and other methods are private, we test them indirectly
            // through the reconstruction process which relies on correct modular arithmetic

            let threshold = 2;
            let participants = 3;
            let shamir =
                AdaptedShamirSSS::new(threshold, participants).unwrap();

            // Use specific values that test modular arithmetic
            let poly = Polynomial::from(vec![Q - 1, 1, Q - 2]); // Values near modulus
            let secret_vector = PolynomialVector::new(vec![poly]);

            let shares = shamir.split_secret(&secret_vector).unwrap();
            let reconstructed =
                shamir.reconstruct_secret(&shares[..threshold]).unwrap();

            assert_eq!(reconstructed, secret_vector);
        }

        #[test]
        fn test_consistency_across_reconstructions() {
            let threshold = 3;
            let participants = 5;
            let shamir =
                AdaptedShamirSSS::new(threshold, participants).unwrap();

            let poly = Polynomial::from(vec![42, 17, 99, 13]);
            let secret_vector = PolynomialVector::new(vec![poly]);

            let shares = shamir.split_secret(&secret_vector).unwrap();

            // Try different combinations of threshold shares
            let combo1 =
                shamir.reconstruct_secret(&shares[0..threshold]).unwrap();
            let combo2 = shamir
                .reconstruct_secret(&shares[1..threshold + 1])
                .unwrap();
            let combo3 = shamir
                .reconstruct_secret(&shares[2..threshold + 2])
                .unwrap();

            assert_eq!(combo1, secret_vector);
            assert_eq!(combo2, secret_vector);
            assert_eq!(combo3, secret_vector);
        }
    }

    // Integration tests
    #[test]
    fn test_full_workflow() {
        let threshold = 5;
        let participants = 9;
        let shamir = AdaptedShamirSSS::new(threshold, participants).unwrap();

        // Create a complex secret
        let poly1 = Polynomial::from(vec![123, 456, 789, 101, 202]);
        let poly2 = Polynomial::from(vec![303, 404, 505, 606, 707]);
        let poly3 = Polynomial::from(vec![808, 909, 111, 222, 333]);
        let secret = PolynomialVector::new(vec![poly1, poly2, poly3]);

        // Split the secret
        let shares = shamir.split_secret(&secret).unwrap();

        // Verify we can reconstruct with any threshold shares
        for start in 0..=(participants - threshold) {
            let selected_shares = &shares[start..start + threshold];
            let reconstructed =
                shamir.reconstruct_secret(selected_shares).unwrap();
            assert_eq!(reconstructed, secret);
        }
    }

    #[test]
    fn test_create_shamir_polynomial() {
        // Test with threshold 2 (minimal case)
        let shamir = AdaptedShamirSSS::new(2, 3).unwrap();
        let secret = 42;
        let poly = shamir.create_shamir_polynomial(secret);

        // Check polynomial has correct length (equal to threshold)
        assert_eq!(poly.len(), 2);
        // Check first coefficient is the secret
        assert_eq!(poly[0], secret);
        // Check second coefficient is in valid range [0, Q)
        assert!(poly[1] >= 0 && poly[1] < Q);

        // Test with threshold 5
        let shamir = AdaptedShamirSSS::new(5, 7).unwrap();
        let secret = 123;
        let poly = shamir.create_shamir_polynomial(secret);

        // Check polynomial has correct length
        assert_eq!(poly.len(), 5);
        // Check first coefficient is the secret
        assert_eq!(poly[0], secret);
        // Check all random coefficients are in valid range
        for i in 1..5 {
            assert!(poly[i] >= 0 && poly[i] < Q);
        }

        // Test with zero secret
        let secret = 0;
        let poly = shamir.create_shamir_polynomial(secret);
        assert_eq!(poly[0], 0);
        assert_eq!(poly.len(), 5);

        // Test with maximum secret value
        let secret = Q - 1;
        let poly = shamir.create_shamir_polynomial(secret);
        assert_eq!(poly[0], secret);
        assert_eq!(poly.len(), 5);

        // Test randomness - create multiple polynomials with same secret
        let secret = 777;
        let poly1 = shamir.create_shamir_polynomial(secret);
        let poly2 = shamir.create_shamir_polynomial(secret);
        let poly3 = shamir.create_shamir_polynomial(secret);

        // All should have same secret
        assert_eq!(poly1[0], secret);
        assert_eq!(poly2[0], secret);
        assert_eq!(poly3[0], secret);

        // But random coefficients should differ (with high probability)
        // Check at least one differs between poly1 and poly2
        let coeffs_differ = (1..5).any(|i| poly1[i] != poly2[i])
            || (1..5).any(|i| poly1[i] != poly3[i]);
        assert!(
            coeffs_differ,
            "Random coefficients should differ between polynomials"
        );

        // Test with large threshold
        let shamir = AdaptedShamirSSS::new(10, 15).unwrap();
        let secret = 999;
        let poly = shamir.create_shamir_polynomial(secret);

        assert_eq!(poly.len(), 10);
        assert_eq!(poly[0], secret);
        for i in 1..10 {
            assert!(poly[i] >= 0 && poly[i] < Q);
        }

        // Test negative secret
        let secret = -42;
        let poly = shamir.create_shamir_polynomial(secret);
        assert_eq!(poly[0], secret);

        // Statistical test for randomness distribution
        // Create many polynomials and check coefficient distribution
        let shamir = AdaptedShamirSSS::new(3, 5).unwrap();
        let mut coeff_sum = 0i64;
        let iterations = 100;

        for _ in 0..iterations {
            let poly = shamir.create_shamir_polynomial(100);
            coeff_sum += poly[1] as i64 + poly[2] as i64;
        }

        // Average should be roughly Q/2 for uniform distribution
        let avg = coeff_sum / (iterations * 2);
        let expected_avg = Q as i64 / 2;
        // Allow 20% deviation
        assert!(
            (avg - expected_avg).abs() < (expected_avg / 5),
            "Random coefficient average {} deviates too much from expected {}",
            avg,
            expected_avg
        );
    }

    #[test]
    fn test_evaluate_polynomial() {
        // Test constant polynomial: f(x) = 42
        let coeffs = vec![42];
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 0), 42);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 1), 42);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 10), 42);

        // Test linear polynomial: f(x) = 3 + 2x
        let coeffs = vec![3, 2];
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 0), 3);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 1), 5);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 2), 7);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 5), 13);

        // Test quadratic polynomial: f(x) = 1 + 2x + 3x^2
        let coeffs = vec![1, 2, 3];
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 0), 1);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 1), 6); // 1 + 2 + 3
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 2), 17); // 1 + 4 + 12
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 3), 34); // 1 + 6 + 27

        // Test with negative coefficients: f(x) = 10 - 5x + 2x^2
        let coeffs = vec![10, -5, 2];
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 0), 10);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 1), 7); // 10 - 5 + 2
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 2), 8); // 10 - 10 + 8

        // Test with zero coefficients: f(x) = 5 + 0x + 3x^2
        let coeffs = vec![5, 0, 3];
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 0), 5);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 1), 8); // 5 + 0 + 3
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 2), 17); // 5 + 0 + 12

        // Test empty polynomial
        let coeffs = vec![];
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 5), 0);

        // Test modular arithmetic with large coefficients
        let large_coeff = Q - 1;
        let coeffs = vec![large_coeff, 2];
        let result = AdaptedShamirSSS::evaluate_polynomial(&coeffs, 1);
        assert_eq!(result, (large_coeff + 2) % Q);

        // Test with large x values
        let coeffs = vec![1, 1, 1]; // f(x) = 1 + x + x^2
        let large_x = 1000;
        let expected = (1 + large_x + large_x * large_x) % Q;
        assert_eq!(
            AdaptedShamirSSS::evaluate_polynomial(&coeffs, large_x),
            expected as i32
        );

        // Verify Horner's method: f(x) = 5 + 3x + 2x^2 + 4x^3
        let coeffs = vec![5, 3, 2, 4];
        let x = 7;
        let manual_result = (5 + 7 * (3 + 7 * (2 + 7 * 4))) % Q;
        assert_eq!(
            AdaptedShamirSSS::evaluate_polynomial(&coeffs, x),
            manual_result as i32
        );

        // Test all-zero polynomial
        let coeffs = vec![0, 0, 0, 0];
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 0), 0);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 100), 0);

        // Test single high-degree term: f(x) = 7x^3
        let coeffs = vec![0, 0, 0, 7];
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 0), 0);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 1), 7);
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 2), 56); // 7 * 8
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 3), 189); // 7 * 27

        // Test Shamir property: evaluation at x=0 returns the secret
        let secret = 42;
        let coeffs = vec![secret, 17, 23, 31];
        assert_eq!(AdaptedShamirSSS::evaluate_polynomial(&coeffs, 0), secret);
    }
}
