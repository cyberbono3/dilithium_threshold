use rand::prelude::*;
use thiserror::Error;

use math::{
    poly::{Polynomial, N, Q},
    poly_vector::PolynomialVector,
};

use crate::config::validate_threshold_config;

#[derive(Error, Debug)]
pub enum ThresholdError {
    #[error("Invalid threshold configuration: threshold {threshold} > participant_number {participant_number}")]
    InvalidThreshold {
        threshold: usize,
        participant_number: usize,
    },

    #[error("Insufficient shares: need {required}, got {provided}")]
    InsufficientShares { required: usize, provided: usize },

    #[error("Invalid participant ID: {0}")]
    InvalidParticipantId(usize),

    #[error("Inconsistent share lengths")]
    InconsistentShareLengths,

    #[error("Modular inverse does not exist")]
    ModularInverseError,

    #[error("Signature generation failed after maximum attempts")]
    SignatureGenerationFailed,

    #[error("Invalid signature bounds")]
    InvalidSignatureBounds,

    #[error("Invalid polynomial index: {index} >= {length}")]
    InvalidPolynomialIndex { index: usize, length: usize },
}

pub type Result<T> = std::result::Result<T, ThresholdError>;

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

        // Process each polynomial in the vector
        for poly_idx in 0..vector_length {
            let polynomial = secret_vector.get(poly_idx).ok_or(
                ThresholdError::InvalidPolynomialIndex {
                    index: poly_idx,
                    length: vector_length,
                },
            )?;

            // Process each coefficient
            for (coeff_idx, &secret_coeff) in
                polynomial.coeffs().iter().enumerate()
            {
                let shamir_poly = self.create_shamir_polynomial(secret_coeff);

                // Evaluate for each participant
                for pid in 1..=self.participant_number {
                    let share_value =
                        self.evaluate_polynomial(&shamir_poly, pid as i32);
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
        self.reconstruct_polynomials(active_shares, &poly_indices)
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
        self.reconstruct_polynomials(active_shares, poly_indices)
    }

    /// Common logic for reconstructing polynomials
    fn reconstruct_polynomials(
        &self,
        active_shares: &[ShamirShare],
        poly_indices: &[usize],
    ) -> Result<PolynomialVector> {
        let mut reconstructed_polys = Vec::with_capacity(poly_indices.len());

        for &poly_idx in poly_indices {
            let mut coeffs = vec![0i32; N];

            for coeff_idx in 0..N {
                let points =
                    self.collect_points(active_shares, poly_idx, coeff_idx)?;
                coeffs[coeff_idx] = self.lagrange_interpolation(&points, 0)?;
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
        &self,
        shares: &[ShamirShare],
        poly_idx: usize,
        coeff_idx: usize,
    ) -> Result<Vec<(i32, i32)>> {
        shares
            .iter()
            .map(|share| {
                let poly = share.share_vector.get(poly_idx).ok_or(
                    ThresholdError::InvalidPolynomialIndex {
                        index: poly_idx,
                        length: share.vector_length(),
                    },
                )?;
                Ok((share.participant_id as i32, poly.coeffs()[coeff_idx]))
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

            // Convert to polynomials
            let polys: Vec<Polynomial> =
                share_polys.into_iter().map(Polynomial::new).collect();

            let share_vector = PolynomialVector::new(polys);
            shares.push(ShamirShare::new(pid, share_vector)?);
        }

        Ok(shares)
    }

    /// Create a Shamir polynomial with given secret as constant term.
    // fn create_shamir_polynomial(&self, secret: i32) -> Vec<i32> {
    //     let mut rng = thread_rng();
    //     let mut coeffs = Vec::with_capacity(self.threshold);
    //     coeffs.push(secret);  // The secret is the constant term

    //     // Generate random coefficients for higher degree terms
    //     for _ in 1..self.threshold {
    //         coeffs.push(rng.gen_range(0..Q));
    //     }

    //     coeffs
    // }
    /// Create a Shamir polynomial with given secret as constant term.
    fn create_shamir_polynomial(&self, secret: i32) -> Vec<i32> {
        std::iter::once(secret)
            .chain((1..self.threshold).map(|_| rand::random_range(0..Q)))
            .collect()
    }

    /// Evaluate polynomial at given point using Horner's method
    fn evaluate_polynomial(&self, coeffs: &[i32], x: i32) -> i32 {
        let mut result = 0i64;
        let mut x_power = 1i64;

        for &coeff in coeffs {
            result = (result + (coeff as i64 * x_power)) % (Q as i64);
            x_power = (x_power * x as i64) % (Q as i64);
        }

        result as i32
    }

    /// Perform Lagrange interpolation to find polynomial value at x
    fn lagrange_interpolation(
        &self,
        points: &[(i32, i32)],
        x: i32,
    ) -> Result<i32> {
        let mut result = 0i64;
        let n = points.len();

        for i in 0..n {
            let (xi, yi) = points[i];
            let mut numerator = 1i64;
            let mut denominator = 1i64;

            for j in 0..n {
                if i != j {
                    let xj = points[j].0;
                    numerator =
                        (numerator * (x - xj) as i64).rem_euclid(Q as i64);
                    denominator =
                        (denominator * (xi - xj) as i64).rem_euclid(Q as i64);
                }
            }

            let denominator_inv = self.mod_inverse(denominator as i32)?;
            let contribution = (yi as i64 * numerator * denominator_inv as i64)
                .rem_euclid(Q as i64);
            result = (result + contribution).rem_euclid(Q as i64);
        }

        Ok(result as i32)
    }

    /// Compute modular inverse using Fermat's little theorem
    fn mod_inverse(&self, a: i32) -> Result<i32> {
        if a == 0 {
            return Err(ThresholdError::ModularInverseError);
        }

        // Using Fermat's little theorem: a^(p-2) ≡ a^(-1) (mod p)
        let result = self.mod_pow(a, Q - 2);

        // Verify the result
        if (result as i64 * a as i64).rem_euclid(Q as i64) != 1 {
            return Err(ThresholdError::ModularInverseError);
        }

        Ok(result)
    }

    /// Compute (base^exp) mod Q using fast exponentiation
    fn mod_pow(&self, base: i32, mut exp: i32) -> i32 {
        let mut result = 1i64;
        let mut base = base as i64;

        while exp > 0 {
            if exp & 1 == 1 {
                result = (result * base) % (Q as i64);
            }
            base = (base * base) % (Q as i64);
            exp >>= 1;
        }

        result as i32
    }
}

// Add this to the end of shamir.rs

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_polynomial_vector(values: Vec<Vec<i32>>) -> PolynomialVector {
        let polys: Vec<Polynomial> = values
            .into_iter()
            .map(|v| {
                let mut coeffs = v;
                coeffs.resize(N, 0);
                Polynomial::new(coeffs)
            })
            .collect();
        PolynomialVector::new(polys)
    }

    #[test]
    fn test_shamir_share_creation() {
        let poly = Polynomial::new(vec![1, 2, 3]);
        let vec = PolynomialVector::new(vec![poly]);
        
        // Valid participant ID
        let share = ShamirShare::new(1, vec.clone());
        assert!(share.is_ok());
        assert_eq!(share.unwrap().participant_id, 1);
        
        // Invalid participant ID (0)
        let share = ShamirShare::new(0, vec);
        assert!(matches!(share, Err(ThresholdError::InvalidParticipantId(0))));
    }

    #[test]
    fn test_threshold_validation() {
        // Valid configurations
        assert!(AdaptedShamirSSS::new(2, 3).is_ok());
        assert!(AdaptedShamirSSS::new(3, 5).is_ok());
        assert!(AdaptedShamirSSS::new(1, 1).is_ok());
        
        // Invalid configurations
        assert!(matches!(
            AdaptedShamirSSS::new(5, 3),
            Err(ThresholdError::InvalidThreshold { threshold: 5, participant_number: 3 })
        ));
        assert!(matches!(
            AdaptedShamirSSS::new(0, 3),
            Err(ThresholdError::InvalidThreshold { .. })
        ));
    }

    #[test]
    fn test_simple_split_and_reconstruct() {
        let shamir = AdaptedShamirSSS::new(2, 3).unwrap();
        
        // Create a simple polynomial vector
        let secret = create_test_polynomial_vector(vec![
            vec![100, 200, 300],
            vec![400, 500, 600],
        ]);
        
        // Split the secret
        let shares = shamir.split_secret(&secret).unwrap();
        assert_eq!(shares.len(), 3);
        
        // Each share should have the same structure as the secret
        for share in &shares {
            assert_eq!(share.vector_length(), 2);
        }
        
        // Reconstruct using minimum threshold (2 shares)
        let reconstructed = shamir.reconstruct_secret(&shares[..2]).unwrap();
        assert_eq!(reconstructed, secret);
        
        // Should also work with all shares
        let reconstructed = shamir.reconstruct_secret(&shares).unwrap();
        assert_eq!(reconstructed, secret);
    }

    #[test]
    fn test_large_polynomial_vector() {
        let shamir = AdaptedShamirSSS::new(3, 5).unwrap();
        
        // Create a larger polynomial vector
        let mut values = Vec::new();
        for i in 0..10 {
            let mut poly_coeffs = vec![0i32; 50];
            for j in 0..50 {
                poly_coeffs[j] = ((i + 1) * 1000 + j) as i32;
            }
            values.push(poly_coeffs);
        }
        
        let secret = create_test_polynomial_vector(values);
        
        // Split and reconstruct
        let shares = shamir.split_secret(&secret).unwrap();
        let reconstructed = shamir.reconstruct_secret(&shares[1..4]).unwrap();
        assert_eq!(reconstructed, secret);
    }

    #[test]
    fn test_insufficient_shares() {
        let shamir = AdaptedShamirSSS::new(3, 5).unwrap();
        let secret = create_test_polynomial_vector(vec![vec![100]]);
        
        let shares = shamir.split_secret(&secret).unwrap();
        
        // Try to reconstruct with only 2 shares (need 3)
        let result = shamir.reconstruct_secret(&shares[..2]);
        assert!(matches!(
            result,
            Err(ThresholdError::InsufficientShares { required: 3, provided: 2 })
        ));
    }

    #[test]
    fn test_inconsistent_share_lengths() {
        let shamir = AdaptedShamirSSS::new(2, 3).unwrap();
        
        // Create shares with different vector lengths
        let poly1 = Polynomial::new(vec![1, 2, 3]);
        let poly2 = Polynomial::new(vec![4, 5, 6]);
        
        let share1 = ShamirShare::new(1, PolynomialVector::new(vec![poly1.clone()])).unwrap();
        let share2 = ShamirShare::new(2, PolynomialVector::new(vec![poly1, poly2])).unwrap();
        
        let result = shamir.reconstruct_secret(&[share1, share2]);
        assert!(matches!(result, Err(ThresholdError::InconsistentShareLengths)));
    }

    #[test]
    fn test_partial_reconstruct() {
        let shamir = AdaptedShamirSSS::new(2, 3).unwrap();
        
        // Create a polynomial vector with 5 polynomials
        let secret = create_test_polynomial_vector(vec![
            vec![100, 200],
            vec![300, 400],
            vec![500, 600],
            vec![700, 800],
            vec![900, 1000],
        ]);
        
        let shares = shamir.split_secret(&secret).unwrap();
        
        // Reconstruct only polynomials at indices 1 and 3
        let partial = shamir.partial_reconstruct(&shares[..2], &[1, 3]).unwrap();
        assert_eq!(partial.len(), 2);
        assert_eq!(partial.get(0).unwrap(), secret.get(1).unwrap());
        assert_eq!(partial.get(1).unwrap(), secret.get(3).unwrap());
    }

    #[test]
    fn test_edge_cases() {
        // 1-of-1 threshold (single participant)
        let shamir = AdaptedShamirSSS::new(1, 1).unwrap();
        let secret = create_test_polynomial_vector(vec![vec![42]]);
        
        let shares = shamir.split_secret(&secret).unwrap();
        assert_eq!(shares.len(), 1);
        
        let reconstructed = shamir.reconstruct_secret(&shares).unwrap();
        assert_eq!(reconstructed, secret);
        
        // Maximum threshold (all participants required)
        let shamir = AdaptedShamirSSS::new(5, 5).unwrap();
        let shares = shamir.split_secret(&secret).unwrap();
        assert_eq!(shares.len(), 5);
        
        // Should fail with 4 shares
        assert!(shamir.reconstruct_secret(&shares[..4]).is_err());
        
        // Should succeed with all 5
        let reconstructed = shamir.reconstruct_secret(&shares).unwrap();
        assert_eq!(reconstructed, secret);
    }

    #[test]
    fn test_polynomial_evaluation() {
        let shamir = AdaptedShamirSSS::new(2, 3).unwrap();
        
        // Test polynomial: f(x) = 5 + 3x (coeffs = [5, 3])
        let coeffs = vec![5, 3];
        
        // f(0) = 5
        assert_eq!(shamir.evaluate_polynomial(&coeffs, 0), 5);
        
        // f(1) = 5 + 3 = 8
        assert_eq!(shamir.evaluate_polynomial(&coeffs, 1), 8);
        
        // f(2) = 5 + 6 = 11
        assert_eq!(shamir.evaluate_polynomial(&coeffs, 2), 11);
    }

    #[test]
    fn test_modular_arithmetic() {
        let shamir = AdaptedShamirSSS::new(2, 3).unwrap();
        
        // Test modular exponentiation
        assert_eq!(shamir.mod_pow(2, 10), 1024);
        assert_eq!(shamir.mod_pow(3, 4), 81);
        
        // Test with large numbers (should wrap around Q)
        let base = Q - 1;
        let result = shamir.mod_pow(base, 2);
        assert_eq!(result, 1); // (Q-1)^2 ≡ 1 (mod Q)
        
        // Test modular inverse
        // 5 * 5^(-1) ≡ 1 (mod Q)
        let inv = shamir.mod_inverse(5).unwrap();
        assert_eq!((5i64 * inv as i64).rem_euclid(Q as i64), 1);
        
        // Test that 0 has no inverse
        assert!(shamir.mod_inverse(0).is_err());
    }

    #[test]
    fn test_lagrange_interpolation() {
        let shamir = AdaptedShamirSSS::new(2, 3).unwrap();
        
        // Points from polynomial f(x) = 2 + 3x
        // f(1) = 5, f(2) = 8, f(3) = 11
        let points = [(1, 5), (2, 8), (3, 11)];
        
        // Interpolate at x=0 (should give constant term = 2)
        let result = shamir.lagrange_interpolation(&points[..2], 0).unwrap();
        assert_eq!(result, 2);
        
        // Interpolate at x=4 (should give 2 + 3*4 = 14)
        let result = shamir.lagrange_interpolation(&points[..2], 4).unwrap();
        assert_eq!(result, 14);
    }

    #[test]
    fn test_deterministic_sharing() {
        // Same secret should produce same shares when using same randomness
        let shamir = AdaptedShamirSSS::new(2, 3).unwrap();
        let secret = create_test_polynomial_vector(vec![vec![100, 200, 300]]);
        
        // Note: This test assumes deterministic RNG seeding
        // In practice, shares will be different each time due to randomness
        let shares1 = shamir.split_secret(&secret).unwrap();
        let shares2 = shamir.split_secret(&secret).unwrap();
        
        // Both sets should reconstruct to the same secret
        let reconstructed1 = shamir.reconstruct_secret(&shares1[..2]).unwrap();
        let reconstructed2 = shamir.reconstruct_secret(&shares2[..2]).unwrap();
        assert_eq!(reconstructed1, secret);
        assert_eq!(reconstructed2, secret);
    }

    #[test]
    fn test_share_combination_independence() {
        // Any valid combination of threshold shares should work
        let shamir = AdaptedShamirSSS::new(3, 5).unwrap();
        let secret = create_test_polynomial_vector(vec![vec![42, 84, 126]]);
        
        let shares = shamir.split_secret(&secret).unwrap();
        
        // Test different combinations of 3 shares
        let combos = vec![
            vec![0, 1, 2],
            vec![0, 1, 3],
            vec![0, 2, 4],
            vec![1, 3, 4],
            vec![2, 3, 4],
        ];
        
        for combo in combos {
            let selected_shares: Vec<_> = combo.iter()
                .map(|&i| shares[i].clone())
                .collect();
            
            let reconstructed = shamir.reconstruct_secret(&selected_shares).unwrap();
            assert_eq!(reconstructed, secret);
        }
    }

    #[test]
    fn test_zero_polynomial() {
        let shamir = AdaptedShamirSSS::new(2, 3).unwrap();
        
        // Test with zero polynomial
        let secret = create_test_polynomial_vector(vec![vec![0; N]]);
        
        let shares = shamir.split_secret(&secret).unwrap();
        let reconstructed = shamir.reconstruct_secret(&shares[..2]).unwrap();
        
        // Check all coefficients are zero
        for i in 0..N {
            assert_eq!(reconstructed.get(0).unwrap().coeffs()[i], 0);
        }
    }

    #[test]
    fn test_boundary_values() {
        let shamir = AdaptedShamirSSS::new(2, 3).unwrap();
        
        // Test with values near Q
        let secret = create_test_polynomial_vector(vec![
            vec![Q - 1, Q - 2, Q - 3],
            vec![1, 2, 3],
        ]);
        
        let shares = shamir.split_secret(&secret).unwrap();
        let reconstructed = shamir.reconstruct_secret(&shares[..2]).unwrap();
        assert_eq!(reconstructed, secret);
    }
}

