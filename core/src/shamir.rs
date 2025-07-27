use rand::prelude::*;

use math::{
    poly::{Polynomial, N, Q},
    poly_vector::PolynomialVector,
};

use crate::config::validate_threshold_config;
use crate::error::{Result, ThresholdError};

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

            // Compute contribution using modular multiplication to avoid overflow
            // contribution = (yi * numerator * denominator_inv) mod Q
            let contribution =
                self.mod_mul_three(yi, numerator as i32, denominator_inv);
            result = (result + contribution as i64).rem_euclid(Q as i64);
        }

        Ok(result as i32)
    }

    /// Multiply three numbers modulo Q without overflow
    fn mod_mul_three(&self, a: i32, b: i32, c: i32) -> i32 {
        // First multiply a * b mod Q
        let ab = self.mod_mul(a, b);
        // Then multiply result by c mod Q
        self.mod_mul(ab, c)
    }

    /// Multiply two numbers modulo Q without overflow
    fn mod_mul(&self, a: i32, b: i32) -> i32 {
        // Convert to i64 to prevent overflow during multiplication
        let a = a as i64;
        let b = b as i64;
        let q = Q as i64;

        // Ensure inputs are in range [0, Q)
        let a = a.rem_euclid(q);
        let b = b.rem_euclid(q);

        // TODO add Montgomery multiplication
        // For Q = 8380417, direct multiplication fits in i64
        ((a * b) % q) as i32
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shamir_split_and_reconstruct() {
        let shamir = AdaptedShamirSSS::new(2, 3).unwrap();

        // Create test polynomial vector
        let poly1 = Polynomial::new(vec![100, 200, 300]);
        let poly2 = Polynomial::new(vec![400, 500, 600]);
        let secret = PolynomialVector::new(vec![poly1, poly2]);

        // Split the secret
        let shares = shamir.split_secret(&secret).unwrap();
        assert_eq!(shares.len(), 3);

        // Reconstruct using first 2 shares
        let reconstructed = shamir.reconstruct_secret(&shares[..2]).unwrap();
        assert_eq!(secret, reconstructed);
    }

    #[test]
    fn test_invalid_threshold() {
        assert!(AdaptedShamirSSS::new(5, 3).is_err());
        assert!(AdaptedShamirSSS::new(0, 3).is_err());
    }

    #[test]
    fn test_insufficient_shares() {
        let shamir = AdaptedShamirSSS::new(3, 5).unwrap();
        let poly = Polynomial::new(vec![100]);
        let secret = PolynomialVector::new(vec![poly]);

        let shares = shamir.split_secret(&secret).unwrap();

        // Try to reconstruct with only 2 shares (need 3)
        assert!(shamir.reconstruct_secret(&shares[..2]).is_err());
    }
}
