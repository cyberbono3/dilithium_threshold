use rand::prelude::*;
use thiserror::Error;

use math::{
    poly::{Polynomial, N, Q},
    poly_vector::PolynomialVector,
};

use crate::config::*;

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
    ///  Initialize the adapted Shamir scheme.
    pub fn new(threshold: usize, participant_number: usize) -> Result<Self> {
        if !validate_threshold_config(threshold, participant_number) {
            return Err(ThresholdError::InvalidThreshold {
                threshold,
                participant_number,
            });
        }

        //   self.participant_ids = list(range(1, participants + 1))
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
        let mut participant_shares: Vec<Vec<Share>> =
            vec![Vec::new(); self.participant_number];

        // Process each polynomial in the vector
        for poly_idx in 0..vector_length {
            // TODO address it
            let polynomial = secret_vector.get(poly_idx).unwrap();

            // Process each coefficient
            for coeff_idx in 0..N {
                let secret_coeff = polynomial.coeffs()[coeff_idx];
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
        let mut shares = Vec::with_capacity(self.participant_number);
        for pid in 1..=self.participant_number {
            let mut share_polys = Vec::with_capacity(vector_length);

            for poly_idx in 0..vector_length {
                let mut coeffs = vec![0i32; N];

                // for (p_idx, c_idx, value) in participant_shares[pid - 1] {
                for share in &participant_shares[pid - 1] {
                    if share.poly_idx == poly_idx {
                        coeffs[share.coeff_idx] = share.value;
                    }
                }

                share_polys.push(Polynomial::new(coeffs));
            }

            let share_vector = PolynomialVector::new(share_polys);
            shares.push(ShamirShare::new(pid, share_vector)?);
        }

        Ok(shares)
    }

    #[cfg(test)]
    /// Partially reconstruct only specified polynomials from the vector.
    // TODO fix code duplications
    pub fn partial_reconstruct(
        &self,
        shares: &[ShamirShare],
        poly_indices: &[usize],
    ) -> Result<PolynomialVector> {
        if shares.len() < self.threshold {
            return Err(ThresholdError::InsufficientShares {
                required: self.threshold,
                provided: shares.len(),
            });
        }

        let active_shares = &shares[..self.threshold];
        let mut reconstructed_polys = Vec::with_capacity(poly_indices.len());

        for poly_index in poly_indices {
            let mut coeffs = vec![0i32; N];
            for coeff_idx in 0..N {
                let points: Vec<(i32, i32)> = active_shares
                    .iter()
                    .map(|share| {
                        (
                            share.participant_id as i32,
                            share
                                .share_vector
                                .get(*poly_index)
                                .unwrap()
                                .coeffs()[coeff_idx],
                        )
                    })
                    .collect();

                coeffs[coeff_idx] = self.lagrange_interpolation(&points, 0)?;
            }
            reconstructed_polys.push(Polynomial::from(coeffs))
        }

        Ok(PolynomialVector::new(reconstructed_polys))
    }

    /// Reconstruct the secret from a sufficient number of shares.
    pub fn reconstruct_secret(
        &self,
        shares: &[ShamirShare],
    ) -> Result<PolynomialVector> {
        if shares.len() < self.threshold {
            return Err(ThresholdError::InsufficientShares {
                required: self.threshold,
                provided: shares.len(),
            });
        }

        let active_shares = &shares[..self.threshold];
        let vector_length = active_shares[0].vector_length();

        // Verify all shares have same length
        if active_shares
            .iter()
            .any(|s| s.vector_length() != vector_length)
        {
            return Err(ThresholdError::InconsistentShareLengths);
        }

        let mut reconstructed_polys = Vec::with_capacity(vector_length);

        for poly_idx in 0..vector_length {
            let mut coeffs = vec![0i32; N];

            for coeff_idx in 0..N {
                let points: Vec<(i32, i32)> = active_shares
                    .iter()
                    .map(|share| {
                        (
                            share.participant_id as i32,
                            share.share_vector.get(poly_idx).unwrap().coeffs()
                                [coeff_idx],
                        )
                    })
                    .collect();

                coeffs[coeff_idx] = self.lagrange_interpolation(&points, 0)?;
            }

            reconstructed_polys.push(Polynomial::new(coeffs));
        }

        Ok(PolynomialVector::new(reconstructed_polys))
    }

    /// Create a Shamir polynomial with given secret as constant term.
    fn create_shamir_polynomial(&self, secret: i32) -> Vec<i32> {
        let mut rng = rand::rng();
        let mut coeffs = vec![secret; self.threshold - 1];

        for _ in 1..self.threshold {
            coeffs.push(rng.random::<u32>().rem_euclid(Q as u32) as i32);
        }

        coeffs
    }

    fn evaluate_polynomial(&self, coeffs: &[i32], x: i32) -> i32 {
        let mut result = 0i64;
        let mut x_power = 1i64;

        for &coeff in coeffs {
            result = (result + (coeff as i64 * x_power)) % (Q as i64);
            x_power = (x_power * x as i64) % (Q as i64);
        }

        result as i32
    }

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
            let contribution =
                ((yi as i64 * numerator * denominator_inv as i64) % (Q as i64))
                    as i64;
            result = (result + contribution).rem_euclid(Q as i64);
        }

        Ok(result as i32)
    }

    /// Compute modular inverse of a modulo m using extended Euclidean algorithm.
    fn mod_inverse(&self, a: i32) -> Result<i32> {
        // Using Fermat's little theorem: a^(p-2) ≡ a^(-1) (mod p)
        let result = self.mod_pow(a, Q - 2);
        if (result as i64 * a as i64).rem_euclid(Q as i64) != 1 {
            return Err(ThresholdError::ModularInverseError);
        }
        Ok(result)
    }

    /// Compute modular inverse of a modulo m using extended Euclidean algorithm.
    fn mod_pow(&self, base: i32, exp: i32) -> i32 {
        let mut result = 1i64;
        let mut base = base as i64;
        let mut exp = exp;

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
