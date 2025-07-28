use rand::prelude::*;
use sha2::{Digest, Sha256};
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};
use std::collections::HashMap;

use crate::{
    config::{
        validate_threshold_config, DilithiumConfig, DEFAULT_SECURITY_LEVEL,
    },
    error::{Result, ThresholdError},
    shamir::{AdaptedShamirSSS, ShamirShare},
    utils::{get_randomness, hash_message},
};

use math::{
    poly::{Polynomial, N, Q},
    poly_vector::PolynomialVector,
};

use crate::dilithium::{Dilithium, DilithiumPublicKey, DilithiumSignature};

/// Represents a participant's share of the threshold key.
///
/// Contains both the Shamir share of the secret key and the
/// participant's identification information.
#[derive(Clone, Debug)]
pub struct ThresholdKeyShare {
    pub participant_id: usize,
    pub s1_share: ShamirShare,
    pub s2_share: ShamirShare,
    pub public_key: DilithiumPublicKey,
}

impl ThresholdKeyShare {
    /// Initialize threshold key share.
    pub fn new(
        participant_id: usize,
        s1_share: ShamirShare,
        s2_share: ShamirShare,
        public_key: DilithiumPublicKey,
    ) -> Self {
        ThresholdKeyShare {
            participant_id,
            s1_share,
            s2_share,
            public_key,
        }
    }
}

impl std::fmt::Display for ThresholdKeyShare {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "ThresholdKeyShare(id={})", self.participant_id)
    }
}

/// Represents a partial signature from one participant.
///
/// Contains the participant's contribution to the threshold signature
/// along with verification information.
#[derive(Clone, Debug)]
pub struct PartialSignature {
    pub participant_id: usize,
    pub z_partial: PolynomialVector,
    pub commitment: PolynomialVector,
    pub challenge: Polynomial,
}

impl PartialSignature {
    /// Initialize partial signature.
    pub fn new(
        participant_id: usize,
        z_partial: PolynomialVector,
        commitment: PolynomialVector,
        challenge: Polynomial,
    ) -> Self {
        Self {
            participant_id,
            z_partial,
            commitment,
            challenge,
        }
    }
}

impl std::fmt::Display for PartialSignature {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "PartialSignature(id={})", self.participant_id)
    }
}

/// Main threshold signature scheme implementation.
///
/// Provides distributed key generation, partial signing, and signature
/// combination functionality while preventing secret leakage.
pub struct ThresholdSignature {
    threshold: usize,
    participants: usize,
    security_level: usize,
    dilithium: Dilithium,
    shamir_s1: AdaptedShamirSSS,
    shamir_s2: AdaptedShamirSSS,
    participant_ids: Vec<usize>,
    config: DilithiumConfig,
}

impl ThresholdSignature {
    /// Initialize threshold signature scheme.
    pub fn new(
        threshold: usize,
        participants: usize,
        security_level: Option<usize>,
    ) -> Result<Self> {
        if !validate_threshold_config(threshold, participants) {
            return Err(ThresholdError::InvalidThreshold {
                threshold,
                participant_number: participants,
            });
        }

        let security_level = security_level.unwrap_or(DEFAULT_SECURITY_LEVEL);
        let config = DilithiumConfig::new(security_level);

        // Initialize underlying schemes
        let dilithium = Dilithium::new(security_level);
        let shamir_s1 = AdaptedShamirSSS::new(threshold, participants)?;
        let shamir_s2 = AdaptedShamirSSS::new(threshold, participants)?;

        // Store participant IDs
        let participant_ids: Vec<usize> = (1..=participants).collect();

        Ok(ThresholdSignature {
            threshold,
            participants,
            security_level,
            dilithium,
            shamir_s1,
            shamir_s2,
            participant_ids,
            config,
        })
    }

    /// Generate threshold keys using distributed key generation.
    ///
    /// This method generates a Dilithium key pair and then splits the
    /// private key into shares using the adapted Shamir scheme.
    pub fn distributed_keygen(
        &self,
        seed: Option<&[u8]>,
    ) -> Result<Vec<ThresholdKeyShare>> {
        // Generate base Dilithium key pair
        let key_pair = self.dilithium.keygen(seed);

        // Split secret vectors using adapted Shamir scheme
        let s1_shares =
            self.shamir_s1.split_secret(&key_pair.private_key.s1)?;
        let s2_shares =
            self.shamir_s2.split_secret(&key_pair.private_key.s2)?;

        let threshold_shares = (0..self.participants)
            .map(|i| {
                ThresholdKeyShare::new(
                    self.participant_ids[i],
                    s1_shares[i].clone(),
                    s2_shares[i].clone(),
                    key_pair.public_key.clone(),
                )
            })
            .collect();

        Ok(threshold_shares)
    }

    /// Create a partial signature using a key share.
    ///
    /// This is the core innovation: each participant can create a partial
    /// signature using only their share, without reconstructing the full
    /// secret key.
    pub fn partial_sign(
        &self,
        message: &[u8],
        key_share: &ThresholdKeyShare,
        randomness: Option<&[u8]>,
    ) -> Result<PartialSignature> {
        let randomness = get_randomness(randomness);

        // Hash message
        let mu = hash_message(message);

        // Generate participant-specific randomness
        let participant_randomness = self.derive_participant_randomness(
            &randomness,
            key_share.participant_id,
        );

        // Sample mask vector y (participant's portion)
        let y_partial = self.sample_partial_y(&participant_randomness);

        // Compute partial commitment w_partial = A * y_partial
        // Note: This is simplified - in practice, we need coordination
        // between participants to compute the full commitment
        let w_partial = Self::compute_partial_commitment(
            &key_share.public_key.m,
            &y_partial,
        );

        // For now, use a simplified challenge generation
        // In practice, this requires coordination between participants
        let challenge = self.generate_partial_challenge(&mu, &w_partial);

        // Compute partial response z_partial = y_partial + c * s1_share
        let c_s1 = key_share.s1_share.share_vector.clone() * challenge;
        let z_partial = c_s1 + y_partial;

        Ok(PartialSignature::new(
            key_share.participant_id,
            z_partial,
            w_partial,
            challenge,
        ))
    }

    /// Combine partial signatures into a complete threshold signature.
    ///
    /// This method reconstructs the full signature from partial signatures
    /// without ever reconstructing the secret key.
    pub fn combine_signatures(
        &self,
        partial_signatures: &[PartialSignature],
        public_key: &DilithiumPublicKey,
    ) -> Result<DilithiumSignature> {
        if partial_signatures.len() < self.threshold {
            return Err(ThresholdError::InsufficientShares {
                required: self.threshold,
                provided: partial_signatures.len(),
            });
        }

        // Verify all partial signatures use the same challenge
        let challenge = &partial_signatures[0].challenge;
        if !partial_signatures
            .iter()
            .all(|ps| &ps.challenge == challenge)
        {
            return Err(ThresholdError::PartialSignatureChallengeMismatch);
        }

        // Use first threshold partial signatures
        let active_partials = &partial_signatures[..self.threshold];

        // Reconstruct z vector using Lagrange interpolation
        let z = self.reconstruct_z_vector(active_partials)?;

        // Reconstruct hint h (simplified for this implementation)
        let h = self.reconstruct_hint(active_partials, public_key)?;

        Ok(DilithiumSignature::new(z, h, *challenge))
    }

    /// Verify a partial signature.
    // pub fn verify_partial_signature(
    //     &self,
    //     message: &[u8],
    //     partial_sig: &PartialSignature,
    //     key_share: &ThresholdKeyShare,
    // ) -> bool {
    //     // Hash message
    //     let mu = hash_message(message);

    //     // Verify challenge consistency
    //     let expected_challenge =
    //         self.generate_partial_challenge(&mu, &partial_sig.commitment);

    //     if partial_sig.challenge != expected_challenge {
    //         return false;
    //     }

    //     // Verify partial signature bounds
    //     self.check_partial_bounds(partial_sig)
    // }

    /// Derive participant-specific randomness.
    fn derive_participant_randomness(
        &self,
        base_randomness: &[u8],
        participant_id: usize,
    ) -> Vec<u8> {
        let mut hasher = Sha256::new();
        Digest::update(&mut hasher, base_randomness);
        Digest::update(&mut hasher, participant_id.to_le_bytes());
        Digest::update(&mut hasher, b"participant_randomness");
        hasher.finalize().to_vec()
    }

    /// Sample partial mask vector y.
    fn sample_partial_y(&self, randomness: &[u8]) -> PolynomialVector {
        let mut polys = Vec::with_capacity(self.config.l);

        for i in 0..self.config.l {
            let mut seed = randomness.to_vec();
            seed.push(i as u8);
            let coeffs = self.sample_gamma1(&seed);
            polys.push(Polynomial::from(coeffs));
        }

        PolynomialVector::new(polys)
    }

    /// Sample coefficients from gamma1 distribution
    fn sample_gamma1(&self, seed: &[u8]) -> Vec<i32> {
        let gamma1 = self.config.gamma1;
        let mut rng = StdRng::from_seed({
            let mut hasher = Sha256::new();
            Digest::update(&mut hasher, seed);
            let hash = hasher.finalize();
            let mut seed_bytes = [0u8; 32];
            seed_bytes.copy_from_slice(&hash);
            seed_bytes
        });

        let mut coeffs = vec![0i32; N];
        for i in 0..N {
            coeffs[i] = rng.gen_range(-gamma1..=gamma1);
        }
        coeffs
    }

    /// Compute partial commitment w_partial.
    fn compute_partial_commitment(
        m: &Vec<Vec<Polynomial>>,
        y_partial: &PolynomialVector,
    ) -> PolynomialVector {
        // Simplified partial commitment computation
        m * y_partial
    }

    /// Generate challenge polynomial (simplified version).
    /// TODO update it
    fn generate_partial_challenge(
        &self,
        mu: &[u8],
        w_partial: &PolynomialVector,
    ) -> Polynomial {
        // Simplified challenge generation
        let mut hasher = Shake256::default();
        hasher.update(mu);

        // Hash w_partial
        for i in 0..w_partial.len() {
            if let Some(poly) = w_partial.get(i) {
                for coeff in poly.coeffs() {
                    hasher.update(&coeff.to_le_bytes());
                }
            }
        }

        let mut reader = hasher.finalize_xof();
        let mut challenge_bytes = vec![0u8; 32];
        reader.read(&mut challenge_bytes);

        // Convert to polynomial with tau non-zero coefficients
        self.sample_challenge(&challenge_bytes)
    }

    /// Sample challenge polynomial
    fn sample_challenge(&self, seed: &[u8]) -> Polynomial {
        let mut rng = StdRng::from_seed({
            let mut hasher = Sha256::new();
            Digest::update(&mut hasher, seed);
            let hash = hasher.finalize();
            let mut seed_bytes = [0u8; 32];
            seed_bytes.copy_from_slice(&hash);
            seed_bytes
        });

        let mut coeffs = vec![0i32; N];
        let mut indices: Vec<usize> = (0..N).collect();
        indices.shuffle(&mut rng);

        // Set tau coefficients to ±1
        for i in 0..self.config.tau {
            coeffs[indices[i]] = if rng.random_bool(0.5) { 1 } else { -1 };
        }

        Polynomial::from(coeffs)
    }

    /// Reconstruct z vector from partial signatures using Lagrange interpolation.
    /// TODO update it
    fn reconstruct_z_vector(
        &self,
        partial_signatures: &[PartialSignature],
    ) -> Result<PolynomialVector> {
        if partial_signatures.is_empty() {
            return Err(ThresholdError::InsufficientShares {
                required: 1,
                provided: 0,
            });
        }

        // Get vector length from first partial signature
        let vector_length = partial_signatures[0].z_partial.len();

        // Reconstruct each polynomial in the vector
        let mut reconstructed_polys = Vec::with_capacity(vector_length);

        for poly_idx in 0..vector_length {
            // Reconstruct each coefficient of this polynomial
            let mut coeffs = vec![0i32; N];

            for coeff_idx in 0..N {
                // Collect points for Lagrange interpolation
                let mut points = Vec::with_capacity(partial_signatures.len());

                for ps in partial_signatures {
                    let x = ps.participant_id as i32;
                    let y = ps
                        .z_partial
                        .get(poly_idx)
                        .ok_or(ThresholdError::InvalidPolynomialIndex {
                            index: poly_idx,
                            length: vector_length,
                        })?
                        .coeffs()[coeff_idx];
                    points.push((x, y));
                }

                // Perform Lagrange interpolation
                let reconstructed_coeff =
                    self.lagrange_interpolation(&points, 0)?;
                coeffs[coeff_idx] = reconstructed_coeff.rem_euclid(Q);
            }

            reconstructed_polys.push(Polynomial::from(coeffs));
        }

        Ok(PolynomialVector::new(reconstructed_polys))
    }

    /// Reconstruct hint vector (simplified implementation).
    /// TODO make onliner
    fn reconstruct_hint(
        &self,
        _partial_signatures: &[PartialSignature],
        _public_key: &DilithiumPublicKey,
    ) -> Result<PolynomialVector> {
        // Simplified hint reconstruction
        // In practice, this would involve more complex coordination
        let mut hint_polys = Vec::with_capacity(self.config.k);
        for _ in 0..self.config.k {
            hint_polys.push(Polynomial::zero());
        }

        Ok(PolynomialVector::new(hint_polys))
    }

    /// Perform Lagrange interpolation.
    fn lagrange_interpolation(
        &self,
        points: &[(i32, i32)],
        x: i32,
    ) -> Result<i32> {
        let mut result = 0i64;
        let n = points.len();

        for i in 0..n {
            let (xi, yi) = points[i];

            // Compute Lagrange basis polynomial L_i(x)
            let mut numerator = 1i64;
            let mut denominator = 1i64;

            for j in 0..n {
                if i != j {
                    let (xj, _) = points[j];
                    numerator =
                        (numerator * (x - xj) as i64).rem_euclid(Q as i64);
                    denominator =
                        (denominator * (xi - xj) as i64).rem_euclid(Q as i64);
                }
            }

            // Compute modular inverse using Fermat's little theorem
            let denominator_inv = self.mod_pow(denominator as i32, Q - 2);

            // Add contribution
            let contribution =
                self.mod_mul_three(yi, numerator as i32, denominator_inv);
            result = (result + contribution as i64).rem_euclid(Q as i64);
        }

        Ok(result as i32)
    }

    /// Compute (base^exp) mod Q using fast exponentiation
    fn mod_pow(&self, base: i32, mut exp: i32) -> i32 {
        let mut result = 1i64;
        let mut base = base as i64;
        let q = Q as i64;

        while exp > 0 {
            if exp & 1 == 1 {
                result = (result * base) % q;
            }
            base = (base * base) % q;
            exp >>= 1;
        }

        result as i32
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
        let a = a as i64;
        let b = b as i64;
        let q = Q as i64;

        let a = a.rem_euclid(q);
        let b = b.rem_euclid(q);

        ((a * b) % q) as i32
    }

    /// Check if partial signature satisfies bound requirements.
    // fn check_partial_bounds(&self, partial_sig: &PartialSignature) -> bool {
    //     let gamma1 = self.dilithium.config.gamma1;
    //     let beta = self.dilithium.config.beta;

    //     partial_sig.z_partial.norm_infinity() < gamma1 - beta
    // }

    /// Get information about the threshold configuration.
    pub fn get_threshold_info(&self) -> HashMap<&'static str, usize> {
        let mut info = HashMap::new();
        info.insert("threshold", self.threshold);
        info.insert("participants", self.participants);
        info.insert("security_level", self.security_level);
        info.insert("min_signers", self.threshold);
        info.insert("max_participants", self.participants);
        info
    }
}

// TODO add comprehensive testing
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_threshold_signature_creation() {
        let threshold_sig = ThresholdSignature::new(3, 5, None).unwrap();
        let info = threshold_sig.get_threshold_info();

        assert_eq!(info.get("threshold"), Some(&3));
        assert_eq!(info.get("participants"), Some(&5));
    }

    #[test]
    fn test_invalid_threshold() {
        assert!(ThresholdSignature::new(6, 5, None).is_err());
        assert!(ThresholdSignature::new(0, 5, None).is_err());
    }
}
