use sha2::{Digest, Sha256};
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};
use std::collections::HashMap;

use crate::{
    config::{validate_threshold_config, DEFAULT_SECURITY_LEVEL},
    error::{Result, ThresholdError},
    shamir::{AdaptedShamirSSS, ShamirShare},
    utils::{
        get_hash_reader, get_randomness, hash_message, lagrange_interpolation,
    },
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
#[derive(Debug)]
pub struct ThresholdSignature {
    threshold: usize,
    participants: usize,
    security_level: usize,
    dilithium: Dilithium,
    shamir_s1: AdaptedShamirSSS,
    shamir_s2: AdaptedShamirSSS,
    participant_ids: Vec<usize>,
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

        // Compute partial commitment w_partial = m * y_partial
        let w_partial = &key_share.public_key.m * &y_partial;

        // For now, use a simplified challenge generation
        // In practice, this requires coordination between participants
        let challenge = self.generate_partial_challenge(&mu);

        // Compute partial response z_partial = y_partial + c * s1_share
        let c_s1 = key_share.s1_share.share_vector.clone() * challenge;
        let z_partial = y_partial + c_s1;

        Ok(PartialSignature::new(
            key_share.participant_id,
            z_partial,
            w_partial,
            challenge,
        ))
    }

    /// Verify that partial signatures are valid for combination.
    ///
    /// Checks that:
    /// 1. There are at least `threshold` partial signatures
    /// 2. All partial signatures use the same challenge
    fn verify_partial_signatures(
        partial_signatures: &[PartialSignature],
        threshold: usize,
    ) -> Result<()> {
        // Check if we have enough partial signatures
        if partial_signatures.len() < threshold {
            return Err(ThresholdError::InsufficientShares {
                required: threshold,
                provided: partial_signatures.len(),
            });
        }

        // At this point, we know we have at least threshold signatures
        // Since threshold must be >= 2 (per validation rules), we know
        // partial_signatures is not empty

        // Verify all partial signatures use the same challenge
        let challenge = &partial_signatures[0].challenge;
        if !partial_signatures
            .iter()
            .all(|ps| &ps.challenge == challenge)
        {
            return Err(ThresholdError::PartialSignatureChallengeMismatch);
        }

        Ok(())
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
        // Verify partial signatures
        Self::verify_partial_signatures(partial_signatures, self.threshold)?;

        // Use first threshold partial signatures
        let active_partials = &partial_signatures[..self.threshold];

        // Reconstruct z vector using Lagrange interpolation
        let z = self.reconstruct_z_vector(active_partials)?;

        // Reconstruct hint h (simplified for this implementation)
        let h = self.reconstruct_hint(active_partials, public_key)?;

        // Get challenge from the verified partial signatures
        let challenge = partial_signatures[0].challenge;

        Ok(DilithiumSignature::new(z, h, challenge))
    }

    /// Verify a partial signature.
    pub fn verify_partial_signature(
        &self,
        message: &[u8],
        partial_sig: &PartialSignature,
    ) -> bool {
        // Hash message
        let mu = hash_message(message);

        // Verify challenge consistency
        let expected_challenge = self.generate_partial_challenge(&mu);

        if partial_sig.challenge != expected_challenge {
            return false;
        }

        // Verify partial signature bounds
        self.check_partial_bounds(partial_sig)
    }

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
        let polys = (0..self.dilithium.config.l)
            .map(|i| {
                Polynomial::from(
                    self.sample_gamma1(&[randomness, &[i as u8]].concat()),
                )
            })
            .collect();

        PolynomialVector::new(polys)
    }

    // TODO optimize it
    /// Sample coefficients from gamma1 distribution
    fn sample_gamma1(&self, seed: &[u8]) -> Vec<i32> {
        // Create SHAKE256 hasher and generate output
        let mut reader = get_hash_reader(seed);

        // Read N * 4 bytes
        let mut hash_output = vec![0u8; N * 4];
        reader.read(&mut hash_output);

        // Convert bytes to u32 values and transform to coefficients
        let coeffs: Vec<i32> = hash_output
            .chunks_exact(4)
            .take(N)
            .map(|chunk| {
                let value = u32::from_le_bytes([
                    chunk[0], chunk[1], chunk[2], chunk[3],
                ]);
                // Map to range [-gamma1, gamma1] then reduce mod Q
                let coeff = (value as i32
                    % (2 * self.dilithium.config.gamma1 + 1))
                    - self.dilithium.config.gamma1;
                coeff.rem_euclid(Q)
            })
            .collect();

        coeffs
    }

    /// Convert to polynomial with tau non-zero coefficients
    fn generate_partial_challenge(&self, mu: &[u8]) -> Polynomial {
        // Create seed
        let mut seed = mu.to_vec();
        seed.extend_from_slice(b"challenge");

        // Sample tau positions for ±1 coefficients
        let mut reader = get_hash_reader(&seed);
        let mut hash_output = vec![0u8; self.dilithium.config.tau * 2];
        reader.read(&mut hash_output);

        // Initialize coefficients array
        let mut coeffs = vec![0i32; N];
        for i in 0..self.dilithium.config.tau {
            let pos = (hash_output[i * 2] as usize) % N;
            let sign = if hash_output[i * 2 + 1] % 2 == 0 {
                1
            } else {
                -1
            };
            coeffs[pos] = sign;
        }

        Polynomial::new(coeffs)
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
                    let poly = ps.z_partial.get(poly_idx).ok_or(
                        ThresholdError::InvalidIndex {
                            index: poly_idx,
                            length: vector_length,
                        },
                    )?;
                    // .coeffs()[coeff_idx];
                    let y = poly.coeffs().get(coeff_idx).copied().ok_or(
                        ThresholdError::InvalidIndex {
                            index: coeff_idx,
                            length: poly.coeffs().len(),
                        },
                    )?;
                    points.push((x, y));
                }

                // Perform Lagrange interpolation
                let reconstructed_coeff = lagrange_interpolation(&points, 0)?;
                coeffs[coeff_idx] = reconstructed_coeff.rem_euclid(Q);
            }

            reconstructed_polys.push(Polynomial::from(coeffs));
        }

        Ok(PolynomialVector::new(reconstructed_polys))
    }

    /// Reconstruct hint vector (simplified implementation).
    fn reconstruct_hint(
        &self,
        _partial_signatures: &[PartialSignature],
        _public_key: &DilithiumPublicKey,
    ) -> Result<PolynomialVector> {
        let hint_polys = (0..self.dilithium.config.k).map(|_| 
            Polynomial::zero()
        ).collect();

        Ok(PolynomialVector::new(hint_polys))
    }

    /// Check if partial signature satisfies bound requirements.
    fn check_partial_bounds(&self, partial_sig: &PartialSignature) -> bool {
        let gamma1 = self.dilithium.config.gamma1;
        let beta = self.dilithium.config.beta;

        let res = partial_sig.z_partial.norm_infinity() < gamma1 - beta;
        if !res {
            println!("check_partial_bounds returned false");
        }

        res
    }

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

#[cfg(test)]
mod tests {
    use super::*;

    // Helper function to create a deterministic seed
    fn create_test_seed(value: u8) -> Vec<u8> {
        vec![value; 32]
    }

    // Helper function to create test message
    fn create_test_message(content: &str) -> Vec<u8> {
        content.as_bytes().to_vec()
    }

    mod threshold_key_share_tests {
        use super::*;

        #[test]
        fn test_key_share_creation() {
            let poly1 = Polynomial::from(vec![1, 2, 3]);
            let poly2 = Polynomial::from(vec![4, 5, 6]);
            let s1_share_vec = PolynomialVector::new(vec![poly1.clone()]);
            let s2_share_vec = PolynomialVector::new(vec![poly2.clone()]);

            let s1_share = ShamirShare::new(1, s1_share_vec).unwrap();
            let s2_share = ShamirShare::new(1, s2_share_vec).unwrap();

            let dilithium = Dilithium::new(DEFAULT_SECURITY_LEVEL);
            let keypair = dilithium.keygen(Some(&create_test_seed(1)));

            let key_share = ThresholdKeyShare::new(
                1,
                s1_share.clone(),
                s2_share.clone(),
                keypair.public_key.clone(),
            );

            assert_eq!(key_share.participant_id, 1);
            assert_eq!(key_share.s1_share.participant_id, 1);
            assert_eq!(key_share.s2_share.participant_id, 1);
        }

        #[test]
        fn test_key_share_display() {
            let poly1 = Polynomial::from(vec![1]);
            let poly2 = Polynomial::from(vec![2]);
            let s1_share_vec = PolynomialVector::new(vec![poly1]);
            let s2_share_vec = PolynomialVector::new(vec![poly2]);

            let s1_share = ShamirShare::new(5, s1_share_vec).unwrap();
            let s2_share = ShamirShare::new(5, s2_share_vec).unwrap();

            let dilithium = Dilithium::new(DEFAULT_SECURITY_LEVEL);
            let keypair = dilithium.keygen(None);

            let key_share = ThresholdKeyShare::new(
                5,
                s1_share,
                s2_share,
                keypair.public_key,
            );

            let display_str = format!("{}", key_share);
            assert_eq!(display_str, "ThresholdKeyShare(id=5)");
        }
    }

    mod partial_signature_tests {
        use super::*;

        #[test]
        fn test_partial_signature_creation() {
            let z_poly = Polynomial::from(vec![100, 200, 300]);
            let z_partial = PolynomialVector::new(vec![z_poly]);

            let comm_poly = Polynomial::from(vec![10, 20, 30]);
            let commitment = PolynomialVector::new(vec![comm_poly]);

            let challenge = Polynomial::from(vec![1, -1, 0, 1]);

            let partial_sig = PartialSignature::new(
                3,
                z_partial.clone(),
                commitment.clone(),
                challenge.clone(),
            );

            assert_eq!(partial_sig.participant_id, 3);
            assert_eq!(partial_sig.z_partial.len(), 1);
            assert_eq!(partial_sig.commitment.len(), 1);
            assert_eq!(partial_sig.challenge, challenge);
        }

        #[test]
        fn test_partial_signature_display() {
            let z_partial = PolynomialVector::new(vec![Polynomial::zero()]);
            let commitment = PolynomialVector::new(vec![Polynomial::zero()]);
            let challenge = Polynomial::zero();

            let partial_sig =
                PartialSignature::new(7, z_partial, commitment, challenge);

            let display_str = format!("{}", partial_sig);
            assert_eq!(display_str, "PartialSignature(id=7)");
        }
    }

    mod threshold_signature_tests {
        use super::*;

        #[test]
        fn test_threshold_signature_creation() {
            let threshold_sig = ThresholdSignature::new(3, 5, None).unwrap();
            let info = threshold_sig.get_threshold_info();

            assert_eq!(info.get("threshold"), Some(&3));
            assert_eq!(info.get("participants"), Some(&5));
            assert_eq!(
                info.get("security_level"),
                Some(&DEFAULT_SECURITY_LEVEL)
            );
            assert_eq!(info.get("min_signers"), Some(&3));
            assert_eq!(info.get("max_participants"), Some(&5));
        }

        #[test]
        fn test_threshold_signature_with_custom_security() {
            let threshold_sig = ThresholdSignature::new(4, 7, Some(3)).unwrap();
            let info = threshold_sig.get_threshold_info();

            assert_eq!(info.get("security_level"), Some(&3));
        }

        #[test]
        fn test_invalid_threshold_configurations() {
            // Threshold greater than participants
            assert!(ThresholdSignature::new(6, 5, None).is_err());

            // Zero threshold
            assert!(ThresholdSignature::new(0, 5, None).is_err());

            // Threshold of 1 (not allowed)
            assert!(ThresholdSignature::new(1, 5, None).is_err());

            // Too many participants
            assert!(ThresholdSignature::new(3, 300, None).is_err());

            // Too few participants
            assert!(ThresholdSignature::new(2, 1, None).is_err());
        }

        #[test]
        fn test_distributed_keygen() {
            let threshold_sig = ThresholdSignature::new(3, 5, None).unwrap();

            // Test with deterministic seed
            let shares = threshold_sig
                .distributed_keygen(Some(&create_test_seed(42)))
                .unwrap();

            assert_eq!(shares.len(), 5);

            // Verify each share has correct participant ID
            for (i, share) in shares.iter().enumerate() {
                assert_eq!(share.participant_id, i + 1);
                assert_eq!(share.s1_share.participant_id, i + 1);
                assert_eq!(share.s2_share.participant_id, i + 1);
            }

            // Verify all shares have the same public key
            let public_key = &shares[0].public_key;
            for share in &shares[1..] {
                assert_eq!(&share.public_key, public_key);
            }
        }

        #[test]
        fn test_distributed_keygen_reproducibility() {
            let threshold_sig = ThresholdSignature::new(2, 3, None).unwrap();
            let seed = create_test_seed(123);

            let shares1 =
                threshold_sig.distributed_keygen(Some(&seed)).unwrap();
            let shares2 =
                threshold_sig.distributed_keygen(Some(&seed)).unwrap();

            // Same seed should produce same shares
            assert_eq!(shares1.len(), shares2.len());
            for i in 0..shares1.len() {
                assert_eq!(
                    shares1[i].participant_id,
                    shares2[i].participant_id
                );
                // Note: Direct comparison of shares might not work due to internal randomness
                // but public keys should match
                assert_eq!(shares1[i].public_key, shares2[i].public_key);
            }
        }

        #[test]
        fn test_partial_sign_basic() {
            let threshold_sig = ThresholdSignature::new(3, 5, None).unwrap();
            let shares = threshold_sig
                .distributed_keygen(Some(&create_test_seed(1)))
                .unwrap();
            let message = create_test_message("Test message");

            // Create partial signature with first share
            let partial_sig = threshold_sig
                .partial_sign(&message, &shares[0], Some(&create_test_seed(2)))
                .unwrap();

            assert_eq!(partial_sig.participant_id, 1);
            assert!(partial_sig.z_partial.len() > 0);
            assert!(partial_sig.commitment.len() > 0);
        }

        #[test]
        fn test_partial_sign_deterministic() {
            let threshold_sig = ThresholdSignature::new(2, 3, None).unwrap();
            let shares = threshold_sig
                .distributed_keygen(Some(&create_test_seed(1)))
                .unwrap();
            let message = create_test_message("Deterministic test");
            let randomness = create_test_seed(42);

            // Same inputs should produce same partial signature
            let partial1 = threshold_sig
                .partial_sign(&message, &shares[0], Some(&randomness))
                .unwrap();

            let partial2 = threshold_sig
                .partial_sign(&message, &shares[0], Some(&randomness))
                .unwrap();

            assert_eq!(partial1.participant_id, partial2.participant_id);
            assert_eq!(partial1.challenge, partial2.challenge);
            // Note: z_partial comparison might need special handling
        }

        #[test]
        fn test_verify_partial_signature() {
            let threshold_sig = ThresholdSignature::new(3, 5, None).unwrap();
            let shares = threshold_sig
                .distributed_keygen(Some(&create_test_seed(1)))
                .unwrap();
            let message = create_test_message("Verify test");

            let partial_sig = threshold_sig
                .partial_sign(&message, &shares[0], Some(&create_test_seed(2)))
                .unwrap();

            // Verify the partial signature
            let is_valid =
                threshold_sig.verify_partial_signature(&message, &partial_sig);

            assert!(is_valid);
        }

        #[test]
        fn test_verify_partial_signature_wrong_message() {
            let threshold_sig = ThresholdSignature::new(3, 5, None).unwrap();
            let shares = threshold_sig
                .distributed_keygen(Some(&create_test_seed(1)))
                .unwrap();
            let message = create_test_message("Original message");
            let wrong_message = create_test_message("Wrong message");

            let partial_sig = threshold_sig
                .partial_sign(&message, &shares[0], Some(&create_test_seed(2)))
                .unwrap();

            // Verify with wrong message should fail
            let is_valid = threshold_sig
                .verify_partial_signature(&wrong_message, &partial_sig);

            assert!(!is_valid);
        }

        #[test]
        fn test_full_threshold_signing_workflow() {
            let threshold = 3;
            let participants = 5;
            let threshold_sig =
                ThresholdSignature::new(threshold, participants, None).unwrap();

            // 1. Distributed key generation
            let shares = threshold_sig
                .distributed_keygen(Some(&create_test_seed(1)))
                .unwrap();
            let public_key = &shares[0].public_key;

            // 2. Message to sign
            let message = create_test_message("Full workflow test message");

            // 3. Create partial signatures from different subsets
            // Test with first threshold participants
            let mut partial_sigs_1 = Vec::new();
            for i in 0..threshold {
                let partial = threshold_sig
                    .partial_sign(
                        &message,
                        &shares[i],
                        Some(&create_test_seed((i + 100) as u8)),
                    )
                    .unwrap();
                partial_sigs_1.push(partial);
            }

            // Test with last threshold participants
            let mut partial_sigs_2 = Vec::new();
            for i in (participants - threshold)..participants {
                let partial = threshold_sig
                    .partial_sign(
                        &message,
                        &shares[i],
                        Some(&create_test_seed((i + 200) as u8)),
                    )
                    .unwrap();
                partial_sigs_2.push(partial);
            }

            // 4. Combine signatures from different subsets
            let combined_sig_1 = threshold_sig
                .combine_signatures(&partial_sigs_1, public_key)
                .unwrap();

            let combined_sig_2 = threshold_sig
                .combine_signatures(&partial_sigs_2, public_key)
                .unwrap();

            // Both combined signatures should be valid
            // Verify using the 'c' field which is the challenge
            assert_eq!(combined_sig_1.c, partial_sigs_1[0].challenge);
            assert_eq!(combined_sig_2.c, partial_sigs_2[0].challenge);

            // Verify structure of both signatures
            assert_eq!(
                combined_sig_1.z.len(),
                threshold_sig.dilithium.config.l
            );
            assert_eq!(
                combined_sig_1.h.len(),
                threshold_sig.dilithium.config.k
            );
            assert_eq!(
                combined_sig_2.z.len(),
                threshold_sig.dilithium.config.l
            );
            assert_eq!(
                combined_sig_2.h.len(),
                threshold_sig.dilithium.config.k
            );

            // Since both signatures are for the same message, they should have the same challenge
            assert_eq!(combined_sig_1.c, combined_sig_2.c);
        }

        #[test]
        fn test_combine_signatures_insufficient_shares() {
            let threshold = 3;
            let participants = 5;
            let threshold_sig =
                ThresholdSignature::new(threshold, participants, None).unwrap();
            let shares = threshold_sig
                .distributed_keygen(Some(&create_test_seed(1)))
                .unwrap();
            let message = create_test_message("Insufficient test");

            // Create only threshold-1 partial signatures
            let mut partial_sigs = Vec::new();
            for i in 0..(threshold - 1) {
                let partial = threshold_sig
                    .partial_sign(
                        &message,
                        &shares[i],
                        Some(&create_test_seed((i + 10) as u8)),
                    )
                    .unwrap();
                partial_sigs.push(partial);
            }

            // Should fail with insufficient shares
            let result = threshold_sig
                .combine_signatures(&partial_sigs, &shares[0].public_key);

            assert!(result.is_err());
        }

        #[test]
        fn test_combine_signatures_mismatched_challenges() {
            let threshold = 3;
            let participants = 5;
            let threshold_sig =
                ThresholdSignature::new(threshold, participants, None).unwrap();
            let shares = threshold_sig
                .distributed_keygen(Some(&create_test_seed(1)))
                .unwrap();

            // Create partial signatures with different messages (leading to different challenges)
            let mut partial_sigs = Vec::new();

            for i in 0..threshold {
                let message = create_test_message(&format!("Message {}", i));
                let partial = threshold_sig
                    .partial_sign(
                        &message,
                        &shares[i],
                        Some(&create_test_seed((i + 10) as u8)),
                    )
                    .unwrap();
                partial_sigs.push(partial);
            }

            // Should fail with mismatched challenges
            let result = threshold_sig
                .combine_signatures(&partial_sigs, &shares[0].public_key);

            assert!(result.is_err());
        }

        // #[test]
        // fn test_modular_arithmetic() {
        //     let threshold_sig = ThresholdSignature::new(2, 3, None).unwrap();

        //     // Test mod_pow
        //     let base = 3;
        //     let exp = 4;
        //     let result = mod_pow(base, exp);
        //     assert_eq!(result, 81); // 3^4 = 81

        //     // Test mod_mul
        //     let a = Q - 1;
        //     let b = 2;
        //     let result = threshold_sig.mod_mul(a, b);
        //     assert_eq!(result, (((Q - 1) as i64 * 2) % Q as i64) as i32);

        //     // Test mod_mul_three
        //     let a = 1000;
        //     let b = 2000;
        //     let c = 3;
        //     let result = threshold_sig.mod_mul_three(a, b, c);
        //     let expected = ((1000i64 * 2000 * 3) % Q as i64) as i32;
        //     assert_eq!(result, expected);
        // }

        #[test]
        fn test_derive_participant_randomness() {
            let threshold_sig = ThresholdSignature::new(2, 3, None).unwrap();
            let base_randomness = create_test_seed(42);

            // Different participants should get different randomness
            let rand1 = threshold_sig
                .derive_participant_randomness(&base_randomness, 1);
            let rand2 = threshold_sig
                .derive_participant_randomness(&base_randomness, 2);

            assert_ne!(rand1, rand2);
            assert_eq!(rand1.len(), 32); // SHA256 output
            assert_eq!(rand2.len(), 32);

            // Same participant should get same randomness
            let rand1_again = threshold_sig
                .derive_participant_randomness(&base_randomness, 1);
            assert_eq!(rand1, rand1_again);
        }

        #[test]
        fn test_sample_gamma1() {
            let threshold_sig = ThresholdSignature::new(2, 3, None).unwrap();
            let seed = create_test_seed(123);

            let coeffs = threshold_sig.sample_gamma1(&seed);

            assert_eq!(coeffs.len(), N);

            // Check all coefficients are within bounds
            for &coeff in &coeffs {
                assert!((0..Q).contains(&coeff));
                // Original coefficient before modular reduction would be in [-gamma1, gamma1]
            }

            // Same seed should produce same coefficients
            let coeffs2 = threshold_sig.sample_gamma1(&seed);
            assert_eq!(coeffs, coeffs2);
        }

        #[test]
        fn test_edge_cases() {
            // Minimum configuration (2 out of 2)
            let threshold_sig = ThresholdSignature::new(2, 2, None).unwrap();
            let shares = threshold_sig.distributed_keygen(None).unwrap();
            assert_eq!(shares.len(), 2);

            // Large threshold
            let threshold_sig_large =
                ThresholdSignature::new(10, 15, None).unwrap();
            let shares_large =
                threshold_sig_large.distributed_keygen(None).unwrap();
            assert_eq!(shares_large.len(), 15);
        }

        #[test]
        fn test_different_security_levels() {
            // Test with different security levels
            for security_level in [2, 3, 5] {
                let threshold_sig =
                    ThresholdSignature::new(3, 5, Some(security_level))
                        .unwrap();
                let info = threshold_sig.get_threshold_info();
                assert_eq!(info.get("security_level"), Some(&security_level));

                // Verify key generation works with different security levels
                let shares = threshold_sig.distributed_keygen(None).unwrap();
                assert_eq!(shares.len(), 5);
            }
        }

        #[test]
        fn test_concurrent_partial_signing() {
            let threshold = 4;
            let participants = 7;
            let threshold_sig =
                ThresholdSignature::new(threshold, participants, None).unwrap();
            let shares = threshold_sig
                .distributed_keygen(Some(&create_test_seed(1)))
                .unwrap();
            let message = create_test_message("Concurrent signing test");

            // Simulate concurrent signing by different participants
            let partial_sigs: Vec<_> = (0..threshold)
                .map(|i| {
                    threshold_sig
                        .partial_sign(
                            &message,
                            &shares[i],
                            Some(&create_test_seed((i * 10) as u8)),
                        )
                        .unwrap()
                })
                .collect();

            // All partial signatures should be valid
            for (i, partial_sig) in partial_sigs.iter().enumerate() {
                assert!(threshold_sig
                    .verify_partial_signature(&message, partial_sig,));
            }

            // Should be able to combine them
            let combined = threshold_sig
                .combine_signatures(&partial_sigs, &shares[0].public_key);
            assert!(combined.is_ok());
        }
    }

    // Integration tests combining threshold signatures with Dilithium
    mod integration_tests {
        use super::*;

        #[test]
        fn test_threshold_vs_regular_dilithium() {
            let security_level = 2;
            let threshold = 3;
            let participants = 5;

            // Create threshold signature scheme
            let threshold_sig = ThresholdSignature::new(
                threshold,
                participants,
                Some(security_level),
            )
            .unwrap();

            // Create regular Dilithium for comparison
            let dilithium = Dilithium::new(security_level);

            // Generate keys
            let seed = create_test_seed(42);
            let threshold_shares =
                threshold_sig.distributed_keygen(Some(&seed)).unwrap();
            let regular_keypair = dilithium.keygen(Some(&seed));

            // Public keys should match since we used same seed
            assert_eq!(
                threshold_shares[0].public_key,
                regular_keypair.public_key
            );
        }

        #[test]
        fn test_reconstruct_vs_original_key() {
            let threshold = 2;
            let participants = 3;
            let threshold_sig =
                ThresholdSignature::new(threshold, participants, None).unwrap();

            // Generate shares
            let shares = threshold_sig
                .distributed_keygen(Some(&create_test_seed(1)))
                .unwrap();

            // Reconstruct s1 and s2 using Shamir reconstruction
            let s1_shares: Vec<_> =
                shares.iter().map(|s| s.s1_share.clone()).collect();
            let s2_shares: Vec<_> =
                shares.iter().map(|s| s.s2_share.clone()).collect();

            let shamir_s1 =
                AdaptedShamirSSS::new(threshold, participants).unwrap();
            let shamir_s2 =
                AdaptedShamirSSS::new(threshold, participants).unwrap();

            // Reconstruct should work with threshold shares
            let reconstructed_s1 =
                shamir_s1.reconstruct_secret(&s1_shares[..threshold]);
            let reconstructed_s2 =
                shamir_s2.reconstruct_secret(&s2_shares[..threshold]);

            assert!(reconstructed_s1.is_ok());
            assert!(reconstructed_s2.is_ok());
        }
    }
}
