use crate::dilithium::error::{ThresholdError, ThresholdResult};
use crate::dilithium::params::{L, N, TAU, validate_threshold_config};
use crate::dilithium::shamir::AdaptedShamirSSS;
use crate::dilithium::sign::DilithiumSignature;
use crate::dilithium::utils::{
    get_hash_reader, get_randomness, hash_message,
    reconstruct_vector_from_points, sample_gamma1,
};
use crate::signature::keypair::{PublicKey, keygen};
use math::{prelude::*, traits::FiniteField};
use num_traits::Zero;
use sha2::{Digest, Sha256};
use sha3::digest::XofReader;

use super::{key_share::ThresholdKeyShare, partial::PartialSignature};

/// Snapshot of key threshold configuration parameters.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ThresholdInfo {
    pub threshold: usize,
    pub participants: usize,
    pub min_signers: usize,
    pub max_participants: usize,
}

/// Main threshold signature scheme implementation.
///
/// Provides distributed key generation, partial signing, and signature
/// combination functionality while preventing secret leakage.
#[derive(Debug)]
pub struct ThresholdSignature {
    threshold: usize,
    participants: usize,
    shamir_s1: AdaptedShamirSSS,
    shamir_s2: AdaptedShamirSSS,
    participant_ids: Vec<usize>,
}

impl ThresholdSignature {
    /// Initialize threshold signature scheme.
    pub fn new(threshold: usize, participants: usize) -> ThresholdResult<Self> {
        if !validate_threshold_config(threshold, participants) {
            return Err(ThresholdError::InvalidThreshold(
                threshold,
                participants,
            ));
        }

        let shamir_s1 = AdaptedShamirSSS::new(threshold, participants)?;
        let shamir_s2 = AdaptedShamirSSS::new(threshold, participants)?;
        let participant_ids: Vec<usize> = (1..=participants).collect();

        Ok(Self {
            threshold,
            participants,
            shamir_s1,
            shamir_s2,
            participant_ids,
        })
    }

    /// Generate threshold keys using distributed key generation.
    ///
    /// This method generates a Dilithium key pair and then splits the
    /// private key into shares using the adapted Shamir scheme.
    pub fn distributed_keygen<FF>(
        &self,
    ) -> ThresholdResult<Vec<ThresholdKeyShare<'static, FF>>>
    where
        FF: FiniteField + From<i64>,
        rand::distributions::Standard: rand::distributions::Distribution<FF>,
    {
        let (pub_key, priv_key) = keygen::<FF>();
        let s1_shares = self.shamir_s1.split_secret(&priv_key.s1)?;
        let s2_shares = self.shamir_s2.split_secret(&priv_key.s2)?;

        let threshold_shares = (0..self.participants)
            .map(|i| {
                ThresholdKeyShare::new(
                    self.participant_ids[i],
                    s1_shares[i].clone(),
                    s2_shares[i].clone(),
                    pub_key.clone(),
                )
            })
            .collect();

        Ok(threshold_shares)
    }

    /// Create a partial signature using a key share.
    ///
    /// Each participant can create a partial signature using only their share,
    /// without reconstructing the full secret key.
    pub fn partial_sign<FF: FiniteField>(
        &self,
        message: &[u8],
        key_share: &ThresholdKeyShare<'static, FF>,
        randomness: Option<&[u8]>,
    ) -> ThresholdResult<PartialSignature<'static, FF>> {
        let randomness = get_randomness(randomness);
        let mu = hash_message(message);

        let participant_randomness = Self::derive_participant_randomness(
            &randomness,
            key_share.participant_id,
        );
        let y_partial = Self::sample_partial_y(&participant_randomness);

        let w_partial = key_share.public_key.a.mul_vector(&y_partial);
        let challenge = Self::generate_partial_challenge(&mu);

        let c_s1: Vec<Polynomial<'static, FF>> = key_share
            .s1_share
            .share_vector
            .iter()
            .map(|p| p.multiply(&challenge))
            .collect();

        let z_partial: Vec<Polynomial<'static, FF>> = y_partial
            .into_iter()
            .zip(c_s1)
            .map(|(p1, p2)| p1 + p2)
            .collect();

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
    pub fn combine_signatures<FF: FiniteField>(
        &self,
        partial_signatures: &[PartialSignature<'static, FF>],
        public_key: &PublicKey<'static, FF>,
    ) -> ThresholdResult<DilithiumSignature<'static, FF>> {
        Self::verify_partial_signatures(partial_signatures, self.threshold)?;

        let active_partials = &partial_signatures[..self.threshold];
        let z = self.reconstruct_z_vector(active_partials)?;
        let h = self.reconstruct_hint(active_partials, public_key)?;
        let challenge = partial_signatures[0].challenge.clone();

        Ok(DilithiumSignature::new(z, h, challenge))
    }

    pub fn verify_partial_signature<FF: FiniteField>(
        &self,
        message: &[u8],
        partial_sig: &PartialSignature<'static, FF>,
    ) -> bool {
        let mu = hash_message(message);
        partial_sig.challenge == Self::generate_partial_challenge(&mu)
    }

    /// Get information about the threshold configuration.
    pub fn get_threshold_info(&self) -> ThresholdInfo {
        ThresholdInfo {
            threshold: self.threshold,
            participants: self.participants,
            min_signers: self.threshold,
            max_participants: self.participants,
        }
    }

    /// Verify that partial signatures are valid for combination.
    fn verify_partial_signatures<FF: FiniteField>(
        partial_signatures: &[PartialSignature<'static, FF>],
        threshold: usize,
    ) -> ThresholdResult<()> {
        if partial_signatures.len() < threshold {
            return Err(ThresholdError::InsufficientShares(
                threshold,
                partial_signatures.len(),
            ));
        }

        let challenge = &partial_signatures[0].challenge;
        if partial_signatures
            .iter()
            .any(|ps| &ps.challenge != challenge)
        {
            return Err(ThresholdError::PartialSignatureChallengeMismatch);
        }

        Ok(())
    }

    /// Derive participant-specific randomness.
    fn derive_participant_randomness(
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
    fn sample_partial_y<FF: FiniteField>(
        randomness: &[u8],
    ) -> Vec<Polynomial<'static, FF>> {
        (0..L)
            .map(|i| sample_gamma1(&[randomness, &[i as u8]].concat()))
            .collect()
    }

    /// Convert to polynomial with tau non-zero coefficients.
    fn generate_partial_challenge<FF: FiniteField>(
        mu: &[u8],
    ) -> Polynomial<'static, FF> {
        let mut seed = Vec::with_capacity(mu.len() + b"challenge".len());
        seed.extend_from_slice(mu);
        seed.extend_from_slice(b"challenge");

        let mut reader = get_hash_reader(&seed);
        let mut buffer = vec![0u8; TAU * 2];
        reader.read(&mut buffer);

        let mut coeffs = vec![0i32; N];
        for chunk in buffer.chunks_exact(2) {
            let position = (chunk[0] as usize) % N;
            let sign = if chunk[1] & 1 == 0 { 1 } else { -1 };
            coeffs[position] = sign;
        }

        poly![coeffs]
    }

    fn reconstruct_z_vector<FF: FiniteField>(
        &self,
        partial_signatures: &[PartialSignature<'static, FF>],
    ) -> ThresholdResult<Vec<Polynomial<'static, FF>>> {
        if partial_signatures.is_empty() {
            return Err(ThresholdError::InsufficientShares(1, 0));
        }

        let active =
            partial_signatures.get(..self.threshold).ok_or_else(|| {
                ThresholdError::InvalidThreshold(
                    self.threshold,
                    partial_signatures.len(),
                )
            })?;

        let vector_length = active
            .first()
            .map(|ps| ps.z_partial.len())
            .unwrap_or_default();

        let indices: Vec<usize> = (0..vector_length).collect();
        reconstruct_vector_from_points::<FF, _>(active, &indices)
    }

    /// Reconstruct hint vector (simplified implementation).
    fn reconstruct_hint<FF: FiniteField>(
        &self,
        _partial_signatures: &[PartialSignature<'static, FF>],
        public_key: &PublicKey<'static, FF>,
    ) -> ThresholdResult<Vec<Polynomial<'static, FF>>> {
        let hint_len = public_key.t.len();
        Ok(vec![Polynomial::<FF>::zero(); hint_len])
    }
}

#[cfg(test)]
impl Default for ThresholdSignature {
    fn default() -> Self {
        Self::new(3, 5).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dilithium::params::{GAMMA1, K, L, N};
    use crate::dilithium::utils::zero_polyvec;
    use math::field_element::FieldElement;

    fn create_test_seed(value: u8) -> Vec<u8> {
        vec![value; 32]
    }

    mod sample_partial_y_tests {
        use super::*;

        fn centered_value(fe: FieldElement) -> i64 {
            let mut v = i64::from(fe);
            let p = FieldElement::P as i64;
            if v > p / 2 {
                v -= p;
            }
            v
        }

        #[test]
        fn produces_l_polynomials_within_bounds() {
            let randomness = vec![0x11u8; 32];

            let polys = ThresholdSignature::sample_partial_y::<FieldElement>(
                &randomness,
            );

            assert_eq!(polys.len(), L);
            for poly in &polys {
                assert_eq!(poly.coefficients().len(), N);
                for coeff in poly.coefficients() {
                    let abs = centered_value(*coeff).abs();
                    assert!(
                        abs <= GAMMA1 as i64,
                        "coefficient {abs} exceeds GAMMA1"
                    );
                }
            }
        }

        #[test]
        fn deterministic_for_same_randomness() {
            let randomness = vec![0x22u8; 32];

            let first = ThresholdSignature::sample_partial_y::<FieldElement>(
                &randomness,
            );
            let second = ThresholdSignature::sample_partial_y::<FieldElement>(
                &randomness,
            );

            assert_eq!(first, second);
        }

        #[test]
        fn different_randomness_changes_output() {
            let randomness_a = vec![0x33u8; 32];
            let randomness_b = vec![0x44u8; 32];

            let polys_a = ThresholdSignature::sample_partial_y::<FieldElement>(
                &randomness_a,
            );
            let polys_b = ThresholdSignature::sample_partial_y::<FieldElement>(
                &randomness_b,
            );

            assert_ne!(polys_a, polys_b);
        }
    }

    mod generate_partial_challenge_tests {
        use super::*;

        fn centered_value(fe: FieldElement) -> i64 {
            let mut v = i64::from(fe);
            let p = FieldElement::P as i64;
            if v > p / 2 {
                v -= p;
            }
            v
        }

        #[test]
        fn deterministic_for_identical_inputs() {
            let mu = b"deterministic-mu";

            let first = ThresholdSignature::generate_partial_challenge::<
                FieldElement,
            >(mu);
            let second = ThresholdSignature::generate_partial_challenge::<
                FieldElement,
            >(mu);

            assert_eq!(first, second);
        }

        #[test]
        fn coefficients_are_sparse_and_small() {
            let mu = b"weight-check";

            let poly = ThresholdSignature::generate_partial_challenge::<
                FieldElement,
            >(mu);

            let mut non_zero = 0usize;
            for &coeff in poly.coefficients() {
                let val = centered_value(coeff);
                assert!(
                    matches!(val, -1 | 0 | 1),
                    "unexpected coefficient value: {val}"
                );
                if val != 0 {
                    non_zero += 1;
                }
            }

            assert!(
                non_zero <= TAU,
                "expected at most {TAU} non-zero coefficients, found {non_zero}"
            );
        }

        #[test]
        fn different_inputs_vary_output() {
            let mu_a = b"challenge-a";
            let mu_b = b"challenge-b";

            let poly_a = ThresholdSignature::generate_partial_challenge::<
                FieldElement,
            >(mu_a);
            let poly_b = ThresholdSignature::generate_partial_challenge::<
                FieldElement,
            >(mu_b);

            assert_ne!(poly_a, poly_b);
        }
    }

    mod reconstruct_z_vector_tests {
        use super::*;

        fn make_partial(
            participant_id: usize,
            poly: &Polynomial<'static, FieldElement>,
        ) -> PartialSignature<'static, FieldElement> {
            PartialSignature::new(
                participant_id,
                vec![poly.clone()],
                vec![Polynomial::<FieldElement>::zero()],
                Polynomial::<FieldElement>::zero(),
            )
        }

        #[test]
        fn reconstructs_constant_polynomials() {
            let threshold_sig = ThresholdSignature::new(2, 3).unwrap();
            let shared_poly: Polynomial<'static, FieldElement> =
                poly![5, -3, 8];

            let partials = vec![
                make_partial(1, &shared_poly),
                make_partial(2, &shared_poly),
            ];

            let reconstructed =
                threshold_sig.reconstruct_z_vector(&partials).unwrap();

            assert_eq!(reconstructed, vec![shared_poly]);
        }

        #[test]
        fn errors_when_not_enough_partials() {
            let threshold_sig = ThresholdSignature::default();
            let poly: Polynomial<'static, FieldElement> = poly![1, 2, 3];
            let partials = vec![make_partial(1, &poly), make_partial(2, &poly)];

            let err =
                threshold_sig.reconstruct_z_vector(&partials).unwrap_err();

            assert!(matches!(
                err,
                ThresholdError::InvalidThreshold(expected, provided)
                if expected == threshold_sig.threshold && provided == partials.len()
            ));
        }

        #[test]
        fn errors_on_empty_input() {
            let threshold_sig = ThresholdSignature::default();
            let err = threshold_sig
                .reconstruct_z_vector::<FieldElement>(&[])
                .unwrap_err();

            assert!(matches!(
                err,
                ThresholdError::InsufficientShares(expected, actual)
                if expected == 1 && actual == 0
            ));
        }
    }

    mod reconstruct_hint_tests {
        use super::*;

        fn zero_partial(id: usize) -> PartialSignature<'static, FieldElement> {
            let z_partial = Vec::from(zero_polyvec::<L, FieldElement>());
            let commitment = Vec::from(zero_polyvec::<K, FieldElement>());
            PartialSignature::new(id, z_partial, commitment, Polynomial::zero())
        }

        #[test]
        fn matches_public_key_dimension() {
            let threshold_sig = ThresholdSignature::default();
            let shares =
                threshold_sig.distributed_keygen::<FieldElement>().unwrap();
            let partials = (0..threshold_sig.threshold)
                .map(|i| zero_partial(i + 1))
                .collect::<Vec<_>>();

            let hints = threshold_sig
                .reconstruct_hint(&partials, &shares[0].public_key)
                .unwrap();

            assert_eq!(hints.len(), shares[0].public_key.t.len());
        }

        #[test]
        fn returns_zero_polynomials() {
            let threshold_sig = ThresholdSignature::default();
            let shares =
                threshold_sig.distributed_keygen::<FieldElement>().unwrap();
            let partials = (0..threshold_sig.threshold)
                .map(|i| zero_partial(i + 1))
                .collect::<Vec<_>>();

            let hints = threshold_sig
                .reconstruct_hint(&partials, &shares[0].public_key)
                .unwrap();

            for poly in hints {
                assert!(poly.coefficients().iter().all(|c| c.is_zero()));
            }
        }
    }

    mod verify_partial_signatures_tests {
        use super::*;

        fn make_partial(id: usize) -> PartialSignature<'static, FieldElement> {
            PartialSignature::new(
                id,
                vec![Polynomial::zero()],
                vec![Polynomial::zero()],
                Polynomial::zero(),
            )
        }

        #[test]
        fn requires_minimum_partials() {
            let err = ThresholdSignature::verify_partial_signatures::<
                FieldElement,
            >(&[], 2)
            .unwrap_err();

            assert!(matches!(
                err,
                ThresholdError::InsufficientShares(expected, actual)
                if expected == 2 && actual == 0
            ));
        }

        #[test]
        fn detects_challenge_mismatch() {
            let a = make_partial(1);
            let mut b = make_partial(2);
            b.challenge = Polynomial::from(vec![fe!(1)]);

            let err = ThresholdSignature::verify_partial_signatures::<
                FieldElement,
            >(&[a, b], 2)
            .unwrap_err();

            assert!(matches!(
                err,
                ThresholdError::PartialSignatureChallengeMismatch
            ));
        }

        #[test]
        fn accepts_matching_partials() {
            let partials = vec![make_partial(1), make_partial(2)];
            assert!(
                ThresholdSignature::verify_partial_signatures::<FieldElement>(
                    &partials, 2,
                )
                .is_ok()
            );
        }
    }

    mod threshold_signature_tests {
        use super::*;

        #[test]
        fn test_threshold_info_matches_configuration() {
            let threshold_sig = ThresholdSignature::new(4, 7).unwrap();
            let info = threshold_sig.get_threshold_info();

            assert_eq!(info.threshold, 4);
            assert_eq!(info.participants, 7);
            assert_eq!(info.min_signers, 4);
            assert_eq!(info.max_participants, 7);
        }

        #[test]
        fn test_invalid_threshold_configurations() {
            assert!(ThresholdSignature::new(6, 5).is_err());
            assert!(ThresholdSignature::new(0, 5).is_err());
            assert!(ThresholdSignature::new(1, 5).is_err());
            assert!(ThresholdSignature::new(3, 300).is_err());
            assert!(ThresholdSignature::new(2, 1).is_err());
        }

        #[test]
        fn test_distributed_keygen() {
            let threshold_sig = ThresholdSignature::default();

            let shares =
                threshold_sig.distributed_keygen::<FieldElement>().unwrap();

            assert_eq!(shares.len(), 5);

            for (i, share) in shares.iter().enumerate() {
                assert_eq!(share.participant_id, i + 1);
                assert_eq!(share.s1_share.participant_id, i + 1);
                assert_eq!(share.s2_share.participant_id, i + 1);
            }

            let public_key = &shares[0].public_key;
            for share in &shares[1..] {
                assert_eq!(&share.public_key, public_key);
            }
        }

        #[test]
        fn test_partial_sign_basic() {
            let threshold_sig = ThresholdSignature::default();
            let shares =
                threshold_sig.distributed_keygen::<FieldElement>().unwrap();
            let message = b"Test message";

            let partial_sig = threshold_sig
                .partial_sign(message, &shares[0], Some(&create_test_seed(2)))
                .unwrap();

            assert_eq!(partial_sig.participant_id, 1);
            assert!(!partial_sig.z_partial.is_empty());
            assert!(!partial_sig.commitment.is_empty());
        }

        #[test]
        fn test_partial_sign_deterministic() {
            let threshold_sig = ThresholdSignature::new(2, 3).unwrap();
            let shares =
                threshold_sig.distributed_keygen::<FieldElement>().unwrap();
            let message = b"Determenistic test";
            let randomness = create_test_seed(42);

            let partial1 = threshold_sig
                .partial_sign(message, &shares[0], Some(&randomness))
                .unwrap();

            let partial2 = threshold_sig
                .partial_sign(message, &shares[0], Some(&randomness))
                .unwrap();

            assert_eq!(partial1.participant_id, partial2.participant_id);
            assert_eq!(partial1.challenge, partial2.challenge);
        }

        #[test]
        fn test_verify_partial_signature() {
            let threshold_sig = ThresholdSignature::new(3, 5).unwrap();
            let shares =
                threshold_sig.distributed_keygen::<FieldElement>().unwrap();
            let message = b"verified_test";

            let partial_sig = threshold_sig
                .partial_sign(message, &shares[0], Some(&create_test_seed(1)))
                .unwrap();

            assert!(
                threshold_sig.verify_partial_signature(message, &partial_sig)
            );
        }

        #[test]
        fn test_verify_partial_signature_wrong_message() {
            let threshold_sig = ThresholdSignature::new(3, 5).unwrap();
            let shares =
                threshold_sig.distributed_keygen::<FieldElement>().unwrap();
            let message = b"Original message";
            let wrong_message = b"Wrong message";

            let partial_sig = threshold_sig
                .partial_sign(message, &shares[0], None)
                .unwrap();

            let is_valid = threshold_sig
                .verify_partial_signature(wrong_message, &partial_sig);

            assert!(!is_valid);
        }
    }
}
