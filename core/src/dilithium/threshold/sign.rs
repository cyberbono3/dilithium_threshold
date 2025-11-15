use crate::basic::keypair::{KeyPair, PublicKey, keygen};
use crate::basic::utils::{make_hints, poly_high, polyvec_sub_scaled};
use crate::dilithium::error::{DilithiumError, DilithiumResult};
use crate::dilithium::params::{K, L, validate_threshold_config};
use crate::dilithium::shamir::AdaptedShamirSSS;
use crate::dilithium::sign::DilithiumSignature;
use crate::dilithium::utils::{
    derive_challenge_polynomial, get_randomness, hash_message,
    reconstruct_vector_from_points, sample_gamma1_vector,
};
use crate::matrix::MatrixMulExt;
use crate::traits::PolyVectorSource;
use math::{prelude::*, traits::FiniteField};
use sha2::{Digest, Sha256};

use super::error::ThresholdError;
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
    /// Build a threshold signing engine after validating the participant configuration.
    pub fn new(threshold: usize, participants: usize) -> DilithiumResult<Self> {
        if !validate_threshold_config(threshold, participants) {
            return Err(DilithiumError::InvalidThreshold(
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

    /// Run distributed key generation, returning per-participant Shamir shares.
    pub fn distributed_keygen<FF>(
        &self,
    ) -> DilithiumResult<Vec<ThresholdKeyShare<'static, FF>>>
    where
        FF: FiniteField + From<i64>,
        rand::distributions::Standard: rand::distributions::Distribution<FF>,
    {
        let KeyPair {
            public: pub_key,
            private: priv_key,
        } = keygen::<FF>()?;
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

    /// Produce a partial signature for `message` using a single participant's share.
    pub fn partial_sign<FF: FiniteField + From<i64>>(
        &self,
        message: &[u8],
        key_share: &ThresholdKeyShare<'static, FF>,
        randomness: Option<&[u8]>,
    ) -> DilithiumResult<PartialSignature<'static, FF>> {
        let randomness = get_randomness(randomness);
        let mu = hash_message(message);

        let participant_randomness = Self::derive_participant_randomness(
            &randomness,
            key_share.participant_id,
        );
        let y_partial = Self::sample_partial_y(&participant_randomness);

        let w_partial = key_share
            .public_key
            .a
            .mul_polynomials(y_partial.as_slice())?;
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

    /// Combine `threshold` partial signatures into a full Dilithium signature.
    pub fn combine_signatures<FF>(
        &self,
        partial_signatures: &[PartialSignature<'static, FF>],
        public_key: &PublicKey<'static, FF>,
    ) -> DilithiumResult<DilithiumSignature<'static, FF>>
    where
        FF: FiniteField + From<i64>,
        i64: From<FF>,
    {
        Self::verify_partial_signatures(partial_signatures, self.threshold)?;

        let active_partials = &partial_signatures[..self.threshold];
        let z = self.reconstruct_z_vector(active_partials)?;
        let az: [Polynomial<'static, FF>; K] =
            public_key.a.mul_polynomials(&z)?;
        let challenge = partial_signatures[0].challenge.clone();
        let az_minus_ct =
            polyvec_sub_scaled::<FF, K>(&az, &challenge, &public_key.t);
        let h = self.reconstruct_hint(active_partials, &az_minus_ct)?;

        Ok(DilithiumSignature::new(z, h, challenge))
    }

    /// Quickly check whether an individual partial signature matches the message.
    pub fn verify_partial_signature<FF: FiniteField + From<i64>>(
        &self,
        message: &[u8],
        partial_sig: &PartialSignature<'static, FF>,
    ) -> bool {
        let mu = hash_message(message);
        partial_sig.challenge == Self::generate_partial_challenge(&mu)
    }

    /// Expose the configuration parameters used to instantiate this instance.
    pub fn get_threshold_info(&self) -> ThresholdInfo {
        ThresholdInfo {
            threshold: self.threshold,
            participants: self.participants,
            min_signers: self.threshold,
            max_participants: self.participants,
        }
    }

    /// Ensure a collection of partial signatures is non-empty and shares the same challenge.
    fn verify_partial_signatures<FF: FiniteField>(
        partial_signatures: &[PartialSignature<'static, FF>],
        threshold: usize,
    ) -> DilithiumResult<()> {
        if partial_signatures.len() < threshold {
            return Err(DilithiumError::InsufficientShares(
                threshold,
                partial_signatures.len(),
            ));
        }

        let challenge = &partial_signatures[0].challenge;
        if partial_signatures
            .iter()
            .any(|ps| &ps.challenge != challenge)
        {
            return Err(
                ThresholdError::PartialSignatureChallengeMismatch.into()
            );
        }

        Ok(())
    }

    /// Expand the global randomness into a participant-specific seed via SHA-256.
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

    /// Draw the `y` mask polynomials used for partial signatures from the shared randomness.
    fn sample_partial_y<FF: FiniteField>(
        randomness: &[u8],
    ) -> Vec<Polynomial<'static, FF>> {
        sample_gamma1_vector(randomness, L)
    }

    /// Derive the sparse challenge polynomial expected by partial proofs.
    fn generate_partial_challenge<FF: FiniteField + From<i64>>(
        mu: &[u8],
    ) -> Polynomial<'static, FF> {
        let mut seed = Vec::with_capacity(mu.len() + b"challenge".len());
        seed.extend_from_slice(mu);
        seed.extend_from_slice(b"challenge");

        derive_challenge_polynomial::<FF>(&seed)
    }

    /// Rebuild the aggregated `z` response vector from the provided partials.
    fn reconstruct_z_vector<FF: FiniteField>(
        &self,
        partial_signatures: &[PartialSignature<'static, FF>],
    ) -> DilithiumResult<[Polynomial<'static, FF>; L]> {
        if partial_signatures.is_empty() {
            return Err(DilithiumError::InsufficientShares(1, 0));
        }

        let active = partial_signatures.get(..self.threshold).ok_or(
            DilithiumError::InvalidThreshold(
                self.threshold,
                partial_signatures.len(),
            ),
        )?;

        let vector_length = active
            .first()
            .map(|ps| ps.z_partial.len())
            .unwrap_or_default();

        if vector_length != L {
            return Err(DilithiumError::InvalidThreshold(L, vector_length));
        }

        let indices: Vec<usize> = (0..L).collect();
        let reconstructed =
            reconstruct_vector_from_points::<FF, _>(active, &indices)?;

        reconstructed.try_into().map_err(|vec: Vec<_>| {
            DilithiumError::InvalidThreshold(L, vec.len())
        })
    }

    fn reconstruct_commitment_vector<FF: FiniteField>(
        &self,
        partial_signatures: &[PartialSignature<'static, FF>],
    ) -> DilithiumResult<[Polynomial<'static, FF>; K]> {
        if partial_signatures.is_empty() {
            return Err(DilithiumError::InsufficientShares(1, 0));
        }

        let vector_length = partial_signatures
            .first()
            .map(|ps| ps.commitment.len())
            .unwrap_or_default();

        if vector_length != K {
            return Err(DilithiumError::InvalidThreshold(K, vector_length));
        }

        let sources: Vec<_> = partial_signatures
            .iter()
            .map(|ps| CommitmentSource {
                participant_id: ps.participant_id,
                commitment: ps.commitment.as_slice(),
            })
            .collect();

        let indices: Vec<usize> = (0..K).collect();
        let reconstructed =
            reconstruct_vector_from_points::<FF, _>(&sources, &indices)?;

        reconstructed.try_into().map_err(|vec: Vec<_>| {
            DilithiumError::InvalidThreshold(K, vec.len())
        })
    }

    fn reconstruct_hint<FF>(
        &self,
        partial_signatures: &[PartialSignature<'static, FF>],
        az_minus_ct: &[Polynomial<'static, FF>; K],
    ) -> DilithiumResult<[Polynomial<'static, FF>; K]>
    where
        FF: FiniteField + From<i64>,
        i64: From<FF>,
    {
        let commitments =
            self.reconstruct_commitment_vector(partial_signatures)?;
        let w1 = std::array::from_fn(|idx| poly_high(&commitments[idx]));
        Ok(make_hints(&w1, az_minus_ct))
    }
}

impl Default for ThresholdSignature {
    fn default() -> Self {
        Self::new(3, 5).expect("default threshold configuration is valid")
    }
}

struct CommitmentSource<'a, FF: FiniteField + 'static> {
    participant_id: usize,
    commitment: &'a [Polynomial<'static, FF>],
}

impl<'a, FF: FiniteField + 'static> PolyVectorSource<FF>
    for CommitmentSource<'a, FF>
{
    fn participant_id(&self) -> usize {
        self.participant_id
    }

    fn polynomials(&self) -> &[Polynomial<'static, FF>] {
        self.commitment
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dilithium::params::{K, L};
    use crate::dilithium::utils::zero_polyvec;
    use math::field_element::FieldElement;
    use num_traits::Zero;

    /// Build a deterministic 32-byte seed filled with the provided value.
    fn create_test_seed(value: u8) -> Vec<u8> {
        vec![value; 32]
    }

    mod reconstruct_z_vector_tests {
        use super::*;

        /// Helper to create a partial signature containing a shared polynomial vector.
        fn make_partial(
            participant_id: usize,
            poly: &Polynomial<'static, FieldElement>,
        ) -> PartialSignature<'static, FieldElement> {
            let mut z_partial = Vec::from(zero_polyvec::<L, FieldElement>());
            if let Some(slot) = z_partial.first_mut() {
                *slot = poly.clone();
            }
            PartialSignature::new(
                participant_id,
                z_partial,
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

            assert_eq!(reconstructed[0], shared_poly);
            for poly in reconstructed.iter().skip(1) {
                assert!(poly.coefficients().iter().all(|c| c.is_zero()));
            }
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
                DilithiumError::InvalidThreshold(expected, provided)
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
                DilithiumError::InsufficientShares(expected, actual)
                if expected == 1 && actual == 0
            ));
        }
    }

    mod reconstruct_hint_tests {
        use super::*;

        /// Build a dummy partial signature with zero-valued polynomials.
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
            let az_minus_ct = zero_polyvec::<K, FieldElement>();

            let hints = threshold_sig
                .reconstruct_hint(&partials, &az_minus_ct)
                .unwrap();

            assert_eq!(hints.len(), shares[0].public_key.t.len());
        }

        #[test]
        fn returns_zero_polynomials() {
            let threshold_sig = ThresholdSignature::default();
            let partials = (0..threshold_sig.threshold)
                .map(|i| zero_partial(i + 1))
                .collect::<Vec<_>>();
            let az_minus_ct = zero_polyvec::<K, FieldElement>();

            let hints = threshold_sig
                .reconstruct_hint(&partials, &az_minus_ct)
                .unwrap();

            for poly in hints.iter() {
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
                DilithiumError::InsufficientShares(expected, actual)
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
                DilithiumError::Threshold(
                    ThresholdError::PartialSignatureChallengeMismatch
                )
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
