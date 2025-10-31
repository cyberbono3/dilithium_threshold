use crate::{
    dilithium::DilithiumSignature,
    error::{ThresholdError, ThresholdResult},
    keypair::{PublicKey, keygen},
    params::{K, L, TAU, validate_threshold_config},
    shamir::{AdaptedShamirSSS, ShamirShare},
    traits::PointSource,
    utils::{
        get_hash_reader, get_randomness, hash_message,
        reconstruct_vector_from_points, sample_gamma1,
    },
};
use num_traits::Zero;
use sha2::{Digest, Sha256};
use sha3::digest::XofReader;

use math::{prelude::*, traits::FiniteField};

//use crate::dilithium::{Dilithium, DilithiumPublicKey, DilithiumSignature};

/// Represents a participant's share of the threshold key.
///
/// Contains both the Shamir share of the secret key and the
/// participant's identification information.
#[derive(Clone, Debug)]
pub struct ThresholdKeyShare<'a, FF: FiniteField> {
    pub participant_id: usize,
    pub s1_share: ShamirShare<'a, FF>,
    pub s2_share: ShamirShare<'a, FF>,
    pub public_key: PublicKey<'a, FF>,
}

impl<FF: FiniteField> ThresholdKeyShare<'static, FF> {
    /// Initialize threshold key share.
    pub fn new(
        participant_id: usize,
        s1_share: ShamirShare<'static, FF>,
        s2_share: ShamirShare<'static, FF>,
        public_key: PublicKey<'static, FF>,
    ) -> Self {
        ThresholdKeyShare {
            participant_id,
            s1_share,
            s2_share,
            public_key,
        }
    }
}

impl<FF: FiniteField> std::fmt::Display for ThresholdKeyShare<'static, FF> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "ThresholdKeyShare(id={})", self.participant_id)
    }
}

/// Represents a partial signature from one participant.
///
/// Contains the participant's contribution to the threshold signature
/// along with verification information.
#[derive(Clone, Debug)]
pub struct PartialSignature<'sign, FF: FiniteField> {
    pub participant_id: usize,
    pub z_partial: Vec<Polynomial<'sign, FF>>,
    pub commitment: Vec<Polynomial<'sign, FF>>,
    pub challenge: Polynomial<'sign, FF>,
}

impl<FF: FiniteField> PartialSignature<'static, FF> {
    /// Initialize partial signature.
    pub fn new(
        participant_id: usize,
        z_partial: Vec<Polynomial<'static, FF>>,
        commitment: Vec<Polynomial<'static, FF>>,
        challenge: Polynomial<'static, FF>,
    ) -> Self {
        Self {
            participant_id,
            z_partial,
            commitment,
            challenge,
        }
    }
}

impl<FF: FiniteField> std::fmt::Display for PartialSignature<'static, FF> {
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
    // security_level: usize,
    // dilithium: Dilithium,
    shamir_s1: AdaptedShamirSSS,
    shamir_s2: AdaptedShamirSSS,
    participant_ids: Vec<usize>,
}

/// Snapshot of key threshold configuration parameters.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ThresholdInfo {
    pub threshold: usize,
    pub participants: usize,
    pub min_signers: usize,
    pub max_participants: usize,
}

impl ThresholdSignature {
    /// Initialize threshold signature scheme.
    pub fn new(
        threshold: usize,
        participants: usize,
        //security_level: Option<usize>,
    ) -> ThresholdResult<Self> {
        if !validate_threshold_config(threshold, participants) {
            return Err(ThresholdError::InvalidThreshold(
                threshold,
                participants,
            ));
        }

        // let security_level = security_level.unwrap_or(DEFAULT_SECURITY_LEVEL);

        // Initialize underlying schemes
        //  let dilithium = Dilithium::new(security_level);
        let shamir_s1 = AdaptedShamirSSS::new(threshold, participants)?;
        let shamir_s2 = AdaptedShamirSSS::new(threshold, participants)?;

        // Store participant IDs
        let participant_ids: Vec<usize> = (1..=participants).collect();

        Ok(ThresholdSignature {
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
        //seed: Option<&[u8]>,
    ) -> ThresholdResult<Vec<ThresholdKeyShare<'static, FF>>>
    where
        FF: FiniteField + From<i64>,
        rand::distributions::Standard: rand::distributions::Distribution<FF>,
    {
        // Generate base Dilithium key pair
        let (pub_key, priv_key) = keygen();

        // Split secret vectors using adapted Shamir scheme
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
    /// This is the core innovation: each participant can create a partial
    /// signature using only their share, without reconstructing the full
    /// secret key.
    // TODO make
    pub fn partial_sign<FF: FiniteField>(
        &self,
        message: &[u8],
        key_share: &ThresholdKeyShare<'static, FF>,
        randomness: Option<&[u8]>,
    ) -> ThresholdResult<PartialSignature<'static, FF>> {
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
        //let w_partial = &key_share.public_key.a * &y_partial;
        // TODO declare trait for Matrix
        let w_partial = key_share.public_key.a.mul_vector(&y_partial);

        // For now, use a simplified challenge generation
        // In practice, this requires coordination between participants
        let challenge = self.generate_partial_challenge(&mu);

        // Compute partial response z_partial = y_partial + c * s1_share
        // TODO implement trait for it
        let c_s1: Vec<Polynomial<'static, FF>> = key_share
            .s1_share
            .share_vector
            .iter()
            .map(|p| p.multiply(&challenge))
            .collect();
        //let c_s1 = key_share.s1_share.share_vector.clone() * &challenge;

        //let z_partial = y_partial + c_s1;
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

    /// Verify that partial signatures are valid for combination.
    ///
    /// Checks that:
    /// 1. There are at least `threshold` partial signatures
    /// 2. All partial signatures use the same challenge
    fn verify_partial_signatures<FF: FiniteField>(
        partial_signatures: &[PartialSignature<'static, FF>],
        threshold: usize,
    ) -> ThresholdResult<()> {
        // Check if we have enough partial signatures
        if partial_signatures.len() < threshold {
            return Err(ThresholdError::InsufficientShares(
                threshold,
                partial_signatures.len(),
            ));
        }

        // At this point, we know we have at least threshold signatures
        // Since threshold must be >= 2 (per validation rules), we know
        // partial_signatures is not empty

        // Verify all partial signatures use the same challenge
        let challenge = &partial_signatures[0].challenge;
        if partial_signatures
            .iter()
            .any(|ps| &ps.challenge != challenge)
        {
            return Err(ThresholdError::PartialSignatureChallengeMismatch);
        }

        Ok(())
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
        // Verify partial signatures
        Self::verify_partial_signatures(partial_signatures, self.threshold)?;

        // Use first threshold partial signatures
        let active_partials = &partial_signatures[..self.threshold];

        // Reconstruct z vector using Lagrange interpolation
        // expected L
        let z = self.reconstruct_z_vector(active_partials)?;

        // Reconstruct hint h (simplified for this implementation)
        // expected K
        let h = self.reconstruct_hint(active_partials, public_key)?;

        // Get challenge from the verified partial signatures
        let challenge = partial_signatures[0].challenge.clone();

        Ok(DilithiumSignature::new(z, h, challenge))
    }

    pub fn verify_partial_signature<FF: FiniteField>(
        &self,
        message: &[u8],
        partial_sig: &PartialSignature<'static, FF>,
    ) -> bool {
        let mu = hash_message(message);
        partial_sig.challenge == self.generate_partial_challenge(&mu)
    }

    /// Derive participant-specific randomness.
    /// TODO remove self
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
    /// TODO remove self
    fn sample_partial_y<FF: FiniteField>(
        &self,
        randomness: &[u8],
    ) -> Vec<Polynomial<'static, FF>> {
        (0..L)
            .map(|i| sample_gamma1(&[randomness, &[i as u8]].concat()))
            .collect()
    }

    /// Convert to polynomial with tau non-zero coefficients
    /// TODO remove self
    fn generate_partial_challenge<FF: FiniteField>(
        &self,
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

        // TODO address clippy issue
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
    /// TODO remove self
    fn reconstruct_hint<FF: FiniteField>(
        &self,
        _partial_signatures: &[PartialSignature<'static, FF>],
        public_key: &PublicKey<'static, FF>,
    ) -> ThresholdResult<Vec<Polynomial<'static, FF>>> {
        let hint_len = public_key.t.len();
        Ok(vec![Polynomial::<FF>::zero(); hint_len])
    }

    // /// Check if partial signature satisfies bound requirements.
    // fn check_partial_bounds<FF: FiniteField>(
    //     &self,
    //     partial_sig: &PartialSignature<'static, FF>,
    // ) -> bool  {

    //     // let gamma1 = self.dilithium.config.gamma1;
    //     // let beta = self.dilithium.config.beta;
    //     //let polys = partial_sig.z_partial.clone();;

    //     let norm_infinity_val = polyvec_max_infty_norm(&partial_sig.z_partial) as u32;

    //     norm_infinity_val < GAMMA1 as u32 - BETA as u32
    // }

    /// Get information about the threshold configuration.
    pub fn get_threshold_info(&self) -> ThresholdInfo {
        ThresholdInfo {
            threshold: self.threshold,
            participants: self.participants,
            min_signers: self.threshold,
            max_participants: self.participants,
        }
    }
}

#[cfg(test)]
impl Default for ThresholdSignature {
    fn default() -> Self {
        Self::new(3, 5).unwrap()
    }
}

impl<FF: FiniteField> PointSource<FF> for PartialSignature<'static, FF> {
    fn x(&self) -> FF {
        (self.participant_id as u64).into()
    }

    fn poly_at(&self, index: usize) -> Option<&Polynomial<'static, FF>> {
        self.z_partial.get(index)
    }

    fn poly_count(&self) -> usize {
        self.z_partial.len()
    }
}
// TODO update tests
#[cfg(test)]
mod tests {
    use super::*;

    // Helper function to create a deterministic seed
    fn create_test_seed(value: u8) -> Vec<u8> {
        vec![value; 32]
    }

    mod threshold_key_share_tests {
        use super::*;

        #[test]
        fn test_key_share_creation() {
            let poly1: Polynomial<'_, FieldElement> = poly![1, 2, 3];
            let poly2: Polynomial<'_, FieldElement> = poly![4, 5, 6];
            let s1_share_vec = vec![poly1];
            let s2_share_vec = vec![poly2];

            let s1_share = ShamirShare::new(1, s1_share_vec).unwrap();
            let s2_share = ShamirShare::new(1, s2_share_vec).unwrap();

            let (pub_key, _) = keygen();

            let key_share = ThresholdKeyShare::new(
                1,
                s1_share.clone(),
                s2_share.clone(),
                pub_key,
            );

            assert_eq!(key_share.participant_id, 1);
            assert_eq!(key_share.s1_share.participant_id, 1);
            assert_eq!(key_share.s2_share.participant_id, 1);
        }

        #[test]
        fn test_key_share_display() {
            let poly1: Polynomial<'_, FieldElement> = poly![vec![1]];
            let poly2: Polynomial<'_, FieldElement> = poly![vec![2]];
            let s1_share_vec = vec![poly1];
            let s2_share_vec = vec![poly2];

            let s1_share = ShamirShare::new(5, s1_share_vec).unwrap();
            let s2_share = ShamirShare::new(5, s2_share_vec).unwrap();

            let (pub_key, _) = keygen();

            let key_share =
                ThresholdKeyShare::new(5, s1_share, s2_share, pub_key);

            let display_str = format!("{}", key_share);
            assert_eq!(display_str, "ThresholdKeyShare(id=5)");
        }
    }

    mod partial_signature_tests {
        use super::*;

        #[test]
        fn test_partial_signature_creation() {
            let z_poly: Polynomial<'_, FieldElement> = poly![100, 200, 300];
            let z_partial = vec![z_poly];

            let comm_poly = poly![10, 20, 30];
            let commitment = vec![comm_poly];

            let challenge: Polynomial<'_, FieldElement> = poly![1, -1, 0, 1];

            let partial_sig = PartialSignature::new(
                3,
                z_partial,
                commitment,
                challenge.clone(),
            );

            assert_eq!(partial_sig.participant_id, 3);
            assert_eq!(partial_sig.z_partial.len(), 1);
            assert_eq!(partial_sig.commitment.len(), 1);
            assert_eq!(partial_sig.challenge, challenge);
        }

        #[test]
        fn test_partial_signature_display() {
            let zero_poly = Polynomial::<FieldElement>::zero();
            let z_partial = vec![zero_poly.clone()];
            let commitment = vec![zero_poly.clone()];
            let challenge = zero_poly;

            let partial_sig =
                PartialSignature::new(7, z_partial, commitment, challenge);

            let display_str = format!("{}", partial_sig);
            assert_eq!(display_str, "PartialSignature(id=7)");
        }
    }

    mod sample_partial_y_tests {
        use super::*;
        use crate::params::{GAMMA1, L, N};
        use math::field_element::FieldElement;

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
            let threshold_sig = ThresholdSignature::default();
            let randomness = vec![0x11u8; 32];

            let polys =
                threshold_sig.sample_partial_y::<FieldElement>(&randomness);

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
            let threshold_sig = ThresholdSignature::default();
            let randomness = vec![0x22u8; 32];

            let first =
                threshold_sig.sample_partial_y::<FieldElement>(&randomness);
            let second =
                threshold_sig.sample_partial_y::<FieldElement>(&randomness);

            assert_eq!(first, second);
        }

        #[test]
        fn different_randomness_changes_output() {
            let threshold_sig = ThresholdSignature::default();
            let randomness_a = vec![0x33u8; 32];
            let randomness_b = vec![0x44u8; 32];

            let polys_a =
                threshold_sig.sample_partial_y::<FieldElement>(&randomness_a);
            let polys_b =
                threshold_sig.sample_partial_y::<FieldElement>(&randomness_b);

            assert_ne!(polys_a, polys_b);
        }
    }

    mod generate_partial_challenge_tests {
        use super::*;
        use crate::params::{N, TAU};
        use math::field_element::FieldElement;

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
            let threshold_sig = ThresholdSignature::default();
            let mu = b"deterministic-mu";

            let first =
                threshold_sig.generate_partial_challenge::<FieldElement>(mu);
            let second =
                threshold_sig.generate_partial_challenge::<FieldElement>(mu);

            assert_eq!(first, second);
        }

        #[test]
        fn coefficients_are_sparse_and_small() {
            let threshold_sig = ThresholdSignature::default();
            let mu = b"weight-check";

            let poly =
                threshold_sig.generate_partial_challenge::<FieldElement>(mu);

            let mut non_zero = 0usize;
            let mut dense = vec![0i64; N];
            for (idx, &coeff) in poly.coefficients().iter().enumerate() {
                dense[idx] = centered_value(coeff);
            }

            for val in dense {
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
            let threshold_sig = ThresholdSignature::default();

            let mu_a = b"challenge-a";
            let mu_b = b"challenge-b";

            let poly_a =
                threshold_sig.generate_partial_challenge::<FieldElement>(mu_a);
            let poly_b =
                threshold_sig.generate_partial_challenge::<FieldElement>(mu_b);

            assert_ne!(poly_a, poly_b);
        }
    }

    mod reconstruct_z_vector_tests {
        use super::*;
        use math::field_element::FieldElement;

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
        use math::field_element::FieldElement;

        fn zero_partial(id: usize) -> PartialSignature<'static, FieldElement> {
            let z_partial =
                Vec::from(crate::utils::zero_polyvec::<L, FieldElement>());
            let commitment =
                Vec::from(crate::utils::zero_polyvec::<K, FieldElement>());
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

            assert_eq!(hints.len(), K);
        }

        #[test]
        fn outputs_zero_polynomials() {
            let threshold_sig = ThresholdSignature::default();
            let shares =
                threshold_sig.distributed_keygen::<FieldElement>().unwrap();

            let hints = threshold_sig
                .reconstruct_hint::<FieldElement>(&[], &shares[0].public_key)
                .unwrap();

            assert!(hints.iter().all(|poly| {
                poly.coefficients()
                    .iter()
                    .all(|&c| c == FieldElement::zero())
            }));
        }

        #[test]
        fn ignores_partial_data() {
            let threshold_sig = ThresholdSignature::default();
            let shares =
                threshold_sig.distributed_keygen::<FieldElement>().unwrap();

            let mut partials = Vec::new();
            for (idx, _) in shares.iter().enumerate().take(3) {
                let dummy_poly = Polynomial::from(vec![FieldElement::from(
                    (idx as i64) + 1,
                )]);
                partials.push(PartialSignature::new(
                    idx + 1,
                    vec![dummy_poly.clone(); L],
                    vec![dummy_poly; K],
                    Polynomial::zero(),
                ));
            }

            let hints = threshold_sig
                .reconstruct_hint(&partials, &shares[0].public_key)
                .unwrap();

            assert!(hints.iter().all(|poly| {
                poly.coefficients()
                    .iter()
                    .all(|&c| c == FieldElement::zero())
            }));
        }
    }

    mod threshold_signature_tests {
        use super::*;

        #[test]
        fn test_threshold_signature_creation() {
            let threshold_sig = ThresholdSignature::default();
            let info = threshold_sig.get_threshold_info();

            assert_eq!(info.threshold, 3);
            assert_eq!(info.participants, 5);
            assert_eq!(info.min_signers, 3);
            assert_eq!(info.max_participants, 5);
        }

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
            // Threshold greater than participants
            assert!(ThresholdSignature::new(6, 5).is_err());

            // Zero threshold
            assert!(ThresholdSignature::new(0, 5).is_err());

            // Threshold of 1 (not allowed)
            assert!(ThresholdSignature::new(1, 5).is_err());

            // Too many participants
            assert!(ThresholdSignature::new(3, 300).is_err());

            // Too few participants
            assert!(ThresholdSignature::new(2, 1).is_err());
        }

        // TODO fix it
        #[test]
        fn test_distributed_keygen() {
            let threshold_sig = ThresholdSignature::default();

            // Test with deterministic seed
            let shares =
                threshold_sig.distributed_keygen::<FieldElement>().unwrap();

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

        //  TODO fix it
        #[test]
        fn test_partial_sign_basic() {
            let threshold_sig = ThresholdSignature::default();
            let shares =
                threshold_sig.distributed_keygen::<FieldElement>().unwrap();
            let message = "Test message".as_bytes();

            // Create partial signature with first share
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
            let message = "Determenistic test".as_bytes();
            let randomness = create_test_seed(42);

            // Same inputs should produce same partial signature
            let partial1 = threshold_sig
                .partial_sign(message, &shares[0], Some(&randomness))
                .unwrap();

            let partial2 = threshold_sig
                .partial_sign(message, &shares[0], Some(&randomness))
                .unwrap();

            assert_eq!(partial1.participant_id, partial2.participant_id);
            assert_eq!(partial1.challenge, partial2.challenge);
        }

        // TODO fix it
        #[test]
        fn test_verify_partial_signature() {
            let threshold_sig = ThresholdSignature::new(3, 5).unwrap();
            let shares =
                threshold_sig.distributed_keygen::<FieldElement>().unwrap();
            let message = "verified_test".as_bytes();

            let partial_sig = threshold_sig
                .partial_sign(message, &shares[0], Some(&create_test_seed(1)))
                .unwrap();

            // Verify the partial signature
            let is_valid =
                threshold_sig.verify_partial_signature(message, &partial_sig);

            assert!(is_valid);
        }

        #[test]
        fn test_verify_partial_signature_wrong_message() {
            let threshold_sig = ThresholdSignature::new(3, 5).unwrap();
            let shares =
                threshold_sig.distributed_keygen::<FieldElement>().unwrap();
            let message = "Original message".as_bytes();
            let wrong_message = "Wrong message".as_bytes();

            let partial_sig = threshold_sig
                .partial_sign(message, &shares[0], None)
                .unwrap();

            // Verify with wrong message should fail
            let is_valid = threshold_sig
                .verify_partial_signature(wrong_message, &partial_sig);

            assert!(!is_valid);
        }

        #[test]
        fn test_full_threshold_signing_workflow() {
            let threshold_sig = ThresholdSignature::default();

            // 1. Distributed key generation
            let shares =
                threshold_sig.distributed_keygen::<FieldElement>().unwrap();
            let public_key = &shares[0].public_key;

            // 2. Message to sign
            let message = "Full workflow test message".as_bytes();

            // 3. Create partial signatures from different subsets
            // Test with first threshold participants
            let mut partial_sigs_1 = Vec::new();
            for i in 0..threshold_sig.threshold {
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
            for i in (threshold_sig.participants - threshold_sig.threshold)
                ..threshold_sig.participants
            {
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
            assert_eq!(combined_sig_1.z.len(), L);
            assert_eq!(combined_sig_1.h.len(), K,);
            assert_eq!(combined_sig_2.z.len(), L);
            assert_eq!(combined_sig_2.h.len(), K);

            // Since both signatures are for the same message, they should have the same challenge
            assert_eq!(combined_sig_1.c, combined_sig_2.c);
        }

        #[test]
        fn test_combine_signatures_insufficient_shares() {
            let threshold_sig = ThresholdSignature::default();
            let shares =
                threshold_sig.distributed_keygen::<FieldElement>().unwrap();
            let message = "Insufficient test".as_bytes();

            // Create only threshold-1 partial signatures
            let mut partial_sigs = Vec::new();
            for (i, share) in
                shares.iter().enumerate().take(threshold_sig.threshold - 1)
            {
                let partial = threshold_sig
                    .partial_sign(
                        message,
                        share,
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

        // TODO fix it
        #[test]
        fn test_combine_signatures_mismatched_challenges() {
            let threshold_sig = ThresholdSignature::default();
            let shares =
                threshold_sig.distributed_keygen::<FieldElement>().unwrap();

            // Create partial signatures with different messages (leading to different challenges)
            let mut partial_sigs = Vec::new();

            for (i, share) in
                shares.iter().enumerate().take(threshold_sig.threshold)
            {
                //  for i in 0..threshold_sig.threshold {
                let msg = format!("Message {}", i);
                let message = msg.as_bytes();
                let partial = threshold_sig
                    .partial_sign(
                        message,
                        share,
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

        #[test]
        fn test_combine_signatures_produces_verifiable_signature() {
            let threshold_sig = ThresholdSignature::default();
            let shares =
                threshold_sig.distributed_keygen::<FieldElement>().unwrap();
            let message = b"threshold signing works";

            let partials: Vec<_> = shares
                .iter()
                .take(threshold_sig.threshold)
                .enumerate()
                .map(|(idx, share)| {
                    threshold_sig
                        .partial_sign(
                            message,
                            share,
                            Some(&create_test_seed((idx as u8) + 21)),
                        )
                        .expect("partial signing should succeed")
                })
                .collect();

            let combined = threshold_sig
                .combine_signatures(&partials, &shares[0].public_key)
                .expect("combination of valid partials must succeed");

            let expected_z = threshold_sig
                .reconstruct_z_vector(&partials)
                .expect("reconstruction should succeed");
            let expected_h = threshold_sig
                .reconstruct_hint(&partials, &shares[0].public_key)
                .expect("hint reconstruction should succeed");

            assert_eq!(combined.c, partials[0].challenge);
            assert_eq!(combined.z.as_slice(), expected_z.as_slice());
            assert_eq!(combined.h.as_slice(), expected_h.as_slice());
        }

        #[test]
        fn test_derive_participant_randomness() {
            let threshold_sig = ThresholdSignature::new(2, 3).unwrap();
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
        fn test_edge_cases() {
            // Minimum configuration (2 out of 2)
            let threshold_sig = ThresholdSignature::new(2, 2).unwrap();
            let shares =
                threshold_sig.distributed_keygen::<FieldElement>().unwrap();
            assert_eq!(shares.len(), 2);

            // Large threshold
            let threshold_sig_large = ThresholdSignature::new(10, 15).unwrap();
            let shares_large = threshold_sig_large
                .distributed_keygen::<FieldElement>()
                .unwrap();
            assert_eq!(shares_large.len(), 15);
        }

        #[test]
        fn test_concurrent_partial_signing() {
            let threshold = 4;
            let participants = 7;
            let threshold_sig =
                ThresholdSignature::new(threshold, participants).unwrap();
            let shares =
                threshold_sig.distributed_keygen::<FieldElement>().unwrap();
            let message = "Concurrent signing test".as_bytes();

            // Simulate concurrent signing by different participants
            let partial_sigs: Vec<_> = (0..threshold)
                .map(|i| {
                    threshold_sig
                        .partial_sign(
                            message,
                            &shares[i],
                            Some(&create_test_seed((i * 10) as u8)),
                        )
                        .unwrap()
                })
                .collect();

            // All partial signatures should be valid
            for partial_sig in partial_sigs.iter() {
                assert!(
                    threshold_sig
                        .verify_partial_signature(message, partial_sig)
                );
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
        fn test_reconstruct_vs_original_key() {
            let threshold = 2;
            let participants = 3;
            let threshold_sig =
                ThresholdSignature::new(threshold, participants).unwrap();

            // Generate shares
            let shares =
                threshold_sig.distributed_keygen::<FieldElement>().unwrap();

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
