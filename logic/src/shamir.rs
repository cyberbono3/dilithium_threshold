use math::{prelude::*, traits::FiniteField};

use crate::{
    error::{ThresholdError, ThresholdResult},
    params::validate_threshold_config,
    traits::PointSource,
    utils::reconstruct_vector_from_points,
};

use rand::Rng;

/// Adapted Shamir's Secret Sharing
#[derive(Clone, Debug, PartialEq)]
pub struct ShamirShare<'a, FF: FiniteField> {
    pub participant_id: usize,
    // pub share_vector: PolynomialVector<'a, FF>,
    pub share_vector: Vec<Polynomial<'a, FF>>,
}

impl<FF: FiniteField> ShamirShare<'static, FF> {
    pub fn new(
        participant_id: usize,
        share_vector: Vec<Polynomial<'static, FF>>,
    ) -> ThresholdResult<Self> {
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
struct Share<FF: FiniteField> {
    poly_idx: usize,
    coeff_idx: usize,
    value: FF,
}

impl<FF: FiniteField> Share<FF> {
    const fn new(poly_idx: usize, coeff_idx: usize, value: FF) -> Self {
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
    pub fn new(
        threshold: usize,
        participant_number: usize,
    ) -> ThresholdResult<Self> {
        if !validate_threshold_config(threshold, participant_number) {
            return Err(ThresholdError::InvalidThreshold(
                threshold,
                participant_number,
            ));
        }

        Ok(AdaptedShamirSSS {
            threshold,
            participant_number,
        })
    }

    /// Split a polynomial vector secret into shares.
    pub fn split_secret<FF>(
        &self,
        secret_vector: &[Polynomial<'static, FF>],
    ) -> ThresholdResult<Vec<ShamirShare<'static, FF>>>
    where
        FF: FiniteField,
        rand::distributions::Standard: rand::distributions::Distribution<FF>,
    {
        let vector_length = secret_vector.len();
        let mut buckets = self.init_share_buckets::<FF>(vector_length);

        for (poly_idx, polynomial) in secret_vector.iter().enumerate() {
            self.distribute_polynomial(poly_idx, polynomial, &mut buckets);
        }

        self.organize_shares(buckets, vector_length)
    }

    fn init_share_buckets<FF>(
        &self,
        vector_length: usize,
    ) -> Vec<Vec<Share<FF>>>
    where
        FF: FiniteField,
    {
        (0..self.participant_number)
            .map(|_| Vec::with_capacity(vector_length * N))
            .collect()
    }

    fn distribute_polynomial<FF>(
        &self,
        poly_idx: usize,
        polynomial: &Polynomial<'static, FF>,
        buckets: &mut [Vec<Share<FF>>],
    ) where
        FF: FiniteField,
        rand::distributions::Standard: rand::distributions::Distribution<FF>,
    {
        for (coeff_idx, secret_coeff) in
            polynomial.coefficients().iter().copied().enumerate()
        {
            let shamir_poly = self.create_shamir_polynomial(&secret_coeff);
            for (pid, bucket) in buckets.iter_mut().enumerate() {
                let point = FF::from((pid + 1) as u64);
                let value = shamir_poly.evaluate_in_same_field(point);
                bucket.push(Share::new(poly_idx, coeff_idx, value));
            }
        }
    }

    /// Reconstruct the secret from a sufficient number of shares.
    // TODO fix code duplciation with partial_reconstruct
    pub fn reconstruct_secret<FF: FiniteField>(
        &self,
        shares: &[ShamirShare<'static, FF>],
    ) -> ThresholdResult<Vec<Polynomial<'static, FF>>> {
        let (active_shares, vector_length) =
            self.select_active_shares(shares)?;
        let poly_indices: Vec<usize> = (0..vector_length).collect();
        self.reconstruct_poly_vector(active_shares, &poly_indices)
    }

    #[cfg(test)]
    /// Partially reconstruct only specified polynomials from the vector.
    pub fn partial_reconstruct<FF: FiniteField>(
        &self,
        shares: &[ShamirShare<'static, FF>],
        poly_indices: &[usize],
    ) -> ThresholdResult<Vec<Polynomial<'static, FF>>> {
        let (active_shares, vector_length) =
            self.select_active_shares(shares)?;
        self.ensure_poly_indices_within_bounds(poly_indices, vector_length)?;
        self.reconstruct_poly_vector(active_shares, poly_indices)
    }

    #[cfg(test)]
    fn ensure_poly_indices_within_bounds(
        &self,
        poly_indices: &[usize],
        vector_length: usize,
    ) -> ThresholdResult<()> {
        if let Some(&invalid_idx) =
            poly_indices.iter().find(|&&idx| idx >= vector_length)
        {
            return Err(ThresholdError::InvalidIndex(
                invalid_idx,
                vector_length,
            ));
        }
        Ok(())
    }

    fn reconstruct_poly_vector<FF: FiniteField>(
        &self,
        shares: &[ShamirShare<'static, FF>],
        poly_indices: &[usize],
    ) -> ThresholdResult<Vec<Polynomial<'static, FF>>> {
        reconstruct_vector_from_points::<FF, _>(shares, poly_indices)
    }

    fn select_active_shares<'a, FF>(
        &self,
        shares: &'a [ShamirShare<'static, FF>],
    ) -> ThresholdResult<(&'a [ShamirShare<'static, FF>], usize)>
    where
        FF: FiniteField,
    {
        if shares.len() < self.threshold {
            return Err(ThresholdError::InsufficientShares(
                self.threshold,
                shares.len(),
            ));
        }

        let active = &shares[..self.threshold];
        let vector_length = self.ensure_consistent_lengths(active)?;
        Ok((active, vector_length))
    }

    fn ensure_consistent_lengths<FF>(
        &self,
        shares: &[ShamirShare<'static, FF>],
    ) -> ThresholdResult<usize>
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
            return Err(ThresholdError::InconsistentShareLengths);
        }

        Ok(expected_length)
    }


    /// Organize participant shares into ShamirShare objects
    fn organize_shares<FF: FiniteField>(
        &self,
        participant_shares: Vec<Vec<Share<FF>>>,
        vector_length: usize,
    ) -> ThresholdResult<Vec<ShamirShare<'static, FF>>> {
        let mut shares = Vec::with_capacity(self.participant_number);

        for (pid_offset, shares_for_pid) in
            participant_shares.into_iter().enumerate()
        {
            let polynomials =
                Self::build_share_polynomials(vector_length, &shares_for_pid);
            shares.push(ShamirShare::new(pid_offset + 1, polynomials)?);
        }

        Ok(shares)
    }

    fn build_share_polynomials<FF: FiniteField>(
        vector_length: usize,
        shares_for_pid: &[Share<FF>],
    ) -> Vec<Polynomial<'static, FF>> {
        let required_lengths =
            Self::collect_required_lengths(vector_length, shares_for_pid);
        let mut buffers =
            Self::allocate_polynomial_buffers::<FF>(&required_lengths);

        for share in shares_for_pid {
            buffers[share.poly_idx][share.coeff_idx] = share.value;
        }

        buffers.into_iter().map(Polynomial::from).collect()
    }

    fn collect_required_lengths<FF: FiniteField>(
        vector_length: usize,
        shares_for_pid: &[Share<FF>],
    ) -> Vec<usize> {
        let mut per_poly_len = vec![0usize; vector_length];
        for share in shares_for_pid {
            let required_len = share.coeff_idx + 1;
            per_poly_len[share.poly_idx] =
                per_poly_len[share.poly_idx].max(required_len);
        }
        per_poly_len
    }

    fn allocate_polynomial_buffers<FF: FiniteField>(
        lengths: &[usize],
    ) -> Vec<Vec<FF>> {
        lengths
            .iter()
            .map(|&len| vec![FF::default(); len])
            .collect()
    }

    fn create_shamir_polynomial<FF>(
        &self,
        secret: &FF,
    ) -> Polynomial<'static, FF>
    where
        FF: FiniteField,
        rand::distributions::Standard: rand::distributions::Distribution<FF>,
    {
        let mut rng = rand::thread_rng();
        let coefficients: Vec<FF> = std::iter::once(*secret)
            .chain((1..self.threshold).map(|_| rng.r#gen()))
            .collect();
        Polynomial::from(coefficients)
    }

    /// Evaluate polynomial at given point using Horner's method
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

impl<FF: FiniteField> PointSource<FF> for ShamirShare<'static, FF> {
    fn x(&self) -> FF {
        // Adjust field name as needed: `id` or `participant_id`.
        let pid: usize = self.participant_id; // or `self.id`
        (pid as u64).into()
    }

    fn poly_at(&self, index: usize) -> Option<&Polynomial<'static, FF>> {
        self.share_vector.get(index)
    }

    fn poly_count(&self) -> usize {
        self.vector_length()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    impl Default for AdaptedShamirSSS {
        fn default() -> Self {
            Self {
                threshold: 3,
                participant_number: 5,
            }
        }
    }

    mod shamir_share_tests {
        use super::*;

        #[test]
        fn test_share_creation() {
            let poly1: Polynomial<'_, FieldElement> = poly![1, 2, 3, 4, 5];
            let poly2: Polynomial<'_, FieldElement> = poly![6, 7, 8, 9, 10];
            let share_vector = vec![poly1, poly2];
            let share = ShamirShare::new(1, share_vector.clone()).unwrap();

            assert_eq!(share.participant_id, 1);
            assert_eq!(share.vector_length(), 2);
            assert_eq!(share.share_vector, share_vector);
        }

        #[test]
        fn test_invalid_participant_id() {
            let poly1: Polynomial<'_, FieldElement> = poly![1, 2, 3, 4, 5];
            let share_vector = vec![poly1];

            // Test that participant ID 0 is rejected
            assert!(ShamirShare::new(0, share_vector.clone()).is_err());
        }

        #[test]
        fn test_share_debug_representation() {
            let poly1: Polynomial<'_, FieldElement> = poly![1, 2, 3];
            let share_vector = vec![poly1];
            let share = ShamirShare::new(1, share_vector).unwrap();

            let debug_str = format!("{:?}", share);
            assert!(debug_str.contains("ShamirShare"));
            assert!(debug_str.contains("participant_id: 1"));
        }
    }

    mod adapted_shamir_sss_tests {
        use num_traits::Zero;

        use super::*;

        fn setup_adapted_shamir(
            threshold: usize,
            participants: usize,
        ) -> ThresholdResult<AdaptedShamirSSS> {
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
            let shamir = AdaptedShamirSSS::default();
            let len = shamir
                .ensure_consistent_lengths::<FieldElement>(&[])
                .expect("empty input should succeed");
            assert_eq!(len, 0);
        }

        #[test]
        fn ensure_consistent_lengths_reports_error() {
            let shamir = AdaptedShamirSSS::default();
            let shares = vec![share_with_length(1, 2), share_with_length(2, 3)];

            let err = shamir
                .ensure_consistent_lengths::<FieldElement>(&shares)
                .expect_err("mismatched lengths must fail");
            assert!(matches!(err, ThresholdError::InconsistentShareLengths));
        }

        #[test]
        fn ensure_consistent_lengths_tracks_length() {
            let shamir = AdaptedShamirSSS::default();
            let shares = vec![share_with_length(1, 4), share_with_length(2, 4)];

            let len = shamir
                .ensure_consistent_lengths::<FieldElement>(&shares)
                .expect("consistent lengths should succeed");
            assert_eq!(len, 4);
        }

        #[test]
        fn ensure_consistent_lengths_single_share() {
            let shamir = AdaptedShamirSSS::default();
            let shares = vec![share_with_length(1, 5)];

            let len = shamir
                .ensure_consistent_lengths::<FieldElement>(&shares)
                .expect("single share should succeed");
            assert_eq!(len, 5);
        }

        #[test]
        fn test_invalid_threshold_config() {
            // Threshold too small
            assert!(matches!(
                setup_adapted_shamir(1, 5),
                Err(ThresholdError::InvalidThreshold(1, 5))
            ));
            assert!(matches!(
                setup_adapted_shamir(0, 5),
                Err(ThresholdError::InvalidThreshold(0, 5))
            ));

            // Threshold larger than participants
            assert!(matches!(
                setup_adapted_shamir(6, 5),
                Err(ThresholdError::InvalidThreshold(6, 5))
            ));

            // Too many participants (assuming max is 255 based on typical limits)
            assert!(matches!(
                setup_adapted_shamir(3, 300),
                Err(ThresholdError::InvalidThreshold(3, 300))
            ));
        }

        // TODO fix it
        #[test]
        fn test_secret_splitting() {
            let shamir = AdaptedShamirSSS::default();

            // Create test secret vector
            let coeffs = fe_vec!(1, 2, 3, 4, 5);
            let secret_poly1: Polynomial<'_, FieldElement> =
                poly!(coeffs.clone());
            assert_eq!(secret_poly1.coefficients().len(), coeffs.len());
            let secret_poly2: Polynomial<'_, FieldElement> =
                poly![10, 20, 30, 40, 50];
            assert_eq!(secret_poly2.coefficients().len(), coeffs.len());
            let secret_vector = poly_vec!(secret_poly1, secret_poly2);

            let shares = shamir.split_secret(&secret_vector).unwrap();

            // Check we get correct number of shares
            assert_eq!(shares.len(), shamir.participant_number);

            // Check each share has correct participant ID
            for (i, share) in shares.iter().enumerate() {
                assert_eq!(share.participant_id, i + 1);
                assert_eq!(share.vector_length(), secret_vector.len());
            }
        }

        #[test]
        fn test_secret_reconstruction() {
            let shamir = AdaptedShamirSSS::default();

            // Create test secret vector
            // TODO declare standadlone function
            let secret_poly1: Polynomial<'_, FieldElement> =
                poly![1, 2, 3, 4, 5];
            let secret_poly2: Polynomial<'_, FieldElement> =
                poly![10, 20, 30, 40, 50];
            let secret_vector = vec![secret_poly1, secret_poly2];

            // Split secret
            let shares = shamir.split_secret(&secret_vector).unwrap();

            // Reconstruct using exactly threshold shares
            let reconstructed = shamir
                .reconstruct_secret(&shares[..shamir.threshold])
                .unwrap();

            // Check reconstruction is correct
            assert_eq!(reconstructed, secret_vector);
        }

        #[test]
        fn test_reconstruction_with_more_shares() {
            let threshold = 3;
            let shamir = setup_adapted_shamir(threshold, 5).unwrap();

            let secret_poly: Polynomial<'_, FieldElement> = poly![42, 17, 99];
            let secret_vector = vec![secret_poly];

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
            let shamir = AdaptedShamirSSS::default();

            let secret_poly: Polynomial<'_, FieldElement> = poly![100, 200];
            let secret_vector = poly_vec!(vec![secret_poly]);

            let shares = shamir.split_secret(&secret_vector).unwrap();

            // Try with threshold - 1 shares
            let result =
                shamir.reconstruct_secret(&shares[..shamir.threshold - 1]);
            assert!(matches!(
                result,
                Err(ThresholdError::InsufficientShares(required, provided))
                    if required == shamir.threshold
                        && provided == shamir.threshold - 1
            ));

            // Try with empty slice
            let empty_shares: &[ShamirShare<FieldElement>] = &[];
            assert!(matches!(
                shamir.reconstruct_secret(empty_shares),
                Err(ThresholdError::InsufficientShares(required, 0))
                    if required == shamir.threshold
            ));
        }

        #[test]
        fn test_partial_reconstruction() {
            let threshold = 3;
            let shamir = setup_adapted_shamir(threshold, 5).unwrap();

            // Create test secret vector
            let secret_poly1: Polynomial<'_, FieldElement> =
                poly![1, 2, 3, 4, 5];
            let secret_poly2: Polynomial<'_, FieldElement> =
                poly![10, 20, 30, 40, 50];
            let secret_vector =
                poly_vec!(secret_poly1.clone(), secret_poly2.clone());

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
            let shamir = setup_adapted_shamir(threshold, 5).unwrap();

            // Test with single polynomial
            let single_poly: Polynomial<'_, FieldElement> = poly![42, 17];
            let single_vector = vec![single_poly];
            let shares = shamir.split_secret(&single_vector).unwrap();
            let reconstructed =
                shamir.reconstruct_secret(&shares[..threshold]).unwrap();
            assert_eq!(reconstructed, single_vector);

            // Test with longer vector
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

            // Create zero polynomial
            let zero_poly: Polynomial<'_, FieldElement> =
                poly![FieldElement::default(); N];
            let zero_vector = vec![zero_poly];

            let shares: Vec<ShamirShare<'static, FieldElement>> =
                shamir.split_secret(&zero_vector).unwrap();
            let reconstructed =
                shamir.reconstruct_secret(&shares[..threshold]).unwrap();

            assert_eq!(reconstructed, zero_vector);

            // Check all coefficients are zero
            let reconstructed_poly = reconstructed.get(0).cloned().unwrap();
            assert!(
                reconstructed_poly
                    .coefficients()
                    .iter()
                    .all(|fe| fe.is_zero())
            );
        }

        #[test]
        fn test_edge_case_thresholds() {
            // Minimum threshold (2 out of 2)
            let shamir_min = AdaptedShamirSSS::new(2, 2).unwrap();
            let secret_poly: Polynomial<'_, FieldElement> = poly![100, 200];
            let secret_vector = vec![secret_poly];

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
            let shamir = setup_adapted_shamir(threshold, 3).unwrap();

            // Use specific values that test modular arithmetic
            let poly: Polynomial<'_, FieldElement> = poly![Q - 1, 1, Q - 2]; // Values near modulus
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

            // Try different combinations of threshold shares
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
        use num_traits::Zero;
        // Test with threshold 2 (minimal case)
        let shamir = AdaptedShamirSSS::new(2, 3).unwrap();
        let secret = FieldElement::new(42u32);
        let poly = shamir.create_shamir_polynomial(&secret);

        // Check polynomial has correct length (equal to threshold)
        assert_eq!(
            poly.coefficients().iter().filter(|&c| !c.is_zero()).count(),
            2
        );
        // Check first coefficient is the secret
        assert_eq!(poly.coefficients()[0], secret);

        // Test with threshold 5
        let shamir = AdaptedShamirSSS::new(5, 7).unwrap();
        let secret = FieldElement::new(123);
        let poly = shamir.create_shamir_polynomial(&secret);

        // Check polynomial has correct length
        assert_eq!(
            poly.coefficients().iter().filter(|&c| !c.is_zero()).count(),
            5
        );
        // Check first coefficient is the secret
        assert_eq!(poly.coefficients()[0], secret);

        // Test with zero secret
        let secret = FieldElement::zero();
        let poly = shamir.create_shamir_polynomial(&secret);
        assert_eq!(poly.coefficients()[0], FieldElement::zero());
        assert_eq!(
            poly.coefficients().iter().filter(|&c| !c.is_zero()).count(),
            4
        );

        // Test randomness - create multiple polynomials with same secret
        let secret = FieldElement::new(777);
        let poly1 = shamir.create_shamir_polynomial(&secret);
        let poly2 = shamir.create_shamir_polynomial(&secret);
        let poly3 = shamir.create_shamir_polynomial(&secret);

        // All should have same secret
        assert_eq!(poly1.coefficients()[0], secret);
        assert_eq!(poly2.coefficients()[0], secret);
        assert_eq!(poly3.coefficients()[0], secret);

        // But random coefficients should differ (with high probability)
        // Check at least one differs between poly1 and poly2
        let coeffs_differ = (1..5)
            .any(|i| poly1.coefficients()[i] != poly2.coefficients()[i])
            || (1..5)
                .any(|i| poly1.coefficients()[i] != poly3.coefficients()[i]);
        assert!(
            coeffs_differ,
            "Random coefficients should differ between polynomials"
        );

        // Test with large threshold
        let shamir = AdaptedShamirSSS::new(10, 15).unwrap();
        let secret = FieldElement::new(999);
        let poly = shamir.create_shamir_polynomial(&secret);

        assert_eq!(
            poly.coefficients().iter().filter(|&c| !c.is_zero()).count(),
            10
        );
        assert_eq!(poly.coefficients()[0], secret);

        // Test negative secret (if FieldElement supports it via From<i32>)
        let secret = FieldElement::from(-42i32);
        let poly = shamir.create_shamir_polynomial(&secret);
        assert_eq!(poly.coefficients()[0], secret);
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
            expected
        );

        // Verify Horner's method: f(x) = 5 + 3x + 2x^2 + 4x^3
        let coeffs = vec![5, 3, 2, 4];
        let x = 7;
        let manual_result = (5 + 7 * (3 + 7 * (2 + 7 * 4))) % Q;
        assert_eq!(
            AdaptedShamirSSS::evaluate_polynomial(&coeffs, x),
            manual_result
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
