use math::{prelude::*, traits::FiniteField};

use crate::{
    error::{Result, ThresholdError},
    params::validate_threshold_config,
};

use rand::Rng;

/// Adapted Shamir's Secret Sharing
#[derive(Clone, Debug, PartialEq)]
pub struct ShamirShare<'a, FF: FiniteField> {
    pub participant_id: usize,
    pub share_vector: PolynomialVector<'a, FF>,
}

impl<FF: FiniteField> ShamirShare<'static, FF> {
    pub fn new(
        participant_id: usize,
        share_vector: PolynomialVector<'static, FF>,
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
struct Share<FF: FiniteField> {
    poly_idx: usize,
    coeff_idx: usize,
    value: FF,
}

impl<FF: FiniteField> Share<FF> {
    pub fn new(poly_idx: usize, coeff_idx: usize, value: FF) -> Self {
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
    pub fn split_secret<FF>(
        &self,
        secret_vector: &PolynomialVector<'static, FF>,
    ) -> Result<Vec<ShamirShare<'static, FF>>>
    where
        FF: FiniteField,
        rand::distributions::Standard: rand::distributions::Distribution<FF>,
    {
        let vector_length = secret_vector.len();
        let mut participant_shares: Vec<Vec<Share<FF>>> = vec![
                Vec::with_capacity(vector_length * N);
                self.participant_number
            ];

        for (poly_idx, poly) in secret_vector.as_slice().iter().enumerate() {
            // Process each coefficient
            // for (coeff_idx, secret_coeff) in
            //     poly.coefficients().iter().enumerate()
            for coeff_idx in 0..N {
                let secret_coeff = poly.coefficients().get(coeff_idx).unwrap();
                let shamir_poly = self.create_shamir_polynomial(secret_coeff);

                // Evaluate for each participant
                for pid in 1..=self.participant_number {
                    let field_element = FF::from(pid as u64);
                    let share_value =
                        shamir_poly.evaluate_in_same_field(field_element);
                    //   Self::evaluate_polynomial(&shamir_poly, pid as i32);
                    let share = Share::new(poly_idx, coeff_idx, share_value);
                    participant_shares[pid - 1].push(share);
                }
            }
        }

        // Organize into shares
        self.organize_shares(participant_shares, vector_length)
    }

    /// Reconstruct the secret from a sufficient number of shares.
    // TODO fix code duplciation with partial_reconstruct
    pub fn reconstruct_secret<FF: FiniteField>(
        &self,
        shares: &[ShamirShare<'static, FF>],
    ) -> Result<PolynomialVector<'static, FF>> {
        let active_shares = &shares[..self.threshold];
        let vector_length = active_shares[0].vector_length();

        self.validate_shares(active_shares, shares, vector_length)?;

        let poly_indices: Vec<usize> = (0..vector_length).collect();
        Self::reconstruct_poly_vector(active_shares, &poly_indices)
    }

    #[cfg(test)]
    /// Partially reconstruct only specified polynomials from the vector.
    pub fn partial_reconstruct<FF: FiniteField>(
        &self,
        shares: &[ShamirShare<'static, FF>],
        poly_indices: &[usize],
    ) -> Result<PolynomialVector<'static, FF>> {
        let active_shares = &shares[..self.threshold];
        let vector_length = active_shares[0].vector_length();
        self.validate_shares(active_shares, shares, vector_length)?;

        Self::reconstruct_poly_vector(active_shares, poly_indices)
    }

    /// Common logic for reconstructing polynomials
    fn reconstruct_poly_vector<FF: FiniteField>(
        active_shares: &[ShamirShare<'static, FF>],
        poly_indices: &[usize],
    ) -> Result<PolynomialVector<'static, FF>> {
        let mut reconstructed_polys: Vec<Polynomial<'static, FF>> =
            Vec::with_capacity(poly_indices.len());

        // TODO fix possible code duplication
        for poly_idx in poly_indices {
            let mut xs = Vec::with_capacity(active_shares.len() * N);
            let mut ys = Vec::with_capacity(active_shares.len() * N);
            for i in 0..N {
                for share in active_shares {
                    let (x, y) = Self::yield_points(share, *poly_idx, i)?;
                    xs.push(x);
                    ys.push(y)
                }
            }
            let poly = Polynomial::lagrange_interpolate(&xs, &ys);
            reconstructed_polys.push(poly);
        }

        Ok(poly_vec![reconstructed_polys])
    }

    /// Validate shares for reconstruction
    fn validate_shares<FF: FiniteField>(
        &self,
        active_shares: &[ShamirShare<'static, FF>],
        shares: &[ShamirShare<'static, FF>],
        vector_length: usize,
    ) -> Result<()> {
        if shares.len() < self.threshold {
            return Err(ThresholdError::InsufficientShares {
                required: self.threshold,
                provided: shares.len(),
            });
        }

        // Check if all shares have the same vector length

        if active_shares
            .iter()
            .any(|s| s.vector_length() != vector_length)
        {
            return Err(ThresholdError::InconsistentShareLengths);
        }

        Ok(())
    }

    fn yield_points<FF: FiniteField>(
        share: &ShamirShare<'static, FF>,
        poly_idx: usize,
        coeff_idx: usize,
    ) -> Result<(FF, FF)> {
        dbg!(
            "yield_points, poly_idx: {}, coeff_idx: {}, length: {}",
            poly_idx,
            coeff_idx,
            share.vector_length()
        );
        let poly = share.share_vector.get(poly_idx).ok_or(
            ThresholdError::InvalidIndex {
                index: poly_idx,
                length: share.vector_length(),
            },
        )?;
        dbg!("share_vector.called");
        dbg!("calling coefficients().get");
        let coeff = poly.coefficients().get(coeff_idx).copied().ok_or(
            ThresholdError::InvalidIndex {
                index: coeff_idx,
                length: share.vector_length(),
            },
        )?;

        dbg!("return");
        Ok((FF::from(share.participant_id as u64), coeff))
    }

    /// Organize participant shares into ShamirShare objects
    fn organize_shares<FF: FiniteField>(
        &self,
        participant_shares: Vec<Vec<Share<FF>>>,
        vector_length: usize,
    ) -> Result<Vec<ShamirShare<'static, FF>>> {
        let mut shares = Vec::with_capacity(self.participant_number);

        for pid in 1..=self.participant_number {
            let mut share_coeffs = vec![vec![FF::default(); N]; vector_length];

            // Fill in the coefficients for each polynomial
            for share in &participant_shares[pid - 1] {
                share_coeffs[share.poly_idx][share.coeff_idx] = share.value;
            }
            //let polys = share_coeffs.into_iter().map(Polynomial::from).collect();
            //let poly: Polynomial<'static, FF> = Polynomial::from(&share_coeffs);
            let poly_vec: Vec<Polynomial<'static, FF>> = share_coeffs
                .into_iter()
                .map(|coeffs| Polynomial::from(coeffs))
                .collect();

            let share_vector = PolynomialVector::new(poly_vec);
            shares.push(ShamirShare::new(pid, share_vector)?);
        }

        Ok(shares)
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
            .chain((1..self.threshold).map(|_| rng.gen()))
            .collect();
        Polynomial::from(coefficients)
    }

    /// Evaluate polynomial at given point using Horner's method
    /// // this is an inplementation
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

    mod shamir_share_tests {
        use super::*;

        #[test]
        fn test_share_creation() {
            let poly1: Polynomial<'_, FieldElement> = poly![1, 2, 3, 4, 5];
            let poly2: Polynomial<'_, FieldElement> = poly![6, 7, 8, 9, 10];
            let share_vector = poly_vec!(poly1, poly2);
            let share = ShamirShare::new(1, share_vector.clone()).unwrap();

            assert_eq!(share.participant_id, 1);
            assert_eq!(share.vector_length(), 2);
            assert_eq!(share.share_vector, share_vector);
        }

        #[test]
        fn test_invalid_participant_id() {
            let poly1: Polynomial<'_, FieldElement> = poly![1, 2, 3, 4, 5];
            let share_vector = poly_vec!(vec![poly1]);

            // Test that participant ID 0 is rejected
            assert!(ShamirShare::new(0, share_vector.clone()).is_err());
        }

        #[test]
        fn test_share_debug_representation() {
            let poly1: Polynomial<'_, FieldElement> = poly![1, 2, 3];
            let share_vector = poly_vec!(vec![poly1]);
            let share = ShamirShare::new(1, share_vector).unwrap();

            let debug_str = format!("{:?}", share);
            assert!(debug_str.contains("ShamirShare"));
            assert!(debug_str.contains("participant_id: 1"));
        }
    }

    mod adapted_shamir_sss_tests {
        use num_traits::Zero;

        use super::*;

        impl Default for AdaptedShamirSSS {
            fn default() -> Self {
                Self {
                    threshold: 3,
                    participant_number: 5,
                }
            }
        }

        fn setup_adapted_shamir(
            threshold: usize,
            participants: usize,
        ) -> Result<AdaptedShamirSSS> {
            AdaptedShamirSSS::new(threshold, participants)
        }

        #[test]
        fn test_shamir_initialization() {
            let adapted_shamir_res = setup_adapted_shamir(3, 5);

            // We can't directly access fields in Rust due to privacy,
            // but we can test that creation succeeds
            assert!(adapted_shamir_res.is_ok());
        }

        #[test]
        fn test_invalid_threshold_config() {
            // Threshold too small
            assert!(matches!(
                setup_adapted_shamir(1, 5),
                Err(ThresholdError::InvalidThreshold {
                    threshold: 1,
                    participant_number: 5
                })
            ));
            assert!(matches!(
                setup_adapted_shamir(0, 5),
                Err(ThresholdError::InvalidThreshold {
                    threshold: 0,
                    participant_number: 5
                })
            ));

            // Threshold larger than participants
            assert!(matches!(
                setup_adapted_shamir(6, 5),
                Err(ThresholdError::InvalidThreshold {
                    threshold: 6,
                    participant_number: 5
                })
            ));

            // Too many participants (assuming max is 255 based on typical limits)
            assert!(matches!(
                setup_adapted_shamir(3, 300),
                Err(ThresholdError::InvalidThreshold {
                    threshold: 3,
                    participant_number: 300
                })
            ));
        }

        // TODO fix it
        #[test]
        fn test_secret_splitting() {
            let participants = 5;
            let shamir = setup_adapted_shamir(3, participants).unwrap();

            // Create test secret vector
            let coeffs = fe_vec!(1, 2, 3, 4, 5);
            let secret_poly1: Polynomial<'_, FieldElement> =
                poly!(coeffs);
            assert_eq!(secret_poly1.coefficients().len(), N);
            let secret_poly2: Polynomial<'_, FieldElement> =
                poly![10, 20, 30, 40, 50];
            assert_eq!(secret_poly2.coefficients().len(), N);
            let secret_vector = poly_vec!(secret_poly1, secret_poly2);

            let shares = shamir.split_secret(&secret_vector).unwrap();

            // Check we get correct number of shares
            assert_eq!(shares.len(), participants);

            // Check each share has correct participant ID
            for (i, share) in shares.iter().enumerate() {
                assert_eq!(share.participant_id, i + 1);
                assert_eq!(share.vector_length(), secret_vector.len());
            }
        }

     
        #[test]
        fn test_secret_reconstruction() {
            let shamir = setup_adapted_shamir(3, 5).unwrap();

            // Create test secret vector
            // TODO declare standadlone function
            let secret_poly1: Polynomial<'_, FieldElement> =
                poly![1, 2, 3, 4, 5];
            let secret_poly2: Polynomial<'_, FieldElement> =
                poly![10, 20, 30, 40, 50];
            let secret_vector = poly_vec!(secret_poly1, secret_poly2);

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
            let secret_vector = poly_vec!(vec![secret_poly]);

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
            let shamir = setup_adapted_shamir(threshold, 5).unwrap();

            let secret_poly: Polynomial<'_, FieldElement> = poly![100, 200];
            let secret_vector = poly_vec!(vec![secret_poly]);

            let shares = shamir.split_secret(&secret_vector).unwrap();

            // Try with threshold - 1 shares
            assert!(shamir
                .reconstruct_secret(&shares[..threshold - 1])
                .is_err());

            // Try with empty slice
            let empty_shares: &[ShamirShare<FieldElement>] = &[];
            assert!(shamir.reconstruct_secret(empty_shares).is_err());
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
            let single_vector: PolynomialVector<'_, FieldElement> =
                poly_vec!(vec![single_poly]);
            let shares = shamir.split_secret(&single_vector).unwrap();
            let reconstructed =
                shamir.reconstruct_secret(&shares[..threshold]).unwrap();
            assert_eq!(reconstructed, single_vector);

            // Test with longer vector
            let poly1: Polynomial<'_, FieldElement> = poly![1, 2, 3];
            let poly2: Polynomial<'_, FieldElement> = poly![10, 20, 30];
            let poly3: Polynomial<'_, FieldElement> = poly![100, 200, 300];
            let long_vector = poly_vec!(poly1, poly2, poly3);
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
            let zero_vector: PolynomialVector<'static, FieldElement> =
                poly_vec!(vec![zero_poly]);

            let shares: Vec<ShamirShare<'static, FieldElement>> =
                shamir.split_secret(&zero_vector).unwrap();
            let reconstructed: PolynomialVector<'static, FieldElement> =
                shamir.reconstruct_secret(&shares[..threshold]).unwrap();

            assert_eq!(reconstructed, zero_vector);

            // Check all coefficients are zero
            let reconstructed_poly = reconstructed.get(0).cloned().unwrap();
            assert!(reconstructed_poly
                .coefficients()
                .iter()
                .all(|fe| fe.is_zero()));
        }

        #[test]
        fn test_edge_case_thresholds() {
            // Minimum threshold (2 out of 2)
            let shamir_min = AdaptedShamirSSS::new(2, 2).unwrap();
            let secret_poly: Polynomial<'_, FieldElement> = poly![100, 200];
            let secret_vector = poly_vec!(vec![secret_poly]);

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
            let secret_vector = poly_vec!(vec![poly]);

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

            let poly: Polynomial<'_, FieldElement> = poly![42, 17, 99, 13];
            let secret_vector = poly_vec!(vec![poly]);

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

    // // Integration tests
    // #[test]
    // fn test_full_workflow() {
    //     let threshold = 5;
    //     let participants = 9;
    //     let shamir = AdaptedShamirSSS::new(threshold, participants).unwrap();

    //     // Create a complex secret
    //     let poly1: Polynomial<'_, FieldElement> =
    //         poly![123, 456, 789, 101, 202];
    //     let poly2: Polynomial<'_, FieldElement> =
    //         poly![303, 404, 505, 606, 707];
    //     let poly3: Polynomial<'_, FieldElement> =
    //         poly![808, 909, 111, 222, 333];
    //     let secret: PolynomialVector<'_, FieldElement> =
    //         poly_vec!(poly1, poly2, poly3);

    //     // Split the secret
    //     let shares = shamir.split_secret(&secret).unwrap();
    //     //dbg!("shares: {:?}", &shares);

    //     // Verify we can reconstruct with any threshold shares
    //     for start in 0..=(participants - threshold) {
    //         let selected_shares = &shares[start..start + threshold];
    //         dbg!("for loop, start", {});
    //         let reconstructed =
    //             shamir.reconstruct_secret(selected_shares).unwrap();
    //         assert_eq!(reconstructed, secret);
    //     }
    // }

    // #[test]
    // fn test_full_workflow() {
    //     let threshold = 5;
    //     let participants = 9;
    //     let shamir = AdaptedShamirSSS::new(threshold, participants).unwrap();

    //     // Create a complex secret with properly sized polynomials
    //     // Each polynomial should have N coefficients (where N is 256)
    //     let mut coeffs1 = vec![0i32; N];
    //     coeffs1[0] = 123;
    //     coeffs1[1] = 456;
    //     coeffs1[2] = 789;
    //     coeffs1[3] = 101;
    //     coeffs1[4] = 202;

    //     let mut coeffs2 = vec![0i32; N];
    //     coeffs2[0] = 303;
    //     coeffs2[1] = 404;
    //     coeffs2[2] = 505;
    //     coeffs2[3] = 606;
    //     coeffs2[4] = 707;

    //     let mut coeffs3 = vec![0i32; N];
    //     coeffs3[0] = 808;
    //     coeffs3[1] = 909;
    //     coeffs3[2] = 111;
    //     coeffs3[3] = 222;
    //     coeffs3[4] = 333;

    //     let poly1: Polynomial<'_, FieldElement> = poly![coeffs1];
    //     let poly2: Polynomial<'_, FieldElement> = poly![coeffs2];
    //     let poly3: Polynomial<'_, FieldElement> = poly![coeffs3];
    //     let secret: PolynomialVector<'_, FieldElement> =
    //         poly_vec!(poly1, poly2, poly3);

    //     // Split the secret
    //     let shares = shamir.split_secret(&secret).unwrap();

    //     // Verify we can reconstruct with any threshold shares
    //     for start in 0..=(participants - threshold) {
    //         let selected_shares = &shares[start..start + threshold];
    //         let reconstructed =
    //             shamir.reconstruct_secret(selected_shares).unwrap();
    //         assert_eq!(reconstructed, secret);
    //     }
    // }

    #[test]
    fn test_create_shamir_polynomial() {
        use num_traits::Zero;
        // Test with threshold 2 (minimal case)
        let shamir = AdaptedShamirSSS::new(2, 3).unwrap();
        let secret = FieldElement::new(42u32);
        let poly = shamir.create_shamir_polynomial(&secret);

        // Check polynomial has correct length (equal to threshold)
        assert_eq!(poly.coefficients().len(), 2);
        // Check first coefficient is the secret
        assert_eq!(poly.coefficients()[0], secret);

        // Test with threshold 5
        let shamir = AdaptedShamirSSS::new(5, 7).unwrap();
        let secret = FieldElement::new(123);
        let poly = shamir.create_shamir_polynomial(&secret);

        // Check polynomial has correct length
        assert_eq!(poly.coefficients().len(), 5);
        // Check first coefficient is the secret
        assert_eq!(poly.coefficients()[0], secret);

        // Test with zero secret
        let secret = FieldElement::zero();
        let poly = shamir.create_shamir_polynomial(&secret);
        assert_eq!(poly.coefficients()[0], FieldElement::zero());
        assert_eq!(poly.coefficients().len(), 5);

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

        assert_eq!(poly.coefficients().len(), 10);
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
            expected as i32
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
