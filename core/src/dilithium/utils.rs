//! Small helpers to keep code DRY and straightforward.

use num_traits::Zero;
use rand::RngCore;

use sha2::digest::core_api::XofReaderCoreWrapper;
use sha3::{
    Shake256, Shake256ReaderCore,
    digest::{ExtendableOutput, Update, XofReader},
};

use crate::dilithium::error::{DilithiumError, DilithiumResult};
use crate::dilithium::params::{GAMMA1, N, TAU};
use crate::dilithium::shamir::error::ShamirError;
use crate::matrix::hash::shake256;
use crate::traits::PointSource;
use math::{poly::Polynomial, traits::FiniteField};

const RANDOMNESS_BYTES: usize = 32;
const HASH_BYTES: usize = 64;

/// Produce a fresh 32-byte array filled with cryptographically secure random data.
pub fn random_bytes() -> [u8; 32] {
    let mut rng = rand::thread_rng();
    let mut tmp = [0u8; 32];
    rng.fill_bytes(&mut tmp);
    tmp
}

/// Either clone the caller-provided randomness or synthesize a fresh buffer.
pub fn get_randomness(randomness: Option<&[u8]>) -> Vec<u8> {
    match randomness {
        Some(bytes) => bytes.to_vec(),
        None => {
            let mut out = vec![0u8; RANDOMNESS_BYTES];
            rand::thread_rng().fill_bytes(&mut out);
            out
        }
    }
}

/// Obtain a SHAKE256 reader over the supplied message bytes.
pub fn get_hash_reader(
    message: &[u8],
) -> XofReaderCoreWrapper<Shake256ReaderCore> {
    let mut hasher = Shake256::default();
    hasher.update(message);
    hasher.finalize_xof()
}

/// Hash a message with SHAKE256, returning a fixed-size digest.
pub fn hash_message(message: &[u8]) -> Vec<u8> {
    let mut reader = get_hash_reader(message);
    let mut output = vec![0u8; HASH_BYTES];
    reader.read(&mut output);
    output
}

/// Convenience wrapper for SHAKE256 that concatenates `seed` with additional parts.
pub fn shake256_squeezed(
    seed: &[u8],
    parts: &[&[u8]],
    out_len: usize,
) -> Vec<u8> {
    let extra: usize = parts.iter().map(|part| part.len()).sum();
    let mut input = Vec::with_capacity(seed.len() + extra);
    input.extend_from_slice(seed);
    for part in parts {
        input.extend_from_slice(part);
    }
    shake256(out_len, &input)
}

/// Sample polynomial with coefficients centered in [-GAMMA1, GAMMA1].
pub fn sample_gamma1<FF: FiniteField>(seed: &[u8]) -> Polynomial<'static, FF> {
    let mut reader = get_hash_reader(seed);
    let mut bytes = [0u8; N * 4];
    reader.read(&mut bytes);

    let coeffs: Vec<FF> = bytes
        .chunks_exact(4)
        .take(N)
        .map(|chunk| {
            let word =
                u32::from_le_bytes(chunk.try_into().expect("chunk length"));
            let sample = word % (2 * GAMMA1 as u32 + 1);
            let centered = sample as i32 - GAMMA1 as i32;
            <i32 as Into<FF>>::into(centered)
        })
        .collect();

    Polynomial::from(coeffs)
}

/// Sample multiple polynomials centered in [-GAMMA1, GAMMA1] using namespaced seeds.
pub fn sample_gamma1_vector<FF: FiniteField>(
    seed: &[u8],
    count: usize,
) -> Vec<Polynomial<'static, FF>> {
    (0..count)
        .map(|index| {
            let mut namespaced = Vec::with_capacity(seed.len() + 1);
            namespaced.extend_from_slice(seed);
            namespaced.push(index as u8);
            sample_gamma1(&namespaced)
        })
        .collect()
}

/// Derive a TAU-sparse challenge polynomial from the provided seed bytes.
pub fn derive_challenge_polynomial<FF: FiniteField + From<i64>>(
    seed: &[u8],
) -> Polynomial<'static, FF> {
    let mut stream = shake256(4 * TAU + 1024, seed);
    let mut used = vec![false; N];
    let mut coeffs = vec![0i64; N];
    let mut filled = 0usize;
    let mut idx = 0usize;

    while filled < TAU {
        if idx + 3 > stream.len() {
            stream.extend_from_slice(&shake256(1024, &stream));
        }

        let position =
            u16::from_le_bytes([stream[idx], stream[idx + 1]]) as usize % N;
        let sign = if stream[idx + 2] & 1 == 1 { 1 } else { -1 };
        idx += 3;

        if !used[position] {
            used[position] = true;
            coeffs[position] = sign;
            filled += 1;
        }
    }

    coeffs.into()
}

/// Reconstruct the requested polynomial indices from a set of point providers.
pub fn reconstruct_vector_from_points<FF, S>(
    items: &[S],
    poly_indices: &[usize],
) -> DilithiumResult<Vec<Polynomial<'static, FF>>>
where
    FF: FiniteField,
    S: PointSource<FF>,
{
    let reference = items
        .first()
        .ok_or(DilithiumError::InsufficientShares(1, 0))?;
    let vector_len = reference.poly_count();
    let xs: Vec<FF> = items.iter().map(|item| item.x()).collect();

    poly_indices
        .iter()
        .map(|&poly_idx| {
            let template = reference
                .poly_at(poly_idx)
                .ok_or(ShamirError::InvalidIndex(poly_idx, vector_len))?;
            let coeff_count = template.coefficients().len();

            let polys = items
                .iter()
                .map(|item| {
                    item.poly_at(poly_idx).ok_or({
                        ShamirError::InvalidIndex(poly_idx, vector_len)
                    })
                })
                .collect::<Result<Vec<_>, ShamirError>>()?;

            let coeffs = (0..coeff_count)
                .map(|coeff_idx| {
                    let ys = polys
                        .iter()
                        .map(|poly| {
                            let coeffs = poly.coefficients();
                            coeffs.get(coeff_idx).copied().ok_or({
                                ShamirError::InvalidIndex(
                                    coeff_idx,
                                    coeffs.len(),
                                )
                            })
                        })
                        .collect::<Result<Vec<_>, ShamirError>>()?;

                    Ok(interpolate_constant_at_zero(
                        xs.as_slice(),
                        ys.as_slice(),
                    ))
                })
                .collect::<DilithiumResult<Vec<_>>>()?;

            Ok(Polynomial::from(coeffs))
        })
        .collect()
}

/// Lagrange interpolate over (xs, ys) and return f(0).
pub fn interpolate_constant_at_zero<FF: FiniteField + Copy + 'static>(
    xs: &[FF],
    ys: &[FF],
) -> FF {
    let f = Polynomial::lagrange_interpolate(xs, ys);
    f.batch_evaluate(&[FF::ZERO])
        .first()
        .copied()
        .unwrap_or(FF::ZERO)
}

#[inline]
/// Create an array of zero-initialised polynomials of compile-time length.
pub fn zero_polyvec<const LEN: usize, FF: FiniteField>()
-> [Polynomial<'static, FF>; LEN] {
    std::array::from_fn(|_| Polynomial::zero())
}

/// Return the maximum ℓ∞ norm across a slice of polynomials (0 if empty).
pub fn polyvec_max_infty_norm<FF: FiniteField>(
    polys: &[Polynomial<'_, FF>],
) -> i64 {
    polys
        .iter()
        .map(|p| p.norm_infinity())
        .max()
        .unwrap_or_default()
        .into()
}

#[cfg(test)]
mod tests {
    use super::GAMMA1;
    use super::*;
    use math::{
        field_element::FieldElement,
        prelude::{N, Polynomial},
    };

    fn centered_i64(value: FieldElement) -> i64 {
        let modulus = FieldElement::P as i64;
        let raw = value.value() as i64;
        if raw > modulus / 2 {
            raw - modulus
        } else {
            raw
        }
    }

    fn coefficients_as_centered_i64(
        poly: &Polynomial<'_, FieldElement>,
    ) -> Vec<i64> {
        poly.coefficients()
            .iter()
            .map(|&coeff| centered_i64(coeff))
            .collect()
    }

    mod get_randomness_tests {
        use super::*;

        #[test]
        fn returns_supplied_randomness() {
            let input = [0xAAu8; 48];
            let result = get_randomness(Some(&input));
            assert_eq!(result, input);
        }

        #[test]
        fn generates_32_bytes_when_none() {
            let random = get_randomness(None);
            assert_eq!(random.len(), 32);
        }

        #[test]
        fn consecutive_none_calls_produce_distinct_vectors() {
            let first = get_randomness(None);
            let second = get_randomness(None);
            assert_eq!(first.len(), 32);
            assert_eq!(second.len(), 32);
            assert_ne!(first, second);
        }
    }

    mod get_hash_reader_tests {
        use super::*;

        #[test]
        fn deterministic_for_identical_input() {
            let mut reader_a = get_hash_reader(b"dilithium");
            let mut reader_b = get_hash_reader(b"dilithium");

            let mut out_a = [0u8; 64];
            let mut out_b = [0u8; 64];
            reader_a.read(&mut out_a);
            reader_b.read(&mut out_b);

            assert_eq!(out_a, out_b);
        }

        #[test]
        fn different_inputs_produce_distinct_streams() {
            let mut reader_a = get_hash_reader(b"message-a");
            let mut reader_b = get_hash_reader(b"message-b");

            let mut out_a = [0u8; 64];
            let mut out_b = [0u8; 64];
            reader_a.read(&mut out_a);
            reader_b.read(&mut out_b);

            assert_ne!(out_a, out_b);
        }

        #[test]
        fn consecutive_reads_concatenate_correctly() {
            let mut reader = get_hash_reader(b"chunked-read");
            let mut first = [0u8; 32];
            let mut second = [0u8; 32];
            reader.read(&mut first);
            reader.read(&mut second);

            let mut combined = [0u8; 64];
            combined[..32].copy_from_slice(&first);
            combined[32..].copy_from_slice(&second);

            let mut single = get_hash_reader(b"chunked-read");
            let mut full = [0u8; 64];
            single.read(&mut full);

            assert_eq!(combined, full);
        }
    }

    mod hash_message_tests {
        use super::*;

        #[test]
        fn produces_fixed_length_output() {
            let digest = hash_message(b"dilithium-threshold");
            assert_eq!(digest.len(), HASH_BYTES);
        }

        #[test]
        fn deterministic_for_same_message() {
            let first = hash_message(b"deterministic");
            let second = hash_message(b"deterministic");
            assert_eq!(first, second);
        }

        #[test]
        fn distinct_messages_hash_differently() {
            let first = hash_message(b"message-one");
            let second = hash_message(b"message-two");
            assert_ne!(first, second);
        }
    }

    mod reconstruct_vector_from_points_tests {
        use super::*;
        use crate::dilithium::error::DilithiumError;
        use crate::dilithium::shamir::error::ShamirError;
        use math::fe;

        #[derive(Clone)]
        struct MockSource {
            x: FieldElement,
            polys: Vec<Polynomial<'static, FieldElement>>,
        }

        impl PointSource<FieldElement> for MockSource {
            fn x(&self) -> FieldElement {
                self.x
            }

            fn poly_at(
                &self,
                index: usize,
            ) -> Option<&Polynomial<'static, FieldElement>> {
                self.polys.get(index)
            }

            fn poly_count(&self) -> usize {
                self.polys.len()
            }
        }

        #[test]
        fn reconstructs_polynomials_from_points() {
            let poly0 = Polynomial::from(vec![fe!(1), fe!(-2), fe!(3)]);
            let poly1 = Polynomial::from(vec![fe!(4)]);

            let sources = vec![
                MockSource {
                    x: fe!(1),
                    polys: vec![poly0.clone(), poly1.clone()],
                },
                MockSource {
                    x: fe!(2),
                    polys: vec![poly0.clone(), poly1.clone()],
                },
            ];

            let reconstructed =
                reconstruct_vector_from_points::<FieldElement, _>(
                    &sources,
                    &[0, 1],
                )
                .expect("reconstruction should succeed");

            assert_eq!(reconstructed.len(), 2);
            assert_eq!(reconstructed[0], poly0);
            assert_eq!(reconstructed[1], poly1);
        }

        #[test]
        fn errors_on_empty_sources() {
            let result = reconstruct_vector_from_points::<FieldElement, _>(
                &Vec::<MockSource>::new(),
                &[0],
            );

            assert!(matches!(
                result,
                Err(DilithiumError::InsufficientShares(required, provided))
                if required == 1 && provided == 0
            ));
        }

        #[test]
        fn errors_when_poly_missing_in_source() {
            let base_poly = Polynomial::from(vec![fe!(7)]);

            let sources = vec![
                MockSource {
                    x: fe!(1),
                    polys: vec![base_poly.clone(), base_poly.clone()],
                },
                MockSource {
                    x: fe!(2),
                    polys: vec![base_poly.clone()],
                },
            ];

            let result = reconstruct_vector_from_points::<FieldElement, _>(
                &sources,
                &[0, 1],
            );

            assert!(matches!(
                result,
                Err(DilithiumError::Shamir(
                    ShamirError::InvalidIndex(index, len)
                ))
                if index == 1 && len == 2
            ));
        }

        #[test]
        fn errors_when_coefficients_missing() {
            let poly_full = Polynomial::from(vec![fe!(1), fe!(2)]);
            let poly_short = Polynomial::from(vec![fe!(1)]);

            let sources = vec![
                MockSource {
                    x: fe!(1),
                    polys: vec![poly_full.clone()],
                },
                MockSource {
                    x: fe!(2),
                    polys: vec![poly_short.clone()],
                },
            ];

            let result = reconstruct_vector_from_points::<FieldElement, _>(
                &sources,
                &[0],
            );

            assert!(matches!(
                result,
                Err(DilithiumError::Shamir(
                    ShamirError::InvalidIndex(index, len)
                ))
                if index == 1 && len == poly_short.coefficients().len()
            ));
        }
    }

    mod interpolate_constant_at_zero_tests {
        use super::*;
        use math::fe;

        #[test]
        fn interpolates_constant_polynomial() {
            let xs = [fe!(1), fe!(2), fe!(3)];
            let ys = [fe!(5), fe!(5), fe!(5)];
            let result = interpolate_constant_at_zero(&xs, &ys);
            assert_eq!(result, fe!(5));
        }

        #[test]
        fn interpolates_linear_polynomial_at_zero() {
            let xs = [fe!(1), fe!(2)];
            let ys = [fe!(3), fe!(5)]; // corresponds to f(x) = 2x + 1
            let result = interpolate_constant_at_zero(&xs, &ys);
            assert_eq!(result, fe!(1));
        }

        #[test]
        fn interpolates_quadratic_polynomial_at_zero() {
            let xs = [fe!(1), fe!(2), fe!(3)];
            let ys = [fe!(2), fe!(5), fe!(10)]; // f(x) = x^2 + 1
            let result = interpolate_constant_at_zero(&xs, &ys);
            assert_eq!(result, fe!(1));
        }
    }

    mod polyvec_max_infty_norm_tests {
        use super::*;
        use math::field_element::FieldElement;

        fn poly(coeffs: &[i64]) -> Polynomial<'static, FieldElement> {
            Polynomial::from(
                coeffs
                    .iter()
                    .map(|&c| FieldElement::from(c))
                    .collect::<Vec<_>>(),
            )
        }

        #[test]
        fn returns_zero_for_empty_slice() {
            let polys: Vec<Polynomial<'static, FieldElement>> = Vec::new();
            assert_eq!(polyvec_max_infty_norm(&polys), 0);
        }

        #[test]
        fn computes_infinity_norm_for_single_poly() {
            let poly = poly(&[3, -7, 5]);
            assert_eq!(polyvec_max_infty_norm(&[poly]), 7);
        }

        #[test]
        fn takes_max_across_polynomials() {
            let polys = [poly(&[1, 2]), poly(&[-8, 3]), poly(&[4, -6])];
            assert_eq!(polyvec_max_infty_norm(&polys), 8);
        }
    }

    #[test]
    fn sample_gamma1_empty_seed_matches_known_prefix() {
        let polynomial = sample_gamma1::<FieldElement>(b"");
        let coeffs = coefficients_as_centered_i64(&polynomial);

        assert_eq!(coeffs.len(), N);

        let expected_prefix =
            [-20913, -23768, 65620, 79162, 37547, 104412, 62380, 51601];
        assert_eq!(&coeffs[..expected_prefix.len()], expected_prefix);
    }

    #[test]
    fn sample_gamma1_is_deterministic_for_reused_seed() {
        let seed = b"dilithium-deterministic-seed";
        let first = sample_gamma1::<FieldElement>(seed);
        let second = sample_gamma1::<FieldElement>(seed);
        assert_eq!(first, second);
    }

    #[test]
    fn sample_gamma1_coefficients_respect_gamma1_bounds() {
        let seed = [0u8; 32];
        let polynomial = sample_gamma1::<FieldElement>(&seed);

        assert!(polynomial.coefficients().len() <= N);

        let bound = GAMMA1;
        let infinity_norm = polynomial.norm_infinity() as i64;
        assert!(
            infinity_norm <= bound,
            "∞-norm {infinity_norm} exceeds bound {bound}"
        );

        for coeff in polynomial.coefficients() {
            let value = centered_i64(*coeff);
            assert!(
                (-bound..=bound).contains(&value),
                "coefficient {value} outside [-{bound}, {bound}]"
            );
        }
    }

    #[test]
    fn sample_gamma1_produces_different_outputs_for_distinct_seeds() {
        let first = sample_gamma1::<FieldElement>(b"seed-one");
        let second = sample_gamma1::<FieldElement>(b"seed-two");
        assert_ne!(first, second);
    }
}
