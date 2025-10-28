//! Small helpers to keep code DRY and straightforward.

use num_traits::Zero;
use rand::RngCore;

use sha2::digest::core_api::XofReaderCoreWrapper;
use sha3::{
    Shake256, Shake256ReaderCore,
    digest::{ExtendableOutput, Update, XofReader},
};

use crate::error::{Result, ThresholdError};
use crate::params::GAMMA1;
use crate::points::PointSource;
use math::{poly::Polynomial, prelude::*, traits::FiniteField};

// Fill byte array of length 32 by random bytes
pub fn random_bytes() -> [u8; 32] {
    let mut rng = rand::thread_rng();
    let mut tmp = [0u8; 32];
    rng.fill_bytes(&mut tmp);
    tmp
}

pub fn get_randomness(randomness: Option<&[u8]>) -> Vec<u8> {
    match randomness {
        Some(r) => r.to_vec(),
        None => {
            let mut rng = rand::thread_rng();
            let mut bytes = vec![0u8; 32];
            rng.fill_bytes(&mut bytes);
            bytes
        }
    }
}

pub fn get_hash_reader(
    message: &[u8],
) -> XofReaderCoreWrapper<Shake256ReaderCore> {
    let mut hasher = Shake256::default();
    hasher.update(message);
    hasher.finalize_xof()
}

pub fn hash_message(message: &[u8]) -> Vec<u8> {
    let mut reader = get_hash_reader(message);
    let mut output = vec![0u8; 64];
    reader.read(&mut output);
    output
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

// Generic reconstruction from any point providers (removes duplication).
pub fn reconstruct_vector_from_points<FF, S>(
    items: &[S],
    poly_indices: &[usize],
) -> Result<Vec<Polynomial<'static, FF>>>
where
    FF: FiniteField,
    S: PointSource<FF>,
{
    if items.is_empty() {
        return Err(ThresholdError::InsufficientShares {
            required: 1,
            provided: 0,
        });
    }

    //assert_eq!(L, poly_indices.len());
    let mut reconstructed = Vec::with_capacity(poly_indices.len());
    let vector_len = items[0].poly_count();

    for &poly_idx in poly_indices {
        let coeff_len = items[0]
            .poly_at(poly_idx)
            .map(|p| p.coefficients().len())
            .ok_or(ThresholdError::InvalidIndex {
                index: poly_idx,
                length: vector_len,
            })?;
        let mut coeffs = vec![FF::ZERO; coeff_len];

        let mut coeffs = vec![FF::ZERO; coeff_len];

        for (coeff_idx, c) in coeffs.iter_mut().enumerate().take(coeff_len) {
            let mut xs = Vec::with_capacity(items.len());
            let mut ys = Vec::with_capacity(items.len());

            for it in items {
                let x = it.x();
                let poly = it.poly_at(poly_idx).ok_or(
                    ThresholdError::InvalidIndex {
                        index: poly_idx,
                        length: vector_len,
                    },
                )?;
                let y = poly.coefficients().get(coeff_idx).copied().ok_or(
                    ThresholdError::InvalidIndex {
                        index: coeff_idx,
                        length: poly.coefficients().len(),
                    },
                )?;
                xs.push(x);
                ys.push(y);
            }

            *c = interpolate_constant_at_zero(&xs, &ys);
        }

        reconstructed.push(Polynomial::from(coeffs));
    }
    // let arr: [Polynomial<'static, FF>; L]  = reconstructed.try_into()
    //     .unwrap_or_else(|v: Vec<Polynomial<'static, FF>>| panic!("Expected a Vec of length {} but it was {}", N, v.len()));
    Ok(reconstructed)
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
pub fn zero_polyvec<const LEN: usize, FF: FiniteField>()
-> [Polynomial<'static, FF>; LEN] {
    std::array::from_fn(|_| Polynomial::zero())
}

/// Return the maximum \ell_\infty norm across a slice of polynomials.
/// Useful for bound checks; returns 0 for an empty slice.
pub fn polyvec_max_infty_norm<FF: FiniteField>(
    polys: &[Polynomial<'_, FF>],
) -> i64
where
    i64: std::convert::From<FF>,
{
    polys
        .iter()
        .map(|p| p.norm_infinity() as i64)
        .max()
        .unwrap_or(0)
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
