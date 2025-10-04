//! Small helpers to keep code DRY and straightforward.

use num_traits::Zero;
use rand::RngCore;

use sha2::digest::core_api::XofReaderCoreWrapper;
use sha3::{
    Shake256, Shake256ReaderCore,
    digest::{ExtendableOutput, Update, XofReader},
};

use crate::error::{Result, ThresholdError};
use crate::params::L;
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
