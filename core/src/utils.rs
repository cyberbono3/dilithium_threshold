use rand::prelude::*;
use sha2::digest::core_api::XofReaderCoreWrapper;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256, Shake256ReaderCore,
};

use crate::error::{Result, ThresholdError};
use crate::points::PointSource;
use math::{prelude::*, traits::FiniteField};

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

/// Lagrange interpolate over (xs, ys) and return f(0).
pub fn interpolate_constant_at_zero<FF: FiniteField + Copy + 'static>(
    xs: &[FF],
    ys: &[FF],
) -> FF {
    let f = Polynomial::lagrange_interpolate(xs, ys);
    f.batch_evaluate(&[FF::ZERO])
        .first()
        .copied()
        .expect("interpolation requires at least one (x, y) pair")
}

/// Generic reconstruction from any point providers (removes duplication).
pub fn reconstruct_vector_from_points<FF, S>(
    items: &[S],
    poly_indices: &[usize],
) -> Result<PolynomialVector<'static, FF>>
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

    let mut reconstructed = Vec::with_capacity(poly_indices.len());
    let vector_len = items[0].poly_count();

    for &poly_idx in poly_indices {
        let coeff_len = items[0]
            .poly_at(poly_idx)
            .map(|p| p.coefficients().len())
            .unwrap_or(0);

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

    Ok(poly_vec!(reconstructed))
}

// TODO icnrease test coeverage 
#[cfg(test)]
mod tests {
    use super::*;

    mod get_randomness_tests {
        use super::*;

        #[test]
        fn test_get_randomness_with_provided_bytes() {
            let input = vec![1, 2, 3, 4, 5];
            let result = get_randomness(Some(&input));
            assert_eq!(result, input);
        }

        #[test]
        fn test_get_randomness_with_empty_slice() {
            let input = vec![];
            let result = get_randomness(Some(&input));
            assert_eq!(result, vec![]);
        }

        #[test]
        fn test_get_randomness_without_input() {
            let result1 = get_randomness(None);
            let result2 = get_randomness(None);

            // Should return 32 bytes
            assert_eq!(result1.len(), 32);
            assert_eq!(result2.len(), 32);

            // Should be different (extremely unlikely to be the same)
            assert_ne!(result1, result2);
        }

        #[test]
        fn test_get_randomness_deterministic_with_seed() {
            let seed = vec![42; 32];
            let result1 = get_randomness(Some(&seed));
            let result2 = get_randomness(Some(&seed));

            // Same input should give same output
            assert_eq!(result1, result2);
            assert_eq!(result1, seed);
        }
    }

    // Integration tests
    mod integration_tests {
        use super::*;

        #[test]
        fn test_randomness_and_hashing_together() {
            // Test that randomness generation and hashing work well together
            let random_bytes = get_randomness(None);
            let hash1 = hash_message(&random_bytes);
            let hash2 = hash_message(&random_bytes);

            // Same input should give same hash
            assert_eq!(hash1, hash2);

            // Different random bytes should give different hash
            let other_random = get_randomness(None);
            let other_hash = hash_message(&other_random);
            assert_ne!(hash1, other_hash);
        }
    }
}
