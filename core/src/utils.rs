use rand::prelude::*;
use sha2::digest::core_api::XofReaderCoreWrapper;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256, Shake256ReaderCore,
};

use crate::error::{Result, ThresholdError};
use math::poly::Q;

pub fn get_randomness(randomness: Option<&[u8]>) -> Vec<u8> {
    match randomness {
        Some(r) => r.to_vec(),
        None => {
            let mut rng = rand::rng();
            let mut bytes = vec![0u8; 32];
            rng.fill_bytes(&mut bytes);
            bytes
        }
    }
}

pub fn get_hash_reader(message: &[u8]) -> XofReaderCoreWrapper<Shake256ReaderCore>{
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

/// Compute (base^exp) mod Q using fast exponentiation
fn mod_pow(base: i32, mut exp: i32) -> i32 {
    let mut result = 1i64;
    let mut base = base as i64;
    let q = Q as i64;

    while exp > 0 {
        if exp & 1 == 1 {
            result = (result * base) % q;
        }
        base = (base * base) % q;
        exp >>= 1;
    }

    result as i32
}

/// Compute modular inverse using Fermat's little theorem
fn mod_inverse(a: i32) -> Result<i32> {
    if a == 0 {
        return Err(ThresholdError::ModularInverseError);
    }

    // Using Fermat's little theorem: a^(p-2) ≡ a^(-1) (mod p)
    let result = mod_pow(a, Q - 2);

    // Verify the result
    if (result as i64 * a as i64).rem_euclid(Q as i64) != 1 {
        return Err(ThresholdError::ModularInverseError);
    }

    Ok(result)
}

fn mod_mul(a: i32, b: i32) -> i32 {
    // Convert to i64 to prevent overflow during multiplication
    let a = a as i64;
    let b = b as i64;
    let q = Q as i64;

    // Ensure inputs are in range [0, Q)
    let a = a.rem_euclid(q);
    let b = b.rem_euclid(q);

    // TODO add Montgomery multiplication
    // For Q = 8380417, direct multiplication fits in i64
    ((a * b) % q) as i32
}

fn mod_mul_three(a: i32, b: i32, c: i32) -> i32 {
    // First multiply a * b mod Q
    let ab = mod_mul(a, b);
    // Then multiply result by c mod Q
    mod_mul(ab, c)
}

/// Perform Lagrange interpolation.
pub fn lagrange_interpolation(points: &[(i32, i32)], x: i32) -> Result<i32> {
    let mut result = 0i64;
    let n = points.len();

    for i in 0..n {
        let (xi, yi) = points[i];

        // Compute Lagrange basis polynomial L_i(x)
        let mut numerator = 1i64;
        let mut denominator = 1i64;

        for j in 0..n {
            if i != j {
                let (xj, _) = points[j];
                numerator = (numerator * (x - xj) as i64).rem_euclid(Q as i64);
                denominator =
                    (denominator * (xi - xj) as i64).rem_euclid(Q as i64);
            }
        }

        // Compute modular inverse using Fermat's little theorem
        let denominator_inv = mod_inverse(denominator as i32)?;
        // Add contribution
        let contribution = mod_mul_three(yi, numerator as i32, denominator_inv);
        result = (result + contribution as i64).rem_euclid(Q as i64);
    }

    Ok(result as i32)
}

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

    mod mod_pow_tests {
        use super::*;

        #[test]
        fn test_mod_pow_basic() {
            assert_eq!(mod_pow(2, 3), 8); // 2^3 = 8
            assert_eq!(mod_pow(3, 4), 81); // 3^4 = 81
            assert_eq!(mod_pow(5, 2), 25); // 5^2 = 25
        }

        #[test]
        fn test_mod_pow_zero_exponent() {
            assert_eq!(mod_pow(12345, 0), 1); // Any number^0 = 1
            assert_eq!(mod_pow(0, 0), 1); // 0^0 = 1 (by convention in modular arithmetic)
        }

        #[test]
        fn test_mod_pow_one_base() {
            assert_eq!(mod_pow(1, 100), 1); // 1^100 = 1
            assert_eq!(mod_pow(1, Q - 1), 1); // 1^(Q-1) = 1
        }

        #[test]
        fn test_mod_pow_large_values() {
            // Test with values that would overflow without modular reduction
            let base = Q - 1;
            let exp = 2;
            let result = mod_pow(base, exp);
            let expected = ((base as i64 * base as i64) % Q as i64) as i32;
            assert_eq!(result, expected);
        }

        #[test]
        fn test_mod_pow_fermats_little_theorem() {
            // For prime Q, a^(Q-1) ≡ 1 (mod Q) for all a not divisible by Q
            assert_eq!(mod_pow(2, Q - 1), 1);
            assert_eq!(mod_pow(3, Q - 1), 1);
            assert_eq!(mod_pow(12345, Q - 1), 1);
        }
    }

    mod mod_inverse_tests {
        use super::*;

        #[test]
        fn test_mod_inverse_basic() {
            // Test that a * a^(-1) ≡ 1 (mod Q)
            let a = 12345;
            let inv = mod_inverse(a).unwrap();
            assert_eq!((a as i64 * inv as i64).rem_euclid(Q as i64), 1);
        }

        #[test]
        fn test_mod_inverse_zero() {
            // Zero has no modular inverse
            assert!(mod_inverse(0).is_err());
        }

        #[test]
        fn test_mod_inverse_one() {
            // 1 is its own inverse
            assert_eq!(mod_inverse(1).unwrap(), 1);
        }

        #[test]
        fn test_mod_inverse_minus_one() {
            // (Q-1) is -1 in modular arithmetic
            let inv = mod_inverse(Q - 1).unwrap();
            assert_eq!(inv, Q - 1); // -1 is its own inverse
        }

        #[test]
        fn test_mod_inverse_various_values() {
            let test_values = vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31];

            for &a in &test_values {
                let inv = mod_inverse(a).unwrap();
                let product = (a as i64 * inv as i64).rem_euclid(Q as i64);
                assert_eq!(product, 1, "Failed for a = {}", a);
            }
        }

        #[test]
        fn test_mod_inverse_large_values() {
            let large_values = vec![Q - 2, Q - 3, Q - 100, Q / 2, Q / 3];

            for &a in &large_values {
                let inv = mod_inverse(a).unwrap();
                let product = (a as i64 * inv as i64).rem_euclid(Q as i64);
                assert_eq!(product, 1, "Failed for a = {}", a);
            }
        }
    }

    mod mod_mul_tests {
        use super::*;

        #[test]
        fn test_mod_mul_basic() {
            assert_eq!(mod_mul(2, 3), 6);
            assert_eq!(mod_mul(5, 7), 35);
            assert_eq!(mod_mul(500000, 500000), 3780473);
        }

        #[test]
        fn test_mod_mul_zero() {
            assert_eq!(mod_mul(0, 12345), 0);
            assert_eq!(mod_mul(12345, 0), 0);
            assert_eq!(mod_mul(0, 0), 0);
        }

        #[test]
        fn test_mod_mul_one() {
            assert_eq!(mod_mul(1, 12345), 12345);
            assert_eq!(mod_mul(12345, 1), 12345);
        }

        #[test]
        fn test_mod_mul_negative_inputs() {
            // mod_mul should handle negative inputs correctly
            let result = mod_mul(-5, 3);
            assert_eq!(result, Q - 15); // -15 mod Q

            let result = mod_mul(5, -3);
            assert_eq!(result, Q - 15); // -15 mod Q

            let result = mod_mul(-5, -3);
            assert_eq!(result, 15); // 15 mod Q
        }

        #[test]
        fn test_mod_mul_large_values() {
            // Test values that would overflow without proper handling
            let a = Q - 1;
            let b = Q - 1;
            let result = mod_mul(a, b);
            // (Q-1)^2 = Q^2 - 2Q + 1 ≡ 1 (mod Q)
            assert_eq!(result, 1);
        }

        #[test]
        fn test_mod_mul_commutative() {
            let test_pairs =
                vec![(2, 3), (100, 200), (Q - 1, 2), (12345, 67890)];

            for (a, b) in test_pairs {
                assert_eq!(mod_mul(a, b), mod_mul(b, a));
            }
        }
    }

    mod mod_mul_three_tests {
        use super::*;

        #[test]
        fn test_mod_mul_three_basic() {
            assert_eq!(mod_mul_three(2, 3, 4), 24);
            assert_eq!(mod_mul_three(5, 6, 7), 210);
        }

        #[test]
        fn test_mod_mul_three_zero() {
            assert_eq!(mod_mul_three(0, 1, 2), 0);
            assert_eq!(mod_mul_three(1, 0, 2), 0);
            assert_eq!(mod_mul_three(1, 2, 0), 0);
        }

        #[test]
        fn test_mod_mul_three_one() {
            assert_eq!(mod_mul_three(1, 1, 1), 1);
            assert_eq!(mod_mul_three(12345, 1, 1), 12345);
        }

        #[test]
        fn test_mod_mul_three_large_values() {
            // Test with values that require modular reduction
            let a = 1000000;
            let b = 2000000;
            let c = 3;
            let result = mod_mul_three(a, b, c);

            // Verify result is in valid range
            assert!(result >= 0 && result < Q);

            // Verify correctness by comparing with direct calculation
            let expected = ((((a as i64 * b as i64) % Q as i64) * c as i64)
                % Q as i64) as i32;
            assert_eq!(result, expected);
        }

        #[test]
        fn test_mod_mul_three_associative() {
            // (a * b) * c should equal a * (b * c) in modular arithmetic
            let a = 12345;
            let b = 67890;
            let c = 11111;

            let result1 = mod_mul_three(a, b, c);
            let result2 = mod_mul(a, mod_mul(b, c));

            assert_eq!(result1, result2);
        }
    }

    mod lagrange_interpolation_comprehensive_test {
        use super::*;

        #[test]
        fn test_lagrange_interpolation_comprehensive() {
            // Test 1: Constant polynomial - f(x) = 42
            {
                let constant = 42;
                let points = vec![(1, constant), (2, constant), (3, constant)];
                assert_eq!(
                    lagrange_interpolation(&points, 0).unwrap(),
                    constant
                );
                assert_eq!(
                    lagrange_interpolation(&points, 5).unwrap(),
                    constant
                );
                assert_eq!(
                    lagrange_interpolation(&points, 100).unwrap(),
                    constant
                );
            }

            // Test 2: Linear polynomial - y = x + 1
            {
                let points = vec![(1, 2), (3, 4)];
                assert_eq!(lagrange_interpolation(&points, 0).unwrap(), 1);
                assert_eq!(lagrange_interpolation(&points, 2).unwrap(), 3);
                assert_eq!(lagrange_interpolation(&points, 5).unwrap(), 6);
            }

            // Test 3: Quadratic polynomial - f(x) = x^2 - 2x + 1
            {
                let points = vec![(1, 0), (2, 1), (3, 4)];
                assert_eq!(lagrange_interpolation(&points, 0).unwrap(), 1);
                assert_eq!(lagrange_interpolation(&points, 4).unwrap(), 9);
            }

            // Test 4: Single point interpolation
            {
                let points = vec![(5, 100)];
                assert_eq!(lagrange_interpolation(&points, 0).unwrap(), 100);
                assert_eq!(lagrange_interpolation(&points, 5).unwrap(), 100);
                assert_eq!(lagrange_interpolation(&points, 10).unwrap(), 100);
            }

            // Test 5: Modular arithmetic with large values
            {
                let points = vec![(1, Q - 1), (2, Q - 2), (3, Q - 3)];
                let result = lagrange_interpolation(&points, 0).unwrap();
                assert!(result >= 0 && result < Q);
            }

            // Test 6: Negative x-values (modular reduction)
            {
                let points = vec![(-1, 10), (0, 5), (1, 2)];
                let result = lagrange_interpolation(&points, 2).unwrap();
                assert!(result >= 0 && result < Q);
            }

            // Test 7: Interpolation at known points
            {
                let points = vec![(1, 100), (2, 200), (3, 150), (4, 250)];
                for &(x, y) in &points {
                    let result = lagrange_interpolation(&points, x).unwrap();
                    assert_eq!(result, y);
                }
            }

            // Test 8: Shamir secret sharing property - secret at f(0)
            {
                let secret = 12345;
                let points = vec![
                    (1, (secret + 100 + 50) % Q),
                    (2, (secret + 200 + 200) % Q),
                    (3, (secret + 300 + 450) % Q),
                ];
                let recovered = lagrange_interpolation(&points, 0).unwrap();
                assert_eq!(recovered, secret);
            }

            // Test 9: Duplicate x-values (shouldn't panic)
            {
                let points = vec![(1, 100), (1, 200), (2, 300)];
                let _ = lagrange_interpolation(&points, 0);
            }

            // Test 10: Many points stress test
            {
                let points: Vec<(i32, i32)> = (1..=10)
                    .map(|x| {
                        let y = ((x * x + 2 * x + 3) % Q) as i32;
                        (x, y)
                    })
                    .collect();
                assert_eq!(lagrange_interpolation(&points, 0).unwrap(), 3);
                assert_eq!(lagrange_interpolation(&points, 11).unwrap(), 146);
            }

            // Test 11: All zeros (zero polynomial)
            {
                let points = vec![(1, 0), (2, 0), (3, 0)];
                assert_eq!(lagrange_interpolation(&points, 0).unwrap(), 0);
                assert_eq!(lagrange_interpolation(&points, 100).unwrap(), 0);
            }

            // Test 12: Random polynomial consistency
            {
                use rand::{prelude::*, rng};
                let mut rng = rng();
                let coeffs: Vec<i32> =
                    (0..5).map(|_| rng.random_range(0..1000)).collect();

                let points: Vec<(i32, i32)> = (1..=5)
                    .map(|x| {
                        let mut y = 0i64;
                        let mut x_power = 1i64;
                        for &coeff in &coeffs {
                            y = (y + coeff as i64 * x_power) % Q as i64;
                            x_power = (x_power * x as i64) % Q as i64;
                        }
                        (x, y as i32)
                    })
                    .collect();

                let interpolated = lagrange_interpolation(&points, 0).unwrap();
                assert_eq!(interpolated, coeffs[0]);
            }

            // Test 13: Modular inverse edge cases
            {
                let points = vec![(1, 100), (Q - 1, 200)];
                let result = lagrange_interpolation(&points, 0).unwrap();
                assert!(result >= 0 && result < Q);
            }

            // Test 14: Comparison with Shamir reconstruction
            {
                let secret = 42;
                let shamir_poly = vec![secret, 10, 5];
                let points: Vec<(i32, i32)> = (1..=3)
                    .map(|x| {
                        let y = shamir_poly[0]
                            + shamir_poly[1] * x
                            + shamir_poly[2] * x * x;
                        (x, y % Q)
                    })
                    .collect();
                let recovered = lagrange_interpolation(&points, 0).unwrap();
                assert_eq!(recovered, secret);
            }

            // Test 15: Numerical stability with large gaps
            {
                let points = vec![
                    (1, 1000),
                    (100, 2000),
                    (1000, 3000),
                    (10000, 4000),
                    (100000, 5000),
                ];
                let result = lagrange_interpolation(&points, 0).unwrap();
                assert!(result >= 0 && result < Q);
            }
        }
    }

    // Integration tests
    mod integration_tests {
        use super::*;

        #[test]
        fn test_mod_operations_consistency() {
            // Test that various mod operations work together correctly
            let a = 12345;
            let b = 67890;

            // Test inverse property: a * a^(-1) ≡ 1 (mod Q)
            let a_inv = mod_inverse(a).unwrap();
            assert_eq!(mod_mul(a, a_inv), 1);

            // Test power and inverse: a^(Q-2) ≡ a^(-1) (mod Q)
            let a_inv_2 = mod_pow(a, Q - 2);
            assert_eq!(a_inv, a_inv_2);

            // Test multiplication chain
            let c = mod_mul(a, b);
            let c_inv = mod_inverse(c).unwrap();
            let result = mod_mul_three(c, c_inv, 1);
            assert_eq!(result, 1);
        }

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
