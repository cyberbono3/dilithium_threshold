use rand::prelude::*;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

use crate::error::{Result, ThresholdError};
use math::poly::Q;

// TODO add testting

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

pub fn hash_message(message: &[u8]) -> Vec<u8> {
    let mut hasher = Shake256::default();
    hasher.update(message);
    let mut reader = hasher.finalize_xof();
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

#[test]
fn test_lagrange_interpolation() {

    // Test simple interpolation
    let points = vec![(1, 2), (2, 3), (3, 4)];
    let result = lagrange_interpolation(&points, 0).unwrap();

    // For points (1,2), (2,3), (3,4), the polynomial is y = x + 1
    // So at x=0, y should be 1
    assert_eq!(result, 1);
}
