//! Shared constants for the Dilithium math primitives.

use crate::field_element::FieldElement;

/// Dilithium prime modulus (alias to the single source of truth).
pub const DILITHIUM_Q: i32 = FieldElement::P as i32;

/// Polynomial degree used across Dilithium (a power of two for NTTs).
pub const DILITHIUM_N: usize = 256;
