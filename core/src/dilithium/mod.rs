//! Core Dilithium building blocks: parameters, utilities, and threshold tooling.

pub mod error;
pub mod params;
pub mod shamir;
pub mod threshold;
pub mod utils;

pub use error::{DilithiumError, DilithiumResult};
pub use threshold::{
    ThresholdError, ThresholdInfo, ThresholdKeyShare, ThresholdSignature,
};
