pub mod error;

mod accumulator;
mod adapted;
mod share;

pub use adapted::AdaptedShamirSSS;
pub use error::{ShamirError, ShamirResult};
pub use share::ShamirShare;
