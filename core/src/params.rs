pub const DEFAULT_SECURITY_LEVEL: usize = 2;
use math::prelude::N;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DilithiumConfig {
    pub k: usize,
    pub l: usize,
    pub eta: u32,
    pub tau: usize,
    pub beta: u32,
    pub gamma1: u32,
    pub gamma2: u32,
    pub d: usize,
}

impl Default for DilithiumConfig {
    fn default() -> Self {
        Self {
            k: 4,
            l: 4,
            eta: 2,
            tau: 39,
            beta: 78,
            gamma1: 131072,
            gamma2: 95232,
            d: 13,
        }
    }
}

impl DilithiumConfig {
    pub fn new(security_level: usize) -> Self {
        match security_level {
            2 => Self::default(),
            3 => Self {
                k: 6,
                l: 5,
                eta: 4,
                tau: 49,
                beta: 196,
                gamma1: 524288,
                gamma2: 261888,
                d: 13,
            },
            5 => Self {
                k: 8,
                l: 7,
                eta: 2,
                tau: 60,
                beta: 120,
                gamma1: 524288,
                gamma2: 261888,
                d: 13,
            },
            _ => panic!("Invalid security level"),
        }
    }
}

/// Validate threshold config
pub fn validate_threshold_config(
    threshold: usize,
    participants: usize,
) -> bool {
    // Threshold must be at least 2 (no single point of failure)
    // Threshold cannot exceed number of participant
    // Reasonable upper limit on participants. Use as a practical limit for the number of participants
    // Participants must be at least 2 (otherwise threshold >= 2 would be impossible)
    (2..=participants).contains(&threshold) && (2..N).contains(&participants)
}
