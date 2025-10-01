//! Small helpers to keep code DRY and straightforward.

use rand::RngCore;

use sha2::digest::core_api::XofReaderCoreWrapper;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256, Shake256ReaderCore,
};


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
