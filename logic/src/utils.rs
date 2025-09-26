//! Small helpers to keep code DRY and straightforward.


use rand::RngCore;


// Fill byte array of length 32 by random bytes
pub fn random_bytes() -> [u8; 32] {
    let mut rng = rand::thread_rng();
    let mut tmp = [0u8; 32];
    rng.fill_bytes(&mut tmp);
    tmp
}
