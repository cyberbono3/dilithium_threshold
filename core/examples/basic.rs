use dilithium_core::basic::sign::DilithiumSignature;
use dilithium_core::basic::{KeyPair, keygen};
use math::field_element::FieldElement;

fn sign_and_verify(
    keypair: &KeyPair<'static, FieldElement>,
    message: &[u8],
) -> DilithiumSignature<'static, FieldElement> {
    let signature = keypair.sign(message).expect("signing should succeed");
    assert!(
        keypair.verify(message, &signature),
        "signature should verify with the keypair"
    );
    signature
}

fn main() {
    let message = b"example Dilithium signing message";
    let keypair =
        keygen::<FieldElement>().expect("key generation should succeed");

    let signature = sign_and_verify(&keypair, message);

    println!(
        "Public key derived from rho seed: {:02X?}",
        keypair.public.rho
    );
    println!(
        "Signature challenge polynomial has {} coefficients",
        signature.c.coefficients().len()
    );
}
