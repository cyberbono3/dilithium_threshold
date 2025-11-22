use dilithium_core::basic::keygen;
use math::field_element::FieldElement;

fn main() {
    let message = b"example Dilithium signing message";
    let keypair =
        keygen::<FieldElement>().expect("key generation should succeed");
    let signature = keypair.sign(message).expect("signing should succeed");
    keypair.verify(message, &signature);

    println!(
        "Public key derived from rho seed: {:02X?}",
        keypair.public.rho
    );
    println!(
        "Signature challenge polynomial has {} coefficients",
        signature.c.coefficients().len()
    );
}
