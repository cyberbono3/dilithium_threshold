use crate::keypair::keygen;
use crate::sign::{sign, verify};

use math::field_element::FieldElement;

#[test]
fn round_trip_sign_verify() {
    let (pub_key, priv_key) = keygen::<FieldElement>();

    let msg = b"The lattice rocks: Dilithium / ML-DSA demo!";
    let sig = sign::<FieldElement>(&priv_key, &pub_key.t, msg);

    assert!(verify::<FieldElement>(&pub_key, msg, &sig));
}

#[test]
fn negative_case_wrong_message() {
    let (pub_key, priv_key) = keygen::<FieldElement>();
    let msg = b"hello world";
    let sig = sign::<FieldElement>(&priv_key, &pub_key.t, msg);

    let wrong = b"hello wurld";
    assert!(!verify::<FieldElement>(&pub_key, wrong, &sig));
}

#[test]
fn empty_and_long_messages_round_trip() {
    let (pub_key, priv_key) = keygen::<FieldElement>();
    let empty = b"";
    let sig_empty = sign::<FieldElement>(&priv_key, &pub_key.t, empty);
    assert!(verify::<FieldElement>(&pub_key, empty, &sig_empty));

    let long = vec![0xABu8; 8192];
    let sig_long = sign::<FieldElement>(&priv_key, &pub_key.t, &long);
    assert!(verify::<FieldElement>(&pub_key, &long, &sig_long));
}
