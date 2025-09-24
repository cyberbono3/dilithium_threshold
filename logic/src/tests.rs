use crate::keypair::keygen;
use crate::sign::{sign, verify};

use math::field_element::FieldElement;

#[test]
fn round_trip_sign_verify() {
    let (pk, sk) = keygen::<FieldElement>();

    let msg = b"The lattice rocks: Dilithium / ML-DSA demo!";
    let sig = sign::<FieldElement>(&sk.a, &sk.s1, &sk.s2, &pk.t, msg);

    assert!(verify::<FieldElement>(&pk.a, &pk.t, msg, &sig));
}

#[test]
fn negative_case_wrong_message() {
    let (pk, sk) = keygen::<FieldElement>();
    let msg = b"hello world";
    let sig = sign::<FieldElement>(&sk.a, &sk.s1, &sk.s2, &pk.t, msg);

    let wrong = b"hello wurld";
    assert!(!verify::<FieldElement>(&pk.a, &pk.t, wrong, &sig));
}

#[test]
fn empty_and_long_messages_round_trip() {
    let (pk, sk) = keygen::<FieldElement>();
    let empty = b"";
    let sig_empty = sign::<FieldElement>(&sk.a, &sk.s1, &sk.s2, &pk.t, empty);
    assert!(verify::<FieldElement>(&pk.a, &pk.t, empty, &sig_empty));

    let long = vec![0xABu8; 8192];
    let sig_long = sign::<FieldElement>(&sk.a, &sk.s1, &sk.s2, &pk.t, &long);
    assert!(verify::<FieldElement>(&pk.a, &pk.t, &long, &sig_long));
}
