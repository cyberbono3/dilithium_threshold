use core::basic::keypair::{self, KeypairSeeds};
use core::basic::{sign, verify};
use criterion::{Criterion, black_box, criterion_group, criterion_main};
use math::field_element::FieldElement;

fn deterministic_keypair() -> (
    core::basic::PublicKey<'static, FieldElement>,
    core::basic::PrivateKey<'static, FieldElement>,
) {
    let seeds = KeypairSeeds::new([0x42; 32], [0x24; 32], [0x18; 32]);
    keypair::keygen_with_seeds::<FieldElement>(seeds)
        .expect("deterministic key generation succeeds")
}

fn bench_sign_only(c: &mut Criterion) {
    let (pub_key, priv_key) = deterministic_keypair();
    let msg = b"hello, sign-only!";

    c.bench_function("sign_only", |b| {
        b.iter(|| {
            let signature =
                sign::<FieldElement>(&priv_key, &pub_key.t, black_box(msg))
                    .expect("signing succeeds");
            black_box(signature);
        });
    });
}

fn bench_sign_and_verify_roundtrip(c: &mut Criterion) {
    let (pub_key, priv_key) = deterministic_keypair();
    let msg = b"hello, sign+verify!";

    c.bench_function("sign_and_verify_roundtrip", |b| {
        b.iter(|| {
            let signature =
                sign::<FieldElement>(&priv_key, &pub_key.t, black_box(msg))
                    .expect("signing succeeds");
            assert!(verify::<FieldElement>(
                &pub_key,
                black_box(msg),
                &signature
            ));
        });
    });
}

criterion_group!(benches, bench_sign_only, bench_sign_and_verify_roundtrip);
criterion_main!(benches);
