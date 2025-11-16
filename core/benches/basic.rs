use criterion::{Criterion, black_box, criterion_group, criterion_main};
use dilithium_core::basic::{KeyPair, keygen};
use math::field_element::FieldElement;

fn fresh_keypair() -> KeyPair<'static, FieldElement> {
    keygen::<FieldElement>().expect("random key generation succeeds")
}

fn bench_sign(c: &mut Criterion) {
    let keypair = fresh_keypair();
    let msg = b"Hello!";

    c.bench_function("sign", move |b| {
        b.iter(|| {
            let signature =
                keypair.sign(black_box(msg)).expect("signing succeeds");
            black_box(signature);
        });
    });
}

fn bench_verify(c: &mut Criterion) {
    let keypair = fresh_keypair();
    let msg = b"Hello";
    let signature = keypair.sign(msg).expect("signing succeeds");
    let pub_key = keypair.public.clone();

    c.bench_function("verify", move |b| {
        b.iter(|| {
            assert!(pub_key.verify(black_box(msg), &signature));
        });
    });
}

criterion_group!(benches, bench_sign, bench_verify);
criterion_main!(benches);
