use criterion::{criterion_group, criterion_main, Criterion};
use num_bigint::BigUint;
use signature_schemes::rsa::RSA;

fn rsa_benchmark(c: &mut Criterion) {
    let rsa = RSA::new(
        BigUint::from(1223u32),
        BigUint::from(1987u32),
        BigUint::from(948047u32),
    );

    let message = BigUint::from(5u32);
    c.bench_function("rsa_benchmark", |b| {
        b.iter(|| {
            let cipher_text = rsa.encryption(&message).unwrap();
            assert_eq!(cipher_text, BigUint::from(915542u32));

            let decrypted = rsa.decryption(&cipher_text).unwrap();
            assert_eq!(decrypted, BigUint::from(5u32));
        });
    });
}

criterion_group!(benches, rsa_benchmark);
criterion_main!(benches);
