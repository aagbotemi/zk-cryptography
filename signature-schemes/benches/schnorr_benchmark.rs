use criterion::{criterion_group, criterion_main, Criterion};
use signature_schemes::{
    interface::SchnorrSigTrait,
    schnorr::{SchnorrPublicKey, SchnorrSig, SchnorrSignature},
};

fn schnorr_benchmark(c: &mut Criterion) {
    let message = b"Schnorr is a signature scheme";

    c.bench_function("schnorr_benchmark", |b| {
        b.iter(|| {
            let (sk, pk) = SchnorrSig::generate_keypair().unwrap();
            let signature = SchnorrSig::sign(&sk, message).unwrap();
            let result = SchnorrSig::verify(&pk, message, &signature).unwrap();

            assert!(result);
        });
    });
}
fn schnorr_batch_verification_benchmark(c: &mut Criterion) {
    let messages: [&[u8]; 7] = [
            b"Heyo, I am Abiodun Awoyemi",
            b"I'm a Blockchain engineer and a ZK engineer.",
            b"I build smart contract on EVM compatible blockchain",
            b"I also work on Bitcoin, contributing to Open source",
            b"I have a collection of cryptographic and ZK tools, find it here: https://github.com/aagbotemi/zk-cryptography",
            b"My tech stack: Rust, Solidity, JavaScript, TypeScript",
            b"Let's create a beautiful world!", 
        ];

    c.bench_function("schnorr_batch_verification_benchmark", |b| {
        b.iter(|| {
            let mut signatures: Vec<SchnorrSignature> = Vec::new();
            let mut public_keys: Vec<SchnorrPublicKey> = Vec::new();

            for i in 0..messages.len() {
                let (sk, pk) = SchnorrSig::generate_keypair().unwrap();
                signatures.push(SchnorrSig::sign(&sk, &messages[i]).unwrap());
                public_keys.push(pk);
            }

            let result = SchnorrSig::batch_verify(&public_keys, &messages, &signatures).unwrap();
            assert!(result);
        });
    });
}

criterion_group!(
    benches,
    schnorr_benchmark,
    schnorr_batch_verification_benchmark
);
criterion_main!(benches);
