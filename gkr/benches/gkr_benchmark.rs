use ark_test_curves::bls12_381::{Bls12_381, Fr};
use circuit::circuit::Circuit;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use gkr::{
    protocol::{GKRProof, GKRProtocol},
    succint_protocol::SuccintGKRProtocol,
    utils::w_mle,
};
use multilinear_kzg::{interface::TrustedSetupInterface, trusted_setup::TrustedSetup};

fn gkr_benchmark(c: &mut Criterion) {
    let circuit = Circuit::random(8);
    let input = black_box(
        (0u8..=255)
            .into_iter()
            .map(|x| Fr::from(x))
            .collect::<Vec<Fr>>(),
    );
    c.bench_function("gkr_benchmark", |b| {
        b.iter(|| {
            let circuit_evaluation = circuit.evaluation(&input);
            let proof: GKRProof<_> = GKRProtocol::prove(&circuit, &circuit_evaluation);
            let verify = GKRProtocol::verify(&circuit, &input, &proof);
            assert!(verify);
        });
    });
}

fn succint_gkr_benchmark(c: &mut Criterion) {
    let circuit = Circuit::random(8);
    let input = black_box(
        (0u8..=255)
            .into_iter()
            .map(|x| Fr::from(x))
            .collect::<Vec<Fr>>(),
    );
    let poly = w_mle(input.clone());

    let points = black_box(
        (0u8..(poly.n_vars as u8))
            .into_iter()
            .map(|x| Fr::from(x))
            .collect::<Vec<Fr>>(),
    );

    c.bench_function("succint_gkr_benchmark", |b| {
        b.iter(|| {
            let circuit_evaluation = circuit.evaluation(&input);

            let tau = TrustedSetup::<Bls12_381>::setup(&points);
            let (commitment, proof) =
                SuccintGKRProtocol::prove(&circuit, &circuit_evaluation, &tau);
            let verify = SuccintGKRProtocol::verify(&circuit, &commitment, &proof, &tau);
            assert!(verify);
        });
    });
}

criterion_group!(benches, gkr_benchmark, succint_gkr_benchmark);
criterion_main!(benches);
