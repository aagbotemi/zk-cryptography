use ark_test_curves::bls12_381::{Bls12_381, Fr};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use kzg::{
    interface::{MultilinearKZGInterface, TrustedSetupInterface},
    multilinear_kzg::{MultilinearKZG, MultilinearKZGProof},
    trusted_setup::TrustedSetup,
};
use polynomial::Multilinear;

fn multilinear_kzg_benchmark(c: &mut Criterion) {
    let poly = black_box(Multilinear::new(
        (0u8..=255)
            .into_iter()
            .map(|x| Fr::from(x))
            .collect::<Vec<Fr>>(),
    ));
    let prover_points = black_box(
        (0u8..poly.n_vars as u8)
            .into_iter()
            .map(|x| Fr::from(x))
            .collect::<Vec<Fr>>(),
    );
    let verifier_points = black_box(
        (0u8..poly.n_vars as u8)
            .into_iter()
            .map(|x| Fr::from(x))
            .collect::<Vec<Fr>>(),
    );

    let tau = black_box(TrustedSetup::<Bls12_381>::setup(&prover_points));

    c.bench_function("multilinear_kzg_benchmark", |b| {
        b.iter(|| {
            let commit = MultilinearKZG::commitment(&poly, &tau);

            let proof: MultilinearKZGProof<Fr, Bls12_381> =
                MultilinearKZG::open(&poly, &verifier_points, &tau);
            let verify_status = MultilinearKZG::verify(&commit, &verifier_points, &proof, &tau);

            assert_eq!(verify_status, true)
        });
    });
}

criterion_group!(benches, multilinear_kzg_benchmark);
criterion_main!(benches);
