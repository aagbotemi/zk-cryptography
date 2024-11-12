use ark_test_curves::bls12_381::{Bls12_381, Fr};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use kzg::{
    interface::UnivariateKZGInterface,
    univariate_kzg::{UnivariateKZG, UnivariateKZGProof},
};
use polynomial::{DenseUnivariatePolynomial, UnivariatePolynomialTrait};

fn univariate_kzg_benchmark(c: &mut Criterion) {
    let poly = black_box(DenseUnivariatePolynomial::new(
        (0u8..=255)
            .into_iter()
            .map(|x| Fr::from(x))
            .collect::<Vec<Fr>>(),
    ));

    let srs = black_box(UnivariateKZG::generate_srs(&Fr::from(10), &4));

    c.bench_function("univariate_kzg_benchmark", |b| {
        b.iter(|| {
            let commit = UnivariateKZG::commitment(&poly, &srs);

            let proof: UnivariateKZGProof<Fr, Bls12_381> =
                UnivariateKZG::open(&poly, Fr::from(2), &srs);
            let verify_status = UnivariateKZG::verify(&commit, &Fr::from(2), &proof, &srs);

            assert_eq!(verify_status, true)
        });
    });
}

criterion_group!(benches, univariate_kzg_benchmark);
criterion_main!(benches);
