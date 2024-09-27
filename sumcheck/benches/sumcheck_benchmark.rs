use ark_test_curves::bls12_381::Fr;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use polynomial::Multilinear;
use sumcheck::sumcheck::Sumcheck;

fn sumcheck_benchmark(c: &mut Criterion) {
    let poly = black_box(Multilinear::new(
        (0u8..=255)
            .into_iter()
            .map(|x| Fr::from(x))
            .collect::<Vec<Fr>>(),
    ));
    c.bench_function("sumcheck_benchmark", |b| {
        b.iter(|| {
            let mut sumcheck = Sumcheck::new(poly.clone());
            sumcheck.poly_sum();
            let proof = sumcheck.prove();
            let verifer = sumcheck.verify(&proof.0);

            assert_eq!(verifer, true);
        });
    });
}

fn sumcheck_benchmark_without_verification(c: &mut Criterion) {
    let poly = black_box(Multilinear::new(
        (0u8..=255)
            .into_iter()
            .map(|x| Fr::from(x))
            .collect::<Vec<Fr>>(),
    ));
    c.bench_function("sumcheck_benchmark_without_verification", |b| {
        b.iter(|| {
            let mut sumcheck = Sumcheck::new(poly.clone());
            sumcheck.poly_sum();
            let proof = sumcheck.prove();
            let verifer = sumcheck.verify(&proof.0);

            assert_eq!(verifer, true);
        });
    });
}

criterion_group!(
    benches,
    sumcheck_benchmark,
    sumcheck_benchmark_without_verification
);
criterion_main!(benches);
