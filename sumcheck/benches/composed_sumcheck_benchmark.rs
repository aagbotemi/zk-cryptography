use ark_test_curves::bls12_381::Fr;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use polynomial::{ComposedMultilinear, Multilinear};
use sumcheck::{composed::composed_sumcheck::ComposedSumcheck, sumcheck::Sumcheck};

fn composed_sumcheck_benchmark(c: &mut Criterion) {
    let poly1 = black_box(Multilinear::new(
        (0u8..=255)
            .into_iter()
            .map(|x| Fr::from(x))
            .collect::<Vec<Fr>>(),
    ));
    let poly2 = black_box(Multilinear::new(
        (0u8..=255)
            .into_iter()
            .map(|x| Fr::from(x))
            .collect::<Vec<Fr>>(),
    ));
    let composed_poly = black_box(ComposedMultilinear::new(vec![poly1, poly2]));

    c.bench_function("composed_sumcheck_benchmark", |b| {
        b.iter(|| {
            let composed_sumcheck = ComposedSumcheck::new(composed_poly.clone());
            let proof = composed_sumcheck.prove();
            let sum = ComposedSumcheck::calculate_poly_sum(&proof.0.poly);
            let verifer = composed_sumcheck.verify(&proof.0, sum);

            assert_eq!(verifer, true);
        });
    });
}

fn composed_sumcheck_benchmark_with_5_polys(c: &mut Criterion) {
    let poly1 = black_box(Multilinear::new(
        (0u8..=255)
            .into_iter()
            .map(|x| Fr::from(x))
            .collect::<Vec<Fr>>(),
    ));
    let poly2 = black_box(Multilinear::new(
        (0u8..=255)
            .into_iter()
            .map(|x| Fr::from(x))
            .collect::<Vec<Fr>>(),
    ));
    let poly3 = black_box(Multilinear::new(
        (0u8..=255)
            .into_iter()
            .map(|x| Fr::from(x))
            .collect::<Vec<Fr>>(),
    ));
    let poly4 = black_box(Multilinear::new(
        (0u8..=255)
            .into_iter()
            .map(|x| Fr::from(x))
            .collect::<Vec<Fr>>(),
    ));
    let poly5 = black_box(Multilinear::new(
        (0u8..=255)
            .into_iter()
            .map(|x| Fr::from(x))
            .collect::<Vec<Fr>>(),
    ));
    let composed_poly = black_box(ComposedMultilinear::new(vec![
        poly1, poly2, poly3, poly4, poly5,
    ]));

    c.bench_function("composed_sumcheck_benchmark_with_5_polys", |b| {
        b.iter(|| {
            let composed_sumcheck = ComposedSumcheck::new(composed_poly.clone());
            let proof = composed_sumcheck.prove();
            let sum = ComposedSumcheck::calculate_poly_sum(&proof.0.poly);
            let verifer = composed_sumcheck.verify(&proof.0, sum);

            assert_eq!(verifer, true);
        });
    });
}

criterion_group!(
    benches,
    composed_sumcheck_benchmark,
    composed_sumcheck_benchmark_with_5_polys
);
criterion_main!(benches);
