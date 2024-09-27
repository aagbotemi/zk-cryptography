use ark_test_curves::bls12_381::Fr;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use polynomial::{ComposedMultilinear, Multilinear};
use sumcheck::composed::multi_composed_sumcheck::{
    MultiComposedSumcheckProver, MultiComposedSumcheckVerifier,
};

fn multi_composed_sumcheck_benchmark(c: &mut Criterion) {
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
    let composed_poly_1 = black_box(ComposedMultilinear::new(vec![poly1, poly2]));
    let composed_poly_2 = black_box(ComposedMultilinear::new(vec![poly3, poly4, poly5]));
    let multi_composed_poly = black_box(vec![composed_poly_1, composed_poly_2]);

    c.bench_function("multi_composed_sumcheck_benchmark", |b| {
        b.iter(|| {
            let sum = MultiComposedSumcheckProver::calculate_poly_sum(&multi_composed_poly);
            let (proof, _) =
                MultiComposedSumcheckProver::prove(&multi_composed_poly, &sum).unwrap();
            let verify =
                MultiComposedSumcheckVerifier::verify(&multi_composed_poly, &proof).unwrap();

            assert!(verify);
        });
    });
}

fn multi_composed_sumcheck_without_verification_benchmark(c: &mut Criterion) {
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
    let composed_poly_1 = black_box(ComposedMultilinear::new(vec![poly1, poly2]));
    let composed_poly_2 = black_box(ComposedMultilinear::new(vec![poly3, poly4, poly5]));
    let multi_composed_poly = black_box(vec![composed_poly_1, composed_poly_2]);

    c.bench_function(
        "multi_composed_sumcheck_without_verification_benchmark",
        |b| {
            b.iter(|| {
                let sum = MultiComposedSumcheckProver::calculate_poly_sum(&multi_composed_poly);
                let (_, _) =
                    MultiComposedSumcheckProver::prove(&multi_composed_poly, &sum).unwrap();
            });
        },
    );
}

fn multi_composed_sumcheck_with_prove_partial_benchmark(c: &mut Criterion) {
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
    let composed_poly_1 = black_box(ComposedMultilinear::new(vec![poly1, poly2]));
    let composed_poly_2 = black_box(ComposedMultilinear::new(vec![poly3, poly4, poly5]));
    let multi_composed_poly = black_box(vec![composed_poly_1, composed_poly_2]);

    c.bench_function(
        "multi_composed_sumcheck_with_prove_partial_benchmark",
        |b| {
            b.iter(|| {
                let sum = MultiComposedSumcheckProver::calculate_poly_sum(&multi_composed_poly);
                let (_, _) =
                    MultiComposedSumcheckProver::prove_partial(&multi_composed_poly, &sum).unwrap();
            });
        },
    );
}

criterion_group!(
    benches,
    multi_composed_sumcheck_benchmark,
    multi_composed_sumcheck_without_verification_benchmark,
    multi_composed_sumcheck_with_prove_partial_benchmark
);
criterion_main!(benches);
