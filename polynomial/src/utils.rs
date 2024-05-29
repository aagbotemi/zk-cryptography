use crate::multilinear::coefficient_form::MultiLinearMonomial;
use ark_ff::PrimeField;

pub fn pick_pairs_with_index<F: PrimeField>(
    terms: &Vec<MultiLinearMonomial<F>>,
) -> Vec<(usize, usize)> {
    let n = terms.len();
    let mut pairs = Vec::with_capacity(n / 2);

    for i in 0..(n / 2) {
        let j = i + n / 2;
        pairs.push((i, j));
    }

    pairs
}

pub fn pick_pairs_with_random_index(
    num_of_evaluations: usize,
    variable_index: usize,
) -> Vec<(usize, usize)> {
    assert!(num_of_evaluations % 2 == 0, "n must be even");
    assert!(
        variable_index < num_of_evaluations / 2,
        "variable_index must be less than n/2"
    );

    let mut result = Vec::new();
    let iters = 1 << variable_index;

    for _ in 0..iters {
        let mut round: Vec<(usize, usize)> = Vec::new();

        for y_1 in 0..((num_of_evaluations / iters) / 2) {
            round.push((
                y_1 + result.len() * 2,
                ((num_of_evaluations / iters) / 2) + y_1 + result.len() * 2,
            ));
        }

        result.extend(round);
    }

    result
}
