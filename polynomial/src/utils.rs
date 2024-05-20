use crate::multilinear::MultiLinearMonomial;
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
