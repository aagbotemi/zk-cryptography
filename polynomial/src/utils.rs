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

pub fn pick_pairs_with_random_n_index(
    num_of_evaluations: usize,
    variable_indices: &[usize],
) -> Vec<(usize, usize)> {
    assert!(num_of_evaluations % 2 == 0, "n must be even");
    assert!(
        variable_indices.len() < num_of_evaluations / 2,
        "variable_index must be less than n/2"
    );

    let mut pairs = Vec::new();

    for i in 0..num_of_evaluations {
        for &variable_index in variable_indices {
            let pair = (i, i ^ (1 << variable_index));
            if !pairs.contains(&pair) && !pairs.contains(&(pair.1, pair.0)) {
                pairs.push(pair);
            }
        }
    }
    pairs
}

pub fn lagrange_basis<F: PrimeField>(points: &[(F, F)], i: usize) -> Vec<F> {
    let mut l_i = vec![F::one()];

    for (j, &(x_j, _)) in points.iter().enumerate() {
        if i != j {
            let mut new_l_i = vec![F::zero(); l_i.len() + 1];
            for (k, &coeff) in l_i.iter().enumerate() {
                new_l_i[k] -= coeff * x_j;
                new_l_i[k + 1] += coeff;
            }
            l_i = new_l_i;
        }
    }

    let denom = points
        .iter()
        .enumerate()
        .filter(|&(j, _)| j != i)
        .fold(F::one(), |acc, (_, &(x_j, _))| acc * (points[i].0 - x_j));
    l_i.into_iter()
        .map(|coeff| coeff * denom.inverse().unwrap())
        .collect()
}

pub fn boolean_hypercube<F: PrimeField>(n: usize) -> Vec<Vec<F>> {
    let mut hypercube = Vec::with_capacity(1 << n);

    for i in 0..(1 << n) {
        let mut vertex = Vec::with_capacity(n);
        for j in (0..n).rev() {
            if (i & (1 << j)) != 0 {
                vertex.push(F::one());
            } else {
                vertex.push(F::zero());
            }
        }
        hypercube.push(vertex);
    }

    hypercube
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Fq;

    #[test]
    fn test_boolean_hypercube() {
        let one: Vec<Vec<Fq>> = boolean_hypercube(1);
        let two: Vec<Vec<Fq>> = boolean_hypercube(2);
        let three: Vec<Vec<Fq>> = boolean_hypercube(3);

        let expected_one = vec![vec![Fq::from(0)], vec![Fq::from(1)]];
        let expected_two = vec![
            vec![Fq::from(0), Fq::from(0)],
            vec![Fq::from(0), Fq::from(1)],
            vec![Fq::from(1), Fq::from(0)],
            vec![Fq::from(1), Fq::from(1)],
        ];
        let expected_three = vec![
            vec![Fq::from(0), Fq::from(0), Fq::from(0)],
            vec![Fq::from(0), Fq::from(0), Fq::from(1)],
            vec![Fq::from(0), Fq::from(1), Fq::from(0)],
            vec![Fq::from(0), Fq::from(1), Fq::from(1)],
            vec![Fq::from(1), Fq::from(0), Fq::from(0)],
            vec![Fq::from(1), Fq::from(0), Fq::from(1)],
            vec![Fq::from(1), Fq::from(1), Fq::from(0)],
            vec![Fq::from(1), Fq::from(1), Fq::from(1)],
        ];

        assert_eq!(one, expected_one);
        assert_eq!(two, expected_two);
        assert_eq!(three, expected_three);
    }
}
