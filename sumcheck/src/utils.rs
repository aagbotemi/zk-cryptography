use ark_ff::{BigInteger, PrimeField};
use polynomial::{interface::MLETrait, MLE};

pub fn boolean_hypercube<F: PrimeField>(n: usize) -> Vec<Vec<F>> {
    let mut hypercube = Vec::new();

    for i in 0..(1 << n) {
        let mut vertex = Vec::new();
        for j in 0..n {
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

pub fn convert_field_to_byte<F: PrimeField>(element: &F) -> Vec<u8> {
    element.into_bigint().to_bytes_be()
}

pub fn skip_first_and_sum_all<F: PrimeField>(current_poly: MLE<F>) -> MLE<F> {
    let rounds = current_poly.n_vars - 1;
    let bh: Vec<Vec<F>> = boolean_hypercube::<F>(rounds);

    let mut bh_sum = MLE::<F>::additive_identity(1);

    for bh_i in bh {
        let mut sum_except_first = current_poly.clone();

        for bh_i_i in bh_i {
            sum_except_first = sum_except_first.partial_evaluation(bh_i_i, 1);
        }
        bh_sum += sum_except_first;
    }

    bh_sum
}
