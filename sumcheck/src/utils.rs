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

pub fn skip_first_and_sum_all<'a, F: PrimeField>(current_poly: MLE<F>) -> MLE<F> {
    let rounds = current_poly.n_vars - 1;
    let bh: Vec<Vec<F>> = boolean_hypercube::<F>(rounds);

    let mut bh_sum = MLE::<F>::additive_identity(1);

    for bh_i in bh {
        let mut sum_except_first = current_poly.clone();

        for bh_i_i in bh_i {
            sum_except_first = sum_except_first.partial_evaluation(&bh_i_i, &1);
        }
        bh_sum += sum_except_first;
    }

    bh_sum
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::MontConfig;
    use ark_ff::{Fp64, MontBackend};

    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    struct FqConfig;
    type Fq = Fp64<MontBackend<FqConfig, 1>>;

    #[test]
    fn test_convert_field_to_byte() {
        let one: Vec<u8> = convert_field_to_byte(&Fq::from(1));
        let hundred: Vec<u8> = convert_field_to_byte(&Fq::from(100));
        let ninety: Vec<u8> = convert_field_to_byte(&Fq::from(90));

        let expected_one = vec![0, 0, 0, 0, 0, 0, 0, 1];
        let expected_hundred = vec![0, 0, 0, 0, 0, 0, 0, 15];
        let incorrect_ninety = vec![0, 0, 0, 0, 0, 0, 0, 10];

        assert_eq!(one, expected_one);
        assert_eq!(hundred, expected_hundred);
        assert_ne!(ninety, incorrect_ninety);
    }

    #[test]
    fn test_boolean_hypercube() {
        let one: Vec<Vec<Fq>> = boolean_hypercube(1);
        let two: Vec<Vec<Fq>> = boolean_hypercube(2);
        let three: Vec<Vec<Fq>> = boolean_hypercube(3);

        let expected_one = vec![vec![Fq::from(0)], vec![Fq::from(1)]];
        let expected_two = vec![
            vec![Fq::from(0), Fq::from(0)],
            vec![Fq::from(1), Fq::from(0)],
            vec![Fq::from(0), Fq::from(1)],
            vec![Fq::from(1), Fq::from(1)],
        ];
        let expected_three = vec![
            vec![Fq::from(0), Fq::from(0), Fq::from(0)],
            vec![Fq::from(1), Fq::from(0), Fq::from(0)],
            vec![Fq::from(0), Fq::from(1), Fq::from(0)],
            vec![Fq::from(1), Fq::from(1), Fq::from(0)],
            vec![Fq::from(0), Fq::from(0), Fq::from(1)],
            vec![Fq::from(1), Fq::from(0), Fq::from(1)],
            vec![Fq::from(0), Fq::from(1), Fq::from(1)],
            vec![Fq::from(1), Fq::from(1), Fq::from(1)],
        ];

        assert_eq!(one, expected_one);
        assert_eq!(two, expected_two);
        assert_eq!(three, expected_three);
    }

    #[test]
    fn test_skip_first_and_sum_all() {
        let poly1 = MLE::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(2),
            Fq::from(2),
            Fq::from(2),
            Fq::from(4),
        ]);
        let poly2 = MLE::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(7),
            Fq::from(3),
            Fq::from(3),
            Fq::from(6),
            Fq::from(11),
        ]);
        let evaluation1 = skip_first_and_sum_all(poly1);
        let evaluation2 = skip_first_and_sum_all(poly2);

        let expected_polynomial1 = MLE::new(vec![Fq::from(2), Fq::from(10)]);
        let expected_polynomial2 = MLE::new(vec![Fq::from(9), Fq::from(6)]);

        assert_eq!(evaluation1, expected_polynomial1);
        assert_eq!(evaluation2, expected_polynomial2);
    }
}
