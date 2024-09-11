use ark_ff::{BigInteger, PrimeField};
use polynomial::{
    interface::MultilinearTrait, utils::boolean_hypercube, ComposedMultilinear, ComposedMultilinearTrait, Multilinear
};

pub fn convert_field_to_byte<F: PrimeField>(element: &F) -> Vec<u8> {
    element.into_bigint().to_bytes_be()
}

pub fn skip_first_and_sum_all<'a, F: PrimeField>(current_poly: Multilinear<F>) -> Multilinear<F> {
    let rounds = current_poly.n_vars - 1;
    let bh: Vec<Vec<F>> = boolean_hypercube::<F>(rounds);

    let mut bh_sum = Multilinear::<F>::additive_identity(1);

    for bh_i in bh {
        let mut sum_except_first = current_poly.clone();

        for bh_i_i in bh_i {
            sum_except_first = sum_except_first.partial_evaluation(&bh_i_i, &1);
        }
        bh_sum += sum_except_first;
    }

    bh_sum
}

pub fn convert_round_poly_to_uni_poly_format<F: PrimeField>(round_poly: &Vec<F>) -> Vec<(F, F)> {
    round_poly
        .iter()
        .enumerate()
        .map(|(i, val)| (F::from(i as u64), *val))
        .collect()
}

pub fn vec_to_bytes<F: PrimeField>(poly: &Vec<F>) -> Vec<u8> {
    let mut bytes = Vec::new();
    for p in poly {
        bytes.extend_from_slice(&p.into_bigint().to_bytes_be());
    }
    bytes
}

pub fn sum_over_boolean_hypercube<F: PrimeField>(poly: &[ComposedMultilinear<F>]) -> F {
    let evaluation: Vec<Vec<F>> = poly.iter().map(|f| f.element_wise_product()).collect();

    (0..evaluation[0].len())
        .map(|i| evaluation.iter().map(|v| v[i]).sum::<F>())
        .sum()
}

pub fn composed_poly_to_bytes<F: PrimeField>(poly: &[ComposedMultilinear<F>]) -> Vec<u8> {
    let mut bytes = Vec::new();
    for p in poly.iter() {
        bytes.extend_from_slice(&p.to_bytes());
    }
    bytes
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::MontConfig;
    use ark_ff::{Fp64, MontBackend};
    use polynomial::ComposedMultilinear;

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
    fn test_skip_first_and_sum_all() {
        let poly1 = Multilinear::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(2),
            Fq::from(2),
            Fq::from(2),
            Fq::from(4),
        ]);
        let poly2 = Multilinear::new(vec![
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

        let expected_polynomial1 = Multilinear::new(vec![Fq::from(2), Fq::from(10)]);
        let expected_polynomial2 = Multilinear::new(vec![Fq::from(9), Fq::from(6)]);

        assert_eq!(evaluation1, expected_polynomial1);
        assert_eq!(evaluation2, expected_polynomial2);
    }

    #[test]
    fn test_convert_round_poly_to_uni_poly_format() {
        let point = vec![Fq::from(1), Fq::from(1), Fq::from(1), Fq::from(1)];
        let uni_poly_points = convert_round_poly_to_uni_poly_format(&point);
        assert_eq!(
            uni_poly_points,
            [
                (Fq::from(0), Fq::from(1)),
                (Fq::from(1), Fq::from(1)),
                (Fq::from(2), Fq::from(1)),
                (Fq::from(3), Fq::from(1))
            ]
        )
    }

    #[test]
    fn test_sum_over_the_boolean_hypercube() {
        let val = vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(5),
            Fq::from(6),
            Fq::from(7),
            Fq::from(8),
        ];

        let poly = ComposedMultilinear::new([Multilinear::new(val)].to_vec());

        let res = sum_over_boolean_hypercube(&[poly]);
        assert!(
            res == Fq::from(36),
            "Incorrect sum over the boolean hypercube"
        );
    }
}
