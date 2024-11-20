use ark_ff::{BigInteger, PrimeField};
use polynomial::{
    interface::MultilinearTrait, utils::boolean_hypercube, ComposedMultilinear,
    ComposedMultilinearTrait, Multilinear,
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
    use ark_test_curves::bls12_381::Fr as Fr_old;
    use field_tracker::Ft;
    use polynomial::ComposedMultilinear;

    type Fr = Ft<4, Fr_old>;

    #[test]
    fn test_convert_field_to_byte() {
        let one: Vec<u8> = convert_field_to_byte(&Fr::from(1));
        let hundred: Vec<u8> = convert_field_to_byte(&Fr::from(100));
        let ninety: Vec<u8> = convert_field_to_byte(&Fr::from(90));

        let expected_one = vec![
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 1,
        ];
        let expected_hundred = vec![
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 100,
        ];
        let incorrect_ninety = vec![
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 10,
        ];

        assert_eq!(one, expected_one);
        assert_eq!(hundred, expected_hundred);
        assert_ne!(ninety, incorrect_ninety);
        // println!("{}", Fr::summary());
    }

    #[test]
    fn test_skip_first_and_sum_all() {
        let poly1 = Multilinear::new(vec![
            Fr::from(0),
            Fr::from(0),
            Fr::from(0),
            Fr::from(2),
            Fr::from(2),
            Fr::from(2),
            Fr::from(2),
            Fr::from(4),
        ]);
        let poly2 = Multilinear::new(vec![
            Fr::from(0),
            Fr::from(0),
            Fr::from(2),
            Fr::from(7),
            Fr::from(3),
            Fr::from(3),
            Fr::from(6),
            Fr::from(11),
        ]);
        let evaluation1 = skip_first_and_sum_all(poly1);
        let evaluation2 = skip_first_and_sum_all(poly2);

        let expected_polynomial1 = Multilinear::new(vec![Fr::from(2), Fr::from(10)]);
        let expected_polynomial2 = Multilinear::new(vec![Fr::from(9), Fr::from(23)]);

        assert_eq!(evaluation1, expected_polynomial1);
        assert_eq!(evaluation2, expected_polynomial2);
        // println!("{}", Fr::summary());
    }

    #[test]
    fn test_convert_round_poly_to_uni_poly_format() {
        let point = vec![Fr::from(1), Fr::from(1), Fr::from(1), Fr::from(1)];
        let uni_poly_points = convert_round_poly_to_uni_poly_format(&point);
        assert_eq!(
            uni_poly_points,
            [
                (Fr::from(0), Fr::from(1)),
                (Fr::from(1), Fr::from(1)),
                (Fr::from(2), Fr::from(1)),
                (Fr::from(3), Fr::from(1))
            ]
        );
        // println!("{}", Fr::summary());
    }

    #[test]
    fn test_sum_over_the_boolean_hypercube() {
        let val = vec![
            Fr::from(1),
            Fr::from(2),
            Fr::from(3),
            Fr::from(4),
            Fr::from(5),
            Fr::from(6),
            Fr::from(7),
            Fr::from(8),
        ];

        let poly = ComposedMultilinear::new([Multilinear::new(val)].to_vec());

        let res = sum_over_boolean_hypercube(&[poly]);
        assert!(
            res == Fr::from(36),
            "Incorrect sum over the boolean hypercube"
        );
        // println!("{}", Fr::summary());
    }
}
