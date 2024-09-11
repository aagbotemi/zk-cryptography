use ark_ff::PrimeField;
use polynomial::{Multilinear, MultilinearTrait};

pub fn get_poly_remainder<F: PrimeField>(
    poly: &Multilinear<F>,
    evaluation_point: &F,
) -> Multilinear<F> {
    poly.partial_evaluation(evaluation_point, &0)
}

pub fn get_poly_quotient<F: PrimeField>(poly: &Multilinear<F>) -> Multilinear<F> {
    let f_1 = poly.partial_evaluation(&F::from(1_u8), &0);
    let f_0 = poly.partial_evaluation(&F::from(0_u8), &0);

    let res = f_1 - f_0;

    res
}

pub fn check_for_zero_and_one<F: PrimeField>(bh: &[F], value: &[F]) -> F {
    assert_eq!(
        bh.len(),
        value.len(),
        "The length of bh and value must be the same"
    );

    bh.iter().zip(value).fold(F::one(), |acc, (&b, &e)| {
        if b.is_zero() {
            acc * (F::one() - e)
        } else {
            acc * e
        }
    })
}

pub fn generate_array_of_points<F: PrimeField>(bh_cube: &[Vec<F>], eval_points: &[F]) -> Vec<F> {
    bh_cube
        .iter()
        .map(|bh| check_for_zero_and_one(bh, eval_points))
        .collect()
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Fr;
    use polynomial::{utils::boolean_hypercube, Multilinear};

    use crate::utils::generate_array_of_points;

    use super::{check_for_zero_and_one, get_poly_quotient, get_poly_remainder};

    #[test]
    fn test_check_for_zero_and_one() {
        let bh_1 = vec![Fr::from(0), Fr::from(0), Fr::from(0)];
        let bh_2 = vec![Fr::from(0), Fr::from(0), Fr::from(1)];
        let bh_3 = vec![Fr::from(0), Fr::from(1), Fr::from(0)];
        let bh_4 = vec![Fr::from(0), Fr::from(1), Fr::from(1)];
        let bh_5 = vec![Fr::from(1), Fr::from(0), Fr::from(0)];
        let bh_6 = vec![Fr::from(1), Fr::from(0), Fr::from(1)];
        let bh_7 = vec![Fr::from(1), Fr::from(1), Fr::from(0)];
        let bh_8 = vec![Fr::from(1), Fr::from(1), Fr::from(1)];
        let value = vec![Fr::from(2), Fr::from(3), Fr::from(4)];

        let checker_1 = check_for_zero_and_one(&bh_1, &value);
        let checker_2 = check_for_zero_and_one(&bh_2, &value);
        let checker_3 = check_for_zero_and_one(&bh_3, &value);
        let checker_4 = check_for_zero_and_one(&bh_4, &value);
        let checker_5 = check_for_zero_and_one(&bh_5, &value);
        let checker_6 = check_for_zero_and_one(&bh_6, &value);
        let checker_7 = check_for_zero_and_one(&bh_7, &value);
        let checker_8 = check_for_zero_and_one(&bh_8, &value);

        assert_eq!(checker_1, Fr::from(-6));
        assert_eq!(checker_2, Fr::from(8));
        assert_eq!(checker_3, Fr::from(9));
        assert_eq!(checker_4, Fr::from(-12));
        assert_eq!(checker_5, Fr::from(12));
        assert_eq!(checker_6, Fr::from(-16));
        assert_eq!(checker_7, Fr::from(-18));
        assert_eq!(checker_8, Fr::from(24));
    }

    #[test]
    fn test_poly_subtraction() {
        let poly1 = Multilinear::new(vec![
            Fr::from(0),
            Fr::from(0),
            Fr::from(0),
            Fr::from(5),
            Fr::from(4),
            Fr::from(4),
            Fr::from(7),
            Fr::from(12),
        ]);
        let poly2 = Multilinear::new(vec![
            Fr::from(0),
            Fr::from(0),
            Fr::from(0),
            Fr::from(2),
            Fr::from(0),
            Fr::from(0),
            Fr::from(1),
            Fr::from(3),
        ]);
        let poly3 = Multilinear::new(vec![
            Fr::from(0),
            Fr::from(0),
            Fr::from(0),
            Fr::from(3),
            Fr::from(4),
            Fr::from(4),
            Fr::from(6),
            Fr::from(9),
        ]);

        let res1 = poly1 - poly2;
        assert_eq!(res1, poly3);
    }

    #[test]
    fn test_get_poly_quotient() {
        let evals = vec![
            Fr::from(0),
            Fr::from(7),
            Fr::from(0),
            Fr::from(5),
            Fr::from(0),
            Fr::from(7),
            Fr::from(4),
            Fr::from(9),
        ];
        let poly1 = Multilinear::new(evals);
        let poly2 = Multilinear::new(vec![Fr::from(0), Fr::from(7), Fr::from(20), Fr::from(25)]);
        let poly3 = Multilinear::new(vec![Fr::from(180), Fr::from(169)]);
        let quotient_x = get_poly_quotient(&poly1);
        let quotient_y = get_poly_quotient(&poly2);
        let quotient_z = get_poly_quotient(&poly3);
        let expected_poly_after_quotient_x =
            Multilinear::new(vec![Fr::from(0), Fr::from(0), Fr::from(4), Fr::from(4)]);
        let expected_poly_after_quotient_y = Multilinear::new(vec![Fr::from(20), Fr::from(18)]);

        assert_eq!(quotient_x, expected_poly_after_quotient_x);
        assert_eq!(quotient_y, expected_poly_after_quotient_y);
        assert_eq!(quotient_z.evaluations[0], Fr::from(-11));
    }

    #[test]
    fn test_get_poly_remainder() {
        let evals = vec![
            Fr::from(0),
            Fr::from(7),
            Fr::from(0),
            Fr::from(5),
            Fr::from(0),
            Fr::from(7),
            Fr::from(4),
            Fr::from(9),
        ];

        let expected_evals_yz = vec![Fr::from(0), Fr::from(7), Fr::from(20), Fr::from(25)];
        let expected_evals_z = vec![Fr::from(180), Fr::from(169)];

        let poly = Multilinear::new(evals);
        let expected_poly_yz = Multilinear::new(expected_evals_yz);
        let expected_poly_z = Multilinear::new(expected_evals_z);

        let remainder_after_x = get_poly_remainder(&poly, &Fr::from(5));
        let remainder_after_y = get_poly_remainder(&expected_poly_yz, &Fr::from(9));
        let remainder_after_z = get_poly_remainder(&expected_poly_z, &Fr::from(6));

        assert_eq!(expected_poly_yz, remainder_after_x);
        assert_eq!(expected_poly_z, remainder_after_y);
        assert_eq!(Fr::from(114), remainder_after_z.evaluations[0]);
    }

    #[test]
    fn test_generate_array_of_points() {
        let expected_evals = vec![
            Fr::from(-6),
            Fr::from(8),
            Fr::from(9),
            Fr::from(-12),
            Fr::from(12),
            Fr::from(-16),
            Fr::from(-18),
            Fr::from(24),
        ];
        let expected_poly = Multilinear::new(expected_evals);

        let bh_cube: Vec<Vec<Fr>> = boolean_hypercube(3);
        let eval_points = vec![Fr::from(2), Fr::from(3_u8), Fr::from(4_u8)];

        let array_of_points = generate_array_of_points(&bh_cube, &eval_points);
        let result_poly = Multilinear::new(array_of_points);

        assert_eq!(expected_poly, result_poly);
    }
}
