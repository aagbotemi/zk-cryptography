use crate::{interface::MultilinearTrait, utils::pick_pairs_with_random_index};
use ark_ff::{BigInteger, PrimeField};
use std::ops::{Add, AddAssign, Mul, Sub, SubAssign};

#[derive(Debug, Clone, PartialEq)]
pub struct Multilinear<F: PrimeField> {
    pub n_vars: usize,
    pub evaluations: Vec<F>,
}

impl<F: PrimeField> Multilinear<F> {
    pub fn new(evaluations: Vec<F>) -> Self {
        let num_evaluations = evaluations.len();
        let n_vars = (num_evaluations as f64).log2() as usize;

        assert_eq!(
            1 << n_vars,
            num_evaluations,
            "Number of evaluations must be a power of 2"
        );

        Self {
            n_vars,
            evaluations,
        }
    }

    pub fn add_distinct(&self, rhs: &Self) -> Self {
        let mut new_evaluations = Vec::new();
        let repeat_sequence = rhs.evaluations.len();

        for i in 0..self.evaluations.len() {
            for j in 0..repeat_sequence {
                new_evaluations.push(self.evaluations[i] + rhs.evaluations[j]);
            }
        }

        Self::new(new_evaluations)
    }

    pub fn mul_distinct(&self, rhs: &Self) -> Self {
        let mut new_evaluations = Vec::new();
        let repeat_sequence = rhs.evaluations.len();

        for i in 0..self.evaluations.len() {
            for j in 0..repeat_sequence {
                new_evaluations.push(self.evaluations[i] * rhs.evaluations[j]);
            }
        }

        Self::new(new_evaluations)
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();

        for evaluation in &self.evaluations {
            bytes.extend(evaluation.into_bigint().to_bytes_be());
        }

        bytes
    }

    pub fn additive_identity(num_vars: usize) -> Self {
        Self::new(vec![F::zero(); 1 << num_vars])
    }

    pub fn split_poly_into_two_and_sum_each_part(&mut self) -> Multilinear<F> {
        let mid_point: usize = self.evaluations.len() / 2;
        let first_half: F = self.evaluations[..mid_point].iter().sum();
        let second_half: F = self.evaluations[mid_point..].iter().sum();

        Self::new(vec![first_half, second_half])
    }

    pub fn is_zero(&self) -> bool {
        self.evaluations.iter().all(|&eval| eval.is_zero())
    }

    pub fn sum_over_the_boolean_hypercube(&self) -> F {
        self.evaluations
            .iter()
            .fold(F::zero(), |acc, val| acc + val)
    }

    pub fn add_to_front(&self, variable_length: &usize) -> Self {
        let new_len = 1 << variable_length;
        let mut res = Vec::with_capacity(new_len);

        for _ in 0..new_len {
            res.extend_from_slice(&self.evaluations);
            res.extend_from_slice(&self.evaluations);
        }

        Self::new(res)
    }

    pub fn add_to_back(&self, variable_length: &usize) -> Self {
        let eval = self.evaluations.clone();
        let input_len = eval.len();
        let output_len = input_len * 2usize.pow(*variable_length as u32);
        let repeat_count = output_len / input_len;

        let res = eval
            .into_iter()
            .flat_map(|num| std::iter::repeat(num).take(repeat_count))
            .collect();

        Self::new(res)
    }

    pub fn duplicate_evaluation(value: &[F]) -> Self {
        let mut res = Vec::with_capacity(2);

        res.extend_from_slice(&value);
        res.extend_from_slice(&value);

        Self::new(res)
    }
}

impl<F: PrimeField> MultilinearTrait<F> for Multilinear<F> {
    fn partial_evaluation(&self, eval_point: &F, variable_index: &usize) -> Self {
        let new_evaluation: &Vec<F> = &self.evaluations;

        let mut result: Vec<F> = Vec::with_capacity(self.evaluations.len() / 2);

        for (i, j) in pick_pairs_with_random_index(self.evaluations.len(), *variable_index) {
            let y1: &F = &new_evaluation[i];
            let y2: &F = &new_evaluation[j];

            // r.y1 + (1-r).y2 straight line formula
            let res_y: F = (*eval_point * y2) + ((F::one() - eval_point) * y1);
            result.push(res_y);
        }

        Self {
            n_vars: self.n_vars - 1,
            evaluations: result,
        }
    }

    fn partial_evaluations(&self, points: &[F], variable_indices: &Vec<usize>) -> Self {
        let mut evaluation = self.clone();

        if points.len() != variable_indices.len() {
            panic!(
                "The length of evaluation_points and variable_indices should be the same: {}, {}",
                points.len(),
                variable_indices.len()
            );
        }

        for i in 0..points.len() {
            evaluation = evaluation.partial_evaluation(&points[i], &variable_indices[i]);
        }

        evaluation
    }

    /// full evaluation of a polynomial - evaluation form
    fn evaluation(&self, evaluation_points: &[F]) -> F {
        assert_eq!(
            evaluation_points.len(),
            self.n_vars,
            "Number of evaluation points must match the number of variables"
        );

        let mut eval_result: Multilinear<F> = self.clone();
        for i in 0..evaluation_points.len() {
            eval_result = eval_result.partial_evaluation(&evaluation_points[i], &0);
        }

        eval_result.evaluations[0]
    }
}

impl<F: PrimeField> Add for Multilinear<F> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let lhs = self.evaluations;
        let mut res = vec![];

        for i in 0..lhs.len() {
            res.push(lhs[i] + rhs.evaluations[i])
        }

        Self {
            n_vars: self.n_vars,
            evaluations: res,
        }
    }
}

impl<F: PrimeField> AddAssign for Multilinear<F> {
    fn add_assign(&mut self, other: Self) {
        // TODO: come up with an algo for handling the case where the number of variables in the two polynomials are not the same
        // if self.n_vars != other.n_vars {
        //     panic!("The number of variables in the two polynomials must be the same");
        // }

        for i in 0..self.evaluations.len() {
            self.evaluations[i] += other.evaluations[i];
        }
    }
}

impl<F: PrimeField> Sub for Multilinear<F> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let lhs = self.evaluations;
        let mut res = vec![];

        for i in 0..lhs.len() {
            res.push(lhs[i] - rhs.evaluations[i])
        }

        Self {
            n_vars: self.n_vars,
            evaluations: res,
        }
    }
}

impl<F: PrimeField> SubAssign for Multilinear<F> {
    fn sub_assign(&mut self, other: Self) {
        for i in 0..self.evaluations.len() {
            self.evaluations[i] -= other.evaluations[i];
        }
    }
}

impl<F: PrimeField> Mul<F> for Multilinear<F> {
    type Output = Self;

    fn mul(self, rhs: F) -> Self::Output {
        let lhs = self.evaluations;
        let mut res = vec![];

        for i in 0..lhs.len() {
            res.push(lhs[i] * rhs)
        }

        Self {
            n_vars: self.n_vars,
            evaluations: res,
        }
    }
}

#[cfg(test)]
mod tests {
    use field_tracker::Ft;

    use crate::interface::MultilinearTrait;
    use crate::multilinear::evaluation_form::Multilinear;
    
    use crate::Fq as Fq_old;

    type Fq = Ft<1, Fq_old>;

    #[test]
    fn test_add_mul_distinct() {
        let polynomial = Multilinear::new(vec![Fq::from(0), Fq::from(0), Fq::from(2), Fq::from(2)]);
        let polynomial2 =
            Multilinear::new(vec![Fq::from(0), Fq::from(3), Fq::from(0), Fq::from(3)]);

        let new_add_poly = Multilinear::new(vec![
            Fq::from(0),
            Fq::from(3),
            Fq::from(0),
            Fq::from(3),
            Fq::from(0),
            Fq::from(3),
            Fq::from(0),
            Fq::from(3),
            Fq::from(2),
            Fq::from(5),
            Fq::from(2),
            Fq::from(5),
            Fq::from(2),
            Fq::from(5),
            Fq::from(2),
            Fq::from(5),
        ]);
        let new_mul_poly = Multilinear::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(6),
            Fq::from(0),
            Fq::from(6),
            Fq::from(0),
            Fq::from(6),
            Fq::from(0),
            Fq::from(6),
        ]);

        let add_distinct = polynomial.add_distinct(&polynomial2);
        let mul_distinct = polynomial.mul_distinct(&polynomial2);

        assert_eq!(add_distinct, new_add_poly);
        assert_eq!(mul_distinct, new_mul_poly);
        println!("{}", Fq::summary());
    }

    #[test]
    fn test_partial_evaluation_1() {
        let evaluations = vec![Fq::from(3), Fq::from(1), Fq::from(2), Fq::from(5)];
        let polynomial = Multilinear::new(evaluations);

        let evaluation_point = Fq::from(5);
        let new_polynomial = polynomial.partial_evaluation(&evaluation_point, &0);

        let expected_polynomial = Multilinear::new(vec![Fq::from(15), Fq::from(4)]);

        assert_eq!(new_polynomial, expected_polynomial);
    }

    #[test]
    fn test_partial_evaluation_2() {
        let evaluations = vec![
            Fq::from(3),
            Fq::from(9),
            Fq::from(7),
            Fq::from(13),
            Fq::from(6),
            Fq::from(12),
            Fq::from(10),
            Fq::from(18),
        ];
        let polynomial = Multilinear::new(evaluations);

        // obtain: f(2,y,z) = 4yz + 4y + 6z + 9 at y = 3, z = 2 = 57
        let new_polynomial_x_1 = polynomial.partial_evaluation(&Fq::from(2), &0);
        // 4yz + 4y + 6z + 9
        let x_1_eval_result = new_polynomial_x_1.evaluation(&vec![Fq::from(3), Fq::from(2)]);
        assert_eq!(x_1_eval_result, Fq::from(57));

        // obtain: f(x,3,z) = 6xz + 3x + 6z + 15 at y = 3, z = 2 = 72
        let new_polynomial_y_1 = polynomial.partial_evaluation(&Fq::from(3), &1);
        // 6xz + 3x + 6z + 15
        let y_1_eval_result = new_polynomial_y_1.evaluation(&vec![Fq::from(3), Fq::from(2)]);
        assert_eq!(y_1_eval_result, Fq::from(72));

        // obtain: f(x,y,1) = 2xy + 3x + 4y + 9  at y = 3, z = 2 = 38
        let new_polynomial_z_1 = polynomial.partial_evaluation(&Fq::from(1), &2);
        // 2xy + 3x + 4y + 9
        let z_1_eval_result = new_polynomial_z_1.evaluation(&vec![Fq::from(3), Fq::from(2)]);
        assert_eq!(z_1_eval_result, Fq::from(38));
        println!("{}", Fq::summary());
    }

    #[test]
    fn test_evaluation_1() {
        let evaluations = vec![Fq::from(3), Fq::from(1), Fq::from(2), Fq::from(5)];
        let polynomial = Multilinear::new(evaluations);

        let points = vec![Fq::from(5), Fq::from(6)];
        let result_polynomial = polynomial.evaluation(&points);

        assert_eq!(result_polynomial, Fq::from(0_u8));
        assert_ne!(result_polynomial, Fq::from(3_u8));

        let evaluations_2 = vec![
            Fq::from(3),
            Fq::from(9),
            Fq::from(7),
            Fq::from(13),
            Fq::from(6),
            Fq::from(12),
            Fq::from(10),
            Fq::from(18),
        ];
        let polynomial_2 = Multilinear::new(evaluations_2);
        let points_2 = vec![Fq::from(2), Fq::from(3), Fq::from(1)];
        let result_polynomial_2 = polynomial_2.evaluation(&points_2);
        assert_eq!(result_polynomial_2, Fq::from(5));
        println!("{}", Fq::summary());
    }

    #[test]
    fn test_evaluation_2() {
        // f(a, b, c) = 2ab + 3bc
        let poly = Multilinear::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(5),
        ]);

        let evaluation_result = poly.evaluation(&[Fq::from(2), Fq::from(3), Fq::from(4)]);
        assert_eq!(evaluation_result, Fq::from(48));
    }

    #[test]
    fn test_split_poly_into_two_and_sum_each_part() {
        let mut poly1 = Multilinear::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(2),
            Fq::from(2),
            Fq::from(2),
            Fq::from(4),
        ]);
        let mut poly2 = Multilinear::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(7),
            Fq::from(3),
            Fq::from(3),
            Fq::from(6),
            Fq::from(11),
        ]);
        let evaluation1 = poly1.split_poly_into_two_and_sum_each_part();
        let evaluation2 = poly2.split_poly_into_two_and_sum_each_part();

        let expected_polynomial1 = Multilinear::new(vec![Fq::from(2), Fq::from(10)]);
        let expected_polynomial2 = Multilinear::new(vec![Fq::from(9), Fq::from(6)]);

        assert_eq!(evaluation1, expected_polynomial1);
        assert_eq!(evaluation2, expected_polynomial2);
        println!("{}", Fq::summary());
    }

    #[test]
    fn test_sum_over_boolean_hypercube() {
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

        let poly = Multilinear::new(val);

        let res = poly.sum_over_the_boolean_hypercube();

        assert!(
            res == Fq::from(36),
            "Incorrect sum over the boolean hypercube"
        );
        println!("{}", Fq::summary());
    }

    #[test]
    fn test_add_to_front_and_back() {
        let poly = Multilinear::new(vec![Fq::from(0), Fq::from(0), Fq::from(4), Fq::from(4)]);
        let poly2 = Multilinear::new(vec![Fq::from(0), Fq::from(4)]);
        let expected_add_to_front_poly = Multilinear::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(4),
            Fq::from(4),
            Fq::from(0),
            Fq::from(0),
            Fq::from(4),
            Fq::from(4),
        ]);
        let expected_add_to_front_poly2 = Multilinear::new(vec![
            Fq::from(0),
            Fq::from(4),
            Fq::from(0),
            Fq::from(4),
            Fq::from(0),
            Fq::from(4),
            Fq::from(0),
            Fq::from(4),
        ]);
        let expected_add_to_back_poly = Multilinear::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(4),
            Fq::from(4),
            Fq::from(4),
            Fq::from(4),
        ]);

        let add_to_front = poly.add_to_front(&0);
        let add_to_front2 = poly2.add_to_front(&1);
        let add_to_back = poly.add_to_back(&1);

        assert_eq!(add_to_front, expected_add_to_front_poly);
        assert_eq!(add_to_front2, expected_add_to_front_poly2);
        assert_eq!(add_to_back, expected_add_to_back_poly);
        println!("{}", Fq::summary());
    }

    #[test]
    fn test_poly_subtraction() {
        let poly1 = Multilinear::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(5),
            Fq::from(4),
            Fq::from(4),
            Fq::from(7),
            Fq::from(12),
        ]);
        let poly2 = Multilinear::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(0),
            Fq::from(0),
            Fq::from(1),
            Fq::from(3),
        ]);
        let poly3 = Multilinear::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
            Fq::from(4),
            Fq::from(4),
            Fq::from(6),
            Fq::from(9),
        ]);

        let res1 = poly1 - poly2;
        assert_eq!(res1, poly3);
        println!("{}", Fq::summary());
    }

    #[test]
    fn test_duplicate_evaluation() {
        let value = vec![Fq::from(11)];
        let duplicate = Multilinear::duplicate_evaluation(&value);

        let expected_poly = Multilinear::new(vec![Fq::from(11), Fq::from(11)]);

        assert_eq!(duplicate, expected_poly);
        println!("{}", Fq::summary());
    }
}
