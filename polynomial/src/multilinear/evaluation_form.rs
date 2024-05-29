use crate::{interface::MLETrait, utils::pick_pairs_with_random_index};
use ark_ff::PrimeField;

#[derive(Debug, Clone, PartialEq)]
struct MLE<F: PrimeField> {
    n_vars: usize,
    evaluations: Vec<F>,
}

impl<F: PrimeField> MLETrait<F> for MLE<F> {
    fn new(evaluations: Vec<F>) -> Self {
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

    fn partial_evaluation(&self, eval_point: F, variable_index: usize) -> Self {
        let new_evaluation: &Vec<F> = &self.evaluations;

        let mut result: Vec<F> = Vec::with_capacity(self.evaluations.len() / 2);

        for (i, j) in pick_pairs_with_random_index(self.evaluations.len(), variable_index) {
            let y1: &F = &new_evaluation[i];
            let y2: &F = &new_evaluation[j];

            // r.y1 + (1-r).y2 straight line formula
            let res_y: F = (eval_point * y2) + ((F::one() - eval_point) * y1);

            result.push(res_y);
        }

        Self {
            n_vars: self.n_vars - 1,
            evaluations: result,
        }
    }

    /// full evaluation of a polynomial - coefficient form
    fn evaluation(&self, evaluation_points: &[F]) -> F {
        assert_eq!(
            evaluation_points.len(),
            self.n_vars,
            "Number of evaluation points must match the number of variables"
        );

        let mut poly_result: F = F::one();
        let mut eval_result: MLE<F> = self.relabel();

        for i in 0..evaluation_points.len() {
            eval_result = eval_result.partial_evaluation(evaluation_points[i], 0);
        }

        poly_result = eval_result.evaluations[0];

        poly_result
    }

    fn relabel(&self) -> Self {
        self.clone()
    }
}

#[cfg(test)]
mod tests {
    use crate::interface::MLETrait;
    use crate::multilinear::evaluation_form::MLE;
    use ark_ff::MontConfig;
    use ark_ff::{Fp64, MontBackend};

    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    struct FqConfig;
    type Fq = Fp64<MontBackend<FqConfig, 1>>;

    #[test]
    fn test_partial_evaluation_1() {
        let evaluations = vec![Fq::from(3), Fq::from(1), Fq::from(2), Fq::from(5)];
        let polynomial = MLE::new(evaluations);

        let evaluation_point = Fq::from(5);
        let new_polynomial = polynomial.partial_evaluation(evaluation_point, 0);

        let expected_polynomial = MLE::new(vec![Fq::from(15), Fq::from(4)]);

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
        let polynomial = MLE::new(evaluations);

        // obtain: f(2,y,z) = 4yz + 4y + 6z + 9 at y = 3, z = 2 = 57
        let new_polynomial_x_1 = polynomial.partial_evaluation(Fq::from(2), 0);
        // 4yz + 4y + 6z + 9
        let x_1_eval_result = new_polynomial_x_1.evaluation(&vec![Fq::from(3), Fq::from(2)]);
        assert_eq!(x_1_eval_result, Fq::from(57));

        // obtain: f(x,3,z) = 6xz + 3x + 6z + 15 at y = 3, z = 2 = 72
        let new_polynomial_y_1 = polynomial.partial_evaluation(Fq::from(3), 1);
        // 6xz + 3x + 6z + 15
        let y_1_eval_result = new_polynomial_y_1.evaluation(&vec![Fq::from(3), Fq::from(2)]);
        assert_eq!(y_1_eval_result, Fq::from(72));

        // obtain: f(x,y,1) = 2xy + 3x + 4y + 9  at y = 3, z = 2 = 38
        let new_polynomial_z_1 = polynomial.partial_evaluation(Fq::from(1), 2);
        // 2xy + 3x + 4y + 9
        let z_1_eval_result = new_polynomial_z_1.evaluation(&vec![Fq::from(3), Fq::from(2)]);
        assert_eq!(z_1_eval_result, Fq::from(38));
    }

    #[test]
    fn test_evaluation_1() {
        let evaluations = vec![Fq::from(3), Fq::from(1), Fq::from(2), Fq::from(5)];
        let polynomial = MLE::new(evaluations);

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
        let polynomial_2 = MLE::new(evaluations_2);
        let points_2 = vec![Fq::from(2), Fq::from(3), Fq::from(1)];
        let result_polynomial_2 = polynomial_2.evaluation(&points_2);
        assert_eq!(result_polynomial_2, Fq::from(5));
    }

    #[test]
    fn test_evaluation_2() {
        // f(a, b, c) = 2ab + 3bc
        let poly = MLE::new(vec![
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
}
