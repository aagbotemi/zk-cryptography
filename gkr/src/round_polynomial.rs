use ark_ff::PrimeField;
use polynomial::{interface::MLETrait, MLE};

#[derive(Debug, Clone, PartialEq)]
pub struct W<F: PrimeField> {
    add_i: MLE<F>,
    mul_i: MLE<F>,
    w_b: MLE<F>,
    w_c: MLE<F>,
    random_r: Vec<F>, // random sampling
}

impl<F: PrimeField> W<F> {
    /// Create a new `W` polynomial.
    pub fn new(add_i: MLE<F>, mul_i: MLE<F>, w_b: MLE<F>, w_c: MLE<F>, random_r: Vec<F>) -> Self {
        Self {
            add_i,
            mul_i,
            w_b,
            w_c,
            random_r,
        }
    }

    pub fn evaluate_of_w(&self, point: &[F]) -> Option<F> {
        let r_b_c = [&self.random_r, point].concat();

        let add_e = self.add_i.evaluation(&r_b_c);
        let mul_e = self.mul_i.evaluation(&r_b_c);
        let w_b = self.w_b.evaluation(&point[0..self.w_b.n_vars]);
        let w_c = self.w_c.evaluation(&point[self.w_c.n_vars..]);

        Some(add_e * (w_b + w_c) + mul_e * (w_b * w_c))
    }

    pub fn partial_evaluations_of_w(&self, points: &[F], variable_indices: &Vec<usize>) -> Self {
        let mut evaluation = self.clone();

        if points.len() != variable_indices.len() {
            panic!(
                "The length of evaluation_points and variable_indices should be the same: {}, {}",
                points.len(),
                variable_indices.len()
            );
        }

        for i in 0..points.len() {
            evaluation = evaluation.partial_evaluation_of_w(&points[i], &variable_indices[i]);
        }

        evaluation
    }

    pub fn partial_evaluation_of_w(&self, point: &F, variable_index: &usize) -> Self {
        W {
            add_i: self.add_i.partial_evaluation(point, variable_index),
            mul_i: self.mul_i.partial_evaluation(point, variable_index),
            w_b: self.w_b.clone(),
            w_c: self.w_c.clone(),
            random_r: self.random_r.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        circuit::{GKRCircuit, GKRCircuitLayer},
        gate::{Gate, GateType},
    };

    use super::*;
    use ark_bls12_381::Fr as F;

    fn generate_add_and_mul_mle() -> (MLE<F>, MLE<F>) {
        let layer_0 = GKRCircuitLayer::new(vec![Gate::new(GateType::Add, [0, 1])]);
        let layer_1 = GKRCircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
        ]);
        let layer_2 = GKRCircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
            Gate::new(GateType::Mul, [4, 5]),
            Gate::new(GateType::Mul, [6, 7]),
        ]);

        let circuit = GKRCircuit::new(vec![layer_0, layer_1, layer_2]);
        let (add_mle, mul_mle) = circuit.add_mul_mle::<F>(1);

        (add_mle, mul_mle)
    }

    #[test]
    fn test_evaluation() {
        let add_i = MLE::<F>::new(vec![
            F::from(4),
            F::from(4),
            F::from(7),
            F::from(7),
            F::from(4),
            F::from(4),
            F::from(7),
            F::from(9),
        ]);
        // f(b) = 4b
        let w_b = MLE::<F>::new(vec![F::from(0), F::from(4)]);
        // f(c) = 3c
        let w_c = MLE::<F>::new(vec![F::from(0), F::from(3)]);
        // f(a,b,c) = 2ab + bc + 3
        let mul_i = MLE::<F>::new(vec![
            F::from(3),
            F::from(3),
            F::from(3),
            F::from(4),
            F::from(3),
            F::from(3),
            F::from(5),
            F::from(6),
        ]);

        let w = W {
            add_i,
            mul_i,
            w_b,
            w_c,
            random_r: vec![F::from(2u32)],
        };

        let expected_evaulation = F::from(1023u32);
        let evaluation = w.evaluate_of_w(&[F::from(3), F::from(1)].to_vec());

        assert_eq!(evaluation, Some(expected_evaulation))
    }

    #[test]
    fn test_partial_evaluation() {
        let add_i = MLE::<F>::new(vec![
            F::from(4),
            F::from(4),
            F::from(7),
            F::from(7),
            F::from(4),
            F::from(4),
            F::from(7),
            F::from(9),
        ]);
        // f(b) = 4b
        let w_b = MLE::<F>::new(vec![F::from(0), F::from(4)]);
        // f(c) = 3c
        let w_c = MLE::<F>::new(vec![F::from(0), F::from(3)]);
        // f(a,b,c) = 2ab + bc + 3
        let mul_i = MLE::<F>::new(vec![
            F::from(3),
            F::from(3),
            F::from(3),
            F::from(4),
            F::from(3),
            F::from(3),
            F::from(5),
            F::from(6),
        ]);
        let random_r = vec![F::from(2u32)];

        let w = W {
            add_i,
            mul_i,
            w_b: w_b.clone(),
            w_c: w_c.clone(),
            random_r: random_r.clone(),
        };

        let expected_partial_evaulation = W {
            add_i: MLE::<F>::new(vec![F::from(4), F::from(13)]),
            mul_i: MLE::<F>::new(vec![F::from(3), F::from(10)]),
            w_b,
            w_c,
            random_r,
        };

        let partial_evaluation =
            w.partial_evaluations_of_w(&[F::from(3), F::from(1)].to_vec(), &[0, 1].to_vec());

        assert_eq!(partial_evaluation, expected_partial_evaulation)
    }
}
