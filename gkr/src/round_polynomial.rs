use ark_ff::PrimeField;
use polynomial::{interface::MLETrait, MLE};

pub struct W<F: PrimeField> {
    add_i: MLE<F>,
    mul_i: MLE<F>,
    w_b: MLE<F>,
    w_c: MLE<F>,
}

impl<F: PrimeField> W<F> {
    /// Create a new `W` polynomial.
    pub fn new(add_i: MLE<F>, mul_i: MLE<F>, w_b: MLE<F>, w_c: MLE<F>) -> Self {
        Self {
            add_i,
            mul_i,
            w_b,
            w_c,
        }
    }

    pub fn evaluate_W(&self, point: &[F], b: &F, c: &F) -> Option<F> {
        let add_e = self.add_i.evaluation(point);
        let mul_e = self.mul_i.evaluation(point);

        let w_b = self.w_b.evaluation(point);
        let w_c = self.w_c.evaluation(point);

        Some(add_e * (w_b + w_c) + mul_e * (w_b * w_c))
    }

    pub fn partial_evaluate_W(&self, point: &[F], variable_index: usize) -> Self {
        W {
            add_i: self
                .add_i
                .partial_evaluation(point[variable_index], variable_index),
            mul_i: self
                .mul_i
                .partial_evaluation(point[variable_index], variable_index),
            w_b: self
                .w_b
                .partial_evaluation(point[variable_index], variable_index),
            w_c: self
                .w_c
                .partial_evaluation(point[variable_index], variable_index),
        }
    }
}
