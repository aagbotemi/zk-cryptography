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

    pub fn evaluate(&self, point: &[F], b: &F, c: &F) -> Option<F> {
        let add_e = self.add_i.evaluation(point);
        let mul_e = self.mul_i.evaluation(point);

        let w_b = self.w_b.evaluation(&[*b]);
        let w_c = self.w_c.evaluation(&[*c]);

        Some(add_e * (w_b + w_c) + mul_e * (w_b * w_c))
    }
}
