use crate::{
    interface::{ComposedMLETrait, MLETrait},
    MLE,
};
use ark_ff::PrimeField;

#[derive(Debug, Clone)]
pub struct ComposedMLE<F: PrimeField> {
    mles: Vec<MLE<F>>,
}

impl<F: PrimeField> ComposedMLE<F> {
    pub fn new(mles: Vec<MLE<F>>) -> Self {
        let n_vars = mles[0].n_vars;
        assert!(mles.iter().all(|p| p.n_vars == n_vars));

        ComposedMLE { mles }
    }

    pub fn n_vars(&self) -> usize {
        self.mles[0].n_vars
    }

    pub fn zero(&self) -> Self {
        Self { mles: vec![] }
    }

    pub fn is_zero(&self) -> bool {
        if self.mles.len() == 0 {
            return true;
        } else {
            if self.mles.iter().all(|p| p.evaluations.is_empty()) {
                return true;
            } else {
                return false;
            }
        }
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();

        for mle in &self.mles {
            bytes.extend_from_slice(&mle.evaluations_to_bytes());
        }

        bytes
    }
}

impl<F: PrimeField> ComposedMLETrait<F> for ComposedMLE<F> {
    fn evaluation(&self, points: &[F]) -> F {
        let mut result = F::one();

        for mle in &self.mles {
            let mle_eval = mle.evaluation(points);
            result *= mle_eval;
        }

        result
    }

    fn partial_evaluation(&self, eval_point: F, variable_index: usize) -> ComposedMLE<F> {
        let mut new_mles = Vec::new();

        for mle in &self.mles {
            new_mles.push(mle.partial_evaluation(eval_point, variable_index));
        }

        ComposedMLE { mles: new_mles }
    }

    fn element_product(&self) -> Vec<F> {
        let length_of_poly = &self.mles[0].evaluations.len();

        (0..*length_of_poly)
            .map(|i| self.mles.iter().map(|v| v.evaluations[i]).product())
            .collect()
    }

    fn element_add(&self) -> Vec<F> {
        let length_of_poly = &self.mles[0].evaluations.len();

        (0..*length_of_poly)
            .map(|i| self.mles.iter().map(|v| v.evaluations[i]).sum())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::MontConfig;
    use ark_ff::{Fp64, MontBackend};

    use super::*;
    use crate::interface::MLETrait;
    use crate::MLE;

    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    struct FqConfig;
    type Fq = Fp64<MontBackend<FqConfig, 1>>;

    #[test]
    fn test_evaluation() {
        let mle1 = MLE::new(vec![Fq::from(0), Fq::from(1), Fq::from(2), Fq::from(3)]);
        let mle2 = MLE::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(1)]);

        let mles = ComposedMLE::new(vec![mle1, mle2]);
        let evaluation = mles.evaluation(&vec![Fq::from(2), Fq::from(3)]);

        assert_eq!(evaluation, Fq::from(42));
    }

    #[test]
    fn test_partial_evaluation() {
        let mle1 = MLE::new(vec![Fq::from(0), Fq::from(1), Fq::from(2), Fq::from(3)]);
        let mle2 = MLE::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(1)]);

        let mles = ComposedMLE::new(vec![mle1, mle2]);
        let partial_evaluation = mles.partial_evaluation(Fq::from(2), 0);

        let evaluation = partial_evaluation.evaluation(&vec![Fq::from(3)]);
        assert_eq!(evaluation, Fq::from(42));
    }

    #[test]
    fn test_element_product() {
        let mle1 = MLE::new(vec![Fq::from(0), Fq::from(1), Fq::from(2), Fq::from(3)]);
        let mle2 = MLE::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(1)]);

        let mles = ComposedMLE::new(vec![mle1, mle2]);
        let element_product = mles.element_product();
        assert_eq!(
            element_product,
            vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(3)]
        );
    }

    #[test]
    fn test_element_add() {
        let mle1 = MLE::new(vec![Fq::from(0), Fq::from(1), Fq::from(2), Fq::from(3)]);
        let mle2 = MLE::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(1)]);

        let mles = ComposedMLE::new(vec![mle1, mle2]);
        let element_product = mles.element_add();
        assert_eq!(
            element_product,
            vec![Fq::from(0), Fq::from(1), Fq::from(2), Fq::from(4)]
        );
    }

    #[test]
    fn test_n_vars() {
        let mle1 = MLE::new(vec![Fq::from(0), Fq::from(1), Fq::from(2), Fq::from(3)]);
        let mle2 = MLE::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(1)]);

        let mles_1 = ComposedMLE::new(vec![mle1, mle2]);
        let n_vars_1 = mles_1.n_vars();
        assert_eq!(n_vars_1, 2);

        let mle3 = MLE::new(vec![Fq::from(0); 32]);
        let mle4 = MLE::new(vec![Fq::from(1); 32]);

        let mles_2 = ComposedMLE::new(vec![mle3, mle4]);
        let n_vars_2 = mles_2.n_vars();
        assert_eq!(n_vars_2, 5);
    }
}
