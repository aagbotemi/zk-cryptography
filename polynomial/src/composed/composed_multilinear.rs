use crate::{
    interface::{ComposedMultilinearTrait, MultilinearTrait},
    Multilinear,
};
use ark_ff::PrimeField;

#[derive(Debug, Clone)]
pub struct ComposedMultilinear<F: PrimeField> {
    polys: Vec<Multilinear<F>>,
}

impl<F: PrimeField> ComposedMultilinear<F> {
    pub fn new(polys: Vec<Multilinear<F>>) -> Self {
        let n_vars = polys[0].n_vars;
        assert!(polys.iter().all(|p| p.n_vars == n_vars));

        ComposedMultilinear { polys }
    }

    pub fn n_vars(&self) -> usize {
        self.polys[0].n_vars
    }

    pub fn zero(&self) -> Self {
        Self { polys: vec![] }
    }

    pub fn is_zero(&self) -> bool {
        if self.polys.len() == 0 {
            return true;
        } else {
            if self.polys.iter().all(|p| p.evaluations.is_empty()) {
                return true;
            } else {
                return false;
            }
        }
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();

        for poly in &self.polys {
            bytes.extend_from_slice(&poly.to_bytes());
        }

        bytes
    }
}

impl<F: PrimeField> MultilinearTrait<F> for ComposedMultilinear<F> {
    fn evaluation(&self, points: &[F]) -> F {
        let mut result = F::one();

        for poly in &self.polys {
            let poly_eval = poly.evaluation(points);
            result *= poly_eval;
        }

        result
    }

    fn partial_evaluation(
        &self,
        evaluation_point: &F,
        variable_index: &usize,
    ) -> ComposedMultilinear<F> {
        let mut new_polys = Vec::new();

        for poly in &self.polys {
            new_polys.push(poly.partial_evaluation(&evaluation_point, &variable_index));
        }

        ComposedMultilinear { polys: new_polys }
    }

    fn partial_evaluations(
        &self,
        evaluation_points: &[F],
        variable_index: &Vec<usize>,
    ) -> ComposedMultilinear<F> {
        let mut new_polys = self.clone();

        if evaluation_points.len() != variable_index.len() {
            panic!(
                "The length of evaluation_points and variable_index should be the same: {}, {}",
                evaluation_points.len(),
                variable_index.len()
            )
        }

        for i in 0..evaluation_points.len() {
            new_polys = new_polys.partial_evaluation(&evaluation_points[i], &variable_index[i])
        }

        new_polys
    }
}

impl<F: PrimeField> ComposedMultilinearTrait<F> for ComposedMultilinear<F> {
    fn max_degree(&self) -> usize {
        self.polys.len()
    }

    fn element_wise_product(&self) -> Vec<F> {
        let length_of_poly = &self.polys[0].evaluations.len();

        (0..*length_of_poly)
            .map(|i| self.polys.iter().map(|v| v.evaluations[i]).product())
            .collect()
    }

    fn element_wise_add(&self) -> Vec<F> {
        let length_of_poly = &self.polys[0].evaluations.len();

        (0..*length_of_poly)
            .map(|i| self.polys.iter().map(|v| v.evaluations[i]).sum())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use field_tracker::Ft;

    use super::*;
    use crate::interface::MultilinearTrait;
    use crate::Fq as Fq_old;
    use crate::Multilinear;

    type Fq = Ft<1, Fq_old>;

    #[test]
    fn test_evaluation() {
        let mle1 = Multilinear::new(vec![Fq::from(0), Fq::from(1), Fq::from(2), Fq::from(3)]);
        let mle2 = Multilinear::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(1)]);

        let polys = ComposedMultilinear::new(vec![mle1, mle2]);
        let evaluation = polys.evaluation(&vec![Fq::from(2), Fq::from(3)]);

        assert_eq!(evaluation, Fq::from(42));
        println!("{}", Fq::summary());
    }

    #[test]
    fn test_partial_evaluation() {
        let mle1 = Multilinear::new(vec![Fq::from(0), Fq::from(1), Fq::from(2), Fq::from(3)]);
        let mle2 = Multilinear::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(1)]);

        let polys = ComposedMultilinear::new(vec![mle1, mle2]);
        let partial_evaluation = polys.partial_evaluation(&Fq::from(2), &0);

        let evaluation = partial_evaluation.evaluation(&vec![Fq::from(3)]);
        assert_eq!(evaluation, Fq::from(42));
        println!("{}", Fq::summary());
    }

    #[test]
    fn test_element_wise_product() {
        let mle1 = Multilinear::new(vec![Fq::from(0), Fq::from(1), Fq::from(2), Fq::from(3)]);
        let mle2 = Multilinear::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(1)]);

        let polys = ComposedMultilinear::new(vec![mle1, mle2]);
        let element_product = polys.element_wise_product();
        assert_eq!(
            element_product,
            vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(3)]
        );
        println!("{}", Fq::summary());
    }

    #[test]
    fn test_element_wise_add() {
        let mle1 = Multilinear::new(vec![Fq::from(0), Fq::from(1), Fq::from(2), Fq::from(3)]);
        let mle2 = Multilinear::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(1)]);

        let polys = ComposedMultilinear::new(vec![mle1, mle2]);
        let element_product = polys.element_wise_add();
        assert_eq!(
            element_product,
            vec![Fq::from(0), Fq::from(1), Fq::from(2), Fq::from(4)]
        );
        println!("{}", Fq::summary());
    }

    #[test]
    fn test_n_vars_and_max_degree() {
        let mle1 = Multilinear::new(vec![Fq::from(0), Fq::from(1), Fq::from(2), Fq::from(3)]);
        let mle2 = Multilinear::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(1)]);

        let mles_1 = ComposedMultilinear::new(vec![mle1, mle2]);
        let n_vars_1 = mles_1.n_vars();
        let max_degree_1 = mles_1.max_degree();
        assert_eq!(n_vars_1, 2);
        assert_eq!(max_degree_1, 2);

        let mle3 = Multilinear::new(vec![Fq::from(0); 32]);
        let mle4 = Multilinear::new(vec![Fq::from(1); 32]);

        let mles_2 = ComposedMultilinear::new(vec![mle3, mle4]);
        let n_vars_2 = mles_2.n_vars();
        let max_degree_2 = mles_2.max_degree();
        assert_eq!(n_vars_2, 5);
        assert_eq!(max_degree_2, 2);
        println!("{}", Fq::summary());
    }
}
