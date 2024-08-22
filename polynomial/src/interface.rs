use crate::UnivariatePolynomial;
use ark_ff::PrimeField;

pub trait UnivariatePolynomialTrait<F: PrimeField> {
    fn new(data: Vec<F>) -> Self;
    fn evaluate(&self, point: F) -> F;
    fn interpolation(points: &[(F, F)]) -> UnivariatePolynomial<F>;
    fn degree(&self) -> F;
}

pub trait MultilinearTrait<F: PrimeField> {
    fn partial_evaluation(&self, eval_point: &F, variable_index: &usize) -> Self;
    fn partial_evaluations(&self, points: &[F], variable_indices: &Vec<usize>) -> Self;
    fn evaluation(&self, evaluation_points: &[F]) -> F;
}

pub trait ComposedMultilinearTrait<F: PrimeField> {
    fn element_wise_product(&self) -> Vec<F>;
    fn element_wise_add(&self) -> Vec<F>;
    fn max_degree(&self) -> usize;
}
