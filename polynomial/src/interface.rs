use crate::{multilinear::coefficient_form::MultiLinearMonomial, UnivariatePolynomial};
use ark_ff::PrimeField;

pub trait UnivariatePolynomialTrait<F: PrimeField>: Clone {
    fn new(data: Vec<F>) -> Self;
    fn evaluate(&self, point: F) -> F;
    fn interpolation(points: &[(F, F)]) -> UnivariatePolynomial<F>;
    fn degree(&self) -> F;
}

pub trait MultiLinearPolynomialTrait<F: PrimeField> {
    fn new(terms: Vec<MultiLinearMonomial<F>>) -> Self;
    fn partial_eval(&self, eval_points: F) -> Self;
    fn evaluation(&self, eval_points: &Vec<F>) -> F;
    fn degree(&self) -> usize;
}

pub trait MLETrait<F: PrimeField> {
    fn new(evaluations: Vec<F>) -> Self;
    fn partial_evaluation(&self, eval_point: F, variable_index: usize) -> Self;
    fn evaluation(&self, evaluation_points: &[F]) -> F;
    fn relabel(&self) -> Self;
}
