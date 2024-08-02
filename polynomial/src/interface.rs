use crate::{multilinear::coefficient_form::MultiLinearMonomial, UnivariatePolynomial};
use ark_ff::PrimeField;

pub trait UnivariatePolynomialTrait<F: PrimeField> {
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
    fn additive_identity(num_vars: usize) -> Self;
    fn to_bytes(&self) -> Vec<u8>;
    fn split_poly_into_two_and_sum_each_part(&mut self) -> Self;
}

pub trait ComposedMLETrait<F: PrimeField> {
    fn evaluation(&self, points: &[F]) -> F;
    fn partial_evaluation(&self, eval_point: F, variable_index: usize) -> Self;
    fn element_wise_product(&self) -> Vec<F>;
    fn element_wise_add(&self) -> Vec<F>;
}
