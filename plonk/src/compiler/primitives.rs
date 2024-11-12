use ark_ff::PrimeField;
use polynomial::{univariate::evaluation::UnivariateEval, DenseUnivariatePolynomial};
use std::collections::HashMap;

pub struct CommonPreprocessedInput<F: PrimeField> {
    pub group_order: u64,
    pub q_l: UnivariateEval<F>,
    // pub q_l: DenseUnivariatePolynomial<F>,
    pub q_r: UnivariateEval<F>,
    // pub q_r: DenseUnivariatePolynomial<F>,
    pub q_m: UnivariateEval<F>,
    // pub q_m: DenseUnivariatePolynomial<F>,
    pub q_o: UnivariateEval<F>,
    // pub q_o: DenseUnivariatePolynomial<F>,
    pub q_c: UnivariateEval<F>,
    // pub q_c: DenseUnivariatePolynomial<F>,
    pub sigma_1: UnivariateEval<F>,
    // pub sigma_1: DenseUnivariatePolynomial<F>,
    pub sigma_2: UnivariateEval<F>,
    // pub sigma_2: DenseUnivariatePolynomial<F>,
    pub sigma_3: UnivariateEval<F>,
    // pub sigma_3: DenseUnivariatePolynomial<F>,
    pub s1_coeff: Option<DenseUnivariatePolynomial<F>>,
    pub s2_coeff: Option<DenseUnivariatePolynomial<F>>,
}

#[derive(Clone)]
pub struct Program<F: PrimeField> {
    pub constraints: Vec<AssemblyEqn<F>>,
    pub group_order: u64,
}

pub struct Witness<F: PrimeField> {
    pub a: UnivariateEval<F>,
    pub b: UnivariateEval<F>,
    pub c: UnivariateEval<F>,
    pub public_poly: UnivariateEval<F>,
}

#[derive(Debug, Clone)]
pub struct GateWire {
    pub left_wire: Option<String>,
    pub right_wire: Option<String>,
    pub output_wire: Option<String>,
}

#[derive(Debug, Clone)]
pub struct Gate<F: PrimeField> {
    pub l: F,
    pub r: F,
    pub m: F,
    pub o: F,
    pub c: F,
}

#[derive(Debug, Clone)]
pub struct AssemblyEqn<F: PrimeField> {
    pub wires: GateWire,
    pub coeffs: HashMap<Option<String>, F>,
}
