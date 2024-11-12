use ark_ff::PrimeField;
use num_complex::Complex64;

use crate::{
    utils::{convert_prime_field_to_f64, fft},
    DenseUnivariatePolynomial, UnivariatePolynomialTrait,
};

#[derive(Clone, PartialEq, Eq, Default, Debug)]
pub struct Domain<F: PrimeField> {
    /// This is a const size of the domain
    pub(crate) size: u64,
    /// This is the generator of the domain, ofter regarded as the root of unity (omega)
    pub(crate) generator: F,
    /// This is the inverse of the group generator
    pub(crate) group_gen_inverse: F,
    /// This is the inverse of the group size
    pub(crate) group_size_inverse: F,
}

#[derive(Debug, Clone)]
pub struct UnivariateEval<F: PrimeField> {
    pub values: Vec<F>,
    pub domain: Domain<F>,
}

impl<F: PrimeField> Domain<F> {
    pub fn new(num_of_coeffs: usize) -> Domain<F> {
        let size = if num_of_coeffs.is_power_of_two() {
            num_of_coeffs
        } else {
            num_of_coeffs.checked_next_power_of_two().unwrap()
        } as u64;

        let generator = F::get_root_of_unity(size).unwrap();
        let group_gen_inverse = generator.inverse().unwrap();
        let group_size_inverse = F::from(size).inverse().unwrap();

        Domain {
            size,
            generator,
            group_gen_inverse,
            group_size_inverse,
        }
    }
}

impl<F: PrimeField> UnivariateEval<F> {
    pub fn new(values: Vec<F>, domain: Domain<F>) -> Self {
        UnivariateEval { values, domain }
    }

    pub fn to_coefficients(&self) -> Vec<F> {
        let evals = self.values.clone();

        let eval_in_complex_form: Vec<Complex64> = evals
            .iter()
            .map(|&x| Complex64::new(convert_prime_field_to_f64(x), 0.0))
            .collect();

        let inverse_fft = fft(&eval_in_complex_form, true);
        inverse_fft
            .iter()
            .take(evals.len())
            .map(|i| F::from(i.re.round() as u64))
            .collect()
    }

    pub fn to_coefficient_poly(&self) -> DenseUnivariatePolynomial<F> {
        let coefficients = self.to_coefficients();
        DenseUnivariatePolynomial::new(coefficients)
    }
}
