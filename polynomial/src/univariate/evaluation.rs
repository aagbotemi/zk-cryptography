use crate::{DenseUnivariatePolynomial, UnivariatePolynomialTrait};
use ark_ff::PrimeField;

use super::domain::Domain;

#[derive(Debug)]
pub struct UnivariateEval<F: PrimeField> {
    /// this is a list of the evaluation of the polynomial
    pub values: Vec<F>,
    /// This is the domian of the polynomal; very important for the FFT and IFFT
    pub domain: Domain<F>,
}

impl<F: PrimeField> UnivariateEval<F> {
    /// This function is used to create a new polynomial from the evaluation form
    pub fn new(values: Vec<F>, domain: Domain<F>) -> Self {
        UnivariateEval { values, domain }
    }

    /// This function is used to create a new polynomial from the evaluation form but also does checks
    pub fn new_checked(values: Vec<F>, domain: Domain<F>) -> Result<Self, &'static str> {
        if values.len() != domain.size() as usize {
            return Err("The size of the values does not match the size of the domain");
        }
        Ok(UnivariateEval { values, domain })
    }

    /// This function performs interpolation on the vaules provided and returns a polynomial
    pub fn interpolate(values: Vec<F>, domain: Domain<F>) -> DenseUnivariatePolynomial<F> {
        let coeffs = domain.ifft(&values);
        DenseUnivariatePolynomial::new(coeffs)
    }

    /// This function is used to convert the coefficient form of the polynomial to the evaluation form
    pub fn from_coefficients(coefficients: Vec<F>) -> Self {
        let mut coeffs = coefficients.clone();
        let domain = Domain::<F>::new(coefficients.len() as usize);
        let evals = domain.fft(&mut coeffs);

        UnivariateEval {
            values: evals,
            domain,
        }
    }

    /// This function is used to convert the evaluation form of the polynomial to the coefficient form
    pub fn to_coefficients(&self) -> Vec<F> {
        let evals = self.values.clone();
        self.domain.ifft(&evals)
    }

    /// This function is used to convert the evaluation form of the polynomial to the coefficient form as a polynomial
    pub fn to_coefficient_poly(&self) -> DenseUnivariatePolynomial<F> {
        let coefficients = self.to_coefficients();
        DenseUnivariatePolynomial::new(coefficients)
    }

    /// This function is used to multiply two polynomials in the evaluation form
    pub fn multiply(
        poly1: &DenseUnivariatePolynomial<F>,
        poly2: &DenseUnivariatePolynomial<F>,
    ) -> DenseUnivariatePolynomial<F> {
        let mut poly1_coeffs = poly1.coefficients.clone();
        let mut poly2_coeffs = poly2.coefficients.clone();

        let length_of_poly_unscaled = poly1_coeffs.len() + poly2_coeffs.len() - 1;
        let length_of_poly = if length_of_poly_unscaled.is_power_of_two() {
            length_of_poly_unscaled
        } else {
            length_of_poly_unscaled.checked_next_power_of_two().unwrap()
        };
        let domain = Domain::<F>::new(length_of_poly);
        poly1_coeffs.resize(length_of_poly, F::ZERO);
        poly2_coeffs.resize(length_of_poly, F::ZERO);

        let poly_1_eval = domain.fft(&poly1_coeffs);
        let poly_2_eval = domain.fft(&poly2_coeffs);

        let mut result = vec![F::ZERO; length_of_poly];
        for i in 0..length_of_poly {
            result[i] = poly_1_eval[i] * poly_2_eval[i];
        }

        let coeff = domain.ifft(&result);
        DenseUnivariatePolynomial::new(coeff[..length_of_poly_unscaled].to_vec())
    }
}
