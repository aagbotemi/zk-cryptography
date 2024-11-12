use ark_ff::PrimeField;
use polynomial::{DenseUnivariatePolynomial, UnivariatePolynomialTrait};

pub fn split_poly_in_3<F: PrimeField>(
    poly: &DenseUnivariatePolynomial<F>,
    group_order: usize,
) -> (
    DenseUnivariatePolynomial<F>,
    DenseUnivariatePolynomial<F>,
    DenseUnivariatePolynomial<F>,
) {
    let poly_low_coeffs = poly.coefficients[0..group_order].to_vec();
    let poly_mid_coeffs = poly.coefficients[group_order..2 * group_order].to_vec();
    let poly_high_coeffs = poly.coefficients[2 * group_order..].to_vec();

    (
        DenseUnivariatePolynomial::new(poly_low_coeffs),
        DenseUnivariatePolynomial::new(poly_mid_coeffs),
        DenseUnivariatePolynomial::new(poly_high_coeffs),
    )
}

pub fn apply_w_to_polynomial<F: PrimeField>(
    poly: &DenseUnivariatePolynomial<F>,
    w: &F,
) -> DenseUnivariatePolynomial<F> {
    let mut result = Vec::new();
    let mut w_power = F::one();

    for coeff in poly.coefficients.iter() {
        result.push(*coeff * w_power);
        w_power *= w;
    }

    DenseUnivariatePolynomial::new(result)
}

pub fn zh_values<F: PrimeField>(group_order: usize) -> Vec<F> {
    let mut zh_values = vec![F::one().neg()];
    for _ in 0..group_order - 1 {
        zh_values.push(F::zero());
    }
    zh_values.push(F::one());
    zh_values
}
