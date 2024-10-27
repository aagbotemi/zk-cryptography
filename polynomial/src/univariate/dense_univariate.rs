use crate::{
    utils::{
        convert_prime_field_to_f64, dense_langrange_basis, fft, remove_trailing_and_redundant_zeros,
    },
    UnivariatePolynomialTrait,
};
use ark_ff::{BigInteger, PrimeField, Zero};
use num_complex::{Complex, Complex64};
use std::{
    fmt::{Display, Formatter, Result},
    ops::{Add, AddAssign, Div, Index, IndexMut, Mul, Neg, Rem, Sub, SubAssign},
};

#[derive(Debug, PartialEq, Clone)]
pub struct DenseUnivariatePolynomial<F: PrimeField> {
    pub coefficients: Vec<F>,
}

impl<F: PrimeField> Index<usize> for DenseUnivariatePolynomial<F> {
    type Output = F;

    fn index(&self, index: usize) -> &Self::Output {
        &self.coefficients[index]
    }
}

impl<F: PrimeField> IndexMut<usize> for DenseUnivariatePolynomial<F> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.coefficients[index]
    }
}

impl<F: PrimeField> DenseUnivariatePolynomial<F> {
    pub fn zero() -> Self {
        DenseUnivariatePolynomial {
            coefficients: vec![],
        }
    }
    pub fn is_zero(&self) -> bool {
        self.coefficients.is_empty()
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();
        for c in self.coefficients.iter() {
            bytes.extend_from_slice(&c.into_bigint().to_bytes_be());
        }
        bytes
    }

    pub fn remove_leading_zeros(&self) -> Self {
        let coefficients = remove_trailing_and_redundant_zeros(&self.coefficients);
        Self { coefficients }
    }

    /// This function creates a new polynomial from a list of coefficients vector
    pub fn from_coefficients_vec(coeffs: Vec<F>) -> Self {
        DenseUnivariatePolynomial {
            coefficients: coeffs,
        }
    }

    /// This function returns an array of coefficients of this polynomial
    pub fn coefficients(&self) -> &[F] {
        &self.coefficients
    }

    pub fn leading_coefficient(&self) -> Option<F> {
        self.coefficients.last().cloned()
    }

    pub fn iter_with_index(&self) -> Vec<(usize, F)> {
        self.coefficients.iter().cloned().enumerate().collect()
    }

    pub fn interpolate(point_ys: Vec<F>, point_xs: Vec<F>) -> Self {
        let langrange_poly_vec = dense_langrange_basis(&point_xs, &point_ys);
        let langrange_poly = langrange_poly_vec
            .iter()
            .fold(DenseUnivariatePolynomial::new(vec![]), |acc, x| {
                acc + x.clone()
            });

        langrange_poly
    }

    /// This function is used for poly division, returning the quotient and remainder
    pub fn divide_with_q_and_r(
        &self,
        divisor: &Self,
    ) -> Option<(DenseUnivariatePolynomial<F>, DenseUnivariatePolynomial<F>)> {
        if self.is_zero() {
            Some((
                DenseUnivariatePolynomial::zero(),
                DenseUnivariatePolynomial::zero(),
            ))
        } else if divisor.is_zero() {
            panic!("Dividing by zero polynomial")
        } else if self.degree() < divisor.degree() {
            Some((DenseUnivariatePolynomial::zero(), self.clone().into()))
        } else {
            // Now we know that self.degree() >= divisor.degree();
            let mut quotient = vec![F::zero(); self.degree() - divisor.degree() + 1];
            let mut remainder: DenseUnivariatePolynomial<F> = self.clone().into();
            // Can unwrap here because we know self is not zero.
            let divisor_leading_inv = divisor.leading_coefficient().unwrap().inverse().unwrap();
            while !remainder.is_zero() && remainder.degree() >= divisor.degree() {
                let cur_q_coeff = *remainder.coefficients.last().unwrap() * divisor_leading_inv;
                let cur_q_degree = remainder.degree() - divisor.degree();
                quotient[cur_q_degree] = cur_q_coeff;

                for (i, div_coeff) in divisor.iter_with_index() {
                    remainder[cur_q_degree + i] -= &(cur_q_coeff * div_coeff);
                }
                while let Some(true) = remainder.coefficients.last().map(|c| c.is_zero()) {
                    remainder.coefficients.pop();
                }
            }
            Some((
                DenseUnivariatePolynomial::from_coefficients_vec(quotient),
                remainder,
            ))
        }
    }

    // (3xy + 2x + 4z + 3) (2xy + 3z + 4)
    // 6x^2y^2 .....+ 25z + 12 // 8

    pub fn fft_mult_poly(
        polya: &DenseUnivariatePolynomial<F>,
        polyb: &DenseUnivariatePolynomial<F>,
    ) -> Self {
        let poly1 = polya.coefficients.clone();
        let poly2 = polyb.coefficients.clone();

        let coefficient_length_of_resultant_poly = poly1.len() + poly2.len() - 1;

        let coefficient_length_of_resultant_poly_pow_of_2 =
            coefficient_length_of_resultant_poly.next_power_of_two();

        let mut poly1_in_complex_form: Vec<Complex64> = poly1
            .iter()
            .map(|&x| Complex64::new(convert_prime_field_to_f64(x), 0.0))
            .collect();
        let mut poly2_in_complex_form: Vec<Complex64> = poly2
            .iter()
            .map(|&x| Complex64::new(convert_prime_field_to_f64(x), 0.0))
            .collect();

        poly1_in_complex_form.resize(
            coefficient_length_of_resultant_poly_pow_of_2,
            Complex64::new(0.0, 0.0),
        );
        poly2_in_complex_form.resize(
            coefficient_length_of_resultant_poly_pow_of_2,
            Complex64::new(0.0, 0.0),
        );

        let fft_poly1 = fft(&poly1_in_complex_form, false);
        let fft_poly2 = fft(&poly2_in_complex_form, false);

        let mut element_wise_product = vec![Complex::zero(); fft_poly1.len()];
        for i in 0..fft_poly1.len() {
            element_wise_product[i] = fft_poly1[i] * fft_poly2[i];
        }

        let inverse_fft = fft(&element_wise_product, true);

        let result: Vec<F> = inverse_fft
            .iter()
            .take(coefficient_length_of_resultant_poly)
            .map(|i| F::from(i.re.round() as u64))
            .collect();

        Self::new(result)
    }
}

impl<F: PrimeField> UnivariatePolynomialTrait<F> for DenseUnivariatePolynomial<F> {
    fn new(coefficients: Vec<F>) -> Self {
        DenseUnivariatePolynomial { coefficients }
    }

    fn evaluate(&self, point: F) -> F {
        // 5 + 2x + 4x^6 at x = 2
        // (5 * 1) + (2 * 2) + (4 * 64)
        // (5 * 1) + (2 * 2) + (4 * 64)
        // 5 + 4 + 256 => 265
        let mut point_evaluation: F = F::zero();
        for (i, coeff) in self.coefficients.iter().enumerate() {
            let power = point.pow([i as u64]);
            point_evaluation += *coeff * power;
        }

        point_evaluation
    }

    /// return the degree of a polynomial
    fn degree(&self) -> usize {
        let coefficients = self.remove_leading_zeros();
        // coefficients.coefficients.len() - 1
        if coefficients.coefficients.is_empty() {
            0
        } else {
            coefficients.coefficients.len() - 1
        }
    }
}

impl<F: PrimeField> Mul for DenseUnivariatePolynomial<F> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        // (3x^2 + 5x + 6)(2x^2 + 4x + 5)
        // (6x^4 + 12x^3 + 15x^2 + 10x^3 + 20x^2 + 25x + 12x^2 + 24x + 30)
        // (6x^4 + 22x^3 + 47x^2  + 49x + 30)

        if self.is_zero() || rhs.is_zero() {
            return DenseUnivariatePolynomial::new(vec![]);
        }

        let degree = self.degree() + rhs.degree();
        let mut result = vec![F::zero(); degree + 1];

        for i in 0..=self.degree() {
            for j in 0..=rhs.degree() {
                result[i + j] += self.coefficients[i] * rhs.coefficients[j];
            }
        }

        DenseUnivariatePolynomial {
            coefficients: result,
        }
    }
}

impl<F: PrimeField> Add for DenseUnivariatePolynomial<F> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let result = if self.degree() >= other.degree() {
            let mut result_coff = Vec::new();

            for i in 0..self.coefficients.len() {
                result_coff
                    .push(self.coefficients[i] + other.coefficients.get(i).unwrap_or(&F::zero()));
            }

            DenseUnivariatePolynomial::new(result_coff)
        } else {
            let mut result_coff = Vec::new();

            for i in 0..other.coefficients.len() {
                result_coff
                    .push(other.coefficients[i] + self.coefficients.get(i).unwrap_or(&F::zero()));
            }

            DenseUnivariatePolynomial::new(result_coff)
        };

        result
    }
}

impl<F: PrimeField> AddAssign for DenseUnivariatePolynomial<F> {
    fn add_assign(&mut self, rhs: Self) {
        if self.degree() >= rhs.degree() {
            for i in 0..self.coefficients.len() {
                self.coefficients[i] += rhs.coefficients.get(i).unwrap_or(&F::zero());
            }
        } else {
            let mut result_coff = self.coefficients.clone();

            for i in 0..rhs.coefficients.len() {
                result_coff
                    .push(rhs.coefficients[i] + self.coefficients.get(i).unwrap_or(&F::zero()));
            }

            self.coefficients = result_coff;
        }
    }
}

impl<F: PrimeField> Sub<F> for DenseUnivariatePolynomial<F> {
    type Output = Self;

    fn sub(self, other: F) -> Self {
        // check for zero polynomials
        if self.is_zero() {
            return DenseUnivariatePolynomial::new(vec![other]);
        }

        let mut sub_coefficients = self.coefficients.clone();
        sub_coefficients[0] -= other;

        DenseUnivariatePolynomial::new(sub_coefficients)
    }
}

impl<F: PrimeField> Sub<F> for &DenseUnivariatePolynomial<F> {
    type Output = DenseUnivariatePolynomial<F>;

    fn sub(self, other: F) -> Self::Output {
        // check for zero polynomials
        if self.is_zero() {
            return DenseUnivariatePolynomial::new(vec![other]);
        }

        let mut sub_coefficients = self.coefficients.clone();
        sub_coefficients[0] -= other;

        DenseUnivariatePolynomial::new(sub_coefficients)
    }
}

impl<F: PrimeField> Sub for DenseUnivariatePolynomial<F> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let negated_other = -other;
        self + negated_other
    }
}

impl<F: PrimeField> SubAssign for DenseUnivariatePolynomial<F> {
    fn sub_assign(&mut self, rhs: Self) {
        *self += -rhs;
    }
}

impl<F: PrimeField> Neg for DenseUnivariatePolynomial<F> {
    type Output = DenseUnivariatePolynomial<F>;

    #[inline]
    fn neg(mut self) -> DenseUnivariatePolynomial<F> {
        self.coefficients.iter_mut().for_each(|coeff| {
            *coeff = -*coeff;
        });

        self
    }
}

// impl<F: PrimeField> SubAssign for DenseUnivariatePolynomial<F> {
//     fn sub_assign(&mut self, rhs: Self) {
//         *self += -rhs;
//     }
// }

impl<F: PrimeField> Div for DenseUnivariatePolynomial<F> {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        self.divide_with_q_and_r(&other).expect("division failed").0
    }
}

impl<F: PrimeField> Rem for DenseUnivariatePolynomial<F> {
    type Output = Self;

    fn rem(self, other: Self) -> Self {
        self.divide_with_q_and_r(&other).expect("division failed").1
    }
}

impl<F: PrimeField> Display for DenseUnivariatePolynomial<F> {
    fn fmt(&self, f: &mut Formatter) -> Result {
        for (i, coeff) in self
            .coefficients
            .iter()
            .enumerate()
            .filter(|(_, c)| !c.is_zero())
        {
            if i == 0 {
                write!(f, "\n{:?}", coeff)?;
            } else if i == 1 {
                write!(f, " + \n{:?} * x", coeff)?;
            } else {
                write!(f, " + \n{:?} * x^{}", coeff, i)?;
            }
        }
        Ok(())
    }
}

mod tests {
    use super::*;
    use crate::utils::generate_random;
    use ark_test_curves::bls12_381::Fr;

    #[test]
    fn test_dense_polynomial_evaluation() {
        // 5 + 2x + 4x^2 at x = 2
        // 5 + 4 + 16
        let data = vec![Fr::from(5), Fr::from(2), Fr::from(4)];
        let polynomial = DenseUnivariatePolynomial::new(data);
        let evaluation = polynomial.evaluate(Fr::from(2));

        assert_eq!(evaluation, Fr::from(25));
    }

    #[test]
    fn test_dense_polynomial_addition() {
        // 5 + 2x + 4x^2
        let polynomial1 =
            DenseUnivariatePolynomial::new(vec![Fr::from(5), Fr::from(2), Fr::from(4)]);
        // 4 + 3x
        let polynomial2 = DenseUnivariatePolynomial::new(vec![Fr::from(4), Fr::from(3)]);

        assert_eq!(
            polynomial1 + polynomial2,
            // 9 + 5x + 4x^2
            DenseUnivariatePolynomial::new(vec![Fr::from(9), Fr::from(5), Fr::from(4)])
        );

        // 5 + 5x^2
        let polynomial1 =
            DenseUnivariatePolynomial::new(vec![Fr::from(5), Fr::from(0), Fr::from(5)]);
        // 2x + 2x^2
        let polynomial2 =
            DenseUnivariatePolynomial::new(vec![Fr::from(0), Fr::from(2), Fr::from(2)]);

        assert_eq!(
            polynomial1 + polynomial2,
            // 5 + 2x + 7x^2
            DenseUnivariatePolynomial::new(vec![Fr::from(5), Fr::from(2), Fr::from(7)])
        );
    }

    #[test]
    fn test_dense_polynomial_multiplication() {
        let polynomial1 =
            DenseUnivariatePolynomial::new(vec![Fr::from(1), Fr::from(3), Fr::from(2)]);
        let polynomial2 = DenseUnivariatePolynomial::new(vec![Fr::from(3), Fr::from(2)]);

        assert_eq!(
            polynomial1 * polynomial2,
            DenseUnivariatePolynomial::new(vec![
                Fr::from(3),
                Fr::from(11),
                Fr::from(12),
                Fr::from(4)
            ])
        );

        // (3x^2 + 5x + 6)
        let polynomial3 =
            DenseUnivariatePolynomial::new(vec![Fr::from(6), Fr::from(5), Fr::from(3)]);
        // (2x^2 + 4x + 5)
        let polynomial4 =
            DenseUnivariatePolynomial::new(vec![Fr::from(5), Fr::from(4), Fr::from(2)]);

        assert_eq!(
            polynomial3 * polynomial4,
            // (6x^4 + 22x^3 + 47x^2  + 49x + 30)
            DenseUnivariatePolynomial::new(vec![
                Fr::from(30_u8),
                Fr::from(49_u8),
                Fr::from(47_u8),
                Fr::from(22_u8),
                Fr::from(6_u8),
            ])
        );
    }

    #[test]
    // #[ignore = "reason"]
    fn test_dense_polynomial_interpolation_1() {
        let point_ys_1 = vec![Fr::from(0), Fr::from(4), Fr::from(16)];
        let point_xs_1 = vec![Fr::from(0), Fr::from(2), Fr::from(4)];

        let poly = DenseUnivariatePolynomial::interpolate(point_ys_1, point_xs_1);
        assert_eq!(
            poly,
            DenseUnivariatePolynomial::new(vec![Fr::from(0), Fr::from(0), Fr::from(1)])
        );
    }

    #[test]
    fn test_dense_polynomial_interpolation_2() {
        let point_ys_1 = vec![Fr::from(5), Fr::from(7), Fr::from(13)];
        let point_xs_1 = vec![Fr::from(0), Fr::from(1), Fr::from(2)];

        let poly = DenseUnivariatePolynomial::interpolate(point_ys_1, point_xs_1);
        assert_eq!(
            poly,
            DenseUnivariatePolynomial::new(vec![Fr::from(5), Fr::from(0), Fr::from(2)])
        );

        // fq_from_vec(vec![0, 1, 3, 4, 5, 8]),
        // fq_from_vec(vec![12, 48, 3150, 11772, 33452, 315020]),
        let point_ys_2 = vec![
            Fr::from(12),
            Fr::from(48),
            Fr::from(3150),
            Fr::from(11772),
            Fr::from(33452),
            Fr::from(315020),
        ];
        let point_xs_2 = vec![
            Fr::from(0),
            Fr::from(1),
            Fr::from(3),
            Fr::from(4),
            Fr::from(5),
            Fr::from(8),
        ];

        let poly = DenseUnivariatePolynomial::interpolate(point_ys_2, point_xs_2);
        let eval = poly.evaluate(Fr::from(3));
        println!("{:?}", eval);
        assert_eq!(
            poly,
            DenseUnivariatePolynomial::new(vec![
                Fr::from(12),
                Fr::from(8),
                Fr::from(1),
                Fr::from(7),
                Fr::from(12),
                Fr::from(8)
            ])
        );

        let point_ys_3 = vec![Fr::from(565), Fr::from(1631), Fr::from(3537), Fr::from(-7)];
        let point_xs_3 = vec![Fr::from(5), Fr::from(7), Fr::from(9), Fr::from(1)];

        let poly = DenseUnivariatePolynomial::interpolate(point_ys_3, point_xs_3);
        assert_eq!(
            poly,
            DenseUnivariatePolynomial::new(vec![
                Fr::from(0),
                Fr::from(-12),
                Fr::from(0),
                Fr::from(5)
            ])
        );
    }

    #[test]
    fn test_dense_polynomial_interpolation_3() {
        let point_ys_1 = vec![Fr::from(0), Fr::from(1)];
        let point_xs_1 = vec![Fr::from(0), Fr::from(2)];
        let interpolation1 = DenseUnivariatePolynomial::interpolate(point_ys_1, point_xs_1);
        // assert_eq!(
        //     interpolation1,
        //     DenseUnivariatePolynomial::new(vec![Fr::from(0), Fr::from(9)])
        // );
        let evaluation1 = interpolation1.evaluate(Fr::from(2));
        assert_eq!(evaluation1, Fr::from(1));

        let point_ys_2 = vec![Fr::from(0), Fr::from(5), Fr::from(14)];
        let point_xs_2 = vec![Fr::from(0), Fr::from(1), Fr::from(2)];
        let interpolation2 = DenseUnivariatePolynomial::interpolate(point_ys_2, point_xs_2);
        assert_eq!(
            interpolation2,
            DenseUnivariatePolynomial::new(vec![Fr::from(0), Fr::from(3), Fr::from(2)])
        );
        let evaluation2 = interpolation2.evaluate(Fr::from(2));
        assert_eq!(evaluation2, Fr::from(14));

        let point_ys_3 = vec![
            Fr::from(6),
            Fr::from(11),
            Fr::from(18),
            Fr::from(27),
            Fr::from(38),
        ];
        let point_xs_3 = vec![
            Fr::from(1),
            Fr::from(2),
            Fr::from(3),
            Fr::from(4),
            Fr::from(5),
        ];

        let interpolation3 = DenseUnivariatePolynomial::interpolate(point_ys_3, point_xs_3);
        assert_eq!(
            interpolation3,
            DenseUnivariatePolynomial::new(vec![
                Fr::from(3),
                Fr::from(2),
                Fr::from(1),
                Fr::from(0),
                Fr::from(0)
            ])
        );
        let evaluation3 = interpolation3.evaluate(Fr::from(2));
        assert_eq!(evaluation3, Fr::from(11));
    }

    #[test]
    fn test_polynomial_degree() {
        let data = vec![Fr::from(2), Fr::from(1), Fr::from(2), Fr::from(4)];

        let polynomial = DenseUnivariatePolynomial::new(data);
        let degree = polynomial.degree();

        assert_eq!(degree, 3);
    }

    #[test]
    fn test_fft_multiplication() {
        // (1 + 2x + 3x^2) * (4 + 5x + 6x^2)
        // 4 * (1 + 2x + 3x^2) + 5x * (1 + 2x + 3x^2) + 6x^2 * (1 + 2x + 3x^2)
        // 4 + 8x + 12x^2 + 5x + 10x^2 + 15x^3 + 6x^2 + 12x^3 + 18x^4
        // 4 + 13x + 28x^2 + 27x^3 + 18x^4
        // [4, 13, 28, 27, 18]
        let poly_a = DenseUnivariatePolynomial::new(vec![Fr::from(1), Fr::from(2), Fr::from(3)]);
        let poly_b = DenseUnivariatePolynomial::new(vec![Fr::from(4), Fr::from(5), Fr::from(6)]);
        let result = DenseUnivariatePolynomial::fft_mult_poly(&poly_a, &poly_b);
        let expected_result = DenseUnivariatePolynomial::new(vec![
            Fr::from(4),
            Fr::from(13),
            Fr::from(28),
            Fr::from(27),
            Fr::from(18),
        ]);

        assert_eq!(result, expected_result);
    }

    #[test]
    fn test_fft_multiplication_2() {
        // (1 + 2x + 3x^2) * (4 + 5x + 6x^2)
        // 4 * (1 + 2x + 3x^2) + 5x * (1 + 2x + 3x^2) + 6x^2 * (1 + 2x + 3x^2)
        // 4 + 8x + 12x^2 + 5x + 10x^2 + 15x^3 + 6x^2 + 12x^3 + 18x^4
        // 4 + 13x + 28x^2 + 27x^3 + 18x^4
        // [4, 13, 28, 27, 18]
        let tt: Vec<Fr> = generate_random(7);
        let poly_a: DenseUnivariatePolynomial<Fr> = DenseUnivariatePolynomial::new(tt);
        println!("poly_a={:?}", poly_a);
        let poly_b: DenseUnivariatePolynomial<Fr> =
            DenseUnivariatePolynomial::new(generate_random(10));
        println!("poly_b={:?}", poly_b);
        let result = DenseUnivariatePolynomial::fft_mult_poly(&poly_a, &poly_b);
        println!("result={:?}", result);
        // let expected_result = DenseUnivariatePolynomial::new(vec![
        //     Fr::from(4),
        //     Fr::from(13),
        //     Fr::from(28),
        //     Fr::from(27),
        //     Fr::from(18),
        // ]);

        // assert_eq!(result, expected_result);
    }
}
