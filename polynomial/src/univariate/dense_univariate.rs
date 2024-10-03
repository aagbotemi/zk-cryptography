use crate::{
    utils::{dense_langrange_basis, get_langrange_basis, remove_trailing_and_redundant_zeros},
    UnivariatePolynomialTrait,
};
use ark_ff::{BigInteger, PrimeField};
use std::{
    fmt::{Display, Formatter, Result},
    ops::{Add, AddAssign, Mul, Neg, Sub, SubAssign},
};

#[derive(Debug, PartialEq, Clone)]
pub struct DenseUnivariatePolynomial<F: PrimeField> {
    pub coefficients: Vec<F>,
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

    pub fn interpolate(point_ys: Vec<F>, point_xs: Vec<F>) -> Self {
        let langrange_poly_vec = dense_langrange_basis(&point_xs, &point_ys);
        let langrange_poly = langrange_poly_vec
            .iter()
            .fold(DenseUnivariatePolynomial::new(vec![]), |acc, x| {
                acc + x.clone()
            });

        langrange_poly
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
    use crate::Fq;

    #[test]
    fn test_dense_polynomial_evaluation() {
        // 5 + 2x + 4x^2 at x = 2
        // 5 + 4 + 16
        let data = vec![Fq::from(5_u8), Fq::from(2_u8), Fq::from(4_u8)];
        let polynomial = DenseUnivariatePolynomial::new(data);
        let evaluation = polynomial.evaluate(Fq::from(2_u8));

        assert_eq!(evaluation, Fq::from(8_u8));
    }

    #[test]
    fn test_dense_polynomial_addition() {
        // 5 + 2x + 4x^2
        let polynomial1 =
            DenseUnivariatePolynomial::new(vec![Fq::from(5_u8), Fq::from(2_u8), Fq::from(4_u8)]);
        // 4 + 3x
        let polynomial2 = DenseUnivariatePolynomial::new(vec![Fq::from(4_u8), Fq::from(3_u8)]);

        assert_eq!(
            polynomial1 + polynomial2,
            // 9 + 5x + 4x^2
            DenseUnivariatePolynomial::new(vec![Fq::from(9_u8), Fq::from(5_u8), Fq::from(4_u8),])
        );

        // 5 + 5x^2
        let polynomial1 =
            DenseUnivariatePolynomial::new(vec![Fq::from(5_u8), Fq::from(0_u8), Fq::from(5_u8)]);
        // 2x + 2x^2
        let polynomial2 =
            DenseUnivariatePolynomial::new(vec![Fq::from(0_u8), Fq::from(2_u8), Fq::from(2_u8)]);

        assert_eq!(
            polynomial1 + polynomial2,
            // 5 + 2x + 7x^2
            DenseUnivariatePolynomial::new(vec![Fq::from(5_u8), Fq::from(2_u8), Fq::from(7_u8),])
        );
    }

    #[test]
    fn test_dense_polynomial_multiplication() {
        let polynomial1 =
            DenseUnivariatePolynomial::new(vec![Fq::from(1), Fq::from(3), Fq::from(2)]);
        let polynomial2 = DenseUnivariatePolynomial::new(vec![Fq::from(3), Fq::from(2)]);

        assert_eq!(
            polynomial1 * polynomial2,
            DenseUnivariatePolynomial::new(vec![
                Fq::from(3),
                Fq::from(11),
                Fq::from(12),
                Fq::from(4)
            ])
        );

        // (3x^2 + 5x + 6)
        let polynomial3 =
            DenseUnivariatePolynomial::new(vec![Fq::from(6), Fq::from(5), Fq::from(3)]);
        // (2x^2 + 4x + 5)
        let polynomial4 =
            DenseUnivariatePolynomial::new(vec![Fq::from(5), Fq::from(4), Fq::from(2)]);

        assert_eq!(
            polynomial3 * polynomial4,
            // (6x^4 + 22x^3 + 47x^2  + 49x + 30)
            DenseUnivariatePolynomial::new(vec![
                Fq::from(30_u8),
                Fq::from(49_u8),
                Fq::from(47_u8),
                Fq::from(22_u8),
                Fq::from(6_u8),
            ])
        );
    }

    #[test]
    // #[ignore = "reason"]
    fn test_dense_polynomial_interpolation_1() {
        let point_ys_1 = vec![Fq::from(0), Fq::from(4), Fq::from(16)];
        let point_xs_1 = vec![Fq::from(0), Fq::from(2), Fq::from(4)];

        let poly = DenseUnivariatePolynomial::interpolate(point_ys_1, point_xs_1);
        assert_eq!(
            poly,
            DenseUnivariatePolynomial::new(vec![Fq::from(0), Fq::from(0), Fq::from(1)])
        );

        let point_xs_2 = vec![Fq::from(2), Fq::from(1), Fq::from(0), Fq::from(4)];
        let point_ys_2 = vec![Fq::from(2), Fq::from(4), Fq::from(1), Fq::from(8)];
        let interpolation = DenseUnivariatePolynomial::interpolate(point_ys_2, point_xs_2);

        let interpolation_check = DenseUnivariatePolynomial::new(vec![
            Fq::from(1),
            Fq::from(9),
            Fq::from(5),
            Fq::from(6),
        ]);
        assert_eq!(interpolation, interpolation_check);

        // to test the evaluation of the polynomial
        let evaluation = interpolation.evaluate(Fq::from(2_u8));
        assert_eq!(evaluation, Fq::from(2_u8));
    }

    #[test]
    fn test_dense_polynomial_interpolation_2() {
        let point_ys_1 = vec![Fq::from(5), Fq::from(7), Fq::from(13)];
        let point_xs_1 = vec![Fq::from(0), Fq::from(1), Fq::from(2)];

        let poly = DenseUnivariatePolynomial::interpolate(point_ys_1, point_xs_1);
        assert_eq!(
            poly,
            DenseUnivariatePolynomial::new(vec![Fq::from(5), Fq::from(0), Fq::from(2)])
        );

        // fq_from_vec(vec![0, 1, 3, 4, 5, 8]),
        // fq_from_vec(vec![12, 48, 3150, 11772, 33452, 315020]),
        let point_ys_2 = vec![
            Fq::from(12),
            Fq::from(48),
            Fq::from(3150),
            Fq::from(11772),
            Fq::from(33452),
            Fq::from(315020),
        ];
        let point_xs_2 = vec![
            Fq::from(0),
            Fq::from(1),
            Fq::from(3),
            Fq::from(4),
            Fq::from(5),
            Fq::from(8),
        ];

        let poly = DenseUnivariatePolynomial::interpolate(point_ys_2, point_xs_2);
        let eval = poly.evaluate(Fq::from(3));
        println!("{:?}", eval);
        assert_eq!(
            poly,
            DenseUnivariatePolynomial::new(vec![
                Fq::from(12),
                Fq::from(8),
                Fq::from(1),
                Fq::from(7),
                Fq::from(12),
                Fq::from(8)
            ])
        );

        let point_ys_3 = vec![Fq::from(565), Fq::from(1631), Fq::from(3537), Fq::from(-7)];
        let point_xs_3 = vec![Fq::from(5), Fq::from(7), Fq::from(9), Fq::from(1)];

        let poly = DenseUnivariatePolynomial::interpolate(point_ys_3, point_xs_3);
        assert_eq!(
            poly,
            DenseUnivariatePolynomial::new(vec![
                Fq::from(0),
                Fq::from(-12),
                Fq::from(0),
                Fq::from(5)
            ])
        );
    }

    #[test]
    fn test_dense_polynomial_interpolation_3() {
        let point_ys_1 = vec![Fq::from(0), Fq::from(1)];
        let point_xs_1 = vec![Fq::from(0), Fq::from(2)];
        let interpolation1 = DenseUnivariatePolynomial::interpolate(point_ys_1, point_xs_1);
        assert_eq!(
            interpolation1,
            DenseUnivariatePolynomial::new(vec![Fq::from(0), Fq::from(9)])
        );
        let evaluation1 = interpolation1.evaluate(Fq::from(2));
        assert_eq!(evaluation1, Fq::from(1));

        let point_ys_2 = vec![Fq::from(0), Fq::from(5), Fq::from(14)];
        let point_xs_2 = vec![Fq::from(0), Fq::from(1), Fq::from(2)];
        let interpolation2 = DenseUnivariatePolynomial::interpolate(point_ys_2, point_xs_2);
        assert_eq!(
            interpolation2,
            DenseUnivariatePolynomial::new(vec![Fq::from(0), Fq::from(3), Fq::from(2)])
        );
        let evaluation2 = interpolation2.evaluate(Fq::from(2));
        assert_eq!(evaluation2, Fq::from(14));

        let point_ys_3 = vec![
            Fq::from(6),
            Fq::from(11),
            Fq::from(18),
            Fq::from(27),
            Fq::from(38),
        ];
        let point_xs_3 = vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(5),
        ];

        let interpolation3 = DenseUnivariatePolynomial::interpolate(point_ys_3, point_xs_3);
        assert_eq!(
            interpolation3,
            DenseUnivariatePolynomial::new(vec![
                Fq::from(3),
                Fq::from(2),
                Fq::from(1),
                Fq::from(0),
                Fq::from(0)
            ])
        );
        let evaluation3 = interpolation3.evaluate(Fq::from(2));
        assert_eq!(evaluation3, Fq::from(11));
    }

    #[test]
    fn test_polynomial_degree() {
        let data = vec![Fq::from(2), Fq::from(1), Fq::from(2), Fq::from(4)];

        let polynomial = DenseUnivariatePolynomial::new(data);
        let degree = polynomial.degree();

        assert_eq!(degree, 3);
    }

    #[test]
    #[ignore = "reason"]
    fn test_fft_multiplication() {
        // 1 * (1 + 2x + 3x^2 ) + x * (1 + 2x + 3x^2 ) + x^2 * (1 + 2x + 3x^2 )
        // 1 + 2x + 3x^2 + x + 2x^2 + 3x^3 + x^2 + 2x^3 + 3x^4
        // 3x^3 + 2x^3 + 3x^4
        // 1 + 3x + 6x^2 + 5x^3 + 3x^4
        let poly_a = vec![Fq::from(1), Fq::from(2), Fq::from(3)]; // Coefficients of A(x)
        let poly_b = vec![Fq::from(1), Fq::from(1), Fq::from(1)]; // Coefficients of B(x)

        // let result = DenseUnivariatePolynomial::multiply_polynomials(&poly_a, &poly_b);

        let poly_1 = DenseUnivariatePolynomial::new(vec![
            Fq::from(1_u8),
            Fq::from(0_u8),
            Fq::from(2_u8),
            Fq::from(1_u8),
            Fq::from(3_u8),
            Fq::from(2_u8),
        ]);

        let poly_2 = DenseUnivariatePolynomial::new(vec![
            Fq::from(1_u8),
            Fq::from(0_u8),
            Fq::from(1_u8),
            Fq::from(1_u8),
            Fq::from(1_u8),
            Fq::from(2_u8),
        ]);

        let result_2 = poly_1 * poly_2;
        // assert_eq!(result, result_2)
    }
}
