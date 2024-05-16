use ark_ff::PrimeField;
use num_bigint::BigUint;
use std::{
    fmt::{Display, Formatter, Result},
    ops::{Add, Mul},
};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Monomial<F: PrimeField> {
    pub coeff: F,
    pub pow: F,
}

#[derive(Debug, Clone, PartialEq)]
pub struct UnivariatePolynomial<F: PrimeField> {
    pub monomial: Vec<Monomial<F>>,
}

pub trait UnivariatePolynomialTrait<F: PrimeField>: Clone {
    fn new(data: Vec<F>) -> Self;
    fn evaluate(&self, point: F) -> F;
    fn interpolation(points: &[(F, F)]) -> UnivariatePolynomial<F>;
    fn degree(&self) -> F;
}

impl<F: PrimeField> UnivariatePolynomialTrait<F> for UnivariatePolynomial<F> {
    fn new(data: Vec<F>) -> Self {
        let mut monomial: Vec<Monomial<F>> = Vec::new();

        for n in (0..data.len()).step_by(2) {
            if n < data.len() - 1 {
                let monomial_value: Monomial<F> = Monomial {
                    coeff: data[n],
                    pow: data[n + 1],
                };
                monomial.push(monomial_value);
            } else if n == data.len() - 1 {
                // Handle the case when the length of data is odd
                let monomial_value: Monomial<F> = Monomial {
                    coeff: data[n],
                    pow: F::zero(),
                };
                monomial.push(monomial_value);
            }
        }

        UnivariatePolynomial { monomial }
    }

    fn evaluate(&self, point: F) -> F {
        let mut point_evaluation: F = F::from(0_u8);

        // 5 + 2x + 4x^6 at x = 2
        // (5 * 1) + (2 * 2) + (4 * 64)
        // 5 + 4 + 256 => 265
        for n in self.monomial.iter() {
            let coefficient = n.coeff;
            let n_pow: <F as PrimeField>::BigInt = n.pow.into();
            let power = point.pow(&n_pow);

            let evaluation = coefficient * power;
            point_evaluation += evaluation;
        }

        point_evaluation
    }

    fn interpolation(points: &[(F, F)]) -> UnivariatePolynomial<F> {
        let mut result_polynomial: UnivariatePolynomial<F> =
            UnivariatePolynomial { monomial: vec![] };
        let zero = F::zero();
        let one = F::one();
        let mut coefficients = vec![zero; points.len()];
        let mut temp_polynomial = Vec::with_capacity(points.len());

        for (i, (x1, y1)) in points.iter().enumerate() {
            temp_polynomial.clear();
            temp_polynomial.push(one);
            let mut denominator = one;

            for (j, (x2, _)) in points.iter().enumerate() {
                if i != j {
                    temp_polynomial.push(zero);
                    for k in (1..temp_polynomial.len()).rev() {
                        let xyz = x2.mul(temp_polynomial[k - 1]);
                        temp_polynomial[k] -= xyz;
                    }
                    denominator *= x1.sub(x2);
                }
            }

            let multiplier = y1.div(denominator);

            for (result_coefficient, temp_coefficient) in
                coefficients.iter_mut().zip(temp_polynomial.iter())
            {
                *result_coefficient += multiplier * temp_coefficient;
            }
        }

        for i in 0..coefficients.len() {
            let coeff = coefficients[i];
            let power = F::from(BigUint::from(i));
            let monomial = Monomial { coeff, pow: power };
            result_polynomial.monomial.push(monomial);
        }

        result_polynomial
    }

    /// return the degree of a polynomial
    fn degree(&self) -> F {
        let mut highest_degree: F = F::from(0_u8);
        for m in self.monomial.iter() {
            if m.pow > highest_degree {
                highest_degree = m.pow;
            }
        }

        highest_degree.try_into().unwrap()
    }
}

impl<F: PrimeField> Mul for UnivariatePolynomial<F> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        // (3x^2 + 5x + 6)(2x^2 + 4x + 5)
        // (6x^4 + 12x^3 + 15x^2 + 10x^3 + 20x^2 + 25x + 12x^2 + 24x + 30)
        // (6x^4 + 22x^3 + 47x^2  + 49x + 30)

        let mut result_monomial: Vec<Monomial<F>> = Vec::new();

        for lhs_mn in &self.monomial {
            for rhs_mn in &rhs.monomial {
                let new_coeff = lhs_mn.coeff * rhs_mn.coeff;
                let new_pow = lhs_mn.pow + rhs_mn.pow;

                let mut found_like_terms = false;
                for res_mn in &mut result_monomial {
                    if res_mn.pow == new_pow {
                        res_mn.coeff += new_coeff;
                        found_like_terms = true;
                        break;
                    }
                }

                if !found_like_terms {
                    result_monomial.push(Monomial {
                        coeff: new_coeff,
                        pow: new_pow,
                    });
                }
            }
        }

        UnivariatePolynomial {
            monomial: result_monomial,
        }
    }
}

impl<F: PrimeField> Add for UnivariatePolynomial<F> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let mut result_monomials = Vec::new();
        let mut lhs_iter = self.monomial.into_iter();
        let mut rhs_iter = rhs.monomial.into_iter();
        let mut lhs_mn = lhs_iter.next();
        let mut rhs_mn = rhs_iter.next();

        while lhs_mn.is_some() || rhs_mn.is_some() {
            match (lhs_mn, rhs_mn) {
                (Some(l), Some(r)) => {
                    if l.pow == r.pow {
                        result_monomials.push(Monomial {
                            coeff: l.coeff + r.coeff,
                            pow: l.pow,
                        });
                        lhs_mn = lhs_iter.next();
                        rhs_mn = rhs_iter.next();
                    } else if l.pow < r.pow {
                        result_monomials.push(l);
                        lhs_mn = lhs_iter.next();
                    } else {
                        result_monomials.push(r);
                        rhs_mn = rhs_iter.next();
                    }
                }
                (Some(l), None) => {
                    result_monomials.push(l);
                    lhs_mn = lhs_iter.next();
                }
                (None, Some(r)) => {
                    result_monomials.push(r);
                    rhs_mn = rhs_iter.next();
                }
                (None, None) => break,
            }
        }

        UnivariatePolynomial {
            monomial: result_monomials,
        }
    }
}

impl<F: PrimeField> Display for UnivariatePolynomial<F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        for (index, mn) in self.monomial.iter().enumerate() {
            if index == 0 {
                write!(f, "{}", mn.coeff)?;
            } else {
                if mn.pow == F::from(0_u8) || mn.coeff == F::from(0_u8) {
                    write!(f, " + {}", mn.coeff)?;
                } else if mn.pow == F::from(1_u8) {
                    write!(f, " + {}x", mn.coeff)?;
                } else {
                    write!(f, " + {}x^{}", mn.coeff, mn.pow)?;
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::univariate::UnivariatePolynomialTrait;

    use super::UnivariatePolynomial;
    use ark_ff::MontConfig;
    use ark_ff::{Fp64, MontBackend};

    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    struct FqConfig;
    type Fq = Fp64<MontBackend<FqConfig, 1>>;

    #[test]
    fn test_polynomial_evaluation() {
        // 5 + 2x + 4x^6 at x = 2
        let data = vec![
            Fq::from(5_u8),
            Fq::from(0_u8),
            Fq::from(2_u8),
            Fq::from(1_u8),
            Fq::from(4_u8),
            Fq::from(6_u8),
        ];
        let polynomial = UnivariatePolynomial::new(data);
        let evaluation = polynomial.evaluate(Fq::from(2_u8));

        assert_eq!(evaluation, Fq::from(10_u8));
    }

    #[test]
    fn test_polynomial_addition() {
        let data = vec![Fq::from(5_u8), Fq::from(0_u8)];
        let data2 = vec![Fq::from(2_u8), Fq::from(1_u8)];

        // 5 + 2x
        // 7x = [7,1]
        let polynomial1 = UnivariatePolynomial::new(data);
        let polynomial2 = UnivariatePolynomial::new(data2);

        assert_eq!(
            polynomial1 + polynomial2,
            UnivariatePolynomial::new(vec![
                Fq::from(5_u8),
                Fq::from(0_u8),
                Fq::from(2_u8),
                Fq::from(1_u8),
            ])
        );

        let data3 = vec![
            Fq::from(5_u8),
            Fq::from(0_u8),
            Fq::from(5_u8),
            Fq::from(2_u8),
        ];
        let data4 = vec![
            Fq::from(2_u8),
            Fq::from(1_u8),
            Fq::from(2_u8),
            Fq::from(2_u8),
        ];

        // 5 + 5x^2  +  2x + 2x^2
        // 5 + 2x + 7x^2 = [5, 0, 2, 1, 7, 2]
        let polynomial1 = UnivariatePolynomial::new(data3);
        let polynomial2 = UnivariatePolynomial::new(data4);

        assert_eq!(
            polynomial1 + polynomial2,
            UnivariatePolynomial::new(vec![
                Fq::from(5_u8),
                Fq::from(0_u8),
                Fq::from(2_u8),
                Fq::from(1_u8),
                Fq::from(7_u8),
                Fq::from(2_u8),
            ])
        );
    }

    #[test]
    fn test_polynomial_multiplication() {
        let data = vec![Fq::from(5_u8), Fq::from(0_u8)];
        let data2 = vec![Fq::from(2_u8), Fq::from(1_u8)];

        // 5 x 2x = 7x = [7,1]
        let polynomial1 = UnivariatePolynomial::new(data);
        let polynomial2 = UnivariatePolynomial::new(data2);

        assert_eq!(
            polynomial1 + polynomial2,
            UnivariatePolynomial::new(vec![
                Fq::from(5_u8),
                Fq::from(0_u8),
                Fq::from(2_u8),
                Fq::from(1_u8),
            ])
        );

        let data3 = vec![
            Fq::from(5_u8),
            Fq::from(0_u8),
            Fq::from(2_u8),
            Fq::from(1_u8),
            Fq::from(4_u8),
            Fq::from(6_u8),
        ];
        let data4 = vec![
            Fq::from(4),
            Fq::from(0),
            Fq::from(3),
            Fq::from(2),
            Fq::from(4),
            Fq::from(5),
        ];

        let polynomial3 = UnivariatePolynomial::new(data3);
        let polynomial4 = UnivariatePolynomial::new(data4);

        // 3 + 15x^2 + 3x^5 + 8x + 6x^3 + 7x^6 + 12x^8 + 16x^11
        assert_eq!(
            polynomial3 * polynomial4,
            UnivariatePolynomial::new(vec![
                Fq::from(3_u8),
                Fq::from(0_u8),
                Fq::from(15_u8),
                Fq::from(2_u8),
                Fq::from(3_u8),
                Fq::from(5_u8),
                Fq::from(8_u8),
                Fq::from(1_u8),
                Fq::from(6_u8),
                Fq::from(3_u8),
                Fq::from(7_u8),
                Fq::from(6_u8),
                Fq::from(12_u8),
                Fq::from(8_u8),
                Fq::from(16_u8),
                Fq::from(11_u8),
            ])
        );
    }

    #[test]
    fn test_polynomial_interpolation() {
        let point_x = vec![Fq::from(2), Fq::from(1), Fq::from(0), Fq::from(4)];
        let point_y = vec![Fq::from(2), Fq::from(4), Fq::from(1), Fq::from(8)];
        let interpolation = UnivariatePolynomial::interpolation(&[
            (Fq::from(1), Fq::from(2)),
            (Fq::from(2), Fq::from(3)),
            (Fq::from(4), Fq::from(11)),
        ]);
        let interpolation_check = UnivariatePolynomial::new(vec![
            Fq::from(1_u8),
            Fq::from(0_u8),
            Fq::from(15_u8),
            Fq::from(1_u8),
            Fq::from(3_u8),
            Fq::from(2_u8),
        ]);
        assert_eq!(interpolation, interpolation_check);

        // to test the evaluation of the polynomial
        let evaluation = interpolation.evaluate(Fq::from(2_u8));
        assert_eq!(evaluation, Fq::from(9_u8));
    }

    #[test]
    fn test_polynomial_degree() {
        let data = vec![Fq::from(2), Fq::from(1), Fq::from(2), Fq::from(4)];

        let polynomial = UnivariatePolynomial::new(data);
        let degree = polynomial.degree();

        assert_eq!(degree, Fq::from(4_u8));
    }
}
