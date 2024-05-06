use ark_ff::PrimeField;
use std::{
    fmt::{Display, Formatter, Result},
    ops::{Add, Mul},
};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Monomial<F: PrimeField> {
    coeff: F,
    pow: F,
}

#[derive(Debug, Clone, PartialEq)]
pub struct UnivariatePolynomial<F>
where
    F: PrimeField,
{
    monomial: Vec<Monomial<F>>,
}

trait UnivariatePolynomialTrait<F: PrimeField>: Clone {
    fn new(data: Vec<F>) -> Self;
    fn evaluate(&self, point: F) -> F;
    fn interpolation(point_x: Vec<F>, point_y: Vec<F>) -> Self;
    fn degree(&self) -> F;
}

impl<F: PrimeField> UnivariatePolynomialTrait<F> for UnivariatePolynomial<F> {
    fn new(data: Vec<F>) -> Self {
        let mut monomial: Vec<Monomial<F>> = Vec::new();

        for n in (0..=data.len()).step_by(2) {
            if n < data.len() {
                let monomial_value: Monomial<F> = Monomial {
                    coeff: data[n],
                    pow: data[n + 1],
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

    fn interpolation(point_x: Vec<F>, point_y: Vec<F>) -> Self {
        // point_x=[2, 1, 0, 4]
        // point_y=[2, 4, 1, 8]
        //
        //     (x-x1)(x-x2)(x-x3)
        // y * ------------------
        //     (x0-x1)(x0-x2)(x0-x3)
        //
        //      (x-(1))(x-0)(x-4)             (x-1)(x-0)(x-4)
        // 2 * ---------------------  =  2 * --------------------------
        //      (2-1)(2-0)(2-4)                 1.2.-2
        //                                     (x^2-x)(x-4)
        //                             =  -2 * ------------
        //                                         -4
        //                                     x^3-5x^2+4x
        //                             =  -2 * -----------
        //                                         -4

        // [(2, 2), (1,4), (0,1), (4,8)]
        let mut lagrange_poly: UnivariatePolynomial<F> = UnivariatePolynomial { monomial: vec![] };

        for (i, &xi) in point_x.iter().enumerate() {
            let mut poly: UnivariatePolynomial<F> = UnivariatePolynomial {
                monomial: vec![Monomial {
                    coeff: F::from(1_u8),
                    pow: F::from(0_u8),
                }],
            };

            for (j, &xj) in point_x.iter().enumerate() {
                if i != j {
                    // Construct (x - xj)
                    let temp_poly: UnivariatePolynomial<F> = UnivariatePolynomial {
                        monomial: vec![
                            Monomial {
                                coeff: -xj,
                                pow: F::from(0_u8),
                            },
                            Monomial {
                                coeff: F::from(1_u8),
                                pow: F::from(1_u8),
                            },
                        ],
                    };

                    poly = poly * temp_poly;
                }
            }

            for monomial in &mut poly.monomial {
                monomial.coeff *= point_y[i];
            }

            lagrange_poly = lagrange_poly + poly;
        }

        lagrange_poly
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
        let interpolation = UnivariatePolynomial::interpolation(point_x, point_y);

        // 9 + 2x + 3x^2 + 15x^3
        assert_eq!(
            interpolation,
            UnivariatePolynomial::new(vec![
                Fq::from(9_u8),
                Fq::from(0_u8),
                Fq::from(2_u8),
                Fq::from(1_u8),
                Fq::from(3_u8),
                Fq::from(2_u8),
                Fq::from(15_u8),
                Fq::from(3_u8),
            ])
        );
    }

    #[test]
    fn test_polynomial_degree() {
        let data = vec![Fq::from(2), Fq::from(1), Fq::from(2), Fq::from(4)];

        let polynomial = UnivariatePolynomial::new(data);
        let degree = polynomial.degree();

        assert_eq!(degree, Fq::from(4_u8));
    }
}
