use crate::{
    interface::UnivariatePolynomialTrait,
    utils::{
        convert_prime_field_to_f64, fft, get_even_indexed_coefficients,
        get_odd_indexed_coefficients, ifft, lagrange_basis,
    },
};
use ark_ff::{BigInteger, PrimeField, Zero};
use num_complex::Complex;
use std::{
    f64::consts::PI,
    fmt::{Display, Formatter, Result},
    ops::{Add, Mul},
};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UnivariateMonomial<F: PrimeField> {
    pub coeff: F,
    pub pow: F,
}

#[derive(Debug, PartialEq)]
pub struct UnivariatePolynomial<F: PrimeField> {
    pub monomial: Vec<UnivariateMonomial<F>>,
}

impl<F: PrimeField> UnivariatePolynomial<F> {
    pub fn zero() -> Self {
        UnivariatePolynomial { monomial: vec![] }
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();
        for p in self.monomial.iter() {
            bytes.extend_from_slice(&p.coeff.into_bigint().to_bytes_be());
            bytes.extend_from_slice(&p.pow.into_bigint().to_bytes_be());
        }
        bytes
    }

    pub fn from_coefficients(&self) -> Vec<F> {
        self.monomial.iter().map(|mn| mn.coeff).collect()
    }

    pub fn multiply_polynomials(poly1: &Vec<F>, poly2: &Vec<F>) -> Vec<F> {
        let mut a = poly1.clone();
        let mut b = poly2.clone();

        if !a.len().is_power_of_two() {
            a.resize(a.len() + 1, F::zero());
        }
        if !b.len().is_power_of_two() {
            b.resize(b.len() + 1, F::zero());
        }

        let fft_a = fft(&a);
        let fft_b = fft(&b);

        let mut fft_result = vec![Complex::zero(); fft_a.len()];
        for i in 0..fft_a.len() {
            fft_result[i] = fft_a[i] * fft_b[i];
        }

        let ifft_result = ifft(&fft_result);

        let result: Vec<_> = ifft_result
            .iter()
            .map(|i| F::from(i.re.round() as u64))
            .collect();

        dbg!(&result);

        result
    }
}

impl<F: PrimeField> UnivariatePolynomialTrait<F> for UnivariatePolynomial<F> {
    fn new(data: Vec<F>) -> Self {
        let mut monomial: Vec<UnivariateMonomial<F>> = Vec::new();

        for n in (0..data.len()).step_by(2) {
            if n < data.len() - 1 {
                let monomial_value: UnivariateMonomial<F> = UnivariateMonomial {
                    coeff: data[n],
                    pow: data[n + 1],
                };
                monomial.push(monomial_value);
            } else if n == data.len() - 1 {
                // Handle the case when the length of data is odd
                let monomial_value: UnivariateMonomial<F> = UnivariateMonomial {
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
        let mut result: Vec<F> = vec![F::zero(); points.len()];

        for (i, &(_, y_i)) in points.iter().enumerate() {
            let l_i: Vec<F> = lagrange_basis(points, i);
            let l_i: Vec<F> = l_i.into_iter().map(|coeff| coeff * y_i).collect();

            for (k, &coeff) in l_i.iter().enumerate() {
                result[k] += coeff;
            }
        }

        let monomial: Vec<UnivariateMonomial<F>> = result
            .into_iter()
            .enumerate()
            .filter(|&(_, coeff)| coeff != F::zero())
            .map(|(pow, coeff)| UnivariateMonomial {
                coeff,
                pow: F::from(pow as u64),
            })
            .collect();

        UnivariatePolynomial { monomial }
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

        let mut result_monomial: Vec<UnivariateMonomial<F>> = Vec::new();

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
                    result_monomial.push(UnivariateMonomial {
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
                        result_monomials.push(UnivariateMonomial {
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

mod tests {
    use std::result;

    use super::*;
    use crate::Fq;

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
        // let point_x = vec![Fq::from(2), Fq::from(1), Fq::from(0), Fq::from(4)];
        // let point_y = vec![Fq::from(2), Fq::from(4), Fq::from(1), Fq::from(8)];
        let interpolation = UnivariatePolynomial::interpolation(&[
            (Fq::from(1), Fq::from(2)),
            (Fq::from(2), Fq::from(3)),
            (Fq::from(4), Fq::from(11)),
        ]);

        let interpolation_check = UnivariatePolynomial::new(vec![
            Fq::from(3_u8),
            Fq::from(0_u8),
            Fq::from(15_u8),
            Fq::from(1_u8),
            Fq::from(1_u8),
            Fq::from(2_u8),
        ]);
        assert_eq!(interpolation, interpolation_check);

        // to test the evaluation of the polynomial
        let evaluation = interpolation.evaluate(Fq::from(2_u8));
        assert_eq!(evaluation, Fq::from(3_u8));
    }

    #[test]
    fn test_polynomial_interpolation_1() {
        // [(0, 0), (1, 2)] // 2x, where x = 2, 4
        let interpolation1 = UnivariatePolynomial::interpolation(&vec![
            (Fq::from(0_u8), Fq::from(0_u8)),
            (Fq::from(1_u8), Fq::from(2_u8)),
        ]);
        let evaluation1 = interpolation1.evaluate(Fq::from(2_u8));
        assert_eq!(evaluation1, Fq::from(4_u8));

        // [(0, 5), (1, 7), (2, 13)] // 2x^2 + 5, where x = 2, 13
        let interpolation2 = UnivariatePolynomial::interpolation(&vec![
            (Fq::from(0_u8), Fq::from(5_u8)),
            (Fq::from(1_u8), Fq::from(7_u8)),
            (Fq::from(2_u8), Fq::from(13_u8)),
        ]);
        let evaluation2 = interpolation2.evaluate(Fq::from(2_u8));
        assert_eq!(evaluation2, Fq::from(13_u8));

        // [(0, 12), (1,48), (3,3150), (4,11772), (5,33452), (8,315020)] // 8x^5 + 12x^4 + 7x^3 + 1x^2 + 8x + 12, where x = 1, 48
        let interpolation3 = UnivariatePolynomial::interpolation(&vec![
            (Fq::from(0_u8), Fq::from(12_u8)),
            (Fq::from(1_u8), Fq::from(48_u8)),
            (Fq::from(3_u8), Fq::from(3150_u16)),
            (Fq::from(4_u8), Fq::from(11772_u16)),
            (Fq::from(5_u8), Fq::from(33452_u16)),
            (Fq::from(8_u8), Fq::from(315020_u32)),
        ]);
        let evaluation3 = interpolation3.evaluate(Fq::from(1_u8));
        assert_eq!(evaluation3, Fq::from(48_u8));

        // [(0,0), (1,5), (2,14)], // 2x2 + 3x, where x = 2, 14
        let interpolation4 = UnivariatePolynomial::interpolation(&vec![
            (Fq::from(0_u8), Fq::from(0_u8)),
            (Fq::from(1_u8), Fq::from(5_u8)),
            (Fq::from(2_u8), Fq::from(14_u8)),
        ]);
        let evaluation4 = interpolation4.evaluate(Fq::from(2_u8));
        assert_eq!(evaluation4, Fq::from(14_u8));

        // [(1, 6), (2, 11), (3, 18), (4, 27), (5, 38)] // x^2 + 2x + 3, where x = 2, 11
        let interpolation5 = UnivariatePolynomial::interpolation(&vec![
            (Fq::from(1), Fq::from(6)),
            (Fq::from(2), Fq::from(11)),
            (Fq::from(3), Fq::from(18)),
            (Fq::from(4), Fq::from(27)),
            (Fq::from(5), Fq::from(38)),
        ]);
        assert_eq!(
            interpolation5,
            UnivariatePolynomial::new(vec![
                Fq::from(3_u8),
                Fq::from(0_u8),
                Fq::from(2_u8),
                Fq::from(1_u8),
                Fq::from(1_u8),
                Fq::from(2_u8),
            ])
        );
        let evaluation5 = interpolation5.evaluate(Fq::from(2_u8));
        assert_eq!(evaluation5, Fq::from(11_u8));
    }

    #[test]
    fn test_polynomial_degree() {
        let data = vec![Fq::from(2), Fq::from(1), Fq::from(2), Fq::from(4)];

        let polynomial = UnivariatePolynomial::new(data);
        let degree = polynomial.degree();

        assert_eq!(degree, Fq::from(4_u8));
    }
}
