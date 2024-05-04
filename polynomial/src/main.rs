use core::num;
use std::{
    fmt::{Display, Error, Formatter, Result},
    ops::{Add, Mul},
};

#[derive(Debug, Clone, Copy)]
pub struct Monomial {
    coeff: i8,
    pow: i8,
}

#[derive(Debug, Clone)]
pub struct Polynomial {
    monomial: Vec<Monomial>,
}

trait PolynomialTrait {
    fn new(data: Vec<i8>) -> Self;
    fn evaluate(&self, point: u8) -> i32;
    fn interpolation(&self, point_x: Vec<i8>, point_y: Vec<i8>) -> Self;
    fn degree(&self) -> usize;
}

impl PolynomialTrait for Polynomial {
    fn new(data: Vec<i8>) -> Self {
        let mut monomial: Vec<Monomial> = Vec::new();

        for n in (0..=data.len()).step_by(2) {
            if n < data.len() {
                let monomial_value: Monomial = Monomial {
                    coeff: data[n],
                    pow: data[n + 1],
                };

                monomial.push(monomial_value);
            }
        }

        Polynomial { monomial }
    }

    fn evaluate(&self, point: u8) -> i32 {
        let mut point_evaluation = 0;

        // 5 + 2x + 4x^6 at x = 2
        // (5 * 1) + (2 * 2) + (4 * 64)
        // 5 + 4 + 256 => 265
        for n in self.monomial.iter() {
            let coefficient = n.coeff;
            let power = point.pow(n.pow as u32);

            let evaluation: i32 = coefficient as i32 * power as i32;
            point_evaluation += evaluation;
        }

        point_evaluation
    }

    fn interpolation(&self, point_x: Vec<i8>, point_y: Vec<i8>) -> Self {
        // point_x=[-2, -1, 0, 4]
        // point_y=[-2, 4, 1, 8]
        //
        //     (x-x1)(x-x2)(x-x3)
        // y * ------------------
        //     (x0-x1)(x0-x2)(x0-x3)
        //
        //      (x-(-1))(x-0)(x-4)             (x-1)(x-0)(x-4)
        // -2 * ---------------------  =  -2 * --------------------------
        //      (-2-(-1))(-2-0)(-2-4)           -1.-2.-6
        //                                     (x^2-x)(x-4)
        //                             =  -2 * ------------
        //                                         -12
        //                                     x^3-5x^2+4x
        //                             =  -2 * -----------
        //                                         -12
        let mut result_polynomial = Polynomial { monomial: vec![] };

        for (i, &x_i) in point_x.iter().enumerate() {
            let mut term = Polynomial {
                monomial: vec![Monomial {
                    coeff: point_y[i],
                    // pow: 1,
                    pow: 0,
                }],
            };
            let mut denominator = 1;
            println!("term=1={:?}", term);

            for (j, &x_j) in point_x.iter().enumerate() {
                if i != j {
                    // Multiply the current term by (x - x_j)
                    term = term
                        * Polynomial {
                            monomial: vec![
                                Monomial { coeff: 1, pow: 1 },
                                Monomial {
                                    coeff: -(x_j),
                                    pow: 0,
                                },
                            ],
                        };

                    denominator *= x_i - x_j;
                }
            }
            
            for monomial in &mut term.monomial {
                monomial.coeff /= denominator;
            }

            result_polynomial = result_polynomial + term;
        }

        result_polynomial
    }

    /// return the degree of a polynomial
    fn degree(&self) -> usize {
        let mut highest_degree = 0;
        for m in self.monomial.iter() {
            if m.pow > highest_degree {
                highest_degree = m.pow;
            }
        }

        highest_degree.try_into().unwrap()
    }
}

impl Mul for Polynomial {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        // (3x^2 + 5x + 6)(2x^2 + 4x + 5)
        // (6x^4 + 12x^3 + 15x^2 + 10x^3 + 20x^2 + 25x + 12x^2 + 24x + 30)
        // (6x^4 + 22x^3 + 47x^2  + 49x + 30)

        let mut result_monomial: Vec<Monomial> = Vec::new();

        for lhs_mn in &self.monomial {
            for rhs_mn in &rhs.monomial {
                let new_coeff = lhs_mn.coeff * rhs_mn.coeff;
                let new_pow = lhs_mn.pow + rhs_mn.pow;

                // Combine like terms
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

        Polynomial {
            monomial: result_monomial,
        }
    }
}

impl Add for Polynomial {
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

        Polynomial {
            monomial: result_monomials,
        }
    }
}

impl Display for Polynomial {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        for (index, mn) in self.monomial.iter().enumerate() {
            if index == 0 {
                write!(f, "{}", mn.coeff)?;
            } else {
                if mn.pow == 0 || mn.coeff == 0 {
                    write!(f, " + {}", mn.coeff)?;
                } else if mn.pow == 1 {
                    write!(f, " + {}x", mn.coeff)?;
                } else {
                    write!(f, " + {}x^{}", mn.coeff, mn.pow)?;
                }
            }
        }
        Ok(())
    }
}

fn main() {
    // 5 + 2x + 4x^6 at x = 2
    let data = vec![5, 0, 2, 1, 4, 6];
    let data_2 = vec![4, 0, 3, 2, 4, 5];
    // (5 + 2x + 4x^6)(4 + 3x^2 + 4x^5)
    // (20 + 15x^2 + 20x^5 + 8x + 6x^3 + 8x^6 + 16x^6 + 12x^8 + 16x^11)
    // (20 + 15x^2 + 20x^5 + 8x + 6x^3 + 24x^6 + 12x^8 + 16x^11)

    /*
    // (3x^2 + 5x + 6)(2x^2 + 4x + 5)
    // (6x^4 + 12x^3 + 15x^2 + 10x^3 + 20x^2 + 25x + 12x^2 + 24x + 30)
    // (6x^4 + 22x^3 + 47x^2  + 49x + 30)
     */
    let polynomial = Polynomial::new(data.clone());
    println!("polynomial={}", polynomial);
    let evaluation = polynomial.evaluate(2);
    // let evaluation = polynomial.degree();
    println!("evaluation={}", evaluation);

    let point_x: Vec<i8> = vec![-2, -1, 0, 4];
    let point_y: Vec<i8> = vec![-2, 4, -1, 8];
    let interpolation = polynomial.interpolation(point_x, point_y);
    println!("interpolation={}", interpolation);

    let polynomial_1 = Polynomial::new(data);
    let polynomial_2 = Polynomial::new(data_2);

    let m = polynomial_1.mul(polynomial_2);
    println!("multiply={}", m);
}
