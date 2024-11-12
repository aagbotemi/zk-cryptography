use crate::utils::pick_pairs_with_index;
use ark_ff::PrimeField;
use std::{
    fmt::{Display, Formatter, Result},
    ops::{Add, Mul},
};

#[derive(Debug, PartialEq)]
pub struct MultiLinearMonomial<F: PrimeField> {
    coefficient: F,
    vars: Vec<bool>,
}

#[derive(Debug, PartialEq)]
pub struct MultiLinearCoefficientPolynomial<F: PrimeField> {
    terms: Vec<MultiLinearMonomial<F>>,
}

impl<F: PrimeField> MultiLinearMonomial<F> {
    pub fn new(coefficient: F, vars: Vec<bool>) -> Self {
        assert!(
            vars.len() > 0,
            "Length of variables must be greater than zero"
        );
        MultiLinearMonomial { coefficient, vars }
    }
}

impl<F: PrimeField> MultiLinearCoefficientPolynomial<F> {
    pub fn new(terms: Vec<MultiLinearMonomial<F>>) -> Self {
        MultiLinearCoefficientPolynomial { terms }
    }

    /// partial evaluation of a polynomial
    pub fn partial_evaluation(&self, eval_points: F) -> Self {
        let mut res: MultiLinearCoefficientPolynomial<F> =
            MultiLinearCoefficientPolynomial { terms: vec![] };

        for (i, j) in pick_pairs_with_index(&self.terms) {
            let y1 = &self.terms[i].coefficient;
            let y2 = &self.terms[j].coefficient;
            // r.y1 + (1-r).y2 straight line formula
            let y = (eval_points * y2) + ((F::one() - eval_points) * y1);

            res.terms.push(MultiLinearMonomial {
                coefficient: y,
                vars: vec![false], // TODO: likely bug here
            })
        }

        res
    }

    /// full evaluation of a polynomial
    pub fn evaluation(&self, eval_points: &Vec<F>) -> F {
        let mut eval_result: F = F::zero();

        for (_, term) in self.terms.iter().enumerate() {
            let mut var_res = F::one();

            for (index_j, ev) in term.vars.iter().enumerate() {
                if *ev == true {
                    var_res *= eval_points[index_j]
                }
            }

            eval_result += term.coefficient * var_res;
        }

        eval_result
    }

    /// degree of polynomial
    pub fn degree(&self) -> usize {
        let mut max_true_count: usize = 0;

        for term in &self.terms {
            let true_count: usize = term.vars.iter().filter(|&var| *var).count();
            max_true_count = max_true_count.max(true_count);
        }

        max_true_count
    }
}

impl<F: PrimeField> Mul for MultiLinearCoefficientPolynomial<F> {
    type Output = Self;
    fn mul(self, _rhs: Self) -> Self {
        MultiLinearCoefficientPolynomial { terms: vec![] }
    }
}

impl<F: PrimeField> Add for MultiLinearCoefficientPolynomial<F> {
    type Output = Self;

    fn add(self, _rhs: Self) -> Self::Output {
        MultiLinearCoefficientPolynomial { terms: vec![] }
    }
}

impl<F: PrimeField> Display for MultiLinearCoefficientPolynomial<F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        let mut first = true;

        for i in 0..self.terms.len() {
            if first {
                first = false;
            } else {
                write!(f, " + ")?;
            }

            write!(f, "{}", self.terms[i].coefficient)?;

            for (j, var) in self.terms[i].vars.iter().enumerate() {
                if *var {
                    write!(f, "{}", (b'a' + j as u8) as char)?;
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use field_tracker::Ft;

    use super::{MultiLinearCoefficientPolynomial, MultiLinearMonomial};
    use ark_test_curves::bls12_381::Fr as Fq_old;

    type Fq = Ft<4, Fq_old>;

    fn init() -> MultiLinearCoefficientPolynomial<Fq> {
        let term1 = MultiLinearMonomial::new(Fq::from(3), vec![false, false]); // 3
        let term2 = MultiLinearMonomial::new(Fq::from(-1), vec![true, false]); // a
        let term3 = MultiLinearMonomial::new(Fq::from(-2), vec![false, true]); // 2b
        let term4 = MultiLinearMonomial::new(Fq::from(5), vec![true, true]); // 5ab

        let polynomial = MultiLinearCoefficientPolynomial::new(vec![term1, term2, term3, term4]);

        polynomial
    }

    #[test]
    fn test_partial_evaluation() {
        let term1 = MultiLinearMonomial::new(Fq::from(3), vec![false, false]); // 3
        let term2 = MultiLinearMonomial::new(Fq::from(1), vec![true, false]); // a
        let term3 = MultiLinearMonomial::new(Fq::from(2), vec![false, true]); // 2b
        let term4 = MultiLinearMonomial::new(Fq::from(5), vec![true, true]); // 5ab
        let polynomial = MultiLinearCoefficientPolynomial::new(vec![term1, term2, term3, term4]);

        let result_term1 = MultiLinearMonomial::new(Fq::from(0), vec![false]);
        let result_term2 = MultiLinearMonomial::new(Fq::from(13), vec![false]);
        let result_polynomial =
            MultiLinearCoefficientPolynomial::new(vec![result_term1, result_term2]);

        let p_eval = polynomial.partial_evaluation(Fq::from(3_u8));
        assert_eq!(p_eval, result_polynomial);
        println!("{}", Fq::summary());
    }

    #[test]
    fn test_full_evaluation_1() {
        let polynomial = init();
        let evaluation = polynomial.evaluation(&vec![Fq::from(5), Fq::from(6)]);
        assert_eq!(evaluation, Fq::from(136));
        println!("{}", Fq::summary());
    }

    #[test]
    fn test_full_evaluation_2() {
        let polynomial = init();
        let evaluation = polynomial.evaluation(&vec![Fq::from(2), Fq::from(3)]);
        assert_eq!(evaluation, Fq::from(25));
        println!("{}", Fq::summary());
    }

    #[test]
    fn test_polynomial_degree() {
        let polynomial = init();
        let degree = polynomial.degree();
        assert_eq!(degree, 2);
        println!("{}", Fq::summary());
    }
}
