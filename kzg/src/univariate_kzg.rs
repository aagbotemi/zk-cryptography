use crate::{interface::UnivariateKZGInterface, trusted_setup::TrustedSetup};
use ark_ec::{pairing::Pairing, Group};
use ark_ff::{Field, PrimeField};
use polynomial::{DenseUnivariatePolynomial, UnivariatePolynomialTrait};
use std::marker::PhantomData;

pub struct UnivariateKZG<P: Pairing> {
    _marker: PhantomData<P>,
}

#[derive(Debug)]
pub struct UnivariateKZGProof<F: PrimeField, P: Pairing> {
    pub evaluation: F,
    pub proof: P::G1,
}

impl<P: Pairing> UnivariateKZGInterface<P> for UnivariateKZG<P> {
    fn generate_srs(tau: &P::ScalarField, max_degree: &usize) -> TrustedSetup<P> {
        let g1 = P::G1::generator();
        let g2 = P::G2::generator();

        let mut powers_of_tau_in_g1 = vec![];
        let mut powers_of_tau_in_g2 = vec![];

        for i in 0..=*max_degree {
            let power_of_tau = tau.pow([i as u64]);
            powers_of_tau_in_g1.push(g1.mul_bigint(power_of_tau.into_bigint()));
            powers_of_tau_in_g2.push(g2.mul_bigint(power_of_tau.into_bigint()));
        }

        TrustedSetup {
            powers_of_tau_in_g1,
            powers_of_tau_in_g2,
        }
    }

    fn commitment(
        poly: &DenseUnivariatePolynomial<P::ScalarField>,
        srs: &TrustedSetup<P>,
    ) -> P::G1 {
        let coefficients: Vec<P::ScalarField> = poly.coefficients.clone();

        assert_eq!(
            srs.powers_of_tau_in_g1.len(),
            coefficients.len(),
            "The length of powers_of_tau_in_g1 and the length of
                the evaluations of the polynomial should tally!"
        );

        let mut commit = P::G1::default();

        for i in 0..coefficients.len() {
            let g1_pow_coeff = srs.powers_of_tau_in_g1[i].mul_bigint(coefficients[i].into_bigint());
            commit += g1_pow_coeff
        }

        commit
    }

    fn open<F: PrimeField>(
        poly_: &DenseUnivariatePolynomial<F>,
        evaluation_points: F,
        srs: &TrustedSetup<P>,
    ) -> UnivariateKZGProof<F, P> {
        let evaluation = poly_.evaluate(evaluation_points);

        let denominator = DenseUnivariatePolynomial::new(vec![-evaluation_points, F::ONE]);

        let numerator = poly_ - evaluation_points;
        let quotient = numerator / denominator;

        let mut proof = P::G1::default();

        for i in 0..quotient.coefficients.len() {
            let g1_pow_coeff =
                srs.powers_of_tau_in_g1[i].mul_bigint(quotient.coefficients[i].into_bigint());
            proof += g1_pow_coeff;
        }

        UnivariateKZGProof { evaluation, proof }
    }

    fn verify<F: PrimeField>(
        commit: &P::G1,
        verifier_point: &P::ScalarField,
        proof: &UnivariateKZGProof<F, P>,
        srs: &TrustedSetup<P>,
    ) -> bool {
        let g1 = P::G1::generator();
        let g2 = P::G2::generator();

        // LHS
        let v = g1.mul_bigint(proof.evaluation.into_bigint());
        let lhs = P::pairing(*commit - v, g2);

        // RHS
        let g2_point = g2.mul_bigint(verifier_point.into_bigint());
        let rhs = P::pairing(proof.proof, &(srs.powers_of_tau_in_g2[1] - g2_point));

        lhs == rhs
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use ark_test_curves::bls12_381::{Bls12_381, Fr};

    #[test]
    fn test_univariate_kzg() {
        let tau = Fr::from(10u64);
        let max_degree = 4 as usize;
        let srs: TrustedSetup<Bls12_381> = UnivariateKZG::generate_srs(&tau, &max_degree);

        let poly = DenseUnivariatePolynomial::new(vec![
            Fr::from(1u64),
            Fr::from(2u64),
            Fr::from(3u64),
            Fr::from(4u64),
            Fr::from(5u64),
        ]);
        let commitment = UnivariateKZG::commitment(&poly, &srs);
        let proof = UnivariateKZG::open(&poly, Fr::from(2u64), &srs);

        let is_valid = UnivariateKZG::verify(&commitment, &Fr::from(2u64), &proof, &srs);

        assert!(is_valid);
    }

    #[test]
    fn test_univariate_kzg_invalid_opening() {
        let tau = Fr::from(10u64);
        let max_degree = 4 as usize;
        let srs: TrustedSetup<Bls12_381> = UnivariateKZG::generate_srs(&tau, &max_degree);

        let poly = DenseUnivariatePolynomial::new(vec![
            Fr::from(1u64),
            Fr::from(2u64),
            Fr::from(3u64),
            Fr::from(4u64),
            Fr::from(5u64),
        ]);
        let commitment = UnivariateKZG::commitment(&poly, &srs);
        let proof = UnivariateKZG::open(&poly, Fr::from(2u64), &srs);
        let is_valid = UnivariateKZG::verify(&commitment, &Fr::from(4u64), &proof, &srs);

        assert_eq!(is_valid, false);
    }
}
