use crate::{
    trusted_setup::TrustedSetup,
    utils::{get_poly_quotient, get_poly_remainder, sum_pairing_results},
};

use ark_bls12_381::{Bls12_381, Config};
use ark_ec::{bls12::G1Projective, pairing::Pairing, Group};
use ark_ff::{PrimeField, Zero};

use polynomial::{Multilinear, MultilinearTrait};

pub struct MultilinearKZG {}

#[derive(Debug)]
pub struct MultilinearKZGProof<P: Pairing> {
    pub evaluation: P::ScalarField,
    pub proofs: Vec<P::G1>,
}

impl MultilinearKZG {
    fn commitment<P: Pairing>(
        poly: &Multilinear<P::ScalarField>,
        powers_of_tau_in_g1: &Vec<P::G1>,
    ) -> P::G1
    where
        P: Pairing,
    {
        let evaluations: Vec<P::ScalarField> = poly.evaluations.clone();

        assert_eq!(
            powers_of_tau_in_g1.len(),
            evaluations.len(),
            "The length of powers_of_tau_in_g1 and the length of
            the evaluations of the polynomial should tally!"
        );

        evaluations
            .iter()
            .zip(powers_of_tau_in_g1.iter())
            .map(|(coefficient, power)| power.mul_bigint(coefficient.into_bigint()))
            .sum()
    }

    fn open<P: Pairing>(
        poly_: &Multilinear<P::ScalarField>,
        evaluation_points: &[P::ScalarField],
        powers_of_tau_in_g1: &Vec<P::G1>,
    ) -> MultilinearKZGProof<P> {
        let evaluation = poly_.evaluation(evaluation_points);

        let mut proofs = vec![];
        let mut poly = poly_.clone();
        let mut final_round_remainder = P::ScalarField::zero();

        for (variable_index, eval_point) in evaluation_points.iter().enumerate() {
            dbg!(&variable_index);
            let mut remainder = Multilinear::additive_identity(variable_index);
            let mut quotient = Multilinear::additive_identity(variable_index);
            let mut blown_poly = Multilinear::additive_identity(variable_index);

            if variable_index != evaluation_points.len() - 1 {
                quotient = get_poly_quotient(&poly);
                dbg!(&quotient);
                remainder = get_poly_remainder(&poly, &eval_point);
                dbg!(&remainder);
                blown_poly = quotient.add_to_front(&(variable_index));
                dbg!(&blown_poly);
            } else {
                quotient = get_poly_quotient(&poly);
                dbg!(&quotient);
                final_round_remainder = poly.evaluation(&[*eval_point]);
                dbg!(&final_round_remainder);

                let duplicate_poly = Multilinear::duplicate_evaluation(&quotient.evaluations);
                dbg!(&duplicate_poly);
                blown_poly = duplicate_poly.add_to_front(&(variable_index - 1));
                dbg!(&blown_poly);
            }

            let proof = MultilinearKZG::commitment::<P>(&blown_poly, &powers_of_tau_in_g1);
            dbg!(&proof);
            poly = remainder;
            dbg!(&poly);
            proofs.push(proof);
            dbg!(&proofs);
        }

        if evaluation != final_round_remainder {
            panic!("Evaluation and final remainder mismatch!");
        }
        dbg!(&evaluation, final_round_remainder);

        MultilinearKZGProof { evaluation, proofs }
    }

    pub fn verify<P: Pairing>(
        commit: &P::G1,
        verifier_points: &[P::ScalarField],
        proof: &MultilinearKZGProof<P>,
        powers_of_tau_in_g2: Vec<P::G2>,
    ) -> bool {
        let g1 = P::G1::generator();
        let g2 = P::G2::generator();
        dbg!(&g1, &g2);

        // LHS
        let v = g1.mul_bigint(proof.evaluation.into_bigint());
        dbg!(&v);
        let lhs = P::pairing(*commit - v, g2);

        // RHS
        let verifier_point_powers_of_tau_in_g2: Vec<P::G2> =
            TrustedSetup::generate_powers_of_tau_in_g2::<P>(verifier_points);
        dbg!(&verifier_point_powers_of_tau_in_g2);

        let rhs = sum_pairing_results::<P>(
            verifier_point_powers_of_tau_in_g2,
            powers_of_tau_in_g2,
            proof.proofs.clone(),
        );

        dbg!(&lhs, &rhs);

        lhs == rhs
    }
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::{Bls12_381, Fr};
    use polynomial::Multilinear;

    use super::MultilinearKZG;
    use crate::{kzg::MultilinearKZGProof, trusted_setup::TrustedSetup};

    #[test]
    fn test_kzg() {
        let prover_points = vec![Fr::from(2), Fr::from(3), Fr::from(4)];
        let verifier_points = vec![Fr::from(5), Fr::from(9), Fr::from(6)];

        let val = vec![
            Fr::from(0),
            Fr::from(7),
            Fr::from(0),
            Fr::from(5),
            Fr::from(0),
            Fr::from(7),
            Fr::from(4),
            Fr::from(9),
        ];
        let poly = Multilinear::new(val);

        let (powers_of_tau_in_g1, powers_of_tau_in_g2) =
            TrustedSetup::setup::<Bls12_381>(&prover_points);

        let commit = MultilinearKZG::commitment::<Bls12_381>(&poly, &powers_of_tau_in_g1);

        let proof: MultilinearKZGProof<Bls12_381> =
            MultilinearKZG::open(&poly, &verifier_points, &powers_of_tau_in_g1);

        let verify_status =
            MultilinearKZG::verify(&commit, &verifier_points, &proof, powers_of_tau_in_g2);

        assert_eq!(verify_status, true)
    }
}
