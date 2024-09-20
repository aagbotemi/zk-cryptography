use ark_ec::{pairing::Pairing, Group};
use ark_ff::{PrimeField, Zero};
use std::marker::PhantomData;

use polynomial::{Multilinear, MultilinearTrait};

use crate::{
    interface::{MultilinearKZGInterface, TrustedSetupInterface},
    trusted_setup::TrustedSetup,
    utils::{get_poly_quotient, get_poly_remainder, sum_pairing_results},
};

pub struct MultilinearKZG<P: Pairing> {
    _marker: PhantomData<P>,
}

#[derive(Debug)]
pub struct MultilinearKZGProof<P: Pairing> {
    pub evaluation: P::ScalarField,
    pub proofs: Vec<P::G1>,
}

impl<P: Pairing> Default for MultilinearKZGProof<P> {
    fn default() -> Self {
        MultilinearKZGProof {
            evaluation: Default::default(),
            proofs: Default::default(),
        }
    }
}

impl<P: Pairing> MultilinearKZGInterface<P> for MultilinearKZG<P> {
    fn commitment(poly: &Multilinear<P::ScalarField>, powers_of_tau_in_g1: &Vec<P::G1>) -> P::G1
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

    fn open(
        poly_: &Multilinear<P::ScalarField>,
        evaluation_points: &[P::ScalarField],
        powers_of_tau_in_g1: &Vec<P::G1>,
    ) -> MultilinearKZGProof<P> {
        let evaluation = poly_.evaluation(evaluation_points);

        let mut proofs = vec![];
        let mut poly = poly_.clone();
        let mut final_round_remainder = P::ScalarField::zero();

        for (variable_index, eval_point) in evaluation_points.iter().enumerate() {
            let mut remainder = Multilinear::additive_identity(variable_index);
            let quotient;
            let blown_poly;

            if variable_index != evaluation_points.len() - 1 {
                quotient = get_poly_quotient(&poly);
                remainder = get_poly_remainder(&poly, &eval_point);
                blown_poly = quotient.add_to_front(&(variable_index));
            } else {
                quotient = get_poly_quotient(&poly);
                final_round_remainder = poly.evaluation(&[*eval_point]);

                let duplicate_poly = Multilinear::duplicate_evaluation(&quotient.evaluations);
                blown_poly = duplicate_poly.add_to_front(&(variable_index - 1));
            }

            let proof = MultilinearKZG::<P>::commitment(&blown_poly, &powers_of_tau_in_g1);
            poly = remainder;
            proofs.push(proof);
        }

        if evaluation != final_round_remainder {
            panic!("Evaluation and final remainder mismatch!");
        }

        MultilinearKZGProof { evaluation, proofs }
    }

    fn verify(
        commit: &P::G1,
        verifier_points: &[P::ScalarField],
        proof: &MultilinearKZGProof<P>,
        powers_of_tau_in_g2: &Vec<P::G2>,
    ) -> bool {
        let g1 = P::G1::generator();
        let g2 = P::G2::generator();

        // LHS
        let v = g1.mul_bigint(proof.evaluation.into_bigint());
        let lhs = P::pairing(*commit - v, g2);

        // RHS
        let verifier_point_powers_of_tau_in_g2: Vec<P::G2> =
            TrustedSetup::<P>::generate_powers_of_tau_in_g2(verifier_points);
        let rhs = sum_pairing_results::<P>(
            powers_of_tau_in_g2,
            &verifier_point_powers_of_tau_in_g2,
            &proof.proofs,
        );

        lhs == rhs
    }
}

#[cfg(test)]
mod tests {
    use ark_test_curves::bls12_381::{Bls12_381, Fr};
    use polynomial::Multilinear;

    use super::MultilinearKZG;
    use crate::{
        interface::TrustedSetupInterface,
        kzg::{MultilinearKZGInterface, MultilinearKZGProof},
        trusted_setup::TrustedSetup,
    };

    #[test]
    fn test_kzg_1() {
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
        let tau = TrustedSetup::<Bls12_381>::setup(&prover_points);
        let commit = MultilinearKZG::<Bls12_381>::commitment(&poly, &tau.powers_of_tau_in_g1);

        let proof: MultilinearKZGProof<Bls12_381> =
            MultilinearKZG::open(&poly, &verifier_points, &tau.powers_of_tau_in_g1);
        let verify_status =
            MultilinearKZG::verify(&commit, &verifier_points, &proof, &tau.powers_of_tau_in_g2);

        assert_eq!(verify_status, true)
    }

    #[test]
    fn test_kzg_2() {
        let prover_points = vec![Fr::from(12), Fr::from(9), Fr::from(28), Fr::from(40)];
        let tampered_prover_points = vec![Fr::from(12), Fr::from(19), Fr::from(28), Fr::from(40)];
        let verifier_points = vec![Fr::from(54), Fr::from(90), Fr::from(76), Fr::from(160)];

        // 4ac + 10bc + 2cd - 12ad
        let value = vec![
            Fr::from(0),
            Fr::from(0),
            Fr::from(0),
            Fr::from(2),
            Fr::from(0),
            Fr::from(0),
            Fr::from(10),
            Fr::from(12),
            Fr::from(0),
            Fr::from(-12),
            Fr::from(4),
            Fr::from(-6),
            Fr::from(0),
            Fr::from(-12),
            Fr::from(14),
            Fr::from(4),
        ];

        let poly = Multilinear::new(value);
        let tau = TrustedSetup::<Bls12_381>::setup(&prover_points);
        let tampered_tau = TrustedSetup::<Bls12_381>::setup(&tampered_prover_points);
        let commit = MultilinearKZG::<Bls12_381>::commitment(&poly, &tau.powers_of_tau_in_g1);

        let proof: MultilinearKZGProof<Bls12_381> =
            MultilinearKZG::open(&poly, &verifier_points, &tau.powers_of_tau_in_g1);
        let verify_status =
            MultilinearKZG::verify(&commit, &verifier_points, &proof, &tau.powers_of_tau_in_g2);
        let tampered_tau_verify_status = MultilinearKZG::verify(
            &commit,
            &verifier_points,
            &proof,
            &tampered_tau.powers_of_tau_in_g2,
        );

        assert_eq!(verify_status, true);
        assert_eq!(tampered_tau_verify_status, false);
    }
}
