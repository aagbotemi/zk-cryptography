use crate::{
    trusted_setup::TrustedSetup,
    utils::{get_poly_quotient, get_poly_remainder, sum_pairing_results},
};

use ark_bls12_381::{Bls12_381, Config};
use ark_ec::{bls12::G1Projective, pairing::Pairing, Group};
use ark_ff::PrimeField;

use polynomial::{Multilinear, MultilinearTrait};

pub struct MultilinearKZG {}

#[derive(Debug)]
pub struct MultilinearKZGProof<F: PrimeField> {
    pub evaluation: F,
    pub proofs: Vec<G1Projective<Config>>,
}

impl MultilinearKZG {
    pub fn commitment<P: Pairing, F: PrimeField>(
        poly: &Multilinear<F>,
        powers_of_tau_in_g1: &Vec<P::G1>,
    ) -> P::G1
    where
        P: Pairing,
        F: PrimeField,
    {
        let evaluations: Vec<F> = poly.evaluations.clone();

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

    pub fn open<P: Pairing, F: PrimeField>(
        poly_: &Multilinear<F>,
        evaluation_points: &[F],
        tau: &TrustedSetup<Bls12_381>,
    ) -> MultilinearKZGProof<F> {
        let evaluation = poly_.evaluation(evaluation_points);
        let mut proofs = vec![];
        let mut poly = poly_.clone();
        let mut final_round_remainder: F = F::zero();

        for (variable_index, eval_point) in evaluation_points.iter().enumerate() {
            let mut remainder: Multilinear<F> = Multilinear::additive_identity(variable_index);
            let mut quotient: Multilinear<F> = Multilinear::additive_identity(variable_index);
            let mut blown_poly: Multilinear<F> = Multilinear::additive_identity(variable_index);

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

            let proof = Self::commitment::<Bls12_381, F>(&blown_poly, &tau.powers_of_tau_in_g1);
            poly = remainder;
            proofs.push(proof);
        }

        if evaluation != final_round_remainder {
            panic!("Evaluation and final remainder mismatch!");
        }

        MultilinearKZGProof { evaluation, proofs }
    }

    pub fn verify<P: Pairing, F: PrimeField>(
        commit: &P::G1,
        verifier_points: &[F],
        proof: &MultilinearKZGProof<F>,
        powers_of_tau_in_g2: Vec<P::G2>,
    ) -> bool {
        let g1 = P::G1::generator();
        let g2 = P::G2::generator();

        // lhs
        let v = g1.mul_bigint(proof.evaluation.into_bigint());
        let lhs = P::pairing(*commit - v, g2);

        // RHS
        let verifier_point_powers_of_tau_in_g2: Vec<P::G2> =
            TrustedSetup::generate_powers_of_tau_in_g2::<F>(verifier_points);

        let g1_power_of_proof = proof
            .proofs
            .iter()
            .map(|point| g1.mul_bigint(point)) // @TODO: fix here
            .collect();

        let rhs = sum_pairing_results(
            verifier_point_powers_of_tau_in_g2,
            powers_of_tau_in_g2,
            g1_power_of_proof,
        );

        lhs == rhs
    }
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::{Bls12_381, Fr};
    use polynomial::Multilinear;

    use super::MultilinearKZG;
    use crate::trusted_setup::TrustedSetup;

    fn test_kzg() {
        let prover_points = vec![Fr::from(2), Fr::from(3_u8), Fr::from(4_u8)];
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
        // let (powers_of_tau_in_g1, powers_of_tau_in_g2) =
        //     TrustedSetup::<Bls12_381>::setup(&prover_points);

        let commit = MultilinearKZG::commitment::<Bls12_381, Fr>(&poly, &tau.powers_of_tau_in_g1);

        let proof = MultilinearKZG::open::<Bls12_381, Fr>(&poly, &verifier_points, &tau);

        let verify_status = MultilinearKZG::verify::<Bls12_381, Fr>(
            &commit,
            &verifier_points,
            &proof,
            tau.powers_of_tau_in_g2,
        );

        assert!(verify_status)
    }
}
