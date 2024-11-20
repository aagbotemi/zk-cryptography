use ark_ec::{pairing::Pairing, Group};
use ark_ff::PrimeField;
use std::fmt::{Debug, Formatter, Result};

use crate::{interface::TrustedSetupInterface, utils::generate_array_of_points};
use polynomial::utils::boolean_hypercube;

#[derive(Clone)]
pub struct TrustedSetup<P: Pairing> {
    pub powers_of_tau_in_g1: Vec<P::G1>,
    pub powers_of_tau_in_g2: Vec<P::G2>,
}

impl<P: Pairing> TrustedSetupInterface<P> for TrustedSetup<P> {
    fn setup<F: PrimeField>(eval_points: &[F]) -> Self {
        let powers_of_tau_in_g1 = Self::generate_powers_of_tau_in_g1(&eval_points);
        let powers_of_tau_in_g2: Vec<P::G2> = Self::generate_powers_of_tau_in_g2(&eval_points);

        TrustedSetup {
            powers_of_tau_in_g1,
            powers_of_tau_in_g2,
        }
    }

    fn generate_powers_of_tau_in_g1<F: PrimeField>(eval_points: &[F]) -> Vec<P::G1> {
        let g1 = P::G1::generator();

        let bh_cube = boolean_hypercube(eval_points.len());
        let array_of_points = generate_array_of_points(&bh_cube, &eval_points);

        array_of_points
            .iter()
            .map(|point| g1.mul_bigint(point.into_bigint()))
            .collect()
    }

    fn generate_powers_of_tau_in_g2<F: PrimeField>(eval_points: &[F]) -> Vec<P::G2> {
        let g2 = P::G2::generator();

        eval_points
            .iter()
            .map(|point| g2.mul_bigint(point.into_bigint()))
            .collect()
    }
}

impl<P: Pairing> Debug for TrustedSetup<P> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        f.debug_struct("TrustedSetup")
            .field("powers_of_tau_in_g1", &self.powers_of_tau_in_g1)
            .field("powers_of_tau_in_g2", &self.powers_of_tau_in_g2)
            .finish()
    }
}
