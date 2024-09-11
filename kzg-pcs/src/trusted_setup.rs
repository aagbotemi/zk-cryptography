use ark_bls12_381::Bls12_381;
use ark_ec::{pairing::Pairing, Group};
use ark_ff::PrimeField;
use std::fmt::{Debug, Formatter, Result};

use crate::utils::generate_array_of_points;
use polynomial::utils::boolean_hypercube;

pub struct TrustedSetup<P: Pairing> {
    pub powers_of_tau_in_g1: Vec<P::G1>,
    pub powers_of_tau_in_g2: Vec<P::G2>,
}

impl<P: Pairing> TrustedSetup<P> {
    pub fn setup<F: PrimeField>(eval_points: &[F]) -> Self
    where
        F: PrimeField,
        P: Pairing,
    {
        let powers_of_tau_in_g1: Vec<P::G1> = Self::generate_powers_of_tau_in_g1(&eval_points);
        let powers_of_tau_in_g2: Vec<P::G2> = Self::generate_powers_of_tau_in_g2(&eval_points);

        TrustedSetup {
            powers_of_tau_in_g1,
            powers_of_tau_in_g2,
        }
    }

    #[inline]
    fn generate_powers_of_tau_in_g1<F: PrimeField>(eval_points: &[F]) -> Vec<P::G1> {
        let g1 = P::G1::generator();

        let bh_cube: Vec<Vec<F>> = boolean_hypercube(eval_points.len());
        let array_of_points: Vec<F> = generate_array_of_points(&bh_cube, &eval_points);

        array_of_points
            .iter()
            .map(|point| g1.mul_bigint(point.into_bigint()))
            .collect()
    }

    #[inline]
    fn generate_powers_of_tau_in_g2<F: PrimeField>(eval_points: &[F]) -> Vec<P::G2> {
        let g2 = P::G2::generator();

        eval_points
            .iter()
            .map(|point| g2.mul_bigint(point.into_bigint()))
            .collect()
    }
}

impl Debug for TrustedSetup<Bls12_381> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        f.debug_struct("TrustedSetup")
            .field("powers_of_tau_in_g1", &self.powers_of_tau_in_g1)
            .field("powers_of_tau_in_g2", &self.powers_of_tau_in_g2)
            .finish()
    }
}
