use ark_ec::{pairing::Pairing, Group};
use ark_ff::PrimeField;
use polynomial::Multilinear;

pub struct MultilinearKZG {}

impl MultilinearKZG {
    pub fn commitment<P: Pairing, F: PrimeField>(
        poly: &Multilinear<F>,
        powers_of_tau_in_g1: &Vec<P::G1>,
    ) -> P::G1
    where
        P: Pairing,
        F: PrimeField,
    {
        let poly: Vec<F> = poly.evaluations.clone();

        assert_eq!(
            powers_of_tau_in_g1.len(),
            poly.len(),
            "The length of powers_of_tau_in_g1 and the length of 
            the evaluations of the polynomial should tally!"
        );

        poly.iter()
            .zip(powers_of_tau_in_g1.iter())
            .map(|(coefficient, power)| power.mul_bigint(coefficient.into_bigint()))
            .sum()
    }
}
