use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use polynomial::Multilinear;

use crate::{kzg::MultilinearKZGProof, trusted_setup::TrustedSetup};

pub trait MultilinearKZGInterface<F: PrimeField, P: Pairing> {
    fn commitment(poly: &Multilinear<F>, tau: &TrustedSetup<P>) -> P::G1;

    fn open(
        poly_: &Multilinear<F>,
        evaluation_points: &[F],
        tau: &TrustedSetup<P>,
    ) -> MultilinearKZGProof<F, P>;

    fn verify(
        commit: &P::G1,
        verifier_points: &[F],
        proof: &MultilinearKZGProof<F, P>,
        tau: &TrustedSetup<P>,
    ) -> bool;
}

pub trait TrustedSetupInterface<P: Pairing> {
    fn setup<F: PrimeField>(eval_points: &[F]) -> Self;

    fn generate_powers_of_tau_in_g1<F: PrimeField>(eval_points: &[F]) -> Vec<P::G1>;

    fn generate_powers_of_tau_in_g2<F: PrimeField>(eval_points: &[F]) -> Vec<P::G2>;
}
