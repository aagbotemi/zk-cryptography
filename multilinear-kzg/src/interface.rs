use ark_ec::pairing::Pairing;
use polynomial::Multilinear;

use crate::kzg::MultilinearKZGProof;

pub trait MultilinearKZGInterface<P: Pairing> {
    fn commitment(poly: &Multilinear<P::ScalarField>, powers_of_tau_in_g1: &Vec<P::G1>) -> P::G1
    where
        P: Pairing;

    fn open(
        poly_: &Multilinear<P::ScalarField>,
        evaluation_points: &[P::ScalarField],
        powers_of_tau_in_g1: &Vec<P::G1>,
    ) -> MultilinearKZGProof<P>;

    fn verify(
        commit: &P::G1,
        verifier_points: &[P::ScalarField],
        proof: &MultilinearKZGProof<P>,
        powers_of_tau_in_g2: Vec<P::G2>,
    ) -> bool;
}

pub trait TrustedSetupInterface<P: Pairing> {
    fn setup(eval_points: &[P::ScalarField]) -> Self
    where
        P: Pairing;

    fn generate_powers_of_tau_in_g1(eval_points: &[P::ScalarField]) -> Vec<P::G1>;

    fn generate_powers_of_tau_in_g2(eval_points: &[P::ScalarField]) -> Vec<P::G2>;
}
