use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use polynomial::{DenseUnivariatePolynomial, Multilinear};

use crate::{
    multilinear_kzg::MultilinearKZGProof, trusted_setup::TrustedSetup,
    univariate_kzg::UnivariateKZGProof,
};

pub trait MultilinearKZGInterface<F: PrimeField, P: Pairing> {
    fn commitment(poly: &Multilinear<F>, srs: &TrustedSetup<P>) -> P::G1;

    fn open(
        poly_: &Multilinear<F>,
        evaluation_points: &[F],
        srs: &TrustedSetup<P>,
    ) -> MultilinearKZGProof<F, P>;

    fn verify(
        commit: &P::G1,
        verifier_points: &[F],
        proof: &MultilinearKZGProof<F, P>,
        srs: &TrustedSetup<P>,
    ) -> bool;
}

pub trait TrustedSetupInterface<P: Pairing> {
    fn setup<F: PrimeField>(eval_points: &[F]) -> Self;

    fn generate_powers_of_tau_in_g1<F: PrimeField>(eval_points: &[F]) -> Vec<P::G1>;

    fn generate_powers_of_tau_in_g2<F: PrimeField>(eval_points: &[F]) -> Vec<P::G2>;
}

pub trait UnivariateKZGInterface<P: Pairing> {
    fn generate_srs(tau: &P::ScalarField, max_degree: &usize) -> TrustedSetup<P>;

    fn commitment(poly: &DenseUnivariatePolynomial<P::ScalarField>, srs: &TrustedSetup<P>)
        -> P::G1;

    fn open<F: PrimeField>(
        poly_: &DenseUnivariatePolynomial<F>,
        evaluation_points: F,
        srs: &TrustedSetup<P>,
    ) -> UnivariateKZGProof<F, P>;

    fn verify<F: PrimeField>(
        commit: &P::G1,
        verifier_point: &P::ScalarField,
        proof: &UnivariateKZGProof<F, P>,
        srs: &TrustedSetup<P>,
    ) -> bool;
}
