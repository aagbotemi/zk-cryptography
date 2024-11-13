use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use polynomial::{DenseUnivariatePolynomial, UnivariatePolynomialTrait};

use super::primitives::{PlonkProof, PlonkRoundTranscript};

pub fn split_poly_in_3<F: PrimeField>(
    poly: &DenseUnivariatePolynomial<F>,
    group_order: usize,
) -> (
    DenseUnivariatePolynomial<F>,
    DenseUnivariatePolynomial<F>,
    DenseUnivariatePolynomial<F>,
) {
    let poly_low_coeffs = poly.coefficients[0..group_order].to_vec();
    let poly_mid_coeffs = poly.coefficients[group_order..2 * group_order].to_vec();
    let poly_high_coeffs = poly.coefficients[2 * group_order..].to_vec();

    (
        DenseUnivariatePolynomial::new(poly_low_coeffs),
        DenseUnivariatePolynomial::new(poly_mid_coeffs),
        DenseUnivariatePolynomial::new(poly_high_coeffs),
    )
}

pub fn apply_w_to_polynomial<F: PrimeField>(
    poly: &DenseUnivariatePolynomial<F>,
    w: &F,
) -> DenseUnivariatePolynomial<F> {
    let mut result = Vec::new();
    let mut w_power = F::one();

    for coeff in poly.coefficients.iter() {
        result.push(*coeff * w_power);
        w_power *= w;
    }

    DenseUnivariatePolynomial::new(result)
}

pub fn zh_values<F: PrimeField>(group_order: usize) -> Vec<F> {
    let mut zh_values = vec![F::one().neg()];
    for _ in 0..group_order - 1 {
        zh_values.push(F::zero());
    }
    zh_values.push(F::one());
    zh_values
}

pub fn l1_values<F: PrimeField>(group_order: u64) -> Vec<F> {
    let mut l1_values = vec![F::one()];
    for _ in 0..group_order - 1 {
        l1_values.push(F::zero());
    }
    l1_values
}

pub fn compute_verifier_challenges<P: Pairing, F: PrimeField>(
    proof: &PlonkProof<P, F>,
) -> (F, F, F, F, F, F) {
    let mut transcript: PlonkRoundTranscript<P> = PlonkRoundTranscript::new();

    // beta and gamma
    let _ = transcript.first_round(proof.a_s, proof.b_s, proof.c_s);
    let challenge: Vec<F> = transcript.challenge_n_round(b"beta_gamma", 2);
    let beta = challenge[0];
    let gamma = challenge[1];

    // alpha
    let _ = transcript.second_round::<F>(proof.accumulator_commitment);
    let alpha: F = transcript.challenge_round(b"alpha");

    // zeta
    let _ = transcript.third_round(proof.t_low, proof.t_mid, proof.t_high);
    let zeta: F = transcript.challenge_round(b"zeta");

    // nu
    let _ = transcript.fourth_round(
        proof.a_s_poly_zeta,
        proof.b_s_poly_zeta,
        proof.c_s_poly_zeta,
        proof.sigma1_poly_zeta,
        proof.sigma2_poly_zeta,
        proof.w_accumulator_poly_zeta,
    );
    let nu: F = transcript.challenge_round(b"nu");

    // mu
    let _ = transcript.fifth_round(proof.w_zeta_commitment, proof.w_zeta_omega_commitment);
    let mu: F = transcript.challenge_round(b"mu");

    (beta, gamma, alpha, zeta, nu, mu)
}
