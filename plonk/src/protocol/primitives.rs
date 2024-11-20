use crate::compiler::primitives::CommonPreprocessedInput;
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use kzg::trusted_setup::TrustedSetup;
use merlin::MerlinTranscript;
use polynomial::DenseUnivariatePolynomial;
use std::marker::PhantomData;

pub struct PlonkProver<F: PrimeField, P: Pairing> {
    pub preprocessed_input: CommonPreprocessedInput<F>,
    pub srs: TrustedSetup<P>,
    pub transcript: PlonkRoundTranscript<P>,
    pub random_number: RandomNumbers<F>,
    pub witness_polys: WitnessPolys<F>,
}

pub struct RandomNumbers<F: PrimeField> {
    pub alpha: F,
    pub beta: F,
    pub gamma: F,
    pub zeta: F,
    pub nu: F,
    pub mu: F,
}

pub struct WitnessPolys<F: PrimeField> {
    pub a_s: DenseUnivariatePolynomial<F>,
    pub b_s: DenseUnivariatePolynomial<F>,
    pub c_s: DenseUnivariatePolynomial<F>,
    pub zh_poly: DenseUnivariatePolynomial<F>,

    pub zh_accumulator_poly: DenseUnivariatePolynomial<F>,
    pub t_low_poly: DenseUnivariatePolynomial<F>,
    pub t_mid_poly: DenseUnivariatePolynomial<F>,
    pub t_high_poly: DenseUnivariatePolynomial<F>,

    pub a_s_poly_zeta: F,
    pub b_s_poly_zeta: F,
    pub c_s_poly_zeta: F,
    pub w_accumulator_poly_zeta: F,
    pub sigma1_poly_zeta: F,
    pub sigma2_poly_zeta: F,

    pub w_zeta_poly: DenseUnivariatePolynomial<F>,
    pub w_zeta_omega_poly: DenseUnivariatePolynomial<F>,
}

pub struct PlonkProof<P: Pairing, F: PrimeField> {
    pub as_commitment: P::G1,
    pub bs_commitment: P::G1,
    pub cs_commitment: P::G1,
    pub accumulator_commitment: P::G1,
    pub t_low: P::G1,
    pub t_mid: P::G1,
    pub t_high: P::G1,
    pub a_s_poly_zeta: F,
    pub b_s_poly_zeta: F,
    pub c_s_poly_zeta: F,
    pub sigma1_poly_zeta: F,
    pub sigma2_poly_zeta: F,
    pub w_accumulator_poly_zeta: F,
    pub w_zeta_commitment: P::G1,
    pub w_zeta_omega_commitment: P::G1,
}

pub struct PlonkRoundTranscript<P: Pairing> {
    pub transcript: MerlinTranscript,
    pub _marker: PhantomData<P>,
}

impl<F: PrimeField> Default for RandomNumbers<F> {
    fn default() -> Self {
        Self {
            alpha: F::ZERO,
            beta: F::ZERO,
            gamma: F::ZERO,
            zeta: F::ZERO,
            nu: F::ZERO,
            mu: F::ZERO,
        }
    }
}

impl<F: PrimeField> Default for WitnessPolys<F> {
    fn default() -> Self {
        Self {
            a_s: DenseUnivariatePolynomial::zero(),
            b_s: DenseUnivariatePolynomial::zero(),
            c_s: DenseUnivariatePolynomial::zero(),
            zh_poly: DenseUnivariatePolynomial::zero(),

            zh_accumulator_poly: DenseUnivariatePolynomial::zero(),
            t_low_poly: DenseUnivariatePolynomial::zero(),
            t_mid_poly: DenseUnivariatePolynomial::zero(),
            t_high_poly: DenseUnivariatePolynomial::zero(),

            a_s_poly_zeta: F::ZERO,
            b_s_poly_zeta: F::ZERO,
            c_s_poly_zeta: F::ZERO,
            w_accumulator_poly_zeta: F::ZERO,
            sigma1_poly_zeta: F::ZERO,
            sigma2_poly_zeta: F::ZERO,

            w_zeta_poly: DenseUnivariatePolynomial::zero(),
            w_zeta_omega_poly: DenseUnivariatePolynomial::zero(),
        }
    }
}
