use std::marker::PhantomData;

use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use merlin::MerlinTranscript;
use polynomial::DenseUnivariatePolynomial;

use super::primitives::PlonkRoundTranscript;

impl<P: Pairing> PlonkRoundTranscript<P> {
    pub fn new() -> Self {
        let transcript = MerlinTranscript::new(b"plonk_protocol");

        Self {
            transcript,
            _marker: PhantomData,
        }
    }

    pub fn first_round(&mut self, a_s: P::G1, b_s: P::G1, c_s: P::G1) {
        self.transcript.append_point::<P>(b"first_round", &a_s);
        self.transcript.append_point::<P>(b"first_round", &b_s);
        self.transcript.append_point::<P>(b"first_round", &c_s);
    }

    pub fn second_round<F: PrimeField>(
        &mut self,
        zh_blinding_accumulator_poly: DenseUnivariatePolynomial<F>,
        accumulator_commitment: P::G1,
    ) {
        self.transcript
            .append_point::<P>(b"second_round", &accumulator_commitment);

        let poly_bytes = zh_blinding_accumulator_poly.to_bytes();
        self.transcript.append_message(b"second_round", &poly_bytes);
    }

    pub fn third_round(&mut self, t_low: P::G1, t_mid: P::G1, t_high: P::G1) {
        self.transcript.append_point::<P>(b"third_round", &t_low);
        self.transcript.append_point::<P>(b"third_round", &t_mid);
        self.transcript.append_point::<P>(b"third_round", &t_high);
    }

    pub fn fourth_round<F: PrimeField>(
        &mut self,
        a_s_poly_zeta: F,
        b_s_poly_zeta: F,
        c_s_poly_zeta: F,
        sigma1_poly_zeta: F,
        sigma2_poly_zeta: F,
        w_accumulator_poly_zeta: F,
    ) {
        self.transcript
            .append_scalar::<F>(b"fourth_round", &a_s_poly_zeta);
        self.transcript
            .append_scalar::<F>(b"fourth_round", &b_s_poly_zeta);
        self.transcript
            .append_scalar::<F>(b"fourth_round", &c_s_poly_zeta);
        self.transcript
            .append_scalar::<F>(b"fourth_round", &sigma1_poly_zeta);
        self.transcript
            .append_scalar::<F>(b"fourth_round", &sigma2_poly_zeta);
        self.transcript
            .append_scalar::<F>(b"fourth_round", &w_accumulator_poly_zeta);
    }

    pub fn fifth_round(&mut self, w_zeta_commitment: P::G1, w_zeta_omega_commitment: P::G1) {
        self.transcript
            .append_point::<P>(b"fifth_round", &w_zeta_commitment);
        self.transcript
            .append_point::<P>(b"fifth_round", &w_zeta_omega_commitment);
    }

    pub fn challenge_n_round<F: PrimeField>(&mut self, label: &[u8], n: usize) -> Vec<F> {
        self.transcript.challenge_n::<F>(label, n)
    }

    pub fn challenge_round<F: PrimeField>(&mut self, label: &[u8]) -> F {
        self.transcript.challenge::<F>(label)
    }
}
