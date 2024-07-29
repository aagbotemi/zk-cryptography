use ark_ff::PrimeField;
use fiat_shamir::{fiat_shamir::FiatShamirTranscript, interface::FiatShamirTranscriptTrait};
use polynomial::UnivariatePolynomial;
use sumcheck::sumcheck::SumcheckProof;

use crate::circuit::{GKRCircuit, GKRCircuitEvaluation};

pub struct GKRProof<'a, F: PrimeField> {
    /// This is the output of the Circuit evaluation
    pub output: Vec<F>,
    /// This is the list of sum check proofs gotten during this protocol
    pub sumcheck_proofs: Vec<SumcheckProof<'a, F>>,
    /// This is the list of q polynomials
    pub q_polynomials: Vec<UnivariatePolynomial<F>>,
}

/// Prove correct circuit evaluation using the GKR protocol
pub fn prove<F: PrimeField>(
    circuit: &GKRCircuit,
    evaluations: &GKRCircuitEvaluation<F>,
) -> GKRProof<'static, F> {
    let mut transcript = FiatShamirTranscript::new();
    let mut sumcheck_proofs: Vec<SumcheckProof<F>> = vec![];
    let mut q_polynomials: Vec<UnivariatePolynomial<F>> = vec![];

    unimplemented!()
}

fn verify<F: PrimeField>(circuit: &GKRCircuit, input: &[F], proof: &GKRProof<F>) -> bool {
    unimplemented!()
}
