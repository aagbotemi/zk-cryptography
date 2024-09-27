use std::marker::PhantomData;

use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use circuit::circuit::Circuit;
use fiat_shamir::{fiat_shamir::FiatShamirTranscript, interface::FiatShamirTranscriptTrait};
use multilinear_kzg::{
    interface::MultilinearKZGInterface,
    kzg::{MultilinearKZG, MultilinearKZGProof},
    trusted_setup::TrustedSetup,
};
use polynomial::{ComposedMultilinear, Multilinear, MultilinearTrait};
use sumcheck::composed::multi_composed_sumcheck::{
    ComposedSumcheckProof, MultiComposedSumcheckProver, MultiComposedSumcheckVerifier,
};

use crate::utils::{
    exponent, generate_layer_one_prove_sumcheck, generate_layer_one_verify_sumcheck, w_mle,
};

#[derive(Debug)]
pub struct SuccintGKRProof<F: PrimeField, P: Pairing> {
    sumcheck_proofs: Vec<ComposedSumcheckProof<F>>,
    wb_s: Vec<F>,
    wc_s: Vec<F>,
    w_0_mle: Multilinear<F>,
    proof_wb_opening: MultilinearKZGProof<F, P>,
    proof_wc_opening: MultilinearKZGProof<F, P>,
}

pub struct SuccintGKRProtocol<F: PrimeField, P: Pairing> {
    _marker: PhantomData<(F, P)>,
}

impl<F: PrimeField, P: Pairing> SuccintGKRProtocol<F, P> {
    /// Prove correct circuit evaluation using the GKR protocol
    pub fn prove(
        circuit: &Circuit,
        circuit_evaluation: &Vec<Vec<F>>,
        tau: &TrustedSetup<P>,
    ) -> (P::G1, SuccintGKRProof<F, P>) {
        let mut transcript = FiatShamirTranscript::new();
        let mut sumcheck_proofs: Vec<ComposedSumcheckProof<F>> = Vec::new();
        let mut wb_s: Vec<F> = Vec::new();
        let mut wc_s: Vec<F> = Vec::new();

        let mut circuit_evaluation_layer_zero_pad = circuit_evaluation[0].clone();
        circuit_evaluation_layer_zero_pad.push(F::zero());

        let w_0_mle = w_mle::<F>(circuit_evaluation_layer_zero_pad.to_vec());
        transcript.commit(&w_0_mle.to_bytes());

        let n_r = transcript.evaluate_n_challenge_into_field(&w_0_mle.n_vars);
        let mut claimed_sum = w_0_mle.evaluation(&n_r);

        let (add_mle_1, mult_mle_1) = circuit.add_mult_mle::<F>(0);
        let w_1_mle = w_mle::<F>(circuit_evaluation[1].to_vec());

        let (claimed, alph, bta, rb, rc) = generate_layer_one_prove_sumcheck(
            &add_mle_1,
            &mult_mle_1,
            &w_1_mle,
            &n_r,
            &claimed_sum,
            &mut transcript,
            &mut sumcheck_proofs,
            &mut wb_s,
            &mut wc_s,
        );

        claimed_sum = claimed;

        let mut alpha: F = alph;
        let mut beta: F = bta;

        let mut r_b: Vec<F> = rb;
        let mut r_c: Vec<F> = rc;

        let mut commitment: P::G1 = Default::default();
        let mut proof_wb_opening: MultilinearKZGProof<F, P> = Default::default();
        let mut proof_wc_opening: MultilinearKZGProof<F, P> = Default::default();

        for layer_index in 2..circuit_evaluation.len() {
            let (add_mle, mult_mle) = circuit.add_mult_mle(layer_index - 1);

            let add_rb_bc = add_mle.partial_evaluations(&r_b, &vec![0; r_b.len()]);
            let mul_rb_bc = mult_mle.partial_evaluations(&r_b, &vec![0; r_b.len()]);

            let add_rc_bc = add_mle.partial_evaluations(&r_c, &vec![0; r_b.len()]);
            let mul_rc_bc = mult_mle.partial_evaluations(&r_c, &vec![0; r_b.len()]);
            let w_i_mle = w_mle::<F>(circuit_evaluation[layer_index].to_vec());

            let wb = w_i_mle.clone();
            let wc = w_i_mle.clone();

            let wb_add_wc = wb.add_distinct(&wc);
            let wb_mul_wc = wb.mul_distinct(&wc);

            // alpha * add(r_b, b, c) + beta * add(r_c, b, c)
            let add_alpha_beta = (add_rb_bc * alpha) + (add_rc_bc * beta);
            // alpha * mul(r_b, b, c) + beta * mult(r_c, b, c)
            let mul_alpha_beta = (mul_rb_bc * alpha) + (mul_rc_bc * beta);

            let fbc_add_alpha_beta = ComposedMultilinear::new(vec![add_alpha_beta, wb_add_wc]);
            let fbc_mul_alpha_beta = ComposedMultilinear::new(vec![mul_alpha_beta, wb_mul_wc]);

            let (sumcheck_proof, challenges) = MultiComposedSumcheckProver::prove_partial(
                &vec![fbc_add_alpha_beta, fbc_mul_alpha_beta],
                &claimed_sum,
            )
            .unwrap();

            transcript.commit(&sumcheck_proof.to_bytes());
            sumcheck_proofs.push(sumcheck_proof);

            let (b, c) = challenges.split_at(&challenges.len() / 2);

            let eval_wb = wb.evaluation(&b);
            let eval_wc = wc.evaluation(&c);

            wb_s.push(eval_wb);
            wc_s.push(eval_wc);

            r_b = b.to_vec();
            r_c = c.to_vec();

            alpha = transcript.evaluate_challenge_into_field::<F>();
            beta = transcript.evaluate_challenge_into_field::<F>();

            if layer_index == circuit_evaluation.len() - 1 {
                let exponent_from_powers_of_tau = exponent(tau.powers_of_tau_in_g1.len());
                let blow_up_var_length = exponent_from_powers_of_tau - w_i_mle.n_vars;

                let poly: Multilinear<F> = w_i_mle.add_to_back(&blow_up_var_length);

                let mut b_clone = b.to_vec();
                let mut c_clone = c.to_vec();

                let padded_zeros_for_b_vec = &vec![F::zero(); &poly.n_vars - b_clone.len()];
                let padded_zeros_for_c_vec = &vec![F::zero(); &poly.n_vars - c_clone.len()];

                b_clone.extend(padded_zeros_for_b_vec);
                c_clone.extend(padded_zeros_for_c_vec);

                commitment = MultilinearKZG::commitment(&poly, &tau);

                proof_wb_opening = MultilinearKZG::open(&poly, &b_clone, &tau);
                proof_wc_opening = MultilinearKZG::open(&poly, &c_clone, &tau);

                claimed_sum = alpha * eval_wb + beta * eval_wc;
            } else {
                claimed_sum = alpha * eval_wb + beta * eval_wc;
            }
        }

        (
            commitment,
            SuccintGKRProof {
                sumcheck_proofs,
                wb_s,
                wc_s,
                w_0_mle,
                proof_wb_opening,
                proof_wc_opening,
            },
        )
    }

    pub fn verify(
        circuit: &Circuit,
        commitment: &P::G1,
        proof: &SuccintGKRProof<F, P>,
        tau: &TrustedSetup<P>,
    ) -> bool {
        if proof.sumcheck_proofs.len() != proof.wb_s.len()
            || proof.sumcheck_proofs.len() != proof.wc_s.len()
        {
            return false;
        }

        let mut transcript = FiatShamirTranscript::new();
        transcript.commit(&proof.w_0_mle.to_bytes());

        let n_r = transcript.evaluate_n_challenge_into_field(&proof.w_0_mle.n_vars);
        let mut claimed_sum = proof.w_0_mle.evaluation(&n_r.clone().as_slice());

        let mut r_b: Vec<F> = vec![];
        let mut r_c: Vec<F> = vec![];
        let mut alpha: F = F::one();
        let mut beta: F = F::one();

        let (add_mle_1, mult_mle_1) = circuit.add_mult_mle(0);
        let (status, r1_sum) = generate_layer_one_verify_sumcheck(
            &add_mle_1,
            &mult_mle_1,
            &proof.sumcheck_proofs[0],
            n_r,
            &claimed_sum,
            &mut transcript,
            &proof.wb_s[0],
            &proof.wc_s[0],
        );

        if !status {
            return false;
        }
        claimed_sum = r1_sum;

        for i in 1..proof.sumcheck_proofs.len() {
            if claimed_sum != proof.sumcheck_proofs[i].sum {
                return false;
            }

            transcript.commit(&proof.sumcheck_proofs[i].to_bytes());

            alpha = transcript.evaluate_challenge_into_field();
            beta = transcript.evaluate_challenge_into_field();

            let verify_subclaim =
                MultiComposedSumcheckVerifier::verify_partial(&proof.sumcheck_proofs[i]).unwrap();

            let (b, c) = verify_subclaim
                .challenges
                .split_at(&verify_subclaim.challenges.len() / 2);

            r_b = b.to_vec();
            r_c = c.to_vec();

            let wb = proof.wb_s[i];
            let wc = proof.wc_s[i];

            claimed_sum = alpha * wb + beta * wc;
        }

        let mut rb_clone = r_b.to_vec();
        let mut rc_clone = r_c.to_vec();

        let length_of_padded_zeros_for_b_vec =
            &vec![F::zero(); &tau.powers_of_tau_in_g2.len() - rb_clone.len()];
        let length_of_padded_zeros_for_c_vec =
            &vec![F::zero(); &tau.powers_of_tau_in_g2.len() - rc_clone.len()];

        rb_clone.extend(length_of_padded_zeros_for_b_vec);
        rc_clone.extend(length_of_padded_zeros_for_c_vec);

        let verify_rb =
            MultilinearKZG::verify(commitment, &rb_clone, &proof.proof_wb_opening, &tau);
        let verify_rc =
            MultilinearKZG::verify(commitment, &rc_clone, &proof.proof_wc_opening, &tau);

        let mut w_mle_rb_input = Default::default();
        let mut w_mle_rc_input = Default::default();

        if verify_rb && verify_rc {
            w_mle_rb_input = proof.proof_wb_opening.evaluation;
            w_mle_rc_input = proof.proof_wc_opening.evaluation;
        }

        let sum = alpha * w_mle_rb_input + beta * w_mle_rc_input;

        if claimed_sum != sum {
            return false;
        }

        true
    }
}

#[cfg(test)]
mod tests {
    use ark_test_curves::bls12_381::{Bls12_381, Fr};
    use circuit::{
        circuit::{Circuit, CircuitLayer},
        gate::{Gate, GateType},
    };
    use multilinear_kzg::{interface::TrustedSetupInterface, trusted_setup::TrustedSetup};

    use crate::succint_protocol::SuccintGKRProtocol;

    #[test]
    fn test_succint_gkr_protocol_1() {
        let layer_0 = CircuitLayer::new(vec![Gate::new(GateType::Mul, [0, 1])]);
        let layer_1 = CircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
        ]);
        let circuit = Circuit::new(vec![layer_0, layer_1]);
        let input = vec![
            Fr::from(2u32),
            Fr::from(3u32),
            Fr::from(4u32),
            Fr::from(5u32),
        ];

        let circuit_evaluation = circuit.evaluation(&input);

        let points = vec![Fr::from(54), Fr::from(90)];
        let tau = TrustedSetup::<Bls12_381>::setup(&points);
        let (commitment, proof) = SuccintGKRProtocol::prove(&circuit, &circuit_evaluation, &tau);
        let verify = SuccintGKRProtocol::verify(&circuit, &commitment, &proof, &tau);

        assert_eq!(verify, true);
    }

    #[test]
    fn test_succint_gkr_protocol_2() {
        let layer_0 = CircuitLayer::new(vec![Gate::new(GateType::Add, [0, 1])]);
        let layer_1 = CircuitLayer::new(vec![
            Gate::new(GateType::Mul, [0, 1]),
            Gate::new(GateType::Add, [2, 3]),
        ]);
        let layer_3 = CircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
            Gate::new(GateType::Mul, [4, 5]),
            Gate::new(GateType::Mul, [6, 7]),
        ]);

        let circuit = Circuit::new(vec![layer_0, layer_1, layer_3]);
        let input = vec![
            Fr::from(4u32),
            Fr::from(3u32),
            Fr::from(7u32),
            Fr::from(6u32),
            Fr::from(6u32),
            Fr::from(1u32),
            Fr::from(4u32),
            Fr::from(2u32),
        ];

        let circuit_evaluation = circuit.evaluation(&input);

        assert_eq!(circuit_evaluation[0][0], Fr::from(308u32));
        let points = vec![Fr::from(54), Fr::from(90), Fr::from(76)];

        let tau = TrustedSetup::<Bls12_381>::setup(&points);

        let (commitment, proof) =
            SuccintGKRProtocol::<Fr, Bls12_381>::prove(&circuit, &circuit_evaluation, &tau);
        let verify =
            SuccintGKRProtocol::<Fr, Bls12_381>::verify(&circuit, &commitment, &proof, &tau);

        assert!(&verify);
    }

    #[test]
    fn test_succint_gkr_protocol_3() {
        let layer_0 = CircuitLayer::new(vec![Gate::new(GateType::Add, [0, 1])]);
        let layer_1 = CircuitLayer::new(vec![
            Gate::new(GateType::Mul, [0, 1]),
            Gate::new(GateType::Add, [2, 3]),
        ]);
        let layer_3 = CircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
            Gate::new(GateType::Mul, [4, 5]),
            Gate::new(GateType::Mul, [6, 7]),
        ]);
        let layer_4 = CircuitLayer::new(vec![
            Gate::new(GateType::Mul, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
            Gate::new(GateType::Mul, [4, 5]),
            Gate::new(GateType::Add, [6, 7]),
            Gate::new(GateType::Mul, [8, 9]),
            Gate::new(GateType::Add, [10, 11]),
            Gate::new(GateType::Mul, [12, 13]),
            Gate::new(GateType::Mul, [14, 15]),
        ]);

        let circuit = Circuit::new(vec![layer_0, layer_1, layer_3, layer_4]);
        let input = vec![
            Fr::from(2u32),
            Fr::from(1u32),
            Fr::from(3u32),
            Fr::from(1u32),
            Fr::from(4u32),
            Fr::from(1u32),
            Fr::from(2u32),
            Fr::from(2u32),
            Fr::from(3u32),
            Fr::from(3u32),
            Fr::from(4u32),
            Fr::from(4u32),
            Fr::from(2u32),
            Fr::from(3u32),
            Fr::from(3u32),
            Fr::from(4u32),
        ];

        let circuit_evaluation = circuit.evaluation(&input);

        assert_eq!(circuit_evaluation[0][0], Fr::from(224u32));
        let points = vec![Fr::from(12), Fr::from(9), Fr::from(28), Fr::from(40)];

        let tau = TrustedSetup::<Bls12_381>::setup(&points);

        let (commitment, proof) = SuccintGKRProtocol::prove(&circuit, &circuit_evaluation, &tau);

        let verify = SuccintGKRProtocol::verify(&circuit, &commitment, &proof, &tau);
        assert!(verify);
    }
}
