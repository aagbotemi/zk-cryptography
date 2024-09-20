use std::marker::PhantomData;

use ark_ec::pairing::Pairing;
use ark_ff::{One, Zero};
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
pub struct SuccintGKRProof<P: Pairing> {
    sumcheck_proofs: Vec<ComposedSumcheckProof<P::ScalarField>>,
    wb_s: Vec<P::ScalarField>,
    wc_s: Vec<P::ScalarField>,
    w_0_mle: Multilinear<P::ScalarField>,
    proof_wb_opening: MultilinearKZGProof<P>,
    proof_wc_opening: MultilinearKZGProof<P>,
}

pub struct SuccintGKRProtocol<P: Pairing> {
    _marker: PhantomData<P>,
}

impl<P: Pairing> SuccintGKRProtocol<P> {
    /// Prove correct circuit evaluation using the GKR protocol
    pub fn prove(
        circuit: &Circuit,
        input: &Vec<P::ScalarField>,
        tau: &TrustedSetup<P>,
    ) -> (P::G1, SuccintGKRProof<P>) {
        let mut transcript = FiatShamirTranscript::new();
        let mut sumcheck_proofs: Vec<ComposedSumcheckProof<P::ScalarField>> = Vec::new();
        let mut wb_s: Vec<P::ScalarField> = Vec::new();
        let mut wc_s: Vec<P::ScalarField> = Vec::new();

        let circuit_evaluation = circuit.evaluation(input);
        let mut circuit_evaluation_layer_zero_pad = circuit_evaluation[0].clone();
        circuit_evaluation_layer_zero_pad.push(P::ScalarField::zero());

        let w_0_mle = w_mle::<P>(circuit_evaluation_layer_zero_pad.to_vec());
        transcript.commit(&w_0_mle.to_bytes());

        let n_r: Vec<P::ScalarField> = transcript.evaluate_n_challenge_into_field(&w_0_mle.n_vars);
        let mut claimed_sum: P::ScalarField = w_0_mle.evaluation(&n_r);

        let (add_mle_1, mult_mle_1) = circuit.add_mult_mle::<P::ScalarField>(0);
        let w_1_mle = w_mle::<P>(circuit_evaluation[1].to_vec());

        let (claimed, alph, bta, rb, rc) = generate_layer_one_prove_sumcheck::<P>(
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

        let mut alpha: P::ScalarField = alph;
        let mut beta: P::ScalarField = bta;

        let mut r_b: Vec<P::ScalarField> = rb;
        let mut r_c: Vec<P::ScalarField> = rc;

        let mut commitment: P::G1 = Default::default();
        let mut proof_wb_opening: MultilinearKZGProof<P> = Default::default();
        let mut proof_wc_opening: MultilinearKZGProof<P> = Default::default();

        for layer_index in 2..circuit_evaluation.len() {
            let (add_mle, mult_mle) = circuit.add_mult_mle::<P::ScalarField>(layer_index - 1);

            let add_rb_bc = add_mle.partial_evaluations(&r_b, &vec![0; r_b.len()]);
            let mul_rb_bc = mult_mle.partial_evaluations(&r_b, &vec![0; r_b.len()]);

            let add_rc_bc = add_mle.partial_evaluations(&r_c, &vec![0; r_b.len()]);
            let mul_rc_bc = mult_mle.partial_evaluations(&r_c, &vec![0; r_b.len()]);
            let w_i_mle = w_mle::<P>(circuit_evaluation[layer_index].to_vec());

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

            alpha = transcript.evaluate_challenge_into_field::<P::ScalarField>();
            beta = transcript.evaluate_challenge_into_field::<P::ScalarField>();

            if layer_index == circuit_evaluation.len() - 1 {
                let exponent_from_powers_of_tau = exponent(tau.powers_of_tau_in_g1.len());
                let blow_up_var_length = exponent_from_powers_of_tau - w_i_mle.n_vars;
                let poly: Multilinear<<P as Pairing>::ScalarField> =
                    w_i_mle.add_to_back(&blow_up_var_length);

                let mut b_clone = b.to_vec();
                let mut c_clone = c.to_vec();

                let padded_zeros_for_b_vec =
                    &vec![P::ScalarField::zero(); &poly.n_vars - b_clone.len()];
                let padded_zeros_for_c_vec =
                    &vec![P::ScalarField::zero(); &poly.n_vars - c_clone.len()];

                b_clone.extend(padded_zeros_for_b_vec);
                c_clone.extend(padded_zeros_for_c_vec);

                commitment = MultilinearKZG::<P>::commitment(&poly, &tau.powers_of_tau_in_g1);
                proof_wb_opening =
                    MultilinearKZG::<P>::open(&poly, &b_clone, &tau.powers_of_tau_in_g1);
                proof_wc_opening =
                    MultilinearKZG::<P>::open(&poly, &c_clone, &tau.powers_of_tau_in_g1);

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
        proof: &SuccintGKRProof<P>,
        tau: &TrustedSetup<P>,
    ) -> bool {
        if proof.sumcheck_proofs.len() != proof.wb_s.len()
            || proof.sumcheck_proofs.len() != proof.wc_s.len()
        {
            return false;
        }

        let mut transcript = FiatShamirTranscript::new();
        transcript.commit(&proof.w_0_mle.to_bytes());

        let n_r: Vec<P::ScalarField> =
            transcript.evaluate_n_challenge_into_field::<P::ScalarField>(&proof.w_0_mle.n_vars);
        let mut claimed_sum = proof.w_0_mle.evaluation(&n_r.clone().as_slice());

        let mut r_b: Vec<P::ScalarField> = vec![];
        let mut r_c: Vec<P::ScalarField> = vec![];
        let mut alpha: P::ScalarField = P::ScalarField::one();
        let mut beta: P::ScalarField = P::ScalarField::one();

        let (add_mle_1, mult_mle_1) = circuit.add_mult_mle::<P::ScalarField>(0);
        let (status, r1_sum) = generate_layer_one_verify_sumcheck::<P>(
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

            alpha = transcript.evaluate_challenge_into_field::<P::ScalarField>();
            beta = transcript.evaluate_challenge_into_field::<P::ScalarField>();

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
            &vec![P::ScalarField::zero(); &tau.powers_of_tau_in_g2.len() - rb_clone.len()];
        let length_of_padded_zeros_for_c_vec =
            &vec![P::ScalarField::zero(); &tau.powers_of_tau_in_g2.len() - rc_clone.len()];

        rb_clone.extend(length_of_padded_zeros_for_b_vec);
        rc_clone.extend(length_of_padded_zeros_for_c_vec);

        let verify_rb = MultilinearKZG::verify(
            commitment,
            &rb_clone,
            &proof.proof_wb_opening,
            &tau.powers_of_tau_in_g2,
        );
        let verify_rc = MultilinearKZG::verify(
            commitment,
            &rc_clone,
            &proof.proof_wc_opening,
            &tau.powers_of_tau_in_g2,
        );

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

    use crate::succint_gkr::SuccintGKRProtocol;

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

        let tau = TrustedSetup::<Bls12_381>::setup(&input);

        let (commitment, proof) = SuccintGKRProtocol::prove(&circuit, &input, &tau);
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

        let evaluation = circuit.evaluation(&input);

        assert_eq!(evaluation[0][0], Fr::from(308u32));

        let tau = TrustedSetup::<Bls12_381>::setup(&input);

        let (commitment, proof) = SuccintGKRProtocol::<Bls12_381>::prove(&circuit, &input, &tau);
        let verify = SuccintGKRProtocol::<Bls12_381>::verify(&circuit, &commitment, &proof, &tau);

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

        let evaluation = circuit.evaluation(&input);

        assert_eq!(evaluation[0][0], Fr::from(224u32));

        let tau = TrustedSetup::<Bls12_381>::setup(&input);

        let (commitment, proof) = SuccintGKRProtocol::prove(&circuit, &input, &tau);

        let verify = SuccintGKRProtocol::verify(&circuit, &commitment, &proof, &tau);
        assert!(verify);
    }
}
