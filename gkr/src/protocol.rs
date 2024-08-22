use ark_ff::PrimeField;
use fiat_shamir::{fiat_shamir::FiatShamirTranscript, interface::FiatShamirTranscriptTrait};
use polynomial::{ComposedMultilinear, Multilinear, MultilinearTrait};
use sumcheck::composed::multi_composed_sumcheck::{
    ComposedSumcheckProof, MultiComposedSumcheckProver, MultiComposedSumcheckVerifier,
};

use crate::{
    circuit::GKRCircuit,
    utils::{generate_layer_one_prove_sumcheck, generate_layer_one_verify_sumcheck, w_mle},
};

pub struct GKRProof<F: PrimeField> {
    sumcheck_proofs: Vec<ComposedSumcheckProof<F>>,
    wb_s: Vec<F>,            // w_mle for layer one onward for rb
    wc_s: Vec<F>,            // w_mle for layer one onward for rc
    w_0_mle: Multilinear<F>, // w_mle for layer
}

pub struct GKRProtocol {}

impl GKRProtocol {
    /// Prove correct circuit evaluation using the GKR protocol
    pub fn prove<'a, F: PrimeField>(circuit: &'a GKRCircuit, input: &'a Vec<F>) -> GKRProof<F> {
        let mut transcript = FiatShamirTranscript::new();
        let mut sumcheck_proofs: Vec<ComposedSumcheckProof<F>> = Vec::new();
        let mut wb_s: Vec<F> = Vec::new();
        let mut wc_s: Vec<F> = Vec::new();

        let circuit_eval = circuit.evaluation(input);
        let mut circuit_eval_layer_zero_pad = circuit_eval.layers[0].clone();
        circuit_eval_layer_zero_pad.push(F::zero());

        let w_0_mle = w_mle(circuit_eval_layer_zero_pad.to_vec());
        transcript.commit(&w_0_mle.to_bytes());

        let n_r: Vec<F> = transcript.evaluate_n_challenge_into_field(&w_0_mle.n_vars);
        let mut claimed_sum: F = w_0_mle.evaluation(&n_r);

        let (add_mle_1, mult_mle_1) = circuit.add_mult_mle::<F>(0);
        let w_1_mle = w_mle(circuit_eval.layers[1].to_vec());

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

        for layer_index in 2..circuit_eval.layers.len() {
            let (add_mle, mult_mle) = circuit.add_mult_mle::<F>(layer_index - 1);

            let add_rb_bc = add_mle.partial_evaluations(&r_b, &vec![0; r_b.len()]);
            let mul_rb_bc = mult_mle.partial_evaluations(&r_b, &vec![0; r_b.len()]);

            let add_rc_bc = add_mle.partial_evaluations(&r_c, &vec![0; r_b.len()]);
            let mul_rc_bc = mult_mle.partial_evaluations(&r_c, &vec![0; r_b.len()]);
            let w_i_mle = w_mle(circuit_eval.layers[layer_index].to_vec());

            let wb = w_i_mle.clone();
            let wc = w_i_mle;

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

            claimed_sum = alpha * eval_wb + beta * eval_wc;
        }

        GKRProof {
            sumcheck_proofs,
            wb_s,
            wc_s,
            w_0_mle,
        }
    }

    pub fn verify<F: PrimeField>(circuit: &GKRCircuit, input: &[F], proof: &GKRProof<F>) -> bool {
        if proof.sumcheck_proofs.len() != proof.wb_s.len()
            || proof.sumcheck_proofs.len() != proof.wc_s.len()
        {
            return false;
        }

        let mut transcript = FiatShamirTranscript::new();
        transcript.commit(&proof.w_0_mle.to_bytes());

        let n_r: Vec<F> = transcript.evaluate_n_challenge_into_field::<F>(&proof.w_0_mle.n_vars);
        let mut claimed_sum = proof.w_0_mle.evaluation(&n_r.clone().as_slice());

        let mut r_b: Vec<F> = vec![];
        let mut r_c: Vec<F> = vec![];
        let mut alpha: F = F::zero();
        let mut beta: F = F::zero();

        let (add_mle_1, mult_mle_1) = circuit.add_mult_mle::<F>(0);
        let (status, sum) = generate_layer_one_verify_sumcheck(
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

        claimed_sum = sum;

        for i in 1..proof.sumcheck_proofs.len() {
            if claimed_sum != proof.sumcheck_proofs[i].sum {
                return false;
            }

            transcript.commit(&proof.sumcheck_proofs[i].to_bytes());

            let verify_subclaim =
                MultiComposedSumcheckVerifier::verify_partial(&proof.sumcheck_proofs[i]).unwrap();
            // println!("verify_subclaim={:?}", verify_subclaim);

            let (b, c) = verify_subclaim
                .challenges
                .split_at(&verify_subclaim.challenges.len() / 2);

            r_b = b.to_vec();
            r_c = c.to_vec();

            let wb = proof.wb_s[i];
            let wc = proof.wc_s[i];

            let alph = transcript.evaluate_challenge_into_field::<F>();
            let bta = transcript.evaluate_challenge_into_field::<F>();

            claimed_sum = alph * wb + bta * wc;

            alpha = alph;
            beta = bta;
        }

        let w_mle_input = w_mle(input.to_vec());

        let w_mle_rb_input = w_mle_input.evaluation(&r_b);
        let w_mle_rc_input = w_mle_input.evaluation(&r_c);

        let sum = alpha * w_mle_rb_input + beta * w_mle_rc_input;

        if claimed_sum != sum {
            return false;
        }

        true
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        circuit::GKRCircuitLayer,
        gate::{Gate, GateType},
    };

    use super::*;

    use ark_ff::MontConfig;
    use ark_ff::{Fp64, MontBackend};

    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    struct FqConfig;
    type Fq = Fp64<MontBackend<FqConfig, 1>>;

    #[test]
    fn test_gkr_protocol_1() {
        let layer_0 = GKRCircuitLayer::new(vec![Gate::new(GateType::Mul, [0, 1])]);
        let layer_1 = GKRCircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
        ]);
        let circuit = GKRCircuit::new(vec![layer_0, layer_1]);
        let input = vec![
            Fq::from(2u32),
            Fq::from(3u32),
            Fq::from(4u32),
            Fq::from(5u32),
        ];

        let proof = GKRProtocol::prove(&circuit, &input);
        let verify = GKRProtocol::verify(&circuit, &input, &proof);

        assert!(verify);
    }

    #[test]
    fn test_gkr_protocol_2() {
        let layer_0 = GKRCircuitLayer::new(vec![Gate::new(GateType::Add, [0, 1])]);
        let layer_1 = GKRCircuitLayer::new(vec![
            Gate::new(GateType::Mul, [0, 1]),
            Gate::new(GateType::Add, [2, 3]),
        ]);
        let layer_3 = GKRCircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
            Gate::new(GateType::Mul, [4, 5]),
            Gate::new(GateType::Mul, [6, 7]),
        ]);
        let layer_4 = GKRCircuitLayer::new(vec![
            Gate::new(GateType::Mul, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
            Gate::new(GateType::Mul, [4, 5]),
            Gate::new(GateType::Add, [6, 7]),
            Gate::new(GateType::Mul, [8, 9]),
            Gate::new(GateType::Add, [10, 11]),
            Gate::new(GateType::Mul, [12, 13]),
            Gate::new(GateType::Mul, [14, 15]),
        ]);

        let circuit = GKRCircuit::new(vec![layer_0, layer_1, layer_3, layer_4]);
        let input = [
            Fq::from(2u32),
            Fq::from(1u32),
            Fq::from(3u32),
            Fq::from(1u32),
            Fq::from(4u32),
            Fq::from(1u32),
            Fq::from(2u32),
            Fq::from(2u32),
            Fq::from(3u32),
            Fq::from(3u32),
            Fq::from(4u32),
            Fq::from(4u32),
            Fq::from(2u32),
            Fq::from(3u32),
            Fq::from(3u32),
            Fq::from(4u32),
        ];

        let evaluation = circuit.evaluation(&input);

        assert_eq!(evaluation.layers[0][0], Fq::from(224u32));

        let proof = GKRProtocol::prove(&circuit, &input.to_vec());

        assert!(GKRProtocol::verify(&circuit, &input, &proof));
    }
}
