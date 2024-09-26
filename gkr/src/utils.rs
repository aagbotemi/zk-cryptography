use ark_ff::PrimeField;
use fiat_shamir::{fiat_shamir::FiatShamirTranscript, interface::FiatShamirTranscriptTrait};
use polynomial::{ComposedMultilinear, Multilinear, MultilinearTrait};
use sumcheck::composed::multi_composed_sumcheck::{
    ComposedSumcheckProof, MultiComposedSumcheckProver, MultiComposedSumcheckVerifier,
};

pub fn w_mle<F: PrimeField>(layer_eval: Vec<F>) -> Multilinear<F> {
    Multilinear::new(layer_eval)
}

pub fn generate_layer_one_prove_sumcheck<F: PrimeField>(
    add_mle: &Multilinear<F>,
    mult_mle: &Multilinear<F>,
    w_1_mle: &Multilinear<F>,
    n_r: &Vec<F>,
    sum: &F,
    transcript: &mut FiatShamirTranscript,
    sumcheck_proofs: &mut Vec<ComposedSumcheckProof<F>>,
    wb_s: &mut Vec<F>,
    wc_s: &mut Vec<F>,
) -> (F, F, F, Vec<F>, Vec<F>) {
    let add_rbc = add_mle.partial_evaluations(&n_r, &vec![0; n_r.len()]);
    let mul_rbc = mult_mle.partial_evaluations(&n_r, &vec![0; n_r.len()]);

    let wb = w_1_mle.clone();
    let wc = w_1_mle;

    let wb_add_wc = wb.add_distinct(&wc);
    let wb_mul_wc = wb.mul_distinct(&wc);

    let add_fbc = ComposedMultilinear::new(vec![add_rbc, wb_add_wc]);
    let mul_fbc = ComposedMultilinear::new(vec![mul_rbc, wb_mul_wc]);

    let (sumcheck_proof, challenges) =
        MultiComposedSumcheckProver::prove_partial(&vec![add_fbc, mul_fbc], &sum).unwrap();
    transcript.commit(&sumcheck_proof.to_bytes());
    sumcheck_proofs.push(sumcheck_proof);

    let (b, c) = challenges.split_at(&challenges.len() / 2);

    let eval_wb = wb.evaluation(b);
    let eval_wc = wc.evaluation(c);
    wb_s.push(eval_wb);
    wc_s.push(eval_wc);

    let alpha = transcript.evaluate_challenge_into_field::<F>();
    let beta = transcript.evaluate_challenge_into_field::<F>();

    let new_claim: F = alpha * eval_wb + beta * eval_wc;

    let claimed_sum = new_claim;
    let rb = b.to_vec();
    let rc = c.to_vec();

    (claimed_sum, alpha, beta, rb, rc)
}

pub fn generate_layer_one_verify_sumcheck<F: PrimeField>(
    add_mle: &Multilinear<F>,
    mult_mle: &Multilinear<F>,
    proof: &ComposedSumcheckProof<F>,
    n_r: Vec<F>,
    sum: &F,
    transcript: &mut FiatShamirTranscript,
    wb: &F,
    wc: &F,
) -> (bool, F) {
    if *sum != proof.sum {
        return (false, F::zero());
    }

    transcript.commit(&proof.to_bytes());

    let verify_subclaim = MultiComposedSumcheckVerifier::verify_partial(proof).unwrap();

    let mut rbc = n_r;
    rbc.extend_from_slice(&verify_subclaim.challenges);

    let add_bc = add_mle.evaluation(&rbc);
    let mul_bc = mult_mle.evaluation(&rbc);

    let fbc_add = add_bc * (*wb + *wc);
    let fbc_mul = mul_bc * (*wb * *wc);

    let fbc_eval = fbc_add + fbc_mul;

    if fbc_eval != verify_subclaim.sum {
        return (false, F::zero());
    }

    let alpha = transcript.evaluate_challenge_into_field::<F>();
    let beta = transcript.evaluate_challenge_into_field::<F>();

    let new_claim: F = alpha * wb + beta * wc;

    (true, new_claim)
}

pub fn exponent(value: usize) -> usize {
    let mut num: usize = value;
    let mut exponent: usize = 0;

    while num > 1 {
        assert_eq!(num % 2, 0, "Value is not a power of 2");
        num /= 2;
        exponent += 1;
    }

    exponent
}
