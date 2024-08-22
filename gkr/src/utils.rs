use ark_ff::PrimeField;
use fiat_shamir::{fiat_shamir::FiatShamirTranscript, interface::FiatShamirTranscriptTrait};
use polynomial::{ComposedMultilinear, Multilinear, MultilinearTrait};
use sumcheck::composed::multi_composed_sumcheck::{
    ComposedSumcheckProof, MultiComposedSumcheckProver, MultiComposedSumcheckVerifier,
};

pub fn size_of_mle_n_var_at_each_layer(layer_index: usize) -> usize {
    if layer_index == 0 {
        return 1 << 3;
    }

    let layer_index_plus_one = layer_index + 1;
    let number_of_variable = layer_index + (2 * layer_index_plus_one);

    1 << number_of_variable
}

pub fn transform_label_to_binary_and_to_decimal(
    layer_index: usize,
    a: usize,
    b: usize,
    c: usize,
) -> usize {
    let a_binary_string: String = binary_string(a, layer_index);
    let b_binary_string: String = binary_string(b, layer_index + 1);
    let c_binary_string: String = binary_string(c, layer_index + 1);

    let combined_binary_string = a_binary_string + &b_binary_string + &c_binary_string;

    usize::from_str_radix(&combined_binary_string, 2).unwrap_or(0)
}

/// Convert a number to a binary string of a given size
pub fn binary_string(index: usize, mut bit_count: usize) -> String {
    if bit_count == 0 {
        bit_count = 1;
    }
    let binary = format!("{:b}", index);
    "0".repeat(bit_count.saturating_sub(binary.len())) + &binary
}

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_size_of_mle_n_var() {
        assert_eq!(size_of_mle_n_var_at_each_layer(0), 8);
        assert_eq!(size_of_mle_n_var_at_each_layer(1), 32);
        assert_eq!(size_of_mle_n_var_at_each_layer(2), 256);
        assert_ne!(size_of_mle_n_var_at_each_layer(2), 128);
        assert_ne!(size_of_mle_n_var_at_each_layer(3), 512);
        assert_eq!(size_of_mle_n_var_at_each_layer(3), 2048);
        assert_eq!(size_of_mle_n_var_at_each_layer(4), 16384);
    }

    #[test]
    fn test_transform_binary_and_to_decimal() {
        // a at layer 0, b & c at layer 1
        assert_eq!(transform_label_to_binary_and_to_decimal(1, 1, 2, 3), 27);
        assert_eq!(transform_label_to_binary_and_to_decimal(2, 1, 2, 3), 83);
    }

    #[test]
    fn test_binary_string() {
        assert_eq!(binary_string(0, 0), "0");
        assert_eq!(binary_string(0, 1), "0");
        assert_eq!(binary_string(0, 2), "00");
        assert_eq!(binary_string(5, 3), "101");
    }
}
