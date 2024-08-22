use ark_ff::PrimeField;
use fiat_shamir::{fiat_shamir::FiatShamirTranscript, interface::FiatShamirTranscriptTrait};
use polynomial::{ComposedMLE, MLETrait, MLE};
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
    // dbg!(&combined_binary_string);



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












    // dbg!(&usize::from_str_radix(&combined_binary, 2).unwrap_or(0));




/// Determines the number of bits needed to represent a number
pub fn bit_count_for_n_elem(size: usize) -> usize {
    // if the size of the array is 2, this will say two binary digits are needed
    // but since array indexing starts at 0 then only 1 binary digit will be needed
    // i.e first element = 0, second element = 1
    // hence we need to subtract 1 from the array size inorder to account for zero indexing
    format!("{:b}", size - 1).len()
}

pub fn transform_label_to_binary_and_to_decimals(a: usize, b: usize, c: usize) -> usize {
    let a_binary: String = format!("{:b}", a);
    let b_binary: String = format!("{:02b}", b);
    let c_binary: String = format!("{:02b}", c);

    let combined_binary = format!("{}{}{}", a_binary, b_binary, c_binary);
    usize::from_str_radix(&combined_binary, 2).unwrap_or(0)
}

pub fn w_mle<F: PrimeField>(layer_eval: Vec<F>) -> MLE<F> {
    MLE::new(layer_eval)
}

pub fn generate_layer_one_prove_sumcheck<F: PrimeField>(
    add_mle: &MLE<F>,
    mult_mle: &MLE<F>,
    w_1_mle: &MLE<F>,
    n_r: &Vec<F>,
    sum: &F,
    transcript: &mut FiatShamirTranscript,
    sumcheck_proofs: &mut Vec<ComposedSumcheckProof<F>>,
    wb_s: &mut Vec<F>,
    wc_s: &mut Vec<F>,
) -> (F, F, F, Vec<F>, Vec<F>) {
    let add_rbc = add_mle.partial_evaluations(&n_r, &vec![0; n_r.len()]);
    let mul_rbc = mult_mle.partial_evaluations(&n_r, &vec![0; n_r.len()]);

    let add_rbc_len = add_rbc.evaluations.len();
    let mul_rbc_len = mul_rbc.evaluations.len();

    let wb = w_1_mle.clone();
    let wc = w_1_mle;

    let wb_add_wc = wb.add_distinct(&wc);
    let wb_mul_wc = wb.mul_distinct(&wc);

    let add_fbc = ComposedMLE::new(vec![add_rbc, wb_add_wc]);
    let mul_fbc = ComposedMLE::new(vec![mul_rbc, wb_mul_wc]);

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
    // dbg!(new_claim);

    let claimed_sum = new_claim;
    let rb = b.to_vec();
    let rc = c.to_vec();

    (claimed_sum, alpha, beta, rb, rc)
}

pub fn generate_layer_one_verify_sumcheck<F: PrimeField>(
    add_mle: &MLE<F>,
    mult_mle: &MLE<F>,
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
    // dbg!(new_claim);

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
        // a at layer 1, b & c at layer 2
        // assert_eq!(transform_label_to_binary_and_to_decimal(1, 1, 2, 3), 27);
        // // a at layer 2, b & c at layer 3
        // assert_eq!(transform_label_to_binary_and_to_decimal(2, 3, 6, 7), 247);
    }

    #[test]
    fn test_binary_string() {
        assert_eq!(binary_string(0, 0), "0");
        assert_eq!(binary_string(0, 1), "0");
        assert_eq!(binary_string(0, 2), "00");
        assert_eq!(binary_string(5, 3), "101");
    }
}
