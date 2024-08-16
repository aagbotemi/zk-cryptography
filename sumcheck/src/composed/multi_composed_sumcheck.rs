use super::composed_sumcheck::ComposedSumcheck;
use crate::utils::{convert_field_to_byte, convert_round_poly_to_uni_poly_format, vec_to_bytes};
use ark_ff::{BigInteger, PrimeField};
use fiat_shamir::{fiat_shamir::FiatShamirTranscript, interface::FiatShamirTranscriptTrait};
use polynomial::{
    interface::ComposedMLETrait, ComposedMLE, UnivariatePolynomial, UnivariatePolynomialTrait,
};

// #[derive(Debug, Clone)]
// pub struct MultiComposedSumcheckProver<F: PrimeField> {
//     pub poly: Vec<ComposedMLE<F>>,
//     pub sum: F,
// }

#[derive(Debug)]
pub struct ComposedSumcheckProof<F: PrimeField> {
    // pub poly: Vec<ComposedMLE<F>>,
    pub round_polys: Vec<UnivariatePolynomial<F>>,
    pub sum: F,
}

#[derive(Debug)]
pub struct SubClaim<F: PrimeField> {
    pub sum: F,
    pub challenges: Vec<F>,
}

impl<F: PrimeField> ComposedSumcheckProof<F> {
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();
        for round_poly in self.round_polys.iter() {
            bytes.extend_from_slice(&round_poly.to_bytes());
        }
        // bytes.extend_from_slice(&self.sum.into_bigint().to_bytes_be());
        bytes
    }
}

pub fn composed_mle_to_bytes<F: PrimeField>(poly: &[ComposedMLE<F>]) -> Vec<u8> {
    let mut bytes = Vec::new();
    for p in poly.iter() {
        bytes.extend_from_slice(&p.to_bytes());
    }
    bytes
}

pub struct MultiComposedSumcheckProver {}
impl MultiComposedSumcheckProver {
    pub fn calculate_poly_sum<F: PrimeField>(poly: &Vec<ComposedMLE<F>>) -> F {
        let mut sum = F::zero();

        for p in poly {
            sum += ComposedSumcheck::calculate_poly_sum(p);
        }

        sum
    }

    pub fn prove<F: PrimeField>(
        poly: &Vec<ComposedMLE<F>>,
        sum: &F,
    ) -> Result<(ComposedSumcheckProof<F>, Vec<F>), &'static str> {
        let mut transcript = FiatShamirTranscript::new();
        transcript.commit(&composed_mle_to_bytes(&poly));
        MultiComposedSumcheckProver::prove_internal(&poly, &sum, &mut transcript)
    }

    pub fn prove_partial<F: PrimeField>(
        poly: &Vec<ComposedMLE<F>>,
        sum: &F,
    ) -> Result<(ComposedSumcheckProof<F>, Vec<F>), &'static str> {
        let mut transcript = FiatShamirTranscript::new();
        MultiComposedSumcheckProver::prove_internal(&poly, &sum, &mut transcript)
    }

    pub fn prove_internal<F: PrimeField>(
        poly: &Vec<ComposedMLE<F>>,
        sum: &F,
        transcript: &mut FiatShamirTranscript,
    ) -> Result<(ComposedSumcheckProof<F>, Vec<F>), &'static str> {
        // append the sum to the transcript
        transcript.commit(&convert_field_to_byte(sum));

        let mut current_poly = poly.clone();
        let mut round_polys = vec![];
        let mut challenges: Vec<F> = vec![];

        for _ in 0..poly[0].n_vars() {
            let mut round_poly = UnivariatePolynomial::zero();

            for p in current_poly.iter() {
                let mut round_i_poly_vec = Vec::new();
                for i in 0..=p.max_degree() {
                    let round: F = p
                        .partial_evaluation(F::from(i as u32), 0)
                        .element_wise_product()
                        .iter()
                        .sum::<F>();

                    round_i_poly_vec.push(round);
                }

                let round_i_poly = UnivariatePolynomial::interpolation(
                    &convert_round_poly_to_uni_poly_format(&round_i_poly_vec),
                );
                round_poly = round_poly + round_i_poly;
            }

            transcript.commit(&round_poly.to_bytes());
            //get the random r
            let random_r: F = transcript.evaluate_challenge_into_field::<F>();

            let mut new_poly = Vec::new();

            for i in 0..current_poly.len() {
                new_poly.push(current_poly[i].partial_evaluation(random_r, 0));
            }

            current_poly = new_poly;

            challenges.push(random_r);
            round_polys.push(round_poly);
        }

        println!("challenges_at_prove={:?}", challenges);

        Ok((
            ComposedSumcheckProof {
                // poly: poly.clone(),
                round_polys,
                sum: *sum,
            },
            challenges,
        ))
    }
}

pub struct MultiComposedSumcheckVerifier {}

impl MultiComposedSumcheckVerifier {
    pub fn verify<F: PrimeField>(
        poly: &Vec<ComposedMLE<F>>,
        proof: &ComposedSumcheckProof<F>,
    ) -> Result<bool, &'static str> {
        let mut transcript = FiatShamirTranscript::new();

        transcript.commit(&composed_mle_to_bytes(&poly));
        let sub_claim = Self::verify_internal(&proof, &mut transcript)?;

        // oracle check
        let mut poly_pe_sum = F::zero();
        for p in poly.iter() {
            poly_pe_sum += p.evaluation(&sub_claim.challenges.as_slice())
        }

        Ok(poly_pe_sum == sub_claim.sum)
    }
    pub fn verify_partial<F: PrimeField>(
        proof: &ComposedSumcheckProof<F>,
    ) -> Result<SubClaim<F>, &'static str> {
        let mut transcript = FiatShamirTranscript::new();
        let sub_claim = Self::verify_internal(&proof, &mut transcript);
        Ok(sub_claim)?
    }

    pub fn verify_internal<F: PrimeField>(
        proof: &ComposedSumcheckProof<F>,
        transcript: &mut FiatShamirTranscript,
    ) -> Result<SubClaim<F>, &'static str> {
        // append the sum to the transcript
        transcript.commit(&convert_field_to_byte(&proof.sum));

        let mut claimed_sum = proof.sum;
        let mut challenges: Vec<F> = vec![];

        for round_poly in proof.round_polys.iter() {
            transcript.commit(&round_poly.to_bytes());
            // genrate the challenge for this round
            let challenge: F = transcript.evaluate_challenge_into_field::<F>();
            challenges.push(challenge);

            let eval_p0_p1 = round_poly.evaluate(F::zero()) + round_poly.evaluate(F::one());
            // println!("eval_p0_p1={:?}, claimed_Sum={:?}", eval_p0_p1, claimed_sum);

            if claimed_sum != eval_p0_p1 {
                return Err("Verification failed");
            }

            // update the sum
            claimed_sum = round_poly.evaluate(challenge);
        }
        println!("challenges_at_verify={:?}", challenges);

        Ok(SubClaim {
            sum: claimed_sum,
            challenges,
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::sumcheck;

    use super::*;
    use ark_ff::MontConfig;
    use ark_ff::{Fp64, MontBackend};
    use polynomial::interface::MLETrait;
    use polynomial::MLE;

    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    struct FqConfig;
    type Fq = Fp64<MontBackend<FqConfig, 1>>;

    #[test]
    fn test_sum_calculation() {
        let mle1 = MLE::new(vec![Fq::from(0), Fq::from(1), Fq::from(2), Fq::from(3)]);
        let mle2 = MLE::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(1)]);
        let composedmle1 = ComposedMLE::new(vec![mle1]);
        let composedmle2 = ComposedMLE::new(vec![mle2]);

        let multi_composed_vec_1 = vec![composedmle1, composedmle2];
        let sum_1 = MultiComposedSumcheckProver::calculate_poly_sum(&multi_composed_vec_1);
        assert_eq!(sum_1, Fq::from(7));

        let mle3 = MLE::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)]);
        let mle4 = MLE::new(vec![Fq::from(0), Fq::from(3), Fq::from(0), Fq::from(3)]);
        let composedmle3 = ComposedMLE::new(vec![mle3]);
        let composedmle4 = ComposedMLE::new(vec![mle4]);

        let multi_composed_vec_2 = vec![composedmle3, composedmle4];
        let sum_2 = MultiComposedSumcheckProver::calculate_poly_sum(&multi_composed_vec_2);
        assert_eq!(sum_2, Fq::from(8));
    }

    #[test]
    fn test_multi_composed_sumcheck_proof() {
        let poly1 = MLE::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)]);
        let poly2 = MLE::new(vec![Fq::from(0), Fq::from(3), Fq::from(0), Fq::from(3)]);

        let composed_1 = ComposedMLE::new(vec![poly1]);
        let composed_2 = ComposedMLE::new(vec![poly2]);

        let multi_composed = vec![composed_1, composed_2];

        let sum = MultiComposedSumcheckProver::calculate_poly_sum(&multi_composed);
        let (proof, _) = MultiComposedSumcheckProver::prove(&multi_composed, &sum).unwrap();
        let verify = MultiComposedSumcheckVerifier::verify(&multi_composed, &proof).unwrap();
        assert!(verify);
    }

    #[test]
    fn test_multi_composed_sumcheck_proof_1() {
        let poly1 = MLE::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)]);
        let poly2 = MLE::new(vec![Fq::from(0), Fq::from(3), Fq::from(0), Fq::from(3)]);

        let composed_1 = ComposedMLE::new(vec![poly1]);
        let composed_2 = ComposedMLE::new(vec![poly2.clone()]);
        let composed_3 = ComposedMLE::new(vec![poly2]);

        let multi_composed = vec![composed_1, composed_2, composed_3];
        let sum = MultiComposedSumcheckProver::calculate_poly_sum(&multi_composed);
        let (proof, _) = MultiComposedSumcheckProver::prove(&multi_composed, &sum).unwrap();
        let verify = MultiComposedSumcheckVerifier::verify(&multi_composed, &proof).unwrap();
        assert!(verify);
    }

    #[test]
    fn test_multi_composed_sum_check_proof_2() {
        let poly1 = MLE::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)]);
        let poly2 = MLE::new(vec![Fq::from(0), Fq::from(3), Fq::from(0), Fq::from(3)]);

        let composed_1 = ComposedMLE::new(vec![poly1.clone(), poly2.clone()]);
        let composed_2 = ComposedMLE::new(vec![poly2.clone(), poly1.clone()]);

        let multi_composed = vec![composed_1, composed_2];
        let sum = MultiComposedSumcheckProver::calculate_poly_sum(&multi_composed);
        let (proof, _) = MultiComposedSumcheckProver::prove(&multi_composed, &sum).unwrap();
        let verify = MultiComposedSumcheckVerifier::verify(&multi_composed, &proof).unwrap();
        assert!(verify);
    }

    #[test]
    fn test_multi_composed_sum_check_proof_2_on_gkr_example() {
        // f(a,b,c) = 2abc + 3b + 4
        let add_i = MLE::<Fq>::new(vec![
            Fq::from(4),
            Fq::from(4),
            Fq::from(7),
            Fq::from(7),
            Fq::from(4),
            Fq::from(4),
            Fq::from(7),
            Fq::from(9),
        ]);
        // f(b) = 4b
        let w_b = MLE::<Fq>::new(vec![Fq::from(0), Fq::from(4)]);
        // f(c) = 3c
        let w_c = MLE::<Fq>::new(vec![Fq::from(0), Fq::from(3)]);
        // f(a,b,c) = 2ab + bc + 3
        let mul_i = MLE::<Fq>::new(vec![
            Fq::from(3),
            Fq::from(3),
            Fq::from(3),
            Fq::from(4),
            Fq::from(3),
            Fq::from(3),
            Fq::from(5),
            Fq::from(6),
        ]);

        let lhs_poly = ComposedMLE::new(vec![
            add_i.partial_evaluation(&Fq::from(2), &0),
            w_b.add_distinct(&w_c),
        ]);
        let rhs_poly = ComposedMLE::new(vec![
            mul_i.partial_evaluation(&Fq::from(2), &0),
            w_b.mul_distinct(&w_c),
        ]);

        let multi_composed = vec![lhs_poly, rhs_poly];
        let sum = MultiComposedSumcheckProver::calculate_poly_sum(&multi_composed);

        let (proof, _) = MultiComposedSumcheckProver::prove(&multi_composed, &sum).unwrap();
        let verify = MultiComposedSumcheckVerifier::verify(&multi_composed, &proof).unwrap();
        assert!(verify);
    }
}
