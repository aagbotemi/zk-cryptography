use super::composed_sumcheck::ComposedSumcheck;
use crate::utils::convert_round_poly_to_uni_poly_format;
use ark_ff::PrimeField;
use fiat_shamir::{fiat_shamir::FiatShamirTranscript, interface::FiatShamirTranscriptTrait};
use polynomial::{
    interface::ComposedMLETrait, ComposedMLE, UnivariatePolynomial, UnivariatePolynomialTrait,
};

#[derive(Debug, Clone)]
pub struct MultiComposedSumcheck<F: PrimeField> {
    pub poly: Vec<ComposedMLE<F>>,
    pub sum: F,
}

#[derive(Debug)]
pub struct SumcheckProof<F: PrimeField> {
    poly: Vec<ComposedMLE<F>>,
    round_polys: Vec<UnivariatePolynomial<F>>,
}

impl<F: PrimeField> MultiComposedSumcheck<F> {
    pub fn new(poly: Vec<ComposedMLE<F>>) -> Self {
        MultiComposedSumcheck {
            poly,
            sum: Default::default(),
        }
    }

    pub fn calculate_poly_sum(poly: &Vec<ComposedMLE<F>>) -> F {
        let mut sum = F::zero();

        for p in poly {
            sum += ComposedSumcheck::calculate_poly_sum(p);
        }

        sum
    }

    pub fn prove(&self) -> (SumcheckProof<F>, Vec<F>) {
        let mut transcript = FiatShamirTranscript::new();

        let mut current_poly = self.poly.clone();
        let mut round_polys = vec![];
        let mut challenges: Vec<F> = vec![];

        for _ in 0..self.poly[0].n_vars() {
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

        (
            SumcheckProof {
                poly: self.poly.clone(),
                round_polys,
            },
            challenges,
        )
    }

    pub fn verify(&self, proof: &SumcheckProof<F>, sum: F) -> bool {
        let mut transcript = FiatShamirTranscript::new();

        let mut claimed_sum = sum;
        let mut challenges: Vec<F> = vec![];

        for round_poly in proof.round_polys.iter() {
            transcript.commit(&round_poly.to_bytes());
            // genrate the challenge for this round
            let challenge: F = transcript.evaluate_challenge_into_field::<F>();
            challenges.push(challenge);

            let eval_p0_p1 = round_poly.evaluate(F::zero()) + round_poly.evaluate(F::one());
            if claimed_sum != eval_p0_p1 {
                return false;
            }

            // update the sum
            claimed_sum = round_poly.evaluate(challenge);
        }

        let mut poly_pe_sum = F::zero();
        for round_poly in proof.poly.iter() {
            poly_pe_sum += round_poly.evaluation(&challenges.as_slice())
        }

        poly_pe_sum == claimed_sum
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
        let sum_1 = MultiComposedSumcheck::calculate_poly_sum(&multi_composed_vec_1);
        assert_eq!(sum_1, Fq::from(7));

        let mle3 = MLE::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)]);
        let mle4 = MLE::new(vec![Fq::from(0), Fq::from(3), Fq::from(0), Fq::from(3)]);
        let composedmle3 = ComposedMLE::new(vec![mle3]);
        let composedmle4 = ComposedMLE::new(vec![mle4]);

        let multi_composed_vec_2 = vec![composedmle3, composedmle4];
        let sum_2 = MultiComposedSumcheck::calculate_poly_sum(&multi_composed_vec_2);
        assert_eq!(sum_2, Fq::from(8));
    }

    #[test]
    fn test_multi_composed_sumcheck_proof() {
        let poly1 = MLE::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)]);
        let poly2 = MLE::new(vec![Fq::from(0), Fq::from(3), Fq::from(0), Fq::from(3)]);

        let composed_1 = ComposedMLE::new(vec![poly1]);
        let composed_2 = ComposedMLE::new(vec![poly2]);

        let multi_composed = vec![composed_1, composed_2];

        let sum = MultiComposedSumcheck::calculate_poly_sum(&multi_composed);
        let sumcheck = MultiComposedSumcheck::new(multi_composed);

        let (proof, _) = sumcheck.prove();
        let verify = sumcheck.verify(&proof, sum);
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
        let sum = MultiComposedSumcheck::calculate_poly_sum(&multi_composed);

        let sumcheck = MultiComposedSumcheck::new(multi_composed);

        let (proof, _) = sumcheck.prove();
        let verify = sumcheck.verify(&proof, sum);
        assert!(verify);
    }

    #[test]
    fn test_multi_composed_sum_check_proof_2() {
        let poly1 = MLE::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)]);
        let poly2 = MLE::new(vec![Fq::from(0), Fq::from(3), Fq::from(0), Fq::from(3)]);

        let composed_1 = ComposedMLE::new(vec![poly1.clone(), poly2.clone()]);
        let composed_2 = ComposedMLE::new(vec![poly2.clone(), poly1.clone()]);

        let multi_composed = vec![composed_1, composed_2];
        let sum = MultiComposedSumcheck::calculate_poly_sum(&multi_composed);

        let sumcheck = MultiComposedSumcheck::new(multi_composed);

        let (proof, _) = sumcheck.prove();
        let verify = sumcheck.verify(&proof, sum);
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
            add_i.partial_evaluation(Fq::from(2), 0),
            w_b.add_distinct(&w_c),
        ]);
        let rhs_poly = ComposedMLE::new(vec![
            mul_i.partial_evaluation(Fq::from(2), 0),
            w_b.mul_distinct(&w_c),
        ]);

        let multi_composed = vec![lhs_poly, rhs_poly];
        let sum = MultiComposedSumcheck::calculate_poly_sum(&multi_composed);

        let sumcheck = MultiComposedSumcheck::new(multi_composed);

        let (proof, _) = sumcheck.prove();
        let verify = sumcheck.verify(&proof, sum);
        assert!(verify);
    }
}
