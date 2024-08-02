use crate::utils::{convert_round_poly_to_uni_poly_format, vec_to_bytes};
use ark_ff::PrimeField;
use fiat_shamir::{fiat_shamir::FiatShamirTranscript, interface::FiatShamirTranscriptTrait};
use polynomial::{
    interface::ComposedMLETrait, ComposedMLE, UnivariatePolynomial, UnivariatePolynomialTrait,
};

#[derive(Debug, Clone)]
pub struct Sumcheck<F: PrimeField> {
    pub poly: ComposedMLE<F>,
    pub sum: F,
}

pub struct SumcheckProof<F: PrimeField> {
    poly: ComposedMLE<F>,
    sum: F,
    round_polys: Vec<Vec<F>>,
}

impl<F: PrimeField> Sumcheck<F> {
    pub fn new(poly: ComposedMLE<F>) -> Self {
        Sumcheck {
            poly,
            sum: Default::default(),
        }
    }

    pub fn calculate_poly_sum(&mut self) {
        self.sum = self.poly.element_wise_product().iter().sum()
    }

    pub fn prove(&self) -> (SumcheckProof<F>, Vec<F>) {
        // send sum as bytes to the transcript
        let mut transcript = FiatShamirTranscript::new();

        let mut current_poly: ComposedMLE<F> = self.poly.clone();
        let mut round_polys: Vec<Vec<F>> = vec![];
        let mut challenges: Vec<F> = vec![];

        for _ in 0..self.poly.n_vars() {
            let mut round_poly: Vec<F> = vec![];
            for i in 0..=current_poly.max_degree() {
                let round: F = current_poly
                    .partial_evaluation(F::from(i as u32), 0)
                    .element_wise_product()
                    .iter()
                    .sum::<F>();

                round_poly.push(round);
            }

            transcript.commit(&vec_to_bytes(&round_poly));
            //get the random r
            let random_r: F = transcript.evaluate_challenge_into_field::<F>();
            challenges.push(random_r);
            round_polys.push(round_poly);

            current_poly = current_poly.partial_evaluation(random_r, 0);
        }

        (
            SumcheckProof {
                poly: self.poly.clone(),
                sum: self.sum,
                round_polys,
            },
            challenges,
        )
    }

    pub fn verify(&self, proof: &SumcheckProof<F>) -> bool {
        let mut transcript = FiatShamirTranscript::new();

        let mut claimed_sum = proof.sum;
        let mut challenges: Vec<F> = vec![];

        for round_poly in proof.round_polys.iter() {
            transcript.commit(&vec_to_bytes(&round_poly));
            // genrate the challenge for this round
            let challenge: F = transcript.evaluate_challenge_into_field::<F>();
            challenges.push(challenge);

            let round_polys_uni: Vec<(F, F)> = convert_round_poly_to_uni_poly_format(&round_poly);
            let uni_poly: UnivariatePolynomial<F> =
                UnivariatePolynomial::interpolation(&round_polys_uni);

            let eval_p0_p1 = uni_poly.evaluate(F::zero()) + uni_poly.evaluate(F::one());
            if claimed_sum != eval_p0_p1 {
                return false;
            }

            // update the sum
            claimed_sum = uni_poly.evaluate(challenge);
        }

        proof.poly.evaluation(challenges.as_slice()) == claimed_sum
    }
}

#[cfg(test)]
mod tests {
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
        let composedmle = ComposedMLE::new(vec![mle1, mle2]);
        let mut prover = Sumcheck::new(composedmle);
        prover.calculate_poly_sum();
        assert_eq!(prover.sum, Fq::from(3));

        let mle1 = MLE::new(vec![Fq::from(3), Fq::from(3), Fq::from(5), Fq::from(5)]);
        let mle2 = MLE::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(1)]);
        let composedmle = ComposedMLE::new(vec![mle1, mle2]);
        let mut prover = Sumcheck::new(composedmle);
        prover.calculate_poly_sum(); // 5
        assert_eq!(prover.sum, Fq::from(5));

        let mle1 = MLE::new(vec![Fq::from(0), Fq::from(1), Fq::from(2), Fq::from(3)]);
        let composedmle = ComposedMLE::new(vec![mle1]);
        let mut prover = Sumcheck::new(composedmle);
        prover.calculate_poly_sum(); // 6
        assert_eq!(prover.sum, Fq::from(6));

        let mle1 = MLE::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(2),
            Fq::from(2),
            Fq::from(2),
            Fq::from(4),
        ]);
        let composedmle = ComposedMLE::new(vec![mle1]);
        let mut prover = Sumcheck::new(composedmle);
        prover.calculate_poly_sum(); // 12
        assert_eq!(prover.sum, Fq::from(12));
    }

    #[test]
    fn test_sum_check_proof() {
        // 2(a^2)b + 3ab
        // (2a + 3)(ab)
        // (2a + 0b + 3)(ab)
        // 00 - 3  | 00 - 0
        // 01 - 3  | 01 - 0
        // 10 - 5  | 10 - 0
        // 11 - 5  | 11 - 1
        let mle1 = MLE::new(vec![Fq::from(3), Fq::from(3), Fq::from(5), Fq::from(5)]);
        let mle2 = MLE::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(1)]);
        let composedmle = ComposedMLE::new(vec![mle1, mle2]);
        let mut sumcheck = Sumcheck::new(composedmle);
        sumcheck.calculate_poly_sum();
        let (proof, _challenges) = &sumcheck.prove();
        let verifer: bool = sumcheck.verify(&proof);
        assert_eq!(verifer, true);
    }

    #[test]
    fn test_sum_check_proof1() {
        let mle = MLE::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(7),
            Fq::from(3),
            Fq::from(3),
            Fq::from(6),
            Fq::from(11),
        ]);
        let composedmle = ComposedMLE::new(vec![mle]);
        let mut sumcheck = Sumcheck::new(composedmle);
        sumcheck.calculate_poly_sum();
        let (proof, _challenges) = &sumcheck.prove();
        let verifer: bool = sumcheck.verify(&proof);
        assert_eq!(verifer, true);
    }

    #[test]
    fn test_sum_check_proof_2() {
        let mle = MLE::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(1),
            Fq::from(1),
            Fq::from(1),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
        ]);
        let composedmle = ComposedMLE::new(vec![mle]);
        let mut sumcheck = Sumcheck::new(composedmle);
        sumcheck.calculate_poly_sum();
        let proof = sumcheck.prove();
        let verifer = sumcheck.verify(&proof.0);

        assert_eq!(verifer, true);
    }

    #[test]
    fn test_sum_check_proof_3() {
        let mle = MLE::new(vec![
            Fq::from(1),
            Fq::from(3),
            Fq::from(5),
            Fq::from(7),
            Fq::from(2),
            Fq::from(4),
            Fq::from(6),
            Fq::from(8),
            Fq::from(3),
            Fq::from(5),
            Fq::from(7),
            Fq::from(9),
            Fq::from(4),
            Fq::from(6),
            Fq::from(8),
            Fq::from(10),
        ]);
        let composedmle = ComposedMLE::new(vec![mle]);
        let mut sumcheck = Sumcheck::new(composedmle);
        sumcheck.calculate_poly_sum();
        let proof = sumcheck.prove();
        let verifer = sumcheck.verify(&proof.0);

        assert_eq!(verifer, true);
    }
}
