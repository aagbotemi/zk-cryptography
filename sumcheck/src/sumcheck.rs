use crate::utils::{convert_field_to_byte, skip_first_and_sum_all};
use ark_ff::PrimeField;
use fiat_shamir::{fiat_shamir::FiatShamirTranscript, interface::FiatShamirTranscriptTrait};
use polynomial::{interface::MLETrait, MLE};

#[derive(Clone, Debug, Default)]
pub struct Sumcheck<F: PrimeField> {
    poly: MLE<F>,
    sum: F,
}

#[derive(Clone, Debug, Default)]
pub struct SumcheckProof<F: PrimeField> {
    poly: MLE<F>,
    sum: F,
    univariate_poly: Vec<MLE<F>>,
}

impl<F: PrimeField> Sumcheck<F> {
    pub fn new(poly: MLE<F>) -> Self {
        Sumcheck {
            poly,
            sum: Default::default(),
        }
    }

    pub fn calculate_poly_sum(&mut self) {
        self.sum = self.poly.evaluations.iter().sum();
    }

    pub fn prove(&mut self) -> (SumcheckProof<F>, Vec<F>) {
        let mut uni_polys = vec![];

        // send sum as bytes to the transcript
        let mut transcript = FiatShamirTranscript::new();
        let poly_sum_bytes = convert_field_to_byte(&self.sum);
        transcript.commit(&poly_sum_bytes);

        let mut challenges: Vec<F> = vec![];
        let mut current_poly = self.poly.relabel();

        for _ in 0..self.poly.n_vars {
            let uni_poly = skip_first_and_sum_all(current_poly.clone());
            uni_polys.push(uni_poly.clone());

            transcript.commit(&uni_poly.evaluations_to_bytes());
            //get the random r
            let random_r: F = transcript.evaluate_challenge_into_field::<F>();
            challenges.push(random_r);

            // update polynomial
            current_poly = current_poly.partial_evaluation(random_r, 0);
        }

        (
            SumcheckProof {
                poly: self.poly.clone(),
                sum: self.sum,
                univariate_poly: uni_polys,
            },
            challenges,
        )
    }

    pub fn verify(&mut self, proof: &SumcheckProof<F>) -> bool {
        // send sum as bytes to the transcript
        let mut transcript = FiatShamirTranscript::new();
        let poly_sum_bytes = convert_field_to_byte(&proof.sum);
        transcript.commit(&poly_sum_bytes);

        let mut claimed_sum = proof.sum;
        let mut challenges: Vec<F> = vec![];

        let univariate_poly = proof.univariate_poly.clone();
        for i in 0..proof.poly.n_vars {
            let uni_poly = &univariate_poly[i];

            // Check if the claimed sum matches the evaluation at 0 and 1
            let eval_p0_p1 =
                uni_poly.evaluation(&vec![F::zero()]) + uni_poly.evaluation(&vec![F::one()]);
            if eval_p0_p1 != claimed_sum {
                return false;
            }

            // Commit the univariate polynomial to the transcript
            transcript.commit(&uni_poly.evaluations_to_bytes());

            // Generate the challenge for this round
            let challenge: F = transcript.evaluate_challenge_into_field::<F>();
            challenges.push(challenge);

            // update the sum
            claimed_sum = uni_poly.evaluation(&vec![challenge]);
        }

        proof.poly.evaluation(challenges.as_slice()) == claimed_sum
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::MontConfig;
    use ark_ff::{Fp64, MontBackend};

    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    struct FqConfig;
    type Fq = Fp64<MontBackend<FqConfig, 1>>;

    #[test]
    fn test_sum_calculation() {
        let poly = MLE::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(2),
            Fq::from(2),
            Fq::from(2),
            Fq::from(4),
        ]);
        let mut prover = Sumcheck::new(poly);
        prover.calculate_poly_sum();
        assert_eq!(prover.sum, Fq::from(12));
    }

    #[test]
    fn test_sum_check_proof() {
        let poly = MLE::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(7),
            Fq::from(3),
            Fq::from(3),
            Fq::from(6),
            Fq::from(11),
        ]);
        let mut sumcheck = Sumcheck::new(poly);
        sumcheck.calculate_poly_sum();
        let proof = sumcheck.prove();
        let verifer = sumcheck.verify(&proof.0);

        assert_eq!(verifer, true);
    }

    #[test]
    fn test_sum_check_proof_2() {
        let poly = MLE::new(vec![
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
        let mut sumcheck = Sumcheck::new(poly);
        sumcheck.calculate_poly_sum();
        let proof = sumcheck.prove();
        let verifer = sumcheck.verify(&proof.0);

        assert_eq!(verifer, true);
    }

    #[test]
    fn test_sum_check_proof_3() {
        let poly = MLE::new(vec![
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
        let mut sumcheck = Sumcheck::new(poly);
        sumcheck.calculate_poly_sum();
        let proof = sumcheck.prove();
        let verifer = sumcheck.verify(&proof.0);

        assert_eq!(verifer, true);
    }
}