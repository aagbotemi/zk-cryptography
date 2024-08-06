use crate::utils::convert_field_to_byte;
use ark_ff::PrimeField;
use fiat_shamir::{fiat_shamir::FiatShamirTranscript, interface::FiatShamirTranscriptTrait};
use polynomial::{interface::MLETrait, MLE};

pub struct Sumcheck<F: PrimeField> {
    poly: MLE<F>,
    sum: F,
}

pub struct SumcheckProof<'a, F: PrimeField> {
    poly: &'a MLE<F>,
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

    pub fn poly_sum(&mut self) {
        self.sum = self.poly.evaluations.iter().sum();
    }

    pub fn prove(&self) -> (SumcheckProof<F>, Vec<F>) {
        let mut uni_polys = vec![];

        // send sum as bytes to the transcript
        let mut transcript = FiatShamirTranscript::new();
        let poly_sum_bytes = convert_field_to_byte(&self.sum);
        transcript.commit(&poly_sum_bytes);

        let mut challenges: Vec<F> = vec![];
        let mut current_poly: MLE<F> = self.poly.clone();

        for _ in 0..self.poly.n_vars {
            let uni_poly = current_poly.split_poly_into_two_and_sum_each_part();
            transcript.commit(&uni_poly.to_bytes());
            uni_polys.push(uni_poly);

            //get the random r
            let random_r = transcript.evaluate_challenge_into_field::<F>();
            challenges.push(random_r);

            // update polynomial
            current_poly = current_poly.partial_evaluation(random_r, 0);
        }

        (
            SumcheckProof {
                poly: &self.poly,
                sum: self.sum,
                univariate_poly: uni_polys,
            },
            challenges,
        )
    }

    pub fn verify(&self, proof: &SumcheckProof<F>) -> bool {
        // send sum as bytes to the transcript
        let mut transcript = FiatShamirTranscript::new();
        let poly_sum_bytes = convert_field_to_byte(&proof.sum);
        transcript.commit(&poly_sum_bytes);

        let mut claimed_sum = proof.sum;
        let mut challenges: Vec<F> = vec![];

        let univariate_poly = &proof.univariate_poly;
        for i in 0..proof.poly.n_vars {
            let uni_poly = &univariate_poly[i];

            // Check if the claimed sum matches the evaluation at 0 and 1
            let eval_p0_p1 =
                uni_poly.evaluation(&vec![F::zero()]) + uni_poly.evaluation(&vec![F::one()]);
            if eval_p0_p1 != claimed_sum {
                return false;
            }

            // Commit the univariate polynomial to the transcript
            transcript.commit(&uni_poly.to_bytes());

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
        prover.poly_sum();
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
        sumcheck.poly_sum();
        let (proof, _challenges) = &sumcheck.prove();
        let verifer: bool = sumcheck.verify(&proof);

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
        sumcheck.poly_sum();
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
        sumcheck.poly_sum();
        let proof = sumcheck.prove();
        let verifer = sumcheck.verify(&proof.0);

        assert_eq!(verifer, true);
    }
}
