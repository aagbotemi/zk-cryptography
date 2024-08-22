use crate::utils::{convert_round_poly_to_uni_poly_format, vec_to_bytes};
use ark_ff::PrimeField;
use fiat_shamir::{fiat_shamir::FiatShamirTranscript, interface::FiatShamirTranscriptTrait};
use polynomial::{
    interface::ComposedMultilinearTrait, ComposedMultilinear, MultilinearTrait, UnivariatePolynomial, UnivariatePolynomialTrait
};

#[derive(Debug, Clone)]
pub struct ComposedSumcheck<F: PrimeField> {
    pub poly: ComposedMultilinear<F>,
    pub sum: F,
}

pub struct ComposedSumcheckProof<F: PrimeField> {
    poly: ComposedMultilinear<F>,
    round_polys: Vec<Vec<F>>,
}

impl<F: PrimeField> ComposedSumcheck<F> {
    pub fn new(poly: ComposedMultilinear<F>) -> Self {
        ComposedSumcheck {
            poly,
            sum: Default::default(),
        }
    }

    pub fn calculate_poly_sum(poly: &ComposedMultilinear<F>) -> F {
        poly.element_wise_product().iter().sum()
    }

    pub fn prove(&self) -> (ComposedSumcheckProof<F>, Vec<F>) {
        let mut transcript = FiatShamirTranscript::new();

        let mut current_poly: ComposedMultilinear<F> = self.poly.clone();
        let mut round_polys: Vec<Vec<F>> = vec![];
        let mut challenges: Vec<F> = vec![];

        for _ in 0..self.poly.n_vars() {
            let mut round_poly: Vec<F> = vec![];
            for i in 0..=current_poly.max_degree() {
                let round: F = current_poly
                    .partial_evaluation(&F::from(i as u32), &0)
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

            current_poly = current_poly.partial_evaluation(&random_r, &0);
        }

        (
            ComposedSumcheckProof {
                poly: self.poly.clone(),
                round_polys,
            },
            challenges,
        )
    }

    pub fn verify(&self, proof: &ComposedSumcheckProof<F>, sum: F) -> bool {
        let mut transcript = FiatShamirTranscript::new();

        let mut claimed_sum = sum;
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
    use polynomial::Multilinear;

    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    struct FqConfig;
    type Fq = Fp64<MontBackend<FqConfig, 1>>;

    #[test]
    fn test_sum_calculation() {
        let poly1 = Multilinear::new(vec![Fq::from(0), Fq::from(1), Fq::from(2), Fq::from(3)]);
        let poly2 = Multilinear::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(1)]);
        let composedpoly1 = ComposedMultilinear::new(vec![poly1, poly2]);
        let sum = ComposedSumcheck::calculate_poly_sum(&composedpoly1);
        assert_eq!(sum, Fq::from(3));

        let poly1 = Multilinear::new(vec![Fq::from(3), Fq::from(3), Fq::from(5), Fq::from(5)]);
        let poly2 = Multilinear::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(1)]);
        let composedpoly2 = ComposedMultilinear::new(vec![poly1, poly2]);
        let sum2 = ComposedSumcheck::calculate_poly_sum(&composedpoly2);
        assert_eq!(sum2, Fq::from(5));

        let poly1 = Multilinear::new(vec![Fq::from(0), Fq::from(1), Fq::from(2), Fq::from(3)]);
        let composedpoly3 = ComposedMultilinear::new(vec![poly1]);
        let sum3 = ComposedSumcheck::calculate_poly_sum(&composedpoly3);
        assert_eq!(sum3, Fq::from(6));

        let poly1 = Multilinear::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(2),
            Fq::from(2),
            Fq::from(2),
            Fq::from(4),
        ]);
        let composedpoly4 = ComposedMultilinear::new(vec![poly1]);
        let sum4 = ComposedSumcheck::calculate_poly_sum(&composedpoly4);
        assert_eq!(sum4, Fq::from(12));
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
        let poly1 = Multilinear::new(vec![Fq::from(3), Fq::from(3), Fq::from(5), Fq::from(5)]);
        let poly2 = Multilinear::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(1)]);
        let composedpoly = ComposedMultilinear::new(vec![poly1, poly2]);
        let sumcheck = ComposedSumcheck::new(composedpoly);
        let (proof, _challenges) = &sumcheck.prove();
        let sum = ComposedSumcheck::calculate_poly_sum(&proof.poly);
        let verifer: bool = sumcheck.verify(&proof, sum);
        assert_eq!(verifer, true);
    }

    #[test]
    fn test_sum_check_proof1() {
        let mle = Multilinear::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(7),
            Fq::from(3),
            Fq::from(3),
            Fq::from(6),
            Fq::from(11),
        ]);
        let composedpoly = ComposedMultilinear::new(vec![mle]);
        let sumcheck = ComposedSumcheck::new(composedpoly);
        let (proof, _challenges) = &sumcheck.prove();
        let sum = ComposedSumcheck::calculate_poly_sum(&proof.poly);
        let verifer: bool = sumcheck.verify(&proof, sum);
        assert_eq!(verifer, true);
    }

    #[test]
    fn test_sum_check_proof_2() {
        let mle = Multilinear::new(vec![
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
        let composedpoly = ComposedMultilinear::new(vec![mle]);
        let sumcheck = ComposedSumcheck::new(composedpoly);
        let proof = sumcheck.prove();
        let sum = ComposedSumcheck::calculate_poly_sum(&proof.0.poly);
        let verifer = sumcheck.verify(&proof.0, sum);

        assert_eq!(verifer, true);
    }

    #[test]
    fn test_sum_check_proof_3() {
        let mle = Multilinear::new(vec![
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
        let composedpoly = ComposedMultilinear::new(vec![mle]);
        let sumcheck = ComposedSumcheck::new(composedpoly);
        let proof = sumcheck.prove();
        let sum = ComposedSumcheck::calculate_poly_sum(&proof.0.poly);
        let verifer = sumcheck.verify(&proof.0, sum);

        assert_eq!(verifer, true);
    }
}
