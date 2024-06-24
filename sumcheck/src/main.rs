use ark_ff::{BigInteger, PrimeField};
use fiat_shamir::{fiat_shamir::FiatShamirTranscript, interface::FiatShamirTranscriptTrait};
use polynomial::{interface::MLETrait, UnivariatePolynomial, MLE};

#[derive(Clone, Debug, Default)]
pub struct SumcheckProver<F: PrimeField> {
    poly: MLE<F>,
    init_poly: MLE<F>,
    // uni_poly: MLE<F>,
    sum: F,
}

#[derive(Clone, Debug, Default)]
pub struct SumcheckProof<F: PrimeField> {
    init_poly: MLE<F>,
    poly: MLE<F>,
    sum: F,
    univariate_poly: Vec<MLE<F>>,
}

#[derive(Debug)]
struct SumCheckVerifier<F: PrimeField> {
    pub poly: MLE<F>,
}

impl<F: PrimeField> SumcheckProver<F> {
    pub fn new(poly: MLE<F>) -> Self {
        SumcheckProver {
            poly,
            init_poly: Default::default(),
            sum: Default::default(),
            // transcript: Default::default(),
        }
    }

    pub fn calculate_poly_sum(&mut self) {
        self.sum = self.poly.evaluations.iter().sum();
    }

    pub fn skip_first_and_sum_all(&mut self) -> MLE<F> {
        let rounds = self.poly.n_vars - 1;
        let bh = boolean_hypercube::<F>(rounds);

        let mut bh_sum = MLE::<F>::additive_identity(1);
        for bh_i in bh {
            let mut sum_except_first = self.poly.clone();
            for bh_i_i in bh_i {
                sum_except_first = sum_except_first.partial_evaluation(bh_i_i, 1);
            }
            bh_sum += sum_except_first;
        }
        // println!("bh_sum={bh_sum:?}");

        // self.init_poly = bh_sum.clone();
        bh_sum
    }

    pub fn prove(&mut self) -> (SumcheckProof<F>, Vec<F>) {
        // let mut uni_polys = MLE::<F>::additive_identity(1);
        let mut uni_polys = vec![];

        // send sum as bytes to the transcript
        let mut transcrpt = FiatShamirTranscript::new();
        let poly_sum_bytes = convert_field_to_byte(&self.sum);
        transcrpt.commit(&poly_sum_bytes);

        // println!("transcrpt-sum={:?}", transcrpt);

        let mut challenges: Vec<F> = vec![];
        let mut current_poly = self.poly.relabel();

        for i in 0..self.poly.n_vars {
            let poly_sum_bytes_1 = self.skip_first_and_sum_all();
            uni_polys.push(poly_sum_bytes_1.clone());
            // println!("poly_sum_bytes_1={poly_sum_bytes_1:?}");

            transcrpt.commit(&poly_sum_bytes_1.evaluations_to_bytes());
            //get the random r
            let random_r: F = transcrpt.evaluate_challenge_into_field::<F>();
            challenges.push(random_r);

            // update polynomial
            current_poly = current_poly.partial_evaluation(random_r, 0);
        }

        (
            SumcheckProof {
                init_poly: self.poly.clone(),
                poly: current_poly,
                sum: self.sum,
                univariate_poly: uni_polys,
            },
            challenges,
        )
    }
}

fn boolean_hypercube<F: PrimeField>(n: usize) -> Vec<Vec<F>> {
    let mut hypercube = Vec::new();

    for i in 0..(1 << n) {
        let mut vertex = Vec::new();
        for j in 0..n {
            if (i & (1 << j)) != 0 {
                vertex.push(F::one());
            } else {
                vertex.push(F::zero());
            }
        }
        hypercube.push(vertex);
    }

    hypercube
}

fn convert_field_to_byte<F: PrimeField>(element: &F) -> Vec<u8> {
    element.into_bigint().to_bytes_be()
}

impl<F: PrimeField> SumCheckVerifier<F> {
    pub fn verify(proof: &SumcheckProof<F>) -> bool {
        // send sum as bytes to the transcript
        let mut transcrpt = FiatShamirTranscript::new();
        let poly_sum_bytes = convert_field_to_byte(&proof.sum);
        transcrpt.commit(&poly_sum_bytes);

        let mut claimed_sum = proof.sum;
        let mut challenges: Vec<F> = vec![];

        let initial_poly = proof.univariate_poly.clone();

        for i in 0..proof.init_poly.n_vars {
            // Verify the univariate polynomial at each round
            let uni_poly = &initial_poly[i];

            // Check if the claimed sum matches the evaluation at 0 and 1
            let eval_p0_p1 =
                uni_poly.evaluation(&vec![F::zero()]) + uni_poly.evaluation(&vec![F::one()]);

            if eval_p0_p1 != proof.sum {
                return false;
            }

            // Commit the univariate polynomial to the transcript
            transcrpt.commit(&uni_poly.evaluations_to_bytes());

            // Generate the challenge for this round
            let challenge: F = transcrpt.evaluate_challenge_into_field::<F>();
            challenges.push(challenge);

            // Update the claimed sum for the next round
            claimed_sum = uni_poly.evaluation(&vec![challenge]);
        }
        
        println!("proof.init_poly={:?}", proof.init_poly.evaluation(&vec![challenges[1]; 3]));
        println!("claimed_sum={:?}", claimed_sum);
        
        proof.init_poly.evaluation(challenges.as_slice()) == claimed_sum
    }
}

fn main() {
    use ark_ff::MontConfig;
    use ark_ff::{Fp64, MontBackend};

    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    struct FqConfig;
    type Fq = Fp64<MontBackend<FqConfig, 1>>;

    // let transcript = FiatShamirTranscript::new();
    // println!("transcript={:?}", transcript);

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
    // println!("polypolypoly={:?}", poly);
    let mut prover = SumcheckProver::new(poly);
    // println!("proverprover={:?}", prover);
    let sum = prover.calculate_poly_sum();
    // let sum = prover.calculate_poly_sum(&vec![Fq::from(5_u8), Fq::from(6_u8), Fq::from(6_u8)]);

    // println!("proverprover={:?}", prover);
    // let mut prover = SumcheckProver::new(poly);
    let proof = prover.prove();
    // println!("profprofprof={:?}", proof);

    let verifier = SumCheckVerifier::verify(&proof.0);
    println!("verifier={:?}", verifier);
}
