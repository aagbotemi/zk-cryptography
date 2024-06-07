use ark_ff::{BigInteger, PrimeField};
use fiat_shamir::{fiat_shamir::FiatShamirTranscript, interface::FiatShamirTranscriptTrait};
use polynomial::{interface::MLETrait, MLE};

#[derive(Clone, Debug, Default)]
pub struct SumcheckProver<F: PrimeField> {
    poly: MLE<F>,
    init_poly: MLE<F>,
    sum: F,
    transcript: FiatShamirTranscript,
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
            transcript: Default::default(),
        }
    }

    pub fn calculate_poly_sum(&mut self) {
        self.sum = self.poly.evaluations.iter().sum();
    }

    pub fn skip_first_and_sum_all(&mut self) -> Vec<u8> {
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
        println!("bh_sum={bh_sum:?}");

        self.init_poly = bh_sum.clone();
        bh_sum.evaluations_to_bytes()
    }

    pub fn prove(&mut self) -> (SumcheckProver<F>, Vec<F>) {
        // send sum as bytes to the transcript
        let mut transcrpt = self.transcript.clone();
        let poly_sum_bytes = convert_field_to_byte(&self.sum);
        transcrpt.commit(&poly_sum_bytes);

        let mut challenges: Vec<F> = vec![];
        let mut current_poly = self.poly.relabel();

        for i in 0..self.poly.n_vars {
            let poly_sum_bytes_1 = self.skip_first_and_sum_all();

            transcrpt.commit(&poly_sum_bytes_1);
            //get the random r
            let random_r: F = transcrpt.evaluate_challenge_into_field::<F>();
            challenges.push(random_r);
            println!("random={random_r:?}");

            // update polynomial
            current_poly = current_poly.partial_evaluation(random_r, 0);
        }

        (
            SumcheckProver {
                poly: current_poly,
                init_poly: self.init_poly.clone(),
                sum: self.sum,
                transcript: transcrpt,
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
    pub fn verify(proof: &SumcheckProver<F>) -> bool {
        // pub fn verify(&self, proof: &SumcheckProver<F>) -> bool {
        // println!("proofproof={:?}", proof);
        let first_sum = proof.init_poly.evaluation(&vec![F::zero()])
            + proof.init_poly.evaluation(&vec![F::one()]);

        if first_sum != proof.sum {
            return false;
        }
        true
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
    // println!("verifier={:?}", verifier);

    // fn boolean_hypercube_2(n: usize) -> Vec<Vec<usize>> {
    //     if n == 0 {
    //         vec![vec![]] // Base case: return a single empty vector
    //     } else {
    //         // Recursively build the hypercube for n - 1
    //         let mut smaller_hypercube = boolean_hypercube_2(n - 1);
    //         // println!("smaller_hypercube={:?}", smaller_hypercube);
    //         let mut hypercube = smaller_hypercube.clone();

    //         // Append 0 to the front of each vector in the smaller hypercube
    //         for vertex in smaller_hypercube.iter_mut() {
    //             vertex.insert(0, 0);
    //         }

    //         // Append 1 to the front of each vector in the cloned smaller hypercube
    //         for vertex in hypercube.iter_mut() {
    //             vertex.insert(0, 1);
    //         }

    //         // Combine both halves to form the hypercube
    //         smaller_hypercube.extend(hypercube);
    //         smaller_hypercube
    //     }
    // }

    // let now_2 = Instant::now();
    // let bh_2 = boolean_hypercube_2(3);
    // println!("Time taken Hypercube 2: {:?}", now_2.elapsed());
    // println!("boolean_hypercube_2_time={:?}", bh_2);

    // let now = Instant::now();
    // let bh = boolean_hypercube::<Fq>(3);
    // println!("Time taken Hypercube 1: {:?}", now.elapsed());
    // println!("boolean_hypercube_1_time={:?}", bh);
}
