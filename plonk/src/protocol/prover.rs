use crate::compiler::{
    primitives::{CommonPreprocessedInput, Witness},
    utils::{root_of_unity, roots_of_unity},
};
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use kzg::{
    interface::UnivariateKZGInterface, trusted_setup::TrustedSetup, univariate_kzg::UnivariateKZG,
};
use polynomial::{
    univariate::{domain::Domain, evaluation::UnivariateEval},
    utils::generate_random_numbers,
    DenseUnivariatePolynomial, UnivariatePolynomialTrait,
};

use super::{
    primitives::{PlonkProof, PlonkProver, PlonkRoundTranscript, RandomNumbers, WitnessPolys},
    utils::{apply_w_to_polynomial, split_poly_in_3, zh_values},
};

impl<F: PrimeField, P: Pairing> PlonkProver<F, P> {
    pub fn new(
        preprocessed_input: CommonPreprocessedInput<F>,
        srs: TrustedSetup<P>,
        transcript: PlonkRoundTranscript<P>,
    ) -> Self {
        PlonkProver {
            preprocessed_input,
            srs,
            transcript,
            random_number: RandomNumbers::default(),
            witness_polys: WitnessPolys::default(),
        }
    }

    pub fn prove(&mut self, witness: &Witness<F>) -> PlonkProof<P, F> {
        // round 1
        let (as_commitment, bs_commitment, cs_commitment) = self.first_round(&witness);
        self.transcript
            .first_round(as_commitment, bs_commitment, cs_commitment);

        // round 2
        let accumulator_commitment = self.second_round(&witness);
        self.transcript.second_round::<F>(accumulator_commitment);

        // round 3
        let zh_accumulator_poly = self.witness_polys.zh_accumulator_poly.clone();
        let (t_low, t_mid, t_high) = self.third_round(&witness, &zh_accumulator_poly);
        self.transcript.third_round(t_low, t_mid, t_high);

        // round 4
        let (
            a_s_poly_zeta,
            b_s_poly_zeta,
            c_s_poly_zeta,
            sigma1_poly_zeta,
            sigma2_poly_zeta,
            w_accumulator_poly_zeta,
        ) = self.fourth_round();
        self.transcript.fourth_round(
            a_s_poly_zeta,
            b_s_poly_zeta,
            c_s_poly_zeta,
            sigma1_poly_zeta,
            sigma2_poly_zeta,
            w_accumulator_poly_zeta,
        );

        // round 5
        let (w_zeta_commitment, w_zeta_omega_commitment) = self.fifth_round(&witness);
        self.transcript
            .fifth_round(w_zeta_commitment, w_zeta_omega_commitment);

        let mu: F = self.transcript.challenge_round(b"mu");
        self.random_number.mu = mu;

        PlonkProof {
            as_commitment,
            bs_commitment,
            cs_commitment,
            accumulator_commitment,
            t_low,
            t_mid,
            t_high,
            a_s_poly_zeta,
            b_s_poly_zeta,
            c_s_poly_zeta,
            sigma1_poly_zeta,
            sigma2_poly_zeta,
            w_accumulator_poly_zeta,
            w_zeta_commitment,
            w_zeta_omega_commitment,
        }
    }

    pub fn first_round(&mut self, witness: &Witness<F>) -> (P::G1, P::G1, P::G1) {
        let rands = generate_random_numbers(6);

        let zh_poly: DenseUnivariatePolynomial<F> =
            DenseUnivariatePolynomial::new(zh_values(self.preprocessed_input.group_order as usize));

        let a_s = DenseUnivariatePolynomial::new(vec![rands[1], rands[0]]) * zh_poly.clone()
            + witness.a.to_coefficient_poly().clone();
        // pub s2_coeff: Option<DenseUnivariatePolynomial<F>>,

        let b_s = DenseUnivariatePolynomial::new(vec![rands[3], rands[2]]) * zh_poly.clone()
            + witness.b.to_coefficient_poly().clone();

        let c_s = DenseUnivariatePolynomial::new(vec![rands[5], rands[4]]) * zh_poly.clone()
            + witness.c.to_coefficient_poly().clone();

        // commit to the polynomials
        let as_commitment = UnivariateKZG::<P>::commitment(&a_s, &self.srs);
        let bs_commitment = UnivariateKZG::<P>::commitment(&b_s, &self.srs);
        let cs_commitment = UnivariateKZG::<P>::commitment(&c_s, &self.srs);

        // dbg!(&as_commitment, &bs_commitment, &cs_commitment);
        self.witness_polys.a_s = a_s;
        self.witness_polys.b_s = b_s;
        self.witness_polys.c_s = c_s;
        (as_commitment, bs_commitment, cs_commitment)
    }

    pub fn second_round(&mut self, witness: &Witness<F>) -> P::G1 {
        let group_order = self.preprocessed_input.group_order as usize;
        let roots_of_unity: Vec<F> = roots_of_unity(group_order as u64);
        let mut accumulator = vec![F::one(); group_order as usize];

        let beta = self.transcript.challenge_round(b"beta");
        let gamma = self.transcript.challenge_round(b"gamma");

        for i in 0..group_order {
            let acc = accumulator[accumulator.len() - 1]
                * ((witness.a.values[i] + (beta * roots_of_unity[i]) + gamma)
                    * (witness.b.values[i] + (beta * F::from(2u8) * roots_of_unity[i]) + gamma)
                    * (witness.c.values[i] + (beta * F::from(3u8) * roots_of_unity[i]) + gamma))
                / ((witness.a.values[i]
                    + (beta + self.preprocessed_input.sigma_1.values[i])
                    + gamma)
                    * (witness.b.values[i]
                        + (beta * self.preprocessed_input.sigma_2.values[i]) * gamma)
                    * witness.c.values[i]
                    + (beta * self.preprocessed_input.sigma_3.values[i]) * gamma);

            accumulator[i] = acc;
        }

        let rands = generate_random_numbers(3);

        let domain: Domain<F> = Domain::new(group_order);
        let accumulator_poly = UnivariateEval::interpolate(accumulator, domain);

        let zh_poly = DenseUnivariatePolynomial::new(zh_values(group_order));

        let zh_blinding_factor = DenseUnivariatePolynomial::new(vec![rands[0], rands[1], rands[2]]);
        let zh_blinding_accumulator_poly =
            accumulator_poly + (zh_blinding_factor * zh_poly.clone());
        let accumulator_commitment =
            UnivariateKZG::<P>::commitment(&zh_blinding_accumulator_poly, &self.srs);

        self.random_number.beta = beta;
        self.random_number.gamma = gamma;

        self.witness_polys.zh_poly = zh_poly;
        self.witness_polys.zh_accumulator_poly = zh_blinding_accumulator_poly.clone();

        accumulator_commitment
    }

    pub fn third_round(
        &mut self,
        witness: &Witness<F>,
        zh_accumulator_poly: &DenseUnivariatePolynomial<F>,
    ) -> (P::G1, P::G1, P::G1) {
        let group_order = self.preprocessed_input.group_order as usize;
        let root_of_unity: F = root_of_unity(group_order as u64);
        let alpha: F = self.transcript.challenge_round(b"alpha");
        let beta = self.random_number.beta;
        let gamma = self.random_number.gamma;

        let zh_poly = DenseUnivariatePolynomial::new(zh_values(group_order));

        // l1 poly
        let mut l1_values = vec![F::zero(); group_order];
        l1_values[0] = F::one();

        let domain = Domain::new(group_order);
        let l1_poly = UnivariateEval::new(l1_values, domain);

        let w_accumulator_poly =
            apply_w_to_polynomial(&zh_accumulator_poly.clone(), &root_of_unity);

        let t_permutation = (((self.witness_polys.a_s.clone()
            * self.witness_polys.b_s.clone()
            * self.preprocessed_input.q_m.to_coefficient_poly())
            + (self.witness_polys.a_s.clone()
                * self.preprocessed_input.q_l.to_coefficient_poly())
            + (self.witness_polys.b_s.clone()
                * self.preprocessed_input.q_r.to_coefficient_poly())
            + (self.witness_polys.c_s.clone()
                * self.preprocessed_input.q_o.to_coefficient_poly())
            + witness.public_poly.to_coefficient_poly()
            + self.preprocessed_input.q_c.to_coefficient_poly())
            / zh_poly.clone())
            + ((((self.witness_polys.a_s.clone()
                + DenseUnivariatePolynomial::new(vec![F::one(), beta, gamma]))
                * (self.witness_polys.b_s.clone()
                    + DenseUnivariatePolynomial::new(vec![
                        F::one(),
                        beta * F::from(2u32),
                        gamma,
                    ]))
                * (self.witness_polys.c_s.clone()
                    + DenseUnivariatePolynomial::new(vec![
                        F::one(),
                        beta * F::from(3u32),
                        gamma,
                    ]))
                * self.witness_polys.zh_accumulator_poly.clone())
                * alpha)
                / zh_poly.clone())
            - ((((self.witness_polys.a_s.clone()
                + (self.preprocessed_input.sigma_1.to_coefficient_poly() * beta)
                + gamma)
                * (self.witness_polys.b_s.clone()
                    + (self.preprocessed_input.sigma_2.to_coefficient_poly() * beta)
                    + gamma)
                * (self.witness_polys.c_s.clone()
                    + (self.preprocessed_input.sigma_3.to_coefficient_poly() * beta)
                    + gamma)
                * w_accumulator_poly.clone())
                * alpha)
                / zh_poly.clone())
            + ((((self.witness_polys.zh_accumulator_poly.clone() - F::ONE)
                * (l1_poly.to_coefficient_poly()))
                * alpha.pow(&[2 as u64]))
                / zh_poly);

        let (t_low, t_mid, t_high) =
            split_poly_in_3(&t_permutation, self.preprocessed_input.group_order as usize);

        //x^n
        let mut x_n_values = vec![F::zero(); group_order as usize + 1];
        x_n_values[group_order as usize] = F::one();

        //x^2n
        let mut x_2n_values = vec![F::zero(); group_order as usize * 2 + 1];
        x_2n_values[group_order as usize * 2] = F::one();

        let rands: Vec<F> = generate_random_numbers(2);

        let t_low_blinding = DenseUnivariatePolynomial::new(x_2n_values.clone()) * rands[0];
        let t_mid_blinding = DenseUnivariatePolynomial::new(x_n_values) * rands[1] - rands[0];
        let t_high_blinding = t_high.clone() + rands[1].neg();

        let t_low_coeff = t_low + t_low_blinding.clone();
        let t_mid_coeff = t_mid + t_mid_blinding.clone();
        let t_high_coeff = t_high + t_high_blinding.clone();

        let t_low_commitment = UnivariateKZG::<P>::commitment(&t_low_coeff, &self.srs);
        let t_mid_commitment = UnivariateKZG::<P>::commitment(&t_mid_coeff, &self.srs);
        let t_high_commitment = UnivariateKZG::<P>::commitment(&t_high_coeff, &self.srs);

        self.random_number.alpha = alpha;
        self.witness_polys.zh_accumulator_poly = w_accumulator_poly.clone();
        self.witness_polys.t_low_poly = t_low_blinding.clone();
        self.witness_polys.t_mid_poly = t_mid_blinding.clone();
        self.witness_polys.t_low_poly = t_low_blinding.clone();

        (t_low_commitment, t_mid_commitment, t_high_commitment)
    }

    pub fn fourth_round(&mut self) -> (F, F, F, F, F, F) {
        let zeta: F = self.transcript.challenge_round(b"zeta");

        let a_s_poly = self.witness_polys.a_s.clone();
        let b_s_poly = self.witness_polys.b_s.clone();
        let c_s_poly = self.witness_polys.c_s.clone();
        let w_accumulator_poly = self.witness_polys.zh_accumulator_poly.clone();
        let sigma1_poly = self.preprocessed_input.sigma_1.to_coefficient_poly();
        let sigma2_poly = self.preprocessed_input.sigma_2.to_coefficient_poly();

        let a_s_poly_zeta = a_s_poly.evaluate(zeta);
        let b_s_poly_zeta = b_s_poly.evaluate(zeta);
        let c_s_poly_zeta = c_s_poly.evaluate(zeta);
        let w_accumulator_poly_zeta = w_accumulator_poly.evaluate(zeta);
        let sigma1_poly_zeta = sigma1_poly.evaluate(zeta);
        let sigma2_poly_zeta = sigma2_poly.evaluate(zeta);

        self.witness_polys.a_s_poly_zeta = a_s_poly_zeta;
        self.witness_polys.b_s_poly_zeta = b_s_poly_zeta;
        self.witness_polys.c_s_poly_zeta = c_s_poly_zeta;
        self.witness_polys.w_accumulator_poly_zeta = w_accumulator_poly_zeta;
        self.witness_polys.sigma1_poly_zeta = sigma1_poly_zeta;
        self.witness_polys.sigma2_poly_zeta = sigma2_poly_zeta;
        self.random_number.zeta = zeta;

        (
            a_s_poly_zeta,
            b_s_poly_zeta,
            c_s_poly_zeta,
            sigma1_poly_zeta,
            sigma2_poly_zeta,
            w_accumulator_poly_zeta,
        )
    }

    pub fn fifth_round(&mut self, witness: &Witness<F>) -> (P::G1, P::G1) {
        let group_order = self.preprocessed_input.group_order as usize;

        let nu: F = self.transcript.challenge_round(b"nu");
        let alpha = self.random_number.alpha;
        let beta = self.random_number.beta;
        let gamma = self.random_number.gamma;
        let zeta = self.random_number.zeta;

        let a_s_poly = self.witness_polys.a_s.clone();
        let b_s_poly = self.witness_polys.b_s.clone();
        let c_s_poly = self.witness_polys.c_s.clone();
        let w_accumulator_poly = self.witness_polys.zh_accumulator_poly.clone();
        let sigma1_poly = self.preprocessed_input.sigma_1.to_coefficient_poly();
        let sigma2_poly = self.preprocessed_input.sigma_2.to_coefficient_poly();

        let a_s_poly_zeta = self.witness_polys.a_s_poly_zeta;
        let b_s_poly_zeta = self.witness_polys.b_s_poly_zeta;
        let c_s_poly_zeta = self.witness_polys.c_s_poly_zeta;
        let w_accumulator_poly_zeta = self.witness_polys.w_accumulator_poly_zeta;
        let sigma1_poly_zeta = self.witness_polys.sigma1_poly_zeta;
        let sigma2_poly_zeta = self.witness_polys.sigma2_poly_zeta;

        let mut l1_values = vec![F::zero(); group_order];
        l1_values[0] = F::one();
        let l1_poly = DenseUnivariatePolynomial::new(l1_values);

        let zh_poly = DenseUnivariatePolynomial::new(zh_values(group_order));
        let root_of_unity: F = root_of_unity(group_order as u64);

        let r_poly =
            ((self.preprocessed_input.q_m.to_coefficient_poly() * a_s_poly_zeta * b_s_poly_zeta)
                + (self.preprocessed_input.q_l.to_coefficient_poly() * a_s_poly_zeta)
                + (self.preprocessed_input.q_r.to_coefficient_poly() * b_s_poly_zeta)
                + (self.preprocessed_input.q_o.to_coefficient_poly() * c_s_poly_zeta)
                + witness.public_poly.to_coefficient_poly().evaluate(zeta)
                + self.preprocessed_input.q_c.to_coefficient_poly())
                + (((w_accumulator_poly.clone()
                    * (a_s_poly_zeta + (beta * zeta) + gamma)
                    * (b_s_poly_zeta + (beta * F::from(2u8) * zeta) + gamma)
                    * (c_s_poly_zeta + (beta * F::from(3u8) * zeta) + gamma))
                    - (((self.preprocessed_input.sigma_3.to_coefficient_poly() * beta)
                        + c_s_poly_zeta
                        + gamma)
                        * (a_s_poly_zeta + (beta * sigma1_poly_zeta) + gamma)
                        * (b_s_poly_zeta + (beta * sigma2_poly_zeta) + gamma)
                        * w_accumulator_poly_zeta))
                    * alpha)
                + (((w_accumulator_poly.clone() - F::ONE) * (l1_poly.evaluate(zeta)))
                    * alpha.pow(&[2 as u64]))
                - ((self.witness_polys.t_low_poly.clone()
                    + (self.witness_polys.t_mid_poly.clone() * zeta.pow(&[group_order as u64]))
                    + (self.witness_polys.t_high_poly.clone()
                        * zeta.pow(&[2 * self.preprocessed_input.group_order])))
                    * zh_poly.evaluate(zeta));

        let x_minus_zeta_poly = DenseUnivariatePolynomial::new(vec![-zeta, F::one()]);

        let w_zeta_poly = (r_poly
            + (a_s_poly.clone() - a_s_poly_zeta) * nu
            + (b_s_poly.clone() - b_s_poly_zeta) * nu.pow(&[2u64])
            + (c_s_poly.clone() - c_s_poly_zeta) * nu.pow(&[3u64])
            + (sigma1_poly.clone() - sigma1_poly_zeta) * nu.pow(&[4u64])
            + (sigma2_poly.clone() - sigma2_poly_zeta) * nu.pow(&[5u64]))
            / x_minus_zeta_poly;

        let x_minus_zeta_omega_poly =
            DenseUnivariatePolynomial::new(vec![(zeta * root_of_unity).neg(), F::one()]);

        let w_zeta_omega_poly =
            (w_accumulator_poly - w_accumulator_poly_zeta) / x_minus_zeta_omega_poly;

        let w_zeta_commitment = UnivariateKZG::<P>::commitment(&w_zeta_poly, &self.srs);
        let w_zeta_omega_commitment = UnivariateKZG::<P>::commitment(&w_zeta_omega_poly, &self.srs);

        self.random_number.nu = nu;

        self.witness_polys.w_zeta_poly = w_zeta_poly;
        self.witness_polys.w_zeta_omega_poly = w_zeta_omega_poly;

        (w_zeta_commitment, w_zeta_omega_commitment)
    }
}
