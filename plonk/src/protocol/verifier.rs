use std::marker::PhantomData;

use ark_ec::{pairing::Pairing, AffineRepr, Group};
use ark_ff::PrimeField;
use kzg::{
    interface::UnivariateKZGInterface, trusted_setup::TrustedSetup, univariate_kzg::UnivariateKZG,
};
use polynomial::{
    univariate::evaluation::UnivariateEval, DenseUnivariatePolynomial, UnivariatePolynomialTrait,
};

use crate::{
    compiler::{primitives::CommonPreprocessedInput, utils::root_of_unity},
    protocol::utils::compute_verifier_challenges,
};

use super::{
    primitives::{PlonkProof, PlonkRoundTranscript},
    utils::l1_values,
};

pub struct Verifier {}

/// This is the verier preprocessed input
pub struct VerifierPreprocessedInput<P: Pairing> {
    pub qm_commitment: P::G1,
    pub ql_commitment: P::G1,
    pub qr_commitment: P::G1,
    pub qo_commitment: P::G1,
    pub qc_commitment: P::G1,
    pub sigma1_commitment: P::G1,
    pub sigma2_commitment: P::G1,
    pub sigma3_commitment: P::G1,
    pub x_2: P::G2,
}

impl<P: Pairing> VerifierPreprocessedInput<P> {
    fn vpi<F: PrimeField>(srs: &TrustedSetup<P>, cpi: &CommonPreprocessedInput<F>) -> Self {
        Self {
            qm_commitment: UnivariateKZG::commitment(&cpi.q_m.to_coefficient_poly(), srs),
            ql_commitment: UnivariateKZG::commitment(&cpi.q_l.to_coefficient_poly(), srs),
            qr_commitment: UnivariateKZG::commitment(&cpi.q_r.to_coefficient_poly(), srs),
            qo_commitment: UnivariateKZG::commitment(&cpi.q_o.to_coefficient_poly(), srs),
            qc_commitment: UnivariateKZG::commitment(&cpi.q_c.to_coefficient_poly(), srs),
            sigma1_commitment: UnivariateKZG::commitment(&cpi.sigma_1.to_coefficient_poly(), srs),
            sigma2_commitment: UnivariateKZG::commitment(&cpi.sigma_2.to_coefficient_poly(), srs),
            sigma3_commitment: UnivariateKZG::commitment(&cpi.sigma_3.to_coefficient_poly(), srs),
            x_2: srs.powers_of_tau_in_g2[1],
        }
    }
}

pub struct PlonkVerifier<P: Pairing, F: PrimeField> {
    pub group_order: u64,
    pub proof: PlonkProof<P, F>,
    pub verifier_preprocessed_input: VerifierPreprocessedInput<P>,
    pub srs: TrustedSetup<P>,
    _marker: PhantomData<F>,
}
impl<P: Pairing, F: PrimeField> PlonkVerifier<P, F> {
    pub fn new(
        group_order: u64,
        proof: PlonkProof<P, F>,
        srs: TrustedSetup<P>,
        verifier_preprocessed_input: VerifierPreprocessedInput<P>,
    ) -> Self {
        Self {
            group_order,
            proof,
            srs,
            verifier_preprocessed_input,
            _marker: PhantomData,
        }
    }

    pub fn verify(&self, public_input_poly: UnivariateEval<F>) -> bool {
        let (beta, gamma, alpha, zeta, nu, mu) = compute_verifier_challenges(&self.proof);

        let group_order = self.group_order;
        let z_h_zeta = zeta.pow(&[group_order]) - F::one();
        let root_of_unity: F = root_of_unity(group_order);

        // let mut l1_values = vec![F::one()];
        // for _ in 0..group_order - 1 {
        //     l1_values.push(F::zero());
        // }
        let l1_poly = DenseUnivariatePolynomial::new(l1_values(group_order));
        let l1_zeta = l1_poly.evaluate(zeta);

        let public_input_poly_at_zeta = public_input_poly.to_coefficient_poly().evaluate(zeta);

        let a_s_zeta = self.proof.a_s_poly_zeta;
        let b_s_zeta = self.proof.b_s_poly_zeta;
        let c_s_zeta = self.proof.c_s_poly_zeta;
        let w_accumulator_poly_zeta = self.proof.w_accumulator_poly_zeta;
        let sigma1_poly_zeta = self.proof.sigma1_poly_zeta;
        let sigma2_poly_zeta = self.proof.sigma2_poly_zeta;

        let r_0 = public_input_poly_at_zeta
            - l1_zeta * alpha.pow(&[2u64])
            - alpha
                * ((a_s_zeta + sigma1_poly_zeta * beta + gamma)
                    * (b_s_zeta + sigma2_poly_zeta * beta + gamma)
                    * (c_s_zeta + gamma)
                    * w_accumulator_poly_zeta);

        let qm = self.verifier_preprocessed_input.qm_commitment;
        let ql = self.verifier_preprocessed_input.ql_commitment;
        let qr = self.verifier_preprocessed_input.qr_commitment;
        let qo = self.verifier_preprocessed_input.qo_commitment;
        let qc = self.verifier_preprocessed_input.qc_commitment;
        let acc = self.proof.accumulator_commitment;
        let sigma3 = self.verifier_preprocessed_input.sigma3_commitment;
        let t_low = self.proof.t_low;
        let t_mid = self.proof.t_mid;
        let t_high = self.proof.t_high;

        let d_1 = (qm.mul_bigint(&(a_s_zeta * b_s_zeta).into_bigint())
            + ql.mul_bigint(&a_s_zeta.into_bigint())
            + qr.mul_bigint(&b_s_zeta.into_bigint())
            + qo.mul_bigint(&c_s_zeta.into_bigint())
            + qc)
            + (acc.mul_bigint(
                ((a_s_zeta + zeta * beta + gamma)
                    * (b_s_zeta + (F::from(2u8) * zeta * beta) + gamma)
                    * (c_s_zeta + (F::from(3u8) * zeta * beta) + gamma)
                    * alpha
                    + l1_zeta * alpha.pow(&[2u64])
                    + mu)
                    .into_bigint(),
            ))
            - (sigma3.mul_bigint(
                ((a_s_zeta + sigma1_poly_zeta * beta + gamma)
                    * (b_s_zeta + sigma2_poly_zeta * beta + gamma)
                    * alpha
                    * beta
                    * w_accumulator_poly_zeta)
                    .into_bigint(),
            ))
            - ((t_low
                + (t_mid.mul_bigint(zeta.pow(&[self.group_order]).into_bigint()))
                + (t_high.mul_bigint(zeta.pow(&[2 * self.group_order]).into_bigint())))
            .mul_bigint(z_h_zeta.into_bigint()));

        let a_s = self.proof.a_s;
        let b_s = self.proof.b_s;
        let c_s = self.proof.c_s;
        let sigma1 = self.verifier_preprocessed_input.sigma1_commitment;
        let sigma2 = self.verifier_preprocessed_input.sigma2_commitment;

        let f_1 = d_1
            + a_s.mul_bigint(nu.into_bigint())
            + b_s.mul_bigint(nu.pow(&[2u64]).into_bigint())
            + c_s.mul_bigint(nu.pow(&[3u64]).into_bigint())
            + sigma1.mul_bigint(nu.pow(&[4u64]).into_bigint())
            + sigma2.mul_bigint(nu.pow(&[5u64]).into_bigint());

        let e_1 = P::G1::generator().mul_bigint(
            (nu * a_s_zeta
                + nu.pow(&[2, 0, 0, 0]) * b_s_zeta
                + nu.pow(&[3, 0, 0, 0]) * c_s_zeta
                + nu.pow(&[4, 0, 0, 0]) * sigma1_poly_zeta
                + nu.pow(&[5, 0, 0, 0]) * sigma2_poly_zeta
                + mu * w_accumulator_poly_zeta
                - r_0)
                .into_bigint(),
        );

        let w_zeta_1 = self.proof.w_zeta_commitment;
        let w_zeta_omega_1 = self.proof.w_zeta_omega_commitment;

        let x_2 = self.verifier_preprocessed_input.x_2;

        let left = P::pairing(
            &(w_zeta_1 + w_zeta_omega_1.mul_bigint(mu.into_bigint())).into(),
            &x_2,
        );

        let right = P::pairing(
            &(w_zeta_1.mul_bigint(zeta.into_bigint())
                + w_zeta_omega_1.mul_bigint((root_of_unity * mu * zeta).into_bigint())
                + f_1
                - e_1)
                .into(),
            P::G2::generator(),
        );

        left == right
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use crate::{
        compiler::primitives::{AssemblyEqn, Program},
        protocol::primitives::PlonkProver,
    };

    use super::*;
    // use crate::{interface::PlonkProverInterface, prover::PlonkProver};
    use ark_test_curves::bls12_381::{Bls12_381, Fr};
    use kzg::{interface::UnivariateKZGInterface, univariate_kzg::UnivariateKZG};
    // use fiat_shamir::FiatShamirTranscript;
    // use kzg_rust::{interface::KZGUnivariateInterface, univariate::UnivariateKZG};
    // use plonk_compiler::{assembly::eq_to_assembly, program::Program};
    // use std::collections::HashMap;

    #[test]
    fn test_plonk_complete_prove_n_verify() {
        let original_constriants = ["e public"];
        let mut assembly_eqns = Vec::new();
        for eq in original_constriants.iter() {
            let assembly_eqn = AssemblyEqn::eq_to_assembly(eq);
            assembly_eqns.push(assembly_eqn);
        }
        let program = Program::new(assembly_eqns, 8);

        let mut variable_assignment = HashMap::new();
        variable_assignment.insert(Some("e".to_string()), Fr::from(3));

        let witness = program.compute_witness(variable_assignment);
        let preprocessed_input = program.common_preprocessed_input();

        let transcript = PlonkRoundTranscript::new();
        let srs: TrustedSetup<Bls12_381> =
            UnivariateKZG::generate_srs(&Fr::from(6), &(program.group_order as usize * 4));
        let verifier_preprocessed_input = VerifierPreprocessedInput::vpi(&srs, &preprocessed_input);
        let mut prover = PlonkProver::new(preprocessed_input, srs.clone(), transcript);
        let proof = prover.prove(&witness);
        let verifer = PlonkVerifier::new(
            program.group_order,
            proof,
            srs.clone(),
            verifier_preprocessed_input,
        );
        let is_valid = verifer.verify(witness.public_poly);
        assert_eq!(is_valid, true);
    }

    // #[test]
    // fn test_plonk_complete_prove_n_verify_1() {
    //     let original_constriants = [
    //         "x public",
    //         "c <== a * b",
    //         "f <== d * e",
    //         "g <== c + f",
    //         "x <== g * y",
    //     ];
    //     let mut assembly_eqns = Vec::new();
    //     for eq in original_constriants.iter() {
    //         let assembly_eqn = eq_to_assembly::<Fr>(eq.to_string());
    //         assembly_eqns.push(assembly_eqn);
    //     }
    //     let program = Program::new(assembly_eqns, 8);

    //     let mut variable_assignment = HashMap::new();
    //     variable_assignment.insert(Some("x".to_string()), Fr::from(258));
    //     variable_assignment.insert(Some("a".to_string()), Fr::from(2));
    //     variable_assignment.insert(Some("b".to_string()), Fr::from(4));
    //     variable_assignment.insert(Some("d".to_string()), Fr::from(5));
    //     variable_assignment.insert(Some("e".to_string()), Fr::from(7));
    //     variable_assignment.insert(Some("y".to_string()), Fr::from(6));

    //     let witness = program.compute_witness_and_public_parameter(variable_assignment);
    //     let circuit_ir = program.common_preproccessed_input();

    //     let transcript = FiatShamirTranscript::new("plonk-protocol".as_bytes().to_vec());
    //     let srs: SRS<Bls12_381> =
    //         UnivariateKZG::generate_srs(&Fr::from(6), &program.group_order as usize * 4);
    //     let mut prover = PlonkProver::new(transcript, circuit_ir.clone(), srs.clone());
    //     let proof = prover.prove(&witness);
    //     let verifer = PlonkVerifier::new(program.group_order, circuit_ir.clone(), srs);
    //     let is_valid = verifer.verify(&proof, witness.pi);
    //     assert_eq!(is_valid, true);
    // }
}
