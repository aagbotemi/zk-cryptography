use std::collections::HashMap;

use ark_test_curves::bls12_381::{Bls12_381, Fr};
use kzg::{
    interface::UnivariateKZGInterface, trusted_setup::TrustedSetup, univariate_kzg::UnivariateKZG,
};
use plonk::{
    compiler::primitives::{AssemblyEqn, Program},
    protocol::{
        primitives::{PlonkProver, PlonkRoundTranscript},
        verifier::{PlonkVerifier, VerifierPreprocessedInput},
    },
};

fn main() {
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
    // dbg!(&preprocessed_input);
    // dbg!(&srs.powers_of_tau_in_g2.len());
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
