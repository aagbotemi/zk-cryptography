use std::collections::HashMap;

use ark_test_curves::bls12_381::{Bls12_381, Fr};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use kzg::{interface::UnivariateKZGInterface, univariate_kzg::UnivariateKZG};
use plonk::{
    compiler::primitives::{AssemblyEqn, Program},
    protocol::{
        primitives::{PlonkProver, PlonkRoundTranscript, VerifierPreprocessedInput},
        verifier::PlonkVerifier,
    },
};

fn plonk_benchmark(c: &mut Criterion) {
    let original_constriants = black_box([
        "x public",
        "c <== a * b",
        "f <== d * e",
        "g <== c + f",
        "x <== g * y",
    ]);

    let assembly_eqns = {
        let mut assembly_eqns = Vec::new();
        for eq in original_constriants.iter() {
            let assembly_eqn = AssemblyEqn::eq_to_assembly(eq);
            assembly_eqns.push(assembly_eqn);
        }
        assembly_eqns
    };

    let program = Program::new(assembly_eqns, 8);

    let variable_assignment = {
        let mut variable_assignment = HashMap::new();
        variable_assignment.insert(Some("x".to_string()), Fr::from(258));
        variable_assignment.insert(Some("a".to_string()), Fr::from(2));
        variable_assignment.insert(Some("b".to_string()), Fr::from(4));
        variable_assignment.insert(Some("d".to_string()), Fr::from(5));
        variable_assignment.insert(Some("e".to_string()), Fr::from(7));
        variable_assignment.insert(Some("y".to_string()), Fr::from(6));
        variable_assignment
    };

    let witness = program.compute_witness_and_public_poly(black_box(variable_assignment));
    let preprocessed_input = program.common_preprocessed_input();
    let srs = UnivariateKZG::generate_srs(&Fr::from(6), &(program.group_order as usize * 4));
    let verifier_preprocessed_input = VerifierPreprocessedInput::vpi(&srs, &preprocessed_input);

    c.bench_function("plonk_benchmark", |b| {
        b.iter(|| {
            let transcript: PlonkRoundTranscript<Bls12_381> = PlonkRoundTranscript::new();
            let mut prover = PlonkProver::new(preprocessed_input.clone(), srs.clone(), transcript);
            let proof = prover.prove(&witness);

            let verifier = PlonkVerifier::new(
                program.group_order,
                proof,
                srs.clone(),
                verifier_preprocessed_input.clone(),
            );
            let is_valid = verifier.verify(witness.public_poly.clone());

            assert!(is_valid);
        });
    });
}

criterion_group!(benches, plonk_benchmark);
criterion_main!(benches);
