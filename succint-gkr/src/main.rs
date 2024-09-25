use ark_test_curves::bls12_381::{Bls12_381, Fr};
use circuit::{
    circuit::{Circuit, CircuitLayer},
    gate::{Gate, GateType},
};
use multilinear_kzg::{interface::TrustedSetupInterface, trusted_setup::TrustedSetup};
use succint_gkr::succint_gkr::SuccintGKRProtocol;

fn main() {
    let layer_0 = CircuitLayer::new(vec![Gate::new(GateType::Add, [0, 1])]);
    let layer_1 = CircuitLayer::new(vec![
        Gate::new(GateType::Mul, [0, 1]),
        Gate::new(GateType::Add, [2, 3]),
    ]);
    let layer_3 = CircuitLayer::new(vec![
        Gate::new(GateType::Add, [0, 1]),
        Gate::new(GateType::Mul, [2, 3]),
        Gate::new(GateType::Mul, [4, 5]),
        Gate::new(GateType::Mul, [6, 7]),
    ]);
    let layer_4 = CircuitLayer::new(vec![
        Gate::new(GateType::Mul, [0, 1]),
        Gate::new(GateType::Mul, [2, 3]),
        Gate::new(GateType::Mul, [4, 5]),
        Gate::new(GateType::Add, [6, 7]),
        Gate::new(GateType::Mul, [8, 9]),
        Gate::new(GateType::Add, [10, 11]),
        Gate::new(GateType::Mul, [12, 13]),
        Gate::new(GateType::Mul, [14, 15]),
    ]);

    let circuit = Circuit::new(vec![layer_0, layer_1, layer_3, layer_4]);
    let input = vec![
        Fr::from(2u32),
        Fr::from(1u32),
        Fr::from(3u32),
        Fr::from(1u32),
        Fr::from(4u32),
        Fr::from(1u32),
        Fr::from(2u32),
        Fr::from(2u32),
        Fr::from(3u32),
        Fr::from(3u32),
        Fr::from(4u32),
        Fr::from(4u32),
        Fr::from(2u32),
        Fr::from(3u32),
        Fr::from(3u32),
        Fr::from(4u32),
    ];

    let evaluation = circuit.evaluation(&input);

    assert_eq!(evaluation[0][0], Fr::from(224u32));

    dbg!("BEGINNING OF TAU IS HERE");
    let tau = TrustedSetup::<Bls12_381>::setup(&input);
    dbg!("TAU IS HERE");

    let (commitment, proof) = SuccintGKRProtocol::prove(&circuit, &input, &tau);
    dbg!("COMMITMENT AND PROOF IS HERE");

    let verify = SuccintGKRProtocol::verify(&circuit, &commitment, &proof, &tau);
    dbg!("VERIFICATION STATUS IS HERE");
    assert!(verify);
}
