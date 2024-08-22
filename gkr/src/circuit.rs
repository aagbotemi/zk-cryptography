use ark_ff::PrimeField;
use polynomial::{interface::MLETrait, MLE};
use std::ops::{Add, Mul};

use crate::{
    gate::{Gate, GateType},
    utils::{
        bit_count_for_n_elem, size_of_mle_n_var_at_each_layer,
        transform_label_to_binary_and_to_decimal,
    },
};

#[derive(Debug)]
pub struct GKRCircuitLayer {
    pub layer: Vec<Gate>,
}

#[derive(Debug)]
pub struct GKRCircuit {
    pub layers: Vec<GKRCircuitLayer>,
}

#[derive(Debug)]
pub struct GKRCircuitEvaluation<F> {
    pub layers: Vec<Vec<F>>,
}

impl GKRCircuitLayer {
    pub fn new(layer: Vec<Gate>) -> Self {
        GKRCircuitLayer { layer }
    }
}

impl<F> GKRCircuitEvaluation<F> {
    pub fn new(layers: Vec<Vec<F>>) -> Self {
        GKRCircuitEvaluation { layers }
    }
}

impl GKRCircuit {
    pub fn new(layers: Vec<GKRCircuitLayer>) -> GKRCircuit {
        Self { layers }
    }
}

impl GKRCircuit {
    pub fn evaluation<F: PrimeField>(&self, input: &[F]) -> GKRCircuitEvaluation<F>
    where
        F: Add<Output = F> + Mul<Output = F> + Copy,
    {
        let mut layers = vec![];
        let mut current_input = input;

        layers.push(input.to_vec());

        for layer in self.layers.iter().rev() {
            let temp_layer: Vec<F> = layer
                .layer
                .iter()
                .map(|e| match e.gate_type {
                    GateType::Add => current_input[e.inputs[0]] + current_input[e.inputs[1]],
                    GateType::Mul => current_input[e.inputs[0]] * current_input[e.inputs[1]],
                })
                .collect();

            layers.push(temp_layer);
            current_input = &layers[layers.len() - 1];
        }

        layers.reverse();
        GKRCircuitEvaluation { layers }
    }

    pub fn add_mult_mle<F: PrimeField>(&self, layer_index: usize) -> (MLE<F>, MLE<F>) {
        // dbg!("constructing for layer = {}", layer_index);
        let layer = &self.layers[layer_index];
        let n_vars = size_of_mle_n_var_at_each_layer(layer_index);

        let mut add_evaluations = vec![F::zero(); n_vars];
        let mut mul_evaluations = vec![F::zero(); n_vars];

        for (gate_index, gate) in layer.layer.iter().enumerate() {
            match gate.gate_type {
                GateType::Add => {
                    let gate_decimal = transform_label_to_binary_and_to_decimal(
                        layer_index,
                        gate_index,
                        gate.inputs[0],
                        gate.inputs[1],
                    );

                    dbg!("add gate_decimal = {}", gate_decimal);
                    add_evaluations[gate_decimal] = F::one()
                }
                GateType::Mul => {
                    let gate_decimal = transform_label_to_binary_and_to_decimal(
                        layer_index,
                        gate_index,
                        gate.inputs[0],
                        gate.inputs[1],
                    );
                    dbg!("mul gate_decimal = {}", gate_decimal);
                    mul_evaluations[gate_decimal] = F::one();
                }
            }
        }

        let add_mle = MLE::new(add_evaluations);
        let mul_mle = MLE::new(mul_evaluations);

        (add_mle, mul_mle)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gate::Gate;
    use ark_test_curves::bls12_381::Fr;
    use polynomial::interface::MLETrait;

    // sample circuit evaluation
    //      100(*)    - layer 0
    //     /     \
    //   5(+)_0   20(*)_1 - layer 1
    //  1/ \2    3/  \4
    //  2   3    4    5
    //
    #[test]
    fn test_circuit_evaluation_1() {
        let layer_0 = GKRCircuitLayer::new(vec![Gate::new(GateType::Mul, [0, 1])]);
        let layer_1 = GKRCircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
        ]);
        let circuit = GKRCircuit::new(vec![layer_0, layer_1]);
        let input = [
            Fr::from(2u32),
            Fr::from(3u32),
            Fr::from(4u32),
            Fr::from(5u32),
        ];
        let evaluation = circuit.evaluation(&input);
        let expected_output = vec![
            vec![Fr::from(100u32)],
            vec![Fr::from(5u32), Fr::from(20u32)],
            vec![
                Fr::from(2u32),
                Fr::from(3u32),
                Fr::from(4u32),
                Fr::from(5u32),
            ],
        ];

        assert_eq!(evaluation.layers, expected_output);
    }

    #[test]
    fn test_circuit_evaluation_2() {
        let layer_0 = GKRCircuitLayer::new(vec![
            Gate::new(GateType::Mul, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
        ]);

        let layer_1 = GKRCircuitLayer::new(vec![
            Gate::new(GateType::Mul, [0, 0]),
            Gate::new(GateType::Mul, [1, 1]),
            Gate::new(GateType::Mul, [1, 2]),
            Gate::new(GateType::Mul, [3, 3]),
        ]);

        let circuit = GKRCircuit::new(vec![layer_0, layer_1]);
        let evaluation = circuit.evaluation(&[
            Fr::from(3u32),
            Fr::from(2u32),
            Fr::from(3u32),
            Fr::from(1u32),
        ]);

        let expected_output = vec![
            vec![Fr::from(36u32), Fr::from(6u32)],
            vec![
                Fr::from(9u32),
                Fr::from(4u32),
                Fr::from(6u32),
                Fr::from(1u32),
            ],
            vec![
                Fr::from(3u32),
                Fr::from(2u32),
                Fr::from(3u32),
                Fr::from(1u32),
            ],
        ];

        assert_eq!(evaluation.layers, expected_output);
    }

    #[test]
    fn test_circuit_evaluation_3() {
        let layer_0 = GKRCircuitLayer::new(vec![Gate::new(GateType::Add, [0, 1])]);

        let layer_1 = GKRCircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
        ]);

        let layer_2 = GKRCircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
            Gate::new(GateType::Mul, [4, 5]),
            Gate::new(GateType::Mul, [6, 7]),
        ]);

        let circuit = GKRCircuit::new(vec![layer_0, layer_1, layer_2]);

        let evaluation = circuit.evaluation(&[
            Fr::from(2u32),
            Fr::from(3u32),
            Fr::from(1u32),
            Fr::from(4u32),
            Fr::from(1u32),
            Fr::from(2u32),
            Fr::from(3u32),
            Fr::from(4u32),
        ]);

        let expected_output = vec![
            vec![Fr::from(33u32)],
            vec![Fr::from(9u32), Fr::from(24u32)],
            vec![
                Fr::from(5u32),
                Fr::from(4u32),
                Fr::from(2u32),
                Fr::from(12u32),
            ],
            vec![
                Fr::from(2u32),
                Fr::from(3u32),
                Fr::from(1u32),
                Fr::from(4u32),
                Fr::from(1u32),
                Fr::from(2u32),
                Fr::from(3u32),
                Fr::from(4u32),
            ],
        ];

        assert_eq!(evaluation.layers, expected_output);
    }

    #[test]
    fn test_get_add_n_mul_mle_layer_0() {
        let layer_0 = GKRCircuitLayer::new(vec![Gate::new(GateType::Add, [0, 1])]);

        let layer_1 = GKRCircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
        ]);

        let layer_2 = GKRCircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
            Gate::new(GateType::Mul, [4, 5]),
            Gate::new(GateType::Mul, [6, 7]),
        ]);

        let circuit = GKRCircuit::new(vec![layer_0, layer_1, layer_2]);

        let (add_mle, mul_mle) = circuit.add_mult_mle::<Fr>(0);

        // there is no mul gate in layer 0, the mul mle should be zero
        assert_eq!(mul_mle.is_zero(), true);
        // there is only one add gate in layer 0, the add mle should be a non-zero value
        assert_eq!(add_mle.is_zero(), false);
        // evaulating the add mle at the correct binary combination should give a one
        assert_eq!(
            add_mle.evaluation(&vec![Fr::from(0u32), Fr::from(0u32), Fr::from(1u32)]),
            Fr::from(1u32)
        );

        // evaulating the add mle at the correct binary combination should give a zero
        assert_eq!(
            add_mle.evaluation(&vec![Fr::from(0u32), Fr::from(0u32), Fr::from(0u32)]),
            Fr::from(0u32)
        );
        assert_eq!(
            add_mle.evaluation(&vec![Fr::from(1u32), Fr::from(0u32), Fr::from(0u32)]),
            Fr::from(0u32)
        );
        assert_eq!(
            add_mle.evaluation(&vec![Fr::from(1u32), Fr::from(0u32), Fr::from(1u32)]),
            Fr::from(0u32)
        );
        assert_eq!(
            add_mle.evaluation(&vec![Fr::from(1u32), Fr::from(1u32), Fr::from(1u32)]),
            Fr::from(0u32)
        );
    }

    #[test]
    fn test_get_add_n_mul_mle_layer_1() {
        let layer_0 = GKRCircuitLayer::new(vec![Gate::new(GateType::Add, [0, 1])]);

        let layer_1 = GKRCircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
        ]);

        let layer_2 = GKRCircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
            Gate::new(GateType::Mul, [4, 5]),
            Gate::new(GateType::Mul, [6, 7]),
        ]);

        let circuit = GKRCircuit::new(vec![layer_0, layer_1, layer_2]);

        let (add_mle, mul_mle) = circuit.add_mult_mle::<Fr>(1);

        // there is one mul gate in layer 0, the mul mle should be non-zero
        assert_eq!(mul_mle.is_zero(), false);
        // there is only one add gate in layer 0, the add mle should be a non-zero value
        assert_eq!(add_mle.is_zero(), false);
        // this num of var for the mle should be 5
        assert_eq!(add_mle.n_vars, 5);
        // this num of var for the mle should be 5
        assert_eq!(mul_mle.n_vars, 5);
        // evaulating the add mle at the correct binary combination should give a one
        assert_eq!(
            add_mle.evaluation(&vec![
                Fr::from(0),
                Fr::from(0),
                Fr::from(0),
                Fr::from(0),
                Fr::from(1)
            ]),
            Fr::from(1u32)
        );
        // evaulating the add mle at the correct binary combination should give a one
        assert_eq!(
            mul_mle.evaluation(&vec![
                Fr::from(1),
                Fr::from(1),
                Fr::from(0),
                Fr::from(1),
                Fr::from(1)
            ]),
            Fr::from(1u32)
        );

        // evaulating the mul mle at the correct binary combination should give a zero
        assert_eq!(
            mul_mle.evaluation(&vec![
                Fr::from(1),
                Fr::from(0),
                Fr::from(0),
                Fr::from(1),
                Fr::from(1)
            ]),
            Fr::from(0u32)
        );
        assert_eq!(
            mul_mle.evaluation(&vec![
                Fr::from(0),
                Fr::from(1),
                Fr::from(0),
                Fr::from(1),
                Fr::from(0)
            ]),
            Fr::from(0u32)
        );
        assert_eq!(
            mul_mle.evaluation(&vec![
                Fr::from(0),
                Fr::from(0),
                Fr::from(0),
                Fr::from(1),
                Fr::from(1)
            ]),
            Fr::from(0u32)
        );
        assert_eq!(
            mul_mle.evaluation(&vec![
                Fr::from(1),
                Fr::from(1),
                Fr::from(1),
                Fr::from(1),
                Fr::from(1)
            ]),
            Fr::from(0u32)
        );

        // evaulating the add mle at the correct binary combination should give a zero
        assert_eq!(
            add_mle.evaluation(&vec![
                Fr::from(1),
                Fr::from(0),
                Fr::from(0),
                Fr::from(1),
                Fr::from(1)
            ]),
            Fr::from(0u32)
        );
        assert_eq!(
            add_mle.evaluation(&vec![
                Fr::from(0),
                Fr::from(1),
                Fr::from(0),
                Fr::from(1),
                Fr::from(0)
            ]),
            Fr::from(0u32)
        );
        assert_eq!(
            add_mle.evaluation(&vec![
                Fr::from(0),
                Fr::from(0),
                Fr::from(0),
                Fr::from(1),
                Fr::from(1)
            ]),
            Fr::from(0u32)
        );
        assert_eq!(
            add_mle.evaluation(&vec![
                Fr::from(1),
                Fr::from(1),
                Fr::from(1),
                Fr::from(1),
                Fr::from(1)
            ]),
            Fr::from(0u32)
        );
    }

    #[test]
    fn test_get_add_n_mul_mle_layer_2() {
        let layer_0 = GKRCircuitLayer::new(vec![Gate::new(GateType::Add, [0, 1])]);

        let layer_1 = GKRCircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
        ]);

        let layer_2 = GKRCircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
            Gate::new(GateType::Mul, [4, 5]),
            Gate::new(GateType::Mul, [6, 7]),
        ]);

        let circuit = GKRCircuit::new(vec![layer_0, layer_1, layer_2]);

        let (add_mle, mul_mle) = circuit.add_mult_mle::<Fr>(2);

        // there is one mul gate in layer 0, the mul mle should be non-zero
        assert_eq!(mul_mle.is_zero(), false);
        // there is only one add gate in layer 0, the add mle should be a non-zero value
        assert_eq!(add_mle.is_zero(), false);
        // this num of var for the mle should be 5
        assert_eq!(add_mle.n_vars, 8);
        // this num of var for the mle should be 5
        assert_eq!(mul_mle.n_vars, 8);

        // evaulating the add mle at the correct binary combination should give a one
        assert_eq!(
            add_mle.evaluation(&vec![
                Fr::from(0),
                Fr::from(0),
                Fr::from(0),
                Fr::from(0),
                Fr::from(0),
                Fr::from(0),
                Fr::from(0),
                Fr::from(1)
            ]),
            Fr::from(1u32)
        );

        // evaulating the mul mle at the correct binary combination should give a one
        assert_eq!(
            mul_mle.evaluation(&vec![
                Fr::from(1),
                Fr::from(0),
                Fr::from(1),
                Fr::from(0),
                Fr::from(0),
                Fr::from(1),
                Fr::from(0),
                Fr::from(1)
            ]),
            Fr::from(1u32)
        );
        assert_eq!(
            mul_mle.evaluation(&vec![
                Fr::from(1),
                Fr::from(1),
                Fr::from(1),
                Fr::from(1),
                Fr::from(0),
                Fr::from(1),
                Fr::from(1),
                Fr::from(1)
            ]),
            Fr::from(1u32)
        );
    }
}
