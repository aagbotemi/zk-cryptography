use crate::{
    gate::{Gate, GateType},
    utils::{size_of_mle_n_var_at_each_layer, transform_label_to_binary_and_to_decimal},
};
use ark_ff::PrimeField;
use polynomial::Multilinear;
use std::ops::{Add, Mul};

#[derive(Debug)]
pub struct CircuitLayer {
    pub layer: Vec<Gate>,
}

#[derive(Debug)]
pub struct Circuit {
    pub layers: Vec<CircuitLayer>,
}

impl CircuitLayer {
    pub fn new(layer: Vec<Gate>) -> Self {
        CircuitLayer { layer }
    }
}

impl Circuit {
    pub fn new(layers: Vec<CircuitLayer>) -> Circuit {
        Self { layers }
    }
}

impl Circuit {
    pub fn evaluation<F: PrimeField>(&self, input: &[F]) -> Vec<Vec<F>>
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
        layers
    }

    pub fn add_mult_mle<F: PrimeField>(
        &self,
        layer_index: usize,
    ) -> (Multilinear<F>, Multilinear<F>) {
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

                    add_evaluations[gate_decimal] = F::one()
                }
                GateType::Mul => {
                    let gate_decimal = transform_label_to_binary_and_to_decimal(
                        layer_index,
                        gate_index,
                        gate.inputs[0],
                        gate.inputs[1],
                    );
                    mul_evaluations[gate_decimal] = F::one();
                }
            }
        }

        let add_mle = Multilinear::new(add_evaluations);
        let mul_mle = Multilinear::new(mul_evaluations);

        (add_mle, mul_mle)
    }

    pub fn random(num_of_layers: usize) -> Self {
        let mut layers = Vec::new();

        for layer_index in 0..num_of_layers {
            let mut layer = Vec::new();
            let number_of_gates = 2usize.pow(layer_index as u32);
            let number_of_inputs = 2usize.pow((layer_index + 1) as u32);

            for gate_index in 0..number_of_gates {
                let input_1 = (gate_index * 2) % number_of_inputs;
                let input_2 = (gate_index * 2 + 1) % number_of_inputs;
                let gate_type = if layer_index % 2 == 0 {
                    GateType::Add
                } else {
                    GateType::Mul
                };
                layer.push(Gate::new(gate_type, [input_1, input_2]));
            }

            layers.push(CircuitLayer::new(layer));
        }

        Circuit::new(layers)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gate::Gate;
    use ark_test_curves::bls12_381::Fr;
    use polynomial::interface::MultilinearTrait;

    // sample circuit evaluation
    //      100(*)    - layer 0
    //     /     \
    //   5(+)_0   20(*)_1 - layer 1
    //  1/ \2    3/  \4
    //  2   3    4    5
    //
    #[test]
    fn test_circuit_evaluation_1() {
        let layer_0 = CircuitLayer::new(vec![Gate::new(GateType::Mul, [0, 1])]);
        let layer_1 = CircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
        ]);
        let circuit = Circuit::new(vec![layer_0, layer_1]);
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

        assert_eq!(evaluation, expected_output);
    }

    #[test]
    fn test_circuit_evaluation_2() {
        let layer_0 = CircuitLayer::new(vec![
            Gate::new(GateType::Mul, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
        ]);

        let layer_1 = CircuitLayer::new(vec![
            Gate::new(GateType::Mul, [0, 0]),
            Gate::new(GateType::Mul, [1, 1]),
            Gate::new(GateType::Mul, [1, 2]),
            Gate::new(GateType::Mul, [3, 3]),
        ]);

        let circuit = Circuit::new(vec![layer_0, layer_1]);
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

        assert_eq!(evaluation, expected_output);
    }

    #[test]
    fn test_circuit_evaluation_3() {
        let layer_0 = CircuitLayer::new(vec![Gate::new(GateType::Add, [0, 1])]);

        let layer_1 = CircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
        ]);

        let layer_2 = CircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
            Gate::new(GateType::Mul, [4, 5]),
            Gate::new(GateType::Mul, [6, 7]),
        ]);

        let circuit = Circuit::new(vec![layer_0, layer_1, layer_2]);

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

        assert_eq!(evaluation, expected_output);
    }

    #[test]
    fn test_get_add_n_mul_mle_layer_0() {
        let layer_0 = CircuitLayer::new(vec![Gate::new(GateType::Add, [0, 1])]);

        let layer_1 = CircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
        ]);

        let layer_2 = CircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
            Gate::new(GateType::Mul, [4, 5]),
            Gate::new(GateType::Mul, [6, 7]),
        ]);

        let circuit = Circuit::new(vec![layer_0, layer_1, layer_2]);

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
        let layer_0 = CircuitLayer::new(vec![Gate::new(GateType::Add, [0, 1])]);

        let layer_1 = CircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
        ]);

        let layer_2 = CircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
            Gate::new(GateType::Mul, [4, 5]),
            Gate::new(GateType::Mul, [6, 7]),
        ]);

        let circuit = Circuit::new(vec![layer_0, layer_1, layer_2]);

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
        let layer_0 = CircuitLayer::new(vec![Gate::new(GateType::Add, [0, 1])]);

        let layer_1 = CircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
        ]);

        let layer_2 = CircuitLayer::new(vec![
            Gate::new(GateType::Add, [0, 1]),
            Gate::new(GateType::Mul, [2, 3]),
            Gate::new(GateType::Mul, [4, 5]),
            Gate::new(GateType::Mul, [6, 7]),
        ]);

        let circuit = Circuit::new(vec![layer_0, layer_1, layer_2]);

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
