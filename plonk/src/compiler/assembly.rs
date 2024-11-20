use ark_ff::PrimeField;
use std::collections::HashMap;

use super::{
    primitives::{AssemblyEqn, Gate, GateWire},
    utils::{evaluate, get_product_key, is_valid_variable_name},
};

impl GateWire {
    pub fn to_vec(&self) -> Vec<Option<String>> {
        vec![
            self.left_wire.clone(),
            self.right_wire.clone(),
            self.output_wire.clone(),
        ]
    }
}

impl<F: PrimeField> AssemblyEqn<F> {
    pub fn left(&self) -> F {
        let value = self.coeffs.get(&self.wires.left_wire);
        match value {
            Some(_) => (*value.unwrap()).neg(),
            None => F::zero(),
        }
    }

    pub fn right(&self) -> F {
        if self.wires.right_wire != self.wires.left_wire {
            let value = self.coeffs.get(&self.wires.right_wire);
            return match value {
                Some(_) => (*value.unwrap()).neg(),
                None => F::zero(),
            };
        }
        F::zero()
    }

    pub fn constant(&self) -> F {
        let value = self.coeffs.get(&None);
        match value {
            Some(_) => (*value.unwrap()).neg(),
            None => F::zero(),
        }
    }
    pub fn output(&self) -> F {
        let value = self.coeffs.get(&Some("$output_coeff".to_owned()));
        match value {
            Some(_) => *value.unwrap(),
            None => F::one(),
        }
    }

    pub fn mul(&self) -> F {
        if !self.wires.to_vec().contains(&None) {
            let value = self.coeffs.get(&get_product_key(
                self.wires.left_wire.clone(),
                self.wires.right_wire.clone(),
            ));
            return match value {
                Some(_) => (*value.unwrap()).neg(),
                None => F::zero(),
            };
        }
        F::zero()
    }

    pub fn gate(&self) -> Gate<F> {
        Gate {
            l: self.left(),
            r: self.right(),
            m: self.mul(),
            o: self.output(),
            c: self.constant(),
        }
    }

    //new eq_to_assembly
    pub fn eq_to_assembly(eq: &str) -> AssemblyEqn<F> {
        let tokens: Vec<&str> = eq.trim().split(" ").collect();
        if tokens[1] == "<==" || tokens[1] == "===" {
            // First token is the output variable
            let mut out = tokens[0];
            // Convert the expression to coefficient map form
            let mut coeffs = evaluate(&tokens[2..].to_vec());
            // Handle the "-x === a * b" case
            if out.chars().nth(0).unwrap() == '-' {
                out = &out[1..];
                coeffs.insert(Some("$output_coeff".to_string()), F::one().neg());
            }
            // Check out variable name validity
            if !is_valid_variable_name(out) {
                panic!("Invalid out variable name: {}", out);
            }
            // Gather list of variables used in the expression
            let mut variables: Vec<&str> = Vec::new();
            for &t in tokens[2..].iter() {
                let var = &t.trim_start_matches("-");
                if is_valid_variable_name(var) && !variables.contains(var) {
                    variables.push(var);
                }
            }
            // Construct the list of allowed coefficients
            let mut allowed_coeffs: Vec<String> =
                variables.iter().map(|&s| s.to_string()).collect();
            allowed_coeffs.extend(vec!["".to_string(), "$output_coeff".to_string()]);

            if variables.is_empty() {
                todo!();
            } else if variables.len() == 1 {
                variables.push(variables[0]);
                let product_key =
                    get_product_key(Some(variables[0].to_owned()), Some(variables[1].to_owned()))
                        .unwrap();
                allowed_coeffs.push(product_key);
            } else if variables.len() == 2 {
                let product_key =
                    get_product_key(Some(variables[0].to_owned()), Some(variables[1].to_owned()))
                        .unwrap();
                allowed_coeffs.push(product_key);
            } else {
                panic!("Max 2 variables, found {}", variables.len());
            }

            // Check that only allowed coefficients are in the coefficient map
            for key_option in coeffs.keys() {
                // Use as_ref to convert Option<String> to Option<&String> so that you can safely access the String reference inside it.
                let key_ref = key_option.as_ref().unwrap();

                // Check if allowed_coeffs contains this reference
                if !allowed_coeffs.contains(key_ref) {
                    panic!("Disallowed multiplication");
                }
            }

            // Return output
            let variables_len = variables.len();
            let mut wires: Vec<Option<&str>> = variables
                .into_iter()
                .map(Some)
                .chain(vec![None; 2 - variables_len])
                .collect();
            wires.push(Some(out));

            return AssemblyEqn {
                wires: GateWire {
                    left_wire: Some(wires[0].unwrap().to_string()),
                    right_wire: Some(wires[1].unwrap().to_string()),
                    output_wire: Some(wires[2].unwrap().to_string()),
                },
                coeffs,
            };
        } else if tokens[1] == "public" {
            let mut coeffs = HashMap::new();
            coeffs.insert(Some(tokens[0].to_string()), F::one().neg());
            coeffs.insert(Some("$output_coeff".to_string()), F::zero());
            coeffs.insert(Some("$public".to_string()), F::one());
            return AssemblyEqn {
                wires: GateWire {
                    left_wire: Some(tokens[0].to_string()),
                    right_wire: None,
                    output_wire: None,
                },
                coeffs,
            };
        } else {
            panic!("Unsupported op: {}", tokens[1]);
        }
    }
}
