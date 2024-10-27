use std::collections::HashMap;

use ark_ff::PrimeField;
use polynomial::DenseUnivariatePolynomial;

use super::{
    assembly::{AssemblyEqn, GateWire},
    domain::Domain,
    utils::{roots_of_unity, Cell, Column},
};

pub struct CommonPreprocessedInput<F: PrimeField> {
    pub group_order: u64,
    pub q_l: UnivariateEval<F>,
    // pub q_l: DenseUnivariatePolynomial<F>,
    pub q_r: UnivariateEval<F>,
    // pub q_r: DenseUnivariatePolynomial<F>,
    pub q_m: UnivariateEval<F>,
    // pub q_m: DenseUnivariatePolynomial<F>,
    pub q_o: UnivariateEval<F>,
    // pub q_o: DenseUnivariatePolynomial<F>,
    pub q_c: UnivariateEval<F>,
    // pub q_c: DenseUnivariatePolynomial<F>,
    pub sigma_1: UnivariateEval<F>,
    // pub sigma_1: DenseUnivariatePolynomial<F>,
    pub sigma_2: UnivariateEval<F>,
    // pub sigma_2: DenseUnivariatePolynomial<F>,
    pub sigma_3: UnivariateEval<F>,
    // pub sigma_3: DenseUnivariatePolynomial<F>,
    pub s1_coeff: Option<DenseUnivariatePolynomial<F>>,
    pub s2_coeff: Option<DenseUnivariatePolynomial<F>>,
}

#[derive(Clone)]
pub struct Program<F: PrimeField> {
    pub constraints: Vec<AssemblyEqn<F>>,
    pub group_order: u64,
}

impl<F: PrimeField> Program<F> {
    pub fn new(constraints: Vec<AssemblyEqn<F>>, group_order: u64) -> Program<F> {
        Program {
            constraints,
            group_order,
        }
    }
    pub fn common_preprocessed_input(&self) -> CommonPreprocessedInput<F> {
        let (q_l, q_r, q_m, q_o, q_c) = self.make_gate_polynomials();
        let (sigma_1, sigma_2, sigma_3) = self.make_s_polynomials();
        CommonPreprocessedInput {
            group_order: self.group_order,
            q_l,
            q_r,
            q_m,
            q_o,
            q_c,
            sigma_1,
            sigma_2,
            sigma_3,
            s1_coeff: None,
            s2_coeff: None,
        }
    }

    pub fn make_gate_polynomials(
        &self,
    ) -> (
        UnivariateEval<F>,
        UnivariateEval<F>,
        UnivariateEval<F>,
        UnivariateEval<F>,
        UnivariateEval<F>,
    ) {
        //make L R M O C gate polynomials
        let mut l = vec![F::zero(); self.group_order as usize];
        let mut r = vec![F::zero(); self.group_order as usize];
        let mut m = vec![F::zero(); self.group_order as usize];
        let mut o = vec![F::zero(); self.group_order as usize];
        let mut c = vec![F::zero(); self.group_order as usize];
        for (i, constraint) in self.constraints.iter().enumerate() {
            let gate = constraint.gate();
            l[i] = gate.l;
            r[i] = gate.r;
            m[i] = gate.m;
            o[i] = gate.o;
            c[i] = gate.c;
        }

        let domain = Domain::new(self.group_order as usize);

        (
            UnivariateEval::new(l, domain.clone()),
            UnivariateEval::new(r, domain.clone()),
            UnivariateEval::new(m, domain.clone()),
            UnivariateEval::new(o, domain.clone()),
            UnivariateEval::new(c, domain.clone()),
        )
    }

    pub fn make_s_polynomials(&self) -> (UnivariateEval<F>, UnivariateEval<F>, UnivariateEval<F>) {
        let mut variable_uses = HashMap::new();
        for (row, constraint) in self.constraints.iter().enumerate() {
            for (column, variable) in constraint.wires.to_vec().into_iter().enumerate() {
                variable_uses.entry(variable).or_insert(vec![]).push(Cell {
                    column: (column + 1).into(),
                    row,
                });
            }
        }

        for row in self.constraints.len()..self.group_order as usize {
            for i in 1..=3 {
                variable_uses.entry(None).or_insert(vec![]).push(Cell {
                    column: i.into(),
                    row,
                })
            }
        }
        let mut s: HashMap<Column, Vec<F>> = HashMap::new();
        s.insert(
            Column::LEFT,
            roots_of_unity(self.group_order)
                .into_iter()
                .map(|element| element)
                .collect(),
        );
        dbg!("HERE");
        s.insert(
            Column::RIGHT,
            roots_of_unity(self.group_order)
                .into_iter()
                .map(|element: F| element * F::from(2u8))
                .collect(),
        );
        s.insert(Column::OUTPUT, vec![F::zero(); self.group_order as usize]);

        // exmaple
        // variable_uses = {"a":[Cell(1,3),Cell(3,4)],"b":[Cell(2,1)]
        for (_, uses) in variable_uses.iter() {
            // _ = "a"
            // uses = [Cell(1,3),Cell(3,4)]
            for (i, cell) in uses.iter().enumerate() {
                let next_i = (i + 1) % uses.len();
                let next_column = uses[next_i].column;
                let next_row = uses[next_i].row;
                if let Some(vec) = s.get_mut(&next_column) {
                    vec[next_row] = cell.label(self.group_order);
                }
            }
        }

        // s1,s2，s3
        let mut s1 = None;
        let mut s2 = None;
        let mut s3 = None;
        for (key, vec) in s.into_iter() {
            let domain = Domain::new(self.group_order as usize);
            match key {
                Column::LEFT => s1 = Some(UnivariateEval::new(vec, domain.clone())),
                Column::RIGHT => s2 = Some(UnivariateEval::new(vec, domain.clone())),
                Column::OUTPUT => s3 = Some(UnivariateEval::new(vec, domain.clone())),
            }
        }
        (s1.unwrap(), s2.unwrap(), s3.unwrap())
    }

    pub fn coeffs(&self) -> Vec<HashMap<Option<String>, F>> {
        let mut coeffs = Vec::new();
        for constraint in self.constraints.iter() {
            coeffs.push(constraint.coeffs.clone());
        }
        coeffs
    }

    pub fn wires(&self) -> Vec<GateWire> {
        let mut wires = Vec::new();
        for constraint in self.constraints.iter() {
            wires.push(constraint.wires.clone());
        }
        return wires;
    }

    pub fn get_public_assignment(&self) -> Vec<Option<String>> {
        let coeffs = self.coeffs();
        let mut out = Vec::new();
        let mut no_more_allowed = false;
        for coeff in coeffs.iter() {
            if coeff.get(&Some("$public".to_string())) != None {
                if no_more_allowed {
                    panic!("Public var declarations must be at the top")
                }
                let mut var_name = Vec::new();
                for (key, _) in coeff.iter() {
                    if key.clone().unwrap().chars().next().unwrap() != '$' {
                        var_name.push(key.clone().unwrap());
                    }
                }

                out.push(Some(var_name.join("")));
            } else {
                no_more_allowed = true;
            }
        }
        out
    }
}

#[derive(Debug, Clone)]
pub struct UnivariateEval<F: PrimeField> {
    /// this is a list of the evaluation of the polynomial
    pub values: Vec<F>,
    /// This is the domian of the polynomal; very important for the FFT and IFFT
    pub domain: Domain<F>,
}

impl<F: PrimeField> UnivariateEval<F> {
    /// This function is used to create a new polynomial from the evaluation form
    pub fn new(values: Vec<F>, domain: Domain<F>) -> Self {
        UnivariateEval { values, domain }
    }
}

#[cfg(test)]
mod test {
    use crate::compiler::{assembly::AssemblyEqn, utils::roots_of_unity};
    use ark_test_curves::bls12_381::Fr;

    use super::Program;

    #[test]
    fn test_make_s_polynomials() {
        //passed
        //L R  O
        //w 2w 3w
        //w^2 2w^2 3w^2

        //a b c
        //a e b
        let original_constriants = ["c <== a * b", "b <== a * e"];
        let mut assembly_eqns = Vec::new();
        for eq in original_constriants.iter() {
            let assembly_eqn: AssemblyEqn<Fr> = AssemblyEqn::eq_to_assembly(eq);
            assembly_eqns.push(assembly_eqn);
        }
        let program = Program::new(assembly_eqns, 8);
        let (s1, s2, s3) = program.make_s_polynomials();

        let unmoved_s1: Vec<_> = roots_of_unity(8);
        let unmoved_s2: Vec<_> = roots_of_unity(8)
            .into_iter()
            .map(|ele: Fr| ele * Fr::from(2))
            .collect();
        let unmoved_s3: Vec<_> = roots_of_unity(8)
            .into_iter()
            .map(|ele: Fr| ele * Fr::from(3))
            .collect();
        assert_eq!(s1.values[0], unmoved_s1[1]);

        assert_eq!(s2.values[0], unmoved_s3[1]);

        // println!("s1:{:?}", s1);
        // println!("s2:{:?}", s2);
        // println!("s3:{:?}", s3);
    }
    #[test]
    fn test_make_gate_polynomials() {
        let original_constriants = ["e public", "c <== a * b", "e <== c * d"];
        let mut assembly_eqns = Vec::new();
        for eq in original_constriants.iter() {
            let assembly_eqn: AssemblyEqn<Fr> = AssemblyEqn::eq_to_assembly(eq);
            assembly_eqns.push(assembly_eqn);
        }
        let program = Program::new(assembly_eqns, 8);
        let (l, r, m, o, c) = program.make_gate_polynomials();
        println!("l:{:?}", l);
        println!("r:{:?}", r);
        println!("m:{:?}", m);
        println!("o:{:?}", o);
        println!("c:{:?}", c);
    }
}
