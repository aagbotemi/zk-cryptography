#[derive(Debug)]
pub enum GateType {
    Add,
    Mul,
}

#[derive(Debug)]
pub struct Gate {
    pub gate_type: GateType,
    pub inputs: [usize; 2],
}

impl Gate {
    pub fn new(gate_type: GateType, inputs: [usize; 2]) -> Self {
        Gate { gate_type, inputs }
    }
}
