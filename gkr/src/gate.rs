pub enum GateType {
    Add,
    Mul,
}

pub struct Gate {
    pub ttype: GateType,
    pub inputs: [usize; 2],
}

impl Gate {
    pub fn new(ttype: GateType, inputs: [usize; 2]) -> Self {
        Gate { ttype, inputs }
    }
}
