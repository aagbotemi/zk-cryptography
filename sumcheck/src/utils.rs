use ark_ff::{BigInteger, PrimeField};

pub fn boolean_hypercube<F: PrimeField>(n: usize) -> Vec<Vec<F>> {
    let mut hypercube = Vec::new();

    for i in 0..(1 << n) {
        let mut vertex = Vec::new();
        for j in 0..n {
            if (i & (1 << j)) != 0 {
                vertex.push(F::one());
            } else {
                vertex.push(F::zero());
            }
        }
        hypercube.push(vertex);
    }

    hypercube
}

pub fn convert_field_to_byte<F: PrimeField>(element: &F) -> Vec<u8> {
    element.into_bigint().to_bytes_be()
}
