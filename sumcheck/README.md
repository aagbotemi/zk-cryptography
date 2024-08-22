# Sumcheck Protocol Implementation
This repository contains an implementation of the Sumcheck Protocol in Rust. The Sumcheck Protocol is a fundamental building block in many zero-knowledge proof systems and interactive proofs.

## Overview
The implementation includes three main components:
1. **Basic Sumcheck:** is used to prove the sum of a multilinear polynomial over a boolean hypercube. It allows a prover to convince a verifier that the sum of the polynomial evaluations over all boolean inputs is equal to a claimed value, without the verifier having to compute all these evaluations.
2. **Composed Sumcheck:** extends the basic protocol to work with more complex polynomial structures. It deals with polynomials that are compositions of simpler polynomials, allowing for more efficient proofs of more intricate relationships.
3. **Multi-Composed Sumcheck:** further extends the protocol to handle multiple composed polynomials simultaneously. This variant is especially useful in complex zero-knowledge proof systems where you need to prove statements about multiple interrelated polynomials. It allows for more efficient proofs by batching the sumcheck process for multiple polynomials.

Each of these variants builds upon the previous one, offering more flexibility and efficiency for different types of polynomial sum proofs. They are crucial in various cryptographic protocols, particularly in the design of efficient zero-knowledge proof systems.

## Features

- Implementation of basic Sumcheck protocol
- Composed Sumcheck for more complex polynomial structures
- Multi-Composed Sumcheck for handling multiple polynomials
- Fiat-Shamir heuristic for non-interactive proofs
- Utility functions for boolean hypercube generation and field conversions

## Main Structures

- `Sumcheck<F: PrimeField>`
- `SumcheckProof<F: PrimeField>`
- `ComposedSumcheck<F: PrimeField>`
- `ComposedSumcheckProof<F: PrimeField>`
- `MultiComposedSumcheckProver`
- `MultiComposedSumcheckVerifier`

## Usage
Here's a basic example of how to use the Sumcheck protocol:
```rs
let poly = Multilinear::new(); // Initialize your polynomial
let mut sumcheck = Sumcheck::new(poly);
sumcheck.poly_sum();
let (proof, challenges) = sumcheck.prove();
let is_valid = sumcheck.verify(&proof);
```

## Contributing
Contributions are welcome! Please feel free to submit a Pull Request.
## License
This project is licensed under the MIT License.
## Disclaimer
This implementation is for educational purposes and has not been audited for production use.