# GKR Protocol Implementation
This repository contains an implementation of the GKR (Goldwasser-Kalai-Rothblum) protocol in Rust. The GKR protocol is an efficient interactive proof system for verifiable computation.

## Overview
The implementation includes the following main components:

- **GKR Protocol:** The core implementation of the GKR protocol.
- **Circuit Representation:** A structure to represent arithmetic circuits.
- **Gate Types:** Support for addition and multiplication gates.
- **Utility Functions:** Helper functions for various computations.

## Features

- **Proof Generation:** Create proofs for correct circuit evaluation.
- **Verification:** Verify the correctness of generated proofs.
- **Circuit Evaluation:** Evaluate arithmetic circuits on given inputs.
- **Multilinear Extension:** Compute multilinear extensions of circuit layers.

## Main Structures

- `GKRProof<F: PrimeField>`
- `GKRProtocol`
- `Circuit`
- `CircuitLayer`
- `Gate`

## Usage
Here's a basic example of how to use the GKR protocol:
```rs
let circuit = Circuit::new(/* ... */);
let input = vec![/* ... */];

// Generate proof
let proof = GKRProtocol::prove(&circuit, &input);

// Verify proof
let is_valid = GKRProtocol::verify(&circuit, &input, &proof);
```

## Implementation Details

- The protocol uses the Fiat-Shamir heuristic for non-interactive proofs.
- Supports arithmetic circuits with addition and multiplication gates.
- Utilizes multilinear polynomials and sumcheck protocols.

## Contributing
Contributions are welcome! Please feel free to submit a Pull Request.
## License
This project is licensed under the MIT License.
## Disclaimer
This implementation is for educational and research purposes. It has not been audited for production use.