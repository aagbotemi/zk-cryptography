# Multilinear KZG Commitment Scheme
This project implements a Multilinear KZG (Kate-Zaverucha-Goldberg) commitment scheme. The implementation includes a trusted setup, commitment generation, opening proofs, and verification.

## Overview
The Multilinear KZG commitment scheme is a cryptographic primitive used in various zero-knowledge proof systems and other cryptographic protocols. This implementation provides the following main components:

1. Trusted Setup
2. Commitment Generation
3. Opening Proofs
4. Verification

## Project Structure
The project is organized into the following modules:

- `interface.rs`: Defines the traits `MultilinearKZGInterface` and `TrustedSetupInterface`.
- `kzg.rs`: Implements the `MultilinearKZG` struct and its methods.
- `trusted_setup.rs`: Implements the `TrustedSetup` struct and its methods.
- `utils.rs`: Contains utility functions used throughout the project.

## Usage
To use this library in your project, you'll need to add it as a dependency and import the necessary modules. Here's a basic example of how to use the Multilinear KZG commitment scheme:

```rust
// Generate trusted setup
let eval_points = vec![Fr::from(2), Fr::from(3), Fr::from(4)];
let tau = TrustedSetup::<Bls12_381>::setup(&eval_points);

// Create a multilinear polynomial
let poly = Multilinear::new(vec![
    Fr::from(0), Fr::from(7), Fr::from(0), Fr::from(5),
    Fr::from(0), Fr::from(7), Fr::from(4), Fr::from(9),
]);

// Generate commitment
let commit = MultilinearKZG::<Bls12_381>::commitment(&poly, &tau.powers_of_tau_in_g1);

// Generate opening proof
let verifier_points = vec![Fr::from(5), Fr::from(9), Fr::from(6)];
let proof = MultilinearKZG::<Bls12_381>::open(&poly, &verifier_points, &tau.powers_of_tau_in_g1);

// Verify the proof
let verify_status = MultilinearKZG::<Bls12_381>::verify(
    &commit,
    &verifier_points,
    &proof,
    tau.powers_of_tau_in_g2,
);

assert!(verify_status, "Verification failed");
```

## Testing
The project includes unit tests for various components. You can run the tests using:
```
cargo test
```