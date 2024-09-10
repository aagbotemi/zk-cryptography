# FiatShamirTranscript
This module implements the Fiat-Shamir heuristic, a cryptographic algorithm that transforms interactive proof systems into non-interactive ones using a hash function.

## Overview
The Fiat-Shamir is used to convert interactive proofs into non-interactive proofs, which is particularly useful in cryptographic protocols to reduce interaction between the prover and verifier.

## Usage
```rs
use ark_ff::PrimeField as F;

// Creating a new transcript
let mut transcript = FiatShamirTranscript::new();

// Committing data to the transcript
transcript.commit(b"some data");

// Generating a challenge
let challenge = transcript.challenge();
println!("Challenge: {:?}", challenge);

// Evaluating the challenge into a field element
let field_element: F = transcript.evaluate_challenge_into_field();
println!("Field Element: {:?}", field_element);

// Evaluating multiple challenges into field elements
let field_elements: Vec<F> = transcript.evaluate_n_challenge_into_field(&3);
println!("Field Elements: {:?}", field_elements);
```
