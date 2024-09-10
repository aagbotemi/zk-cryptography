# Shamir's Secret Sharing
This module implements Shamir's Secret Sharing scheme, a cryptographic algorithm that allows a secret to be divided into parts, giving each participant its own unique part.
## Overview
Shamir's Secret Sharing is used to secure a secret in a distributed way, particularly useful when you want to reduce the risk of losing the secret or protect against unauthorized access.

## Usage
```rs
use ark_ff::PrimeField as F;

// Creating shares
let secret = F::from(1234); // Your secret
let threshold = 3; // Minimum shares needed to reconstruct
let total_shares = 5; // Total number of shares to create
let shares = create_shares(secret, threshold, total_shares);

// Reconstructing the secret
let reconstructed_secret = reconstruct_secret(&shares[0..3], F::from(0));
assert_eq!(secret, reconstructed_secret);
```
