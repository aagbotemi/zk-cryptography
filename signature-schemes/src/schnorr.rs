use ark_bls12_381::{Fr as ScalarField, G1Affine, G1Projective};
use ark_ec::{CurveGroup, Group};
use ark_ff::UniformRand;
use rand::rngs::ThreadRng;
use rand::thread_rng;
use rayon::prelude::*;
use std::ops::Mul;

use crate::{
    interface::{SchnorrError, SchnorrSigTrait},
    utils::hash_message_and_point,
};

pub struct SchnorrSig;

#[derive(Debug, Clone)]
pub struct SchnorrPublicKey(pub G1Affine);

#[derive(Debug)]
pub struct SchnorrPrivateKey(pub ScalarField);

#[derive(Debug)]
pub struct SchnorrSignature {
    pub r: G1Affine,
    pub sig: ScalarField,
}

impl SchnorrSigTrait for SchnorrSig {
    /// Generate keypair
    fn generate_keypair() -> Result<(SchnorrPrivateKey, SchnorrPublicKey), SchnorrError> {
        let mut rng = thread_rng();
        let private_key = ScalarField::rand(&mut rng);
        let public_key = (G1Projective::generator() * private_key).into_affine();

        if !public_key.is_on_curve() || !public_key.is_in_correct_subgroup_assuming_on_curve() {
            return Err(SchnorrError::InvalidPublicKey(
                "Invalid public key".to_owned(),
            ));
        }

        Ok((SchnorrPrivateKey(private_key), SchnorrPublicKey(public_key)))
    }

    /// Sign a message
    fn sign(
        private_key: &SchnorrPrivateKey,
        message: &[u8],
    ) -> Result<SchnorrSignature, SchnorrError> {
        let mut rng: ThreadRng = thread_rng();

        // a random nonce
        let nonce = ScalarField::rand(&mut rng);
        // compute r = g^nonce
        let r = (G1Projective::generator() * nonce).into_affine();
        // hash message and r to get c
        let c = hash_message_and_point(message, &r)?;
        // signature
        let sig = nonce + c * private_key.0;

        Ok(SchnorrSignature { r, sig })
    }

    /// Verify a signature
    fn verify(
        public_key: &SchnorrPublicKey,
        message: &[u8],
        signature: &SchnorrSignature,
    ) -> Result<bool, SchnorrError> {
        if !public_key.0.is_on_curve() || !public_key.0.is_in_correct_subgroup_assuming_on_curve() {
            return Err(SchnorrError::InvalidPublicKey(
                "Invalid public key".to_owned(),
            ));
        }

        let c = hash_message_and_point(message, &signature.r)?;
        let lhs = G1Projective::generator() * signature.sig;
        let rhs = signature.r + public_key.0.mul(c);

        Ok(lhs == rhs)
    }

    /// Batch verification of signatures
    fn batch_verify(
        public_keys: &Vec<SchnorrPublicKey>,
        messages: &[&[u8]],
        signatures: &[SchnorrSignature],
    ) -> Result<bool, SchnorrError> {
        assert_eq!(public_keys.len(), messages.len(), "Length Mismatch");
        assert_eq!(public_keys.len(), signatures.len(), "Length Mismatch");

        // Perform batch verification
        let verification_results: Vec<Result<bool, SchnorrError>> = public_keys
            .par_iter()
            .zip(messages.par_iter())
            .zip(signatures.par_iter())
            .map(|((pk, msg), sig)| SchnorrSig::verify(pk, msg, sig))
            .collect();

        // Check if all verifications passed
        verification_results
            .into_iter()
            .all(|r| r == Ok(true))
            .then_some(true)
            .ok_or(SchnorrError::InvalidSignature(
                "Signature is Invalid".to_owned(),
            ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sign_verify() {
        let (sk, pk) = SchnorrSig::generate_keypair().unwrap();
        let message = b"Schnorr is a signature scheme";

        let signature = SchnorrSig::sign(&sk, message).unwrap();
        let result = SchnorrSig::verify(&pk, message, &signature).unwrap();

        assert!(result);
    }

    #[test]
    fn test_batch_verify() {
        let messages: [&[u8]; 7] = [
            b"Heyo, I am Abiodun Awoyemi",
            b"I'm a Blockchain engineer and a ZK engineer.",
            b"I build smart contract on EVM compatible blockchain",
            b"I also work on Bitcoin, contributing to Open source",
            b"I have a collection of cryptographic and ZK tools, find it here: https://github.com/aagbotemi/zk-cryptography",
            b"My tech stack: Rust, Solidity, JavaScript, TypeScript",
            b"Let's create a beautiful world!", 
        ];

        let mut signatures: Vec<SchnorrSignature> = Vec::new();
        let mut public_keys: Vec<SchnorrPublicKey> = Vec::new();

        for i in 0..messages.len() {
            let (sk, pk) = SchnorrSig::generate_keypair().unwrap();
            signatures.push(SchnorrSig::sign(&sk, &messages[i]).unwrap());
            public_keys.push(pk);
        }

        let result = SchnorrSig::batch_verify(&public_keys, &messages, &signatures).unwrap();
        assert!(result);
    }

    #[test]
    fn test_batch_verify_tampered_message() {
        let messages: [&[u8]; 2] = [
            b"Heyo, I am Abiodun Awoyemi",
            b"I'm a Blockchain engineer and a ZK engineer.",
        ];

        let mut signatures: Vec<SchnorrSignature> = Vec::new();
        let mut public_keys: Vec<SchnorrPublicKey> = Vec::new();

        for i in 0..messages.len() {
            let (sk, pk) = SchnorrSig::generate_keypair().unwrap();
            signatures.push(SchnorrSig::sign(&sk, &messages[i]).unwrap());
            public_keys.push(pk);
        }

        let tampered_messages: [&[u8]; 2] = [
            b"Heyo, I am Abiodun Awoyemi",
            b"I'm a Blockchain engineer & a ZK engineer.",
        ];

        let result = SchnorrSig::batch_verify(&public_keys, &tampered_messages, &signatures);
        assert_eq!(
            result,
            Err(SchnorrError::InvalidSignature(
                "Signature is Invalid".to_owned()
            ))
        );
    }
}
