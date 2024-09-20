use ark_bls12_381::{Fr as ScalarField, G1Affine};
use ark_ff::Field;
use ark_serialize::CanonicalSerialize;
use blake2::{
    digest::{consts::U64, generic_array::GenericArray},
    Blake2b, Digest,
};
use num_bigint::BigUint;
use num_traits::ToPrimitive;

use crate::interface::{RSAError, SchnorrError};

pub fn hash_message_and_point(
    message: &[u8],
    point: &G1Affine,
) -> Result<ScalarField, SchnorrError> {
    let mut hasher = Blake2b::new();
    hasher.update(message);

    let mut serialized_point = Vec::new();
    point
        .serialize_compressed(&mut serialized_point)
        .map_err(|_| SchnorrError::SerializationError("Serialization failed".to_owned()))?;

    hasher.update(&serialized_point);

    for i in 0..100 {
        let hash_result: GenericArray<u8, U64> = hasher.finalize_reset();
        hasher.update(&[i as u8]);

        if let Some(scalar) = ScalarField::from_random_bytes(&hash_result) {
            return Ok(scalar);
        }
    }

    Err(SchnorrError::ScalarConversionError(
        "Failed to convert bytes to scalar field after multiple attempts".to_owned(),
    ))
}

/// Convert BigUint to i64 safely, handling conversion errors
pub fn convert_biguint_i64(n: &BigUint) -> Result<i64, RSAError> {
    n.to_i64()
        .ok_or_else(|| RSAError::ConversionFailure("Failed to convert BigUint to i64".into()))
}

/// Euler's totient function, Φ(n) = (p-1)(q-1)
pub fn euler_totient(p: &BigUint, q: &BigUint) -> BigUint {
    let one_val = BigUint::from(1u32);

    (p - &one_val) * (q - &one_val)
}

/// Euclidean algorithm to check if gcd(a, b) = 1
pub fn euclidean_algorithm(a: &BigUint, b: &BigUint) -> Result<bool, RSAError> {
    let mut x = a.clone();
    let mut y = b.clone();
    while y != BigUint::from(0u32) {
        let temp = y.clone();
        y = x % y;
        x = temp;
    }

    if x != BigUint::from(1u32) {
        return Err(RSAError::OperationFailure("GCD is not 1".into()));
    }
    Ok(true)
}
