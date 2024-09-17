use ark_bls12_381::{Fr as ScalarField, G1Affine};
use ark_ff::Field;
use ark_serialize::CanonicalSerialize;
use blake2::{
    digest::{consts::U64, generic_array::GenericArray},
    Blake2b, Digest,
};

use crate::interface::SchnorrError;

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
