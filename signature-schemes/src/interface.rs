use crate::schnorr::{SchnorrPrivateKey, SchnorrPublicKey, SchnorrSignature};

#[derive(Debug, PartialEq)]
pub enum SchnorrError {
    SerializationError(String),
    ScalarConversionError(String),
    InvalidPublicKey(String),
    VerificationFailed(String),
    InvalidSignature(String),
}

#[derive(Debug)]
pub enum RSAError {
    BadArgument(String),
    OperationFailure(String),
    ConversionFailure(String),
}

pub trait SchnorrSigTrait {
    fn generate_keypair() -> Result<(SchnorrPrivateKey, SchnorrPublicKey), SchnorrError>;

    fn sign(
        private_key: &SchnorrPrivateKey,
        message: &[u8],
    ) -> Result<SchnorrSignature, SchnorrError>;

    fn verify(
        public_key: &SchnorrPublicKey,
        message: &[u8],
        signature: &SchnorrSignature,
    ) -> Result<bool, SchnorrError>;

    fn batch_verify(
        public_keys: &Vec<SchnorrPublicKey>,
        messages: &[&[u8]],
        signatures: &[SchnorrSignature],
    ) -> Result<bool, SchnorrError>;
}
