use ark_ff::PrimeField;
use ark_test_curves::pairing::Pairing;
use sha2::{Digest, Sha256};

#[derive(Debug)]
pub struct MerlinTranscript {
    hasher: Sha256,
    buffer: Vec<u8>,
}

impl MerlinTranscript {
    pub fn new(label: &[u8]) -> Self {
        let mut hasher = Sha256::new();
        hasher.update(b"Merlin Transcript");
        hasher.update(label);

        Self {
            hasher,
            buffer: Vec::new(),
        }
    }

    pub fn append_message(&mut self, label: &[u8], message: &[u8]) {
        self.hasher.update(label);
        let message_length = message.len().to_le_bytes();
        self.hasher.update(message_length);
        self.hasher.update(message);
    }

    pub fn append_scalar<F: PrimeField>(&mut self, label: &[u8], scalar: &F) {
        self.buffer.clear();
        scalar.serialize_compressed(&mut self.buffer).unwrap();
        let message = self.buffer.clone();
        self.append_message(label, &message);
    }

    pub fn append_point<P: Pairing>(&mut self, label: &[u8], point: &P::G1) {
        let stringed_message = point.to_string();
        let message = stringed_message.as_bytes();
        self.append_message(label, message);
    }

    pub fn challenge<F: PrimeField>(&mut self, label: &[u8]) -> F {
        let challenge = self.hasher.finalize_reset();
        self.hasher.update(&label);
        let challenge_bytes: [u8; 32] = challenge.into();

        F::from_be_bytes_mod_order(&challenge_bytes)
    }

    pub fn challenge_n<F: PrimeField>(&mut self, label: &[u8], n: usize) -> Vec<F> {
        let mut result = Vec::new();

        for _ in 0..n {
            // for i in 0..n {
            // let tt = format!("{:?} {}", label, i);
            // result.push(self.challenge(tt.as_bytes()));
            result.push(self.challenge(label));
        }
        result
    }
}

impl Clone for MerlinTranscript {
    fn clone(&self) -> Self {
        Self {
            hasher: self.hasher.clone(),
            buffer: self.buffer.clone(),
        }
    }
}

impl Default for MerlinTranscript {
    fn default() -> Self {
        Self::new(b"default")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::Zero;
    use ark_test_curves::bls12_381::Fr;

    #[test]
    fn test_transcript() {
        let mut transcript = MerlinTranscript::new(b"test_protocol");

        transcript.append_message(b"public_input", b"hello, world");
        transcript.append_scalar::<Fr>(&b"secret_scalar"[..], &Fr::from(42u64));

        let challenge = transcript.challenge::<Fr>(b"challenge");
        assert_ne!(challenge, Fr::zero());
    }
}
