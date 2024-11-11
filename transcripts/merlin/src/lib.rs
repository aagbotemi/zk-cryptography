use ark_ff::PrimeField;
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

    pub fn challenge<F: PrimeField>(&mut self, label: &[u8]) -> F {
        self.hasher.update(label);
        let challenge_bytes = self.hasher.finalize_reset();

        let hasher = Sha256::new();
        self.hasher.update(challenge_bytes);
        let final_bytes = hasher.finalize();

        F::from_random_bytes(&final_bytes).expect("Failed to generate challenge")
    }

    pub fn challenge_n<F: PrimeField>(&mut self, label: &[u8], n: usize) -> Vec<F> {
        let mut result = Vec::new();

        for _ in 0..n {
            result.push(self.challenge(&label));
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
