use crate::interface::FiatShamirTranscriptTrait;
use ark_ff::PrimeField;
use sha2::{Digest, Sha256};

#[derive(Debug, Default)]
pub struct FiatShamirTranscript {
    hasher: Sha256,
}

impl FiatShamirTranscriptTrait for FiatShamirTranscript {
    fn new() -> Self {
        Self {
            hasher: Sha256::new(),
        }
    }

    fn commit(&mut self, new_data: &[u8]) {
        self.hasher.update(new_data);
    }

    fn challenge(&mut self) -> [u8; 32] {
        let response = self.hasher.finalize_reset();
        self.hasher.update(&response);
        response.into()
    }

    fn evaluate_challenge_into_field<F: PrimeField>(&mut self) -> F {
        F::from_be_bytes_mod_order(&self.challenge())
    }

    fn evaluate_n_challenge_into_field<F: PrimeField>(&mut self, n: &usize) -> Vec<F> {
        let mut result = Vec::new();

        for _ in 0..*n {
            result.push(F::from_be_bytes_mod_order(&self.challenge()));
        }

        result
    }
}
