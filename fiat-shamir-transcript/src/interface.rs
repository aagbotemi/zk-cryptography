use ark_ff::PrimeField;

pub trait FiatShamirTranscriptTrait {
    fn new() -> Self;
    fn commit(&mut self, new_data: &[u8]);
    fn challenge(&mut self) -> [u8; 32];
    fn evaluate_challenge_into_field<F: PrimeField>(&mut self) -> F;
    fn evaluate_n_challenge_into_field<F: PrimeField>(&mut self, n: &usize) -> Vec<F>;
}
