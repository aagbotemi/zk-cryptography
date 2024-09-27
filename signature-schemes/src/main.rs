use signature_schemes::{
    interface::SchnorrSigTrait,
    schnorr::{SchnorrPublicKey, SchnorrSig, SchnorrSignature},
};

fn main() {
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
