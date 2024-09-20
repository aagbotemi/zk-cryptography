use modinverse::modinverse;
use num_bigint::BigUint;
use num_traits::FromPrimitive;

use crate::{
    interface::RSAError,
    utils::{convert_biguint_i64, euclidean_algorithm, euler_totient},
};

pub struct RSA {
    p: BigUint,
    q: BigUint,
    pub_key: BigUint,
}

impl RSA {
    pub fn new(p: BigUint, q: BigUint, pub_key: BigUint) -> Self {
        Self { p, q, pub_key }
    }

    /// RSA Encryption: c = m^e mod n
    pub fn encryption(&self, m: &BigUint) -> Result<BigUint, RSAError> {
        let n = &self.p * &self.q;
        let phi_n = euler_totient(&self.p, &self.q);
        let pub_ver_exp = euclidean_algorithm(&self.pub_key, &phi_n).unwrap();

        if pub_ver_exp {
            Ok(m.modpow(&self.pub_key, &n))
        } else {
            Err(RSAError::OperationFailure(
                "Public key exponent does not satisfy conditions".into(),
            ))
        }
    }

    /// RSA Decryption: m = c^d mod n
    pub fn decryption(&self, cipher_text: &BigUint) -> Result<BigUint, RSAError> {
        let n = &self.p * &self.q;
        let phi_n = euler_totient(&self.p, &self.q);

        let d = modinverse(
            convert_biguint_i64(&self.pub_key).unwrap(),
            convert_biguint_i64(&phi_n).unwrap(),
        )
        .ok_or_else(|| RSAError::OperationFailure("Failed to find modular inverse".into()))?;

        let d_big = BigUint::from_i64(d)
            .ok_or_else(|| RSAError::ConversionFailure("Failed to convert d to BigUint".into()))?;
        Ok(cipher_text.modpow(&d_big, &n))
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_encryption_decryption() {
        let rsa = RSA::new(
            BigUint::from(13u32),
            BigUint::from(17u32),
            BigUint::from(35u32),
        );

        let message = BigUint::from(5u32);
        let cipher_text = rsa.encryption(&message).unwrap();
        assert_eq!(cipher_text, BigUint::from(125u32));

        let decrypted = rsa.decryption(&cipher_text).unwrap();
        assert_eq!(decrypted, BigUint::from(5u32));
    }

    #[test]
    fn test_large_prime() {
        let rsa = RSA::new(
            BigUint::from(1223u32),
            BigUint::from(1987u32),
            BigUint::from(948047u32),
        );

        let message = BigUint::from(5u32);
        let cipher_text = rsa.encryption(&message).unwrap();
        assert_eq!(cipher_text, BigUint::from(915542u32));

        let decrypted = rsa.decryption(&cipher_text).unwrap();
        assert_eq!(decrypted, BigUint::from(5u32));
    }
}
