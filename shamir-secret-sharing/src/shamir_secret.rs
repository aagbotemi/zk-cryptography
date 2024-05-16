use ark_ff::PrimeField;
use polynomial::*;
use rand::thread_rng;

pub fn create_shares<F: PrimeField>(
    secret: F,
    threshold: usize,
    total_shares: usize,
) -> Vec<(F, F)> {
    let mut rng = thread_rng();

    let mut secret_shares: Vec<(F, F)> = Vec::with_capacity(threshold);

    for i in 0..threshold {
        let index = F::from(i as u64);
        if i == 0 {
            secret_shares.push((index, secret));
        } else {
            secret_shares.push((index, F::rand(&mut rng)));
        }
    }

    let polynom: UnivariatePolynomial<F> = UnivariatePolynomial::interpolation(&secret_shares);

    let mut shares: Vec<(F, F)> = Vec::with_capacity(threshold);
    for i in 1..=total_shares {
        let evaluation = UnivariatePolynomial::evaluate(&polynom, F::from(i as u64));
        shares.push((F::from(i as u64), evaluation));
    }

    shares
}

pub fn reconstruct_secret<F: PrimeField>(shares: &[(F, F)], point: F) -> F {
    let polynom: UnivariatePolynomial<F> = UnivariatePolynomial::interpolation(shares);
    let evaluation: F = UnivariatePolynomial::evaluate(&polynom, point);
    evaluation
}

#[cfg(test)]
mod tests {
    use ark_ff::MontConfig;
    use ark_ff::{Fp64, MontBackend};

    use crate::shamir_secret::{create_shares, reconstruct_secret};

    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    struct FqConfig;
    type Fq = Fp64<MontBackend<FqConfig, 1>>;

    #[test]
    fn test_create_shares_and_reconstruct_secret() {
        let secret = Fq::from(20);
        let threshold = 3;
        let total_shares = 5;

        let shares = create_shares(secret, threshold, total_shares);
        let picked_points = shares[..threshold].to_vec();

        let reconstruct_at = Fq::from(0);
        let reconstructed_secret = reconstruct_secret(&picked_points, reconstruct_at);

        assert_eq!(secret, reconstructed_secret)
    }

    #[test]
    fn test_create_shares_and_reconstruct_secret_2() {
        let secret = Fq::from(102057);
        let threshold = 17; // threshold should be equal to or less than the modulus
        let total_shares = 40;

        let shares = create_shares(secret, threshold, total_shares);
        let picked_points = shares[..threshold].to_vec();

        let reconstruct_at = Fq::from(0);
        let reconstructed_secret = reconstruct_secret(&picked_points, reconstruct_at);

        assert_eq!(secret, reconstructed_secret)
    }

    #[test]
    fn test_create_shares_and_reconstruct_secret_fail_with_points_above_and_below_threshold() {
        let secret = Fq::from(20);
        let threshold = 3;
        let total_shares = 5;

        let shares = create_shares(secret, threshold, total_shares);

        let picked_points_above_threshold = shares[..4].to_vec();
        let picked_points_below_threshold = shares[..2].to_vec();

        let reconstruct_at = Fq::from(0);
        let reconstructed_secret_above_threshold =
            reconstruct_secret(&picked_points_above_threshold, reconstruct_at);
        let reconstructed_secret_below_threshold =
            reconstruct_secret(&picked_points_below_threshold, reconstruct_at);

        assert_ne!(secret, reconstructed_secret_above_threshold);
        assert_ne!(secret, reconstructed_secret_below_threshold);
    }
}
