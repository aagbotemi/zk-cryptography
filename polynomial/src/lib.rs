pub mod composed;
pub mod interface;
pub mod multilinear;
pub mod univariate;
pub mod utils;

pub use composed::composed_multilinear::ComposedMultilinear;
pub use interface::{ComposedMultilinearTrait, MultilinearTrait, UnivariatePolynomialTrait};
pub use multilinear::evaluation_form::Multilinear;
pub use univariate::{
    dense_univariate::DenseUnivariatePolynomial,
    sparse_univariate::{SparseUnivariatePolynomial, UnivariateMonomial},
};

use ark_ff::MontConfig;
use ark_ff::{Fp64, MontBackend};

#[derive(MontConfig)]
#[modulus = "17"]
#[generator = "3"]
pub struct FqConfig;
pub type Fq = Fp64<MontBackend<FqConfig, 1>>;
