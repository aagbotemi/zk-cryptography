pub mod interface;
mod multilinear;
mod univariate;
pub mod utils;

pub use interface::UnivariatePolynomialTrait;
pub use univariate::{UnivariateMonomial, UnivariatePolynomial};
