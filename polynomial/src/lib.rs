mod composed;
pub mod interface;
mod multilinear;
mod univariate;
pub mod utils;

pub use composed::multilinear_poly::ComposedMLE;
pub use interface::UnivariatePolynomialTrait;
pub use multilinear::evaluation_form::MLE;
pub use univariate::{UnivariateMonomial, UnivariatePolynomial};
