mod composed;
pub mod interface;
mod multilinear;
mod univariate;
pub mod utils;

pub use composed::composed_multilinear::ComposedMultilinear;
pub use interface::{ComposedMultilinearTrait, MultilinearTrait, UnivariatePolynomialTrait};
pub use multilinear::evaluation_form::Multilinear;
pub use univariate::univariate::{UnivariateMonomial, UnivariatePolynomial};
