mod composed;
pub mod interface;
mod multilinear;
mod univariate;
pub mod utils;

pub use composed::multilinear::ComposedMLE;
pub use interface::{
    ComposedMLETrait, MLETrait, MultiLinearPolynomialTrait, UnivariatePolynomialTrait,
};
pub use multilinear::evaluation_form::MLE;
pub use univariate::univariate::{UnivariateMonomial, UnivariatePolynomial};
