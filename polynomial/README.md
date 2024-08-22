# Polynomial Library
This library provides implementations for various types of polynomials, including univariate polynomials, multilinear polynomials, and composed multilinear polynomials. It is designed to work with prime fields and offers efficient operations for cryptographic applications.

## Features

1. **Univariate Polynomials:** Representation and operations for single-variable polynomials.
2. **Multilinear Polynomials:** Support for polynomials with multiple variables, each of degree at most 1.
3. **Composed Multilinear Polynomials:** Handling of polynomials composed of multiple multilinear polynomials.


Prime Field Operations: All polynomial operations are compatible with prime fields.

## Main Structures
- `UnivariatePolynomial<F: PrimeField>`
- `MultiLinearCoefficientPolynomial<F: PrimeField>`
- `Multilinear<F: PrimeField>`
- `ComposedMultilinear<F: PrimeField>`

## Implementation Details
#### Univariate Polynomial

- Represented as a vector of UnivariateMonomial<F> structures, each containing a coefficient and power.
- Supports operations like evaluation, interpolation, addition, and multiplication.
- Implements efficient degree calculation and polynomial arithmetic.

#### Multilinear Coefficient Polynomial

- Represented as MultiLinearCoefficientPolynomial<F> containing MultiLinearMonomial<F> terms.
- Each monomial has a coefficient and a vector of boolean variables.
Supports partial evaluation and full evaluation over multiple variables.
- Implements degree calculation based on the number of true variables in monomials.

#### Multilinear Polynomial (Evaluation Form)

- Represented as Multilinear<F> with a vector of evaluations over the boolean hypercube.
- Supports efficient partial evaluation and full evaluation.
- Implements operations like addition, scalar multiplication, and distinct polynomial multiplication.

#### Composed Multilinear Polynomial

- Represented as ComposedMultilinear<F> containing a vector of Multilinear<F> polynomials.
- Supports element-wise product and addition of the constituent polynomials.
- Implements partial and full evaluation over the composed structure.


## Usage
```rs
use use ark_ff::PrimeField as F;

// Univariate Polynomial
let poly = UnivariatePolynomial::new(vec![F::from(1), F::from(2), F::from(3)]); // 1 + 2x + 3x^2
let evaluation = poly.evaluate(F::from(2)); // Evaluate at x = 2
let interpolated = UnivariatePolynomial::interpolation(&[(F::from(0), F::from(1)), (F::from(1), F::from(6))]); // Interpolate (0,1) and (1,6)

// Multilinear Polynomial (Coefficient Form)
let term1 = MultiLinearMonomial::new(F::from(2), vec![true, false]); // 2x
let term2 = MultiLinearMonomial::new(F::from(3), vec![true, true]); // 3xy
let multi_coeff_poly = MultiLinearCoefficientPolynomial::new(vec![term1, term2]);
let eval_result = multi_coeff_poly.evaluation(&vec![F::from(2), F::from(3)]); // Evaluate at x=2, y=3

// Multilinear Polynomial (Evaluation Form)
let multi_poly = Multilinear::new(vec![F::from(1), F::from(2), F::from(3), F::from(4)]); // 2-variable polynomial
let partial_eval = multi_poly.partial_evaluation(&F::from(2), &0); // Partial evaluation with x=2

// Composed Multilinear Polynomial
let poly1 = Multilinear::new(vec![F::from(1), F::from(2)]);
let poly2 = Multilinear::new(vec![F::from(3), F::from(4)]);
let composed_poly = ComposedMultilinear::new(vec![poly1, poly2]);
let element_wise_product = composed_poly.element_wise_product();
let evaluation = composed_poly.evaluation(&[F::from(2)]); // Evaluate at x=2
```

## Contributing
Contributions are welcome! Please feel free to submit a Pull Request.
## License
This project is licensed under the MIT License.
## Disclaimer
This library is for educational and research purposes. It has not been audited for production use in cryptographic applications. 