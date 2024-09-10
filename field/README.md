# Field
This module implements a finite field element, which is a fundamental component in many cryptographic algorithms, including Elliptic Curve Cryptography (ECC).

## Overview
A finite field (or Galois field) is a field that contains a finite number of elements. Finite fields are widely used in cryptography, coding theory, and other areas of computer science and mathematics.

## Usage
```rs
// Creating a new field element
let modulus = 17; // Example modulus
let value = 5; // Example value
let field_element = Field::new(value, modulus);

// Performing field operations
let another_value = 3;
let another_field_element = Field::new(another_value, modulus);

// Addition
let sum = field_element + another_field_element;
println!("Sum: {:?}", sum);

// Subtraction
let difference = field_element - another_field_element;
println!("Difference: {:?}", difference);

// Multiplication
let product = field_element * another_field_element;
println!("Product: {:?}", product);

// Division
let quotient = field_element / another_field_element;
println!("Quotient: {:?}", quotient);

// Inverse
let inverse = field_element.inverse().unwrap();
println!("Inverse: {:?}", inverse);

// Power
let exponent = 3;
let power = field_element.pow(exponent);
println!("Power: {:?}", power);

// Square root
let sqrt = field_element.sqrt().unwrap();
println!("Square Root: {:?}", sqrt);

// Zero and One
let zero = field_element.zero();
let one = field_element.one();
println!("Zero: {:?}", zero);
println!("One: {:?}", one);
```
