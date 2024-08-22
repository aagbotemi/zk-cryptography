# EllipticCurve
This module implements a short Weierstrass form of Elliptic Curve Cryptography (ECC), which is widely used in cryptographic protocols for secure communication.

## Overview
Elliptic Curve Cryptography (ECC) is a public-key cryptography approach based on the algebraic structure of elliptic curves over finite fields. It provides similar security to traditional public-key cryptography but with smaller key sizes.


## Usage
```rs
use field::field::{Field, FieldTrait};

// Define the curve parameters
let a = Field::new(2, 17); // Example field element
let b = Field::new(3, 17); // Example field element
let curve = EllipticCurve::new(a, b);

// Create a point on the curve
let x = Field::new(5, 17); // Example field element
let y = Field::new(1, 17); // Example field element
let point = ECPoint::new(x, y, curve);

// Check if the point is on the curve
assert!(curve.is_on_curve(&point));

// Perform point addition
let point_b = ECPoint::new(Field::new(6, 17), Field::new(3, 17), curve);
let result = curve.add(&point, &point_b).unwrap();
println!("Result of addition: {:?}", result);

// Perform point doubling
let doubled_point = curve.double(&point).unwrap();
println!("Result of doubling: {:?}", doubled_point);

// Perform scalar multiplication
let scalar = 3;
let multiplied_point = curve.scalar_multiplication(&point, scalar).unwrap();
println!("Result of scalar multiplication: {:?}", multiplied_point);
```


## License
This project is licensed under the MIT License.

## Disclaimer
This implementation is for educational purposes. For use in production systems, please ensure it meets your specific security requirements and has been properly audited.