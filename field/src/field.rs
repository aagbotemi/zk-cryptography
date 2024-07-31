use crate::utils::mod_pow;
use std::ops::{Add, Div, Mul, Sub};

pub trait FieldTrait {
    fn modulus(&self) -> usize;
    fn inverse(&self) -> Option<Field>;
    fn pow(&self, exponent: usize) -> Field;
    fn sqrt(&self) -> Option<Field>;
    fn zero(&self) -> Self;
    fn one(&self) -> Self;
}

#[derive(Debug, Clone, Copy, Default)]
pub struct Field {
    value: usize,
    modulus: usize,
}

impl Field {
    pub fn new(value: usize, modulus: usize) -> Self {
        assert!(modulus > 1, "Modulus should be greater than 1");
        Self {
            value: value % modulus,
            modulus,
        }
    }
}

impl FieldTrait for Field {
    fn modulus(&self) -> usize {
        self.modulus
    }

    fn inverse(&self) -> Option<Field> {
        for i in 1..self.modulus {
            if (self.value * i) % self.modulus == 1 {
                return Some(Field::new(i, self.modulus));
            }
        }
        None
    }

    fn pow(&self, exponent: usize) -> Field {
        assert!(exponent > 0, "Exponent must be non-negative");
        if exponent == 0 {
            return Field::new(1, self.modulus);
        }
        let result_value = mod_pow(self.value, exponent, self.modulus);
        Field::new(result_value, self.modulus)
    }

    fn sqrt(&self) -> Option<Self> {
        if self.value <= 0 {
            return None;
        }

        let result_value = (self.value as f64).sqrt() as usize;

        Some(Field::new(result_value, self.modulus))
    }

    fn zero(&self) -> Self {
        Field::new(0, self.modulus())
    }

    fn one(&self) -> Self {
        Field::new(1, self.modulus())
    }
}

impl Add for Field {
    type Output = Self;
    fn add(self, other: Field) -> Self {
        assert_eq!(
            self.modulus, other.modulus,
            "Add Operation should be within the same Field"
        );
        Field {
            value: (self.value + other.value) % self.modulus,
            modulus: self.modulus,
        }
    }
}

impl Sub for Field {
    type Output = Self;

    fn sub(self, other: Field) -> Self {
        assert_eq!(
            self.modulus, other.modulus,
            "Sub Operation should be within the same Field"
        );

        let value = if self.value >= other.value {
            self.value - other.value
        } else {
            self.modulus - (other.value - self.value) % self.modulus
        };

        Field {
            value: value % self.modulus,
            modulus: self.modulus,
        }
    }
}

impl Mul for Field {
    type Output = Self;
    fn mul(self, other: Field) -> Self {
        assert_eq!(
            self.modulus, other.modulus,
            "Mul Operation should be within the same Field"
        );

        Field {
            value: (self.value * other.value) % self.modulus,
            modulus: self.modulus,
        }
    }
}

impl Div for Field {
    type Output = Self;
    fn div(self, other: Field) -> Self {
        assert_eq!(
            self.modulus, other.modulus,
            "Div Operation should be within the same Field"
        );
        assert_ne!(other.value, 0, "Division by zero");

        let inverse = other.inverse().expect("No multiplicative inverse exists");
        self * inverse
    }
}

impl PartialEq for Field {
    fn eq(&self, other: &Self) -> bool {
        assert_eq!(
            self.modulus, other.modulus,
            "You can only compare between same Field"
        );
        self.value == other.value
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sqrt_and_pow() {
        // square root
        let field_1 = Field::new(28, 6);
        let sqrt_result = field_1.sqrt();
        let expected_sqrt_result = Some(Field::new(2, 6));
        assert_eq!(sqrt_result, expected_sqrt_result);

        // raise to pow
        let field_2 = Field::new(2, 9);
        let pow_result = field_2.pow(4);
        let expected_pow_result = Field::new(7, 9);
        assert_eq!(pow_result, expected_pow_result);
    }

    #[test]
    fn test_zero_and_one() {
        let field = Field::new(28, 6);

        let zero = field.zero();
        assert_eq!(zero, Field::new(0, 6));

        let one = field.one();
        assert_eq!(one, Field::new(1, 6));
    }

    #[test]
    fn test_add_mul_sub_div_eq_and_modulus() {
        // addition
        let field_1 = Field::new(15, 10);
        let field_2 = Field::new(12, 10);
        let expected_add_field = Field::new(7, 10);
        assert_eq!(field_1 + field_2, expected_add_field);

        // modulus
        let modulus_1 = expected_add_field.modulus();
        assert_eq!(modulus_1, 10);

        // multiplication
        let field_3 = Field::new(15, 10);
        let field_4 = Field::new(3, 10);
        let expected_mul_field = Field::new(5, 10);
        assert_eq!(field_3 * field_4, expected_mul_field);

        // subtraction
        let field_5 = Field::new(10, 3);
        let field_6 = Field::new(2, 3);
        let expected_sub_field = Field::new(8, 3);
        assert_eq!(field_5 - field_6, expected_sub_field);

        // modulus
        let modulus_2 = expected_sub_field.modulus();
        assert_eq!(modulus_2, 3);

        // division
        let field_7 = Field::new(6, 3);
        let field_8 = Field::new(2, 3);
        let expected_div_field = Field::new(0, 3);
        assert_eq!(field_7 / field_8, expected_div_field);

        // equality
        let field_1 = Field::new(6, 3);
        let field_2 = Field::new(6, 3);
        assert_eq!(field_1, field_2);

        // equality
        let field_3 = Field::new(9, 3);
        let field_4 = Field::new(5, 3);
        assert_ne!(field_3, field_4)
    }

    #[test]
    #[should_panic(expected = "Add Operation should be within the same Field")]
    fn test_add_fails() {
        let field_1 = Field::new(15, 10);
        let field_2 = Field::new(12, 9);
        let expected_field = Field::new(7, 10);
        assert_eq!(field_1 + field_2, expected_field);
    }

    #[test]
    #[should_panic(expected = "Mul Operation should be within the same Field")]
    fn test_mul_fails() {
        let field_1 = Field::new(15, 10);
        let field_2 = Field::new(3, 7);
        let expected_field = Field::new(5, 10);
        assert_eq!(field_1 * field_2, expected_field)
    }
}
