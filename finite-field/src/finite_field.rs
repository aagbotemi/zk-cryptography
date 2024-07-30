use crate::utils::mod_pow;
use std::ops::{Add, Div, Mul, Sub};

trait FiniteFieldTrait {
    fn modulus(&self) -> usize;
    fn inverse(&self) -> Option<FiniteField>;
    fn pow(&self, exponent: usize) -> FiniteField;
    fn sqrt(&self) -> Option<FiniteField>;
}

#[derive(Debug)]
pub struct FiniteField {
    value: usize,
    modulus: usize,
}

impl FiniteField {
    fn new(value: usize, modulus: usize) -> Self {
        Self {
            value: value % modulus,
            modulus,
        }
    }
}

impl FiniteFieldTrait for FiniteField {
    fn modulus(&self) -> usize {
        self.modulus
    }

    fn inverse(&self) -> Option<FiniteField> {
        for i in 1..self.modulus {
            if (self.value * i) % self.modulus == 1 {
                return Some(FiniteField::new(i, self.modulus));
            }
        }
        None
    }

    fn pow(&self, exponent: usize) -> FiniteField {
        assert!(exponent > 0, "Exponent must be non-negative");
        if exponent == 0 {
            return FiniteField::new(1, self.modulus);
        }
        let result_value = mod_pow(self.value, exponent, self.modulus);
        FiniteField::new(result_value, self.modulus)
    }

    fn sqrt(&self) -> Option<Self> {
        if self.value <= 0 {
            return None;
        }

        let result_value = (self.value as f64).sqrt() as usize;

        Some(FiniteField::new(result_value, self.modulus))
    }
}

impl Add for FiniteField {
    type Output = Self;
    fn add(self, other: FiniteField) -> Self {
        assert_eq!(
            self.modulus, other.modulus,
            "Add Operation should be within the same FiniteField"
        );
        FiniteField {
            value: (self.value + other.value) % self.modulus,
            modulus: self.modulus,
        }
    }
}

impl Sub for FiniteField {
    type Output = Self;

    fn sub(self, other: FiniteField) -> Self {
        assert_eq!(
            self.modulus, other.modulus,
            "Sub Operation should be within the same FiniteField"
        );

        let value = if self.value >= other.value {
            self.value - other.value
        } else {
            self.modulus - (other.value - self.value) % self.modulus
        };

        FiniteField {
            value: value % self.modulus,
            modulus: self.modulus,
        }
    }
}

impl Mul for FiniteField {
    type Output = Self;
    fn mul(self, other: FiniteField) -> Self {
        assert_eq!(
            self.modulus, other.modulus,
            "Mul Operation should be within the same FiniteField"
        );

        FiniteField {
            value: (self.value * other.value) % self.modulus,
            modulus: self.modulus,
        }
    }
}

impl Div for FiniteField {
    type Output = Self;
    fn div(self, other: FiniteField) -> Self {
        assert_eq!(
            self.modulus, other.modulus,
            "Div Operation should be within the same FiniteField"
        );
        assert_ne!(other.value, 0, "Division by zero");

        let inverse = other.inverse().expect("No multiplicative inverse exists");
        self * inverse
    }
}

impl PartialEq for FiniteField {
    fn eq(&self, other: &Self) -> bool {
        assert_eq!(
            self.modulus, other.modulus,
            "You can only compare between same FiniteField"
        );
        self.value == other.value
    }
}

#[cfg(test)]
mod tests {
    use crate::finite_field::FiniteFieldTrait;
    use super::*;

    #[test]
    fn test_sqrt_and_pow() {
        // square root
        let field_1 = FiniteField::new(28, 6);
        let sqrt_result = field_1.sqrt();
        let expected_sqrt_result = Some(FiniteField::new(2, 6));
        assert_eq!(sqrt_result, expected_sqrt_result);

        // raise to pow
        let field_2 = FiniteField::new(2, 9);
        let pow_result = field_2.pow(4);
        let expected_pow_result = FiniteField::new(7, 9);
        assert_eq!(pow_result, expected_pow_result);
    }

    #[test]
    fn test_add_mul_sub_div_eq_and_modulus() {
        // addition
        let field_1 = FiniteField::new(15, 10);
        let field_2 = FiniteField::new(12, 10);
        let expected_add_field = FiniteField::new(7, 10);
        assert_eq!(field_1 + field_2, expected_add_field);

        // modulus
        let modulus_1 = expected_add_field.modulus();
        assert_eq!(modulus_1, 10);

        // multiplication
        let field_3 = FiniteField::new(15, 10);
        let field_4 = FiniteField::new(3, 10);
        let expected_mul_field = FiniteField::new(5, 10);
        assert_eq!(field_3 * field_4, expected_mul_field);

        // subtraction
        let field_5 = FiniteField::new(10, 3);
        let field_6 = FiniteField::new(2, 3);
        let expected_sub_field = FiniteField::new(8, 3);
        assert_eq!(field_5 - field_6, expected_sub_field);

        // modulus
        let modulus_2 = expected_sub_field.modulus();
        assert_eq!(modulus_2, 3);

        // division
        let field_7 = FiniteField::new(6, 3);
        let field_8 = FiniteField::new(2, 3);
        let expected_div_field = FiniteField::new(0, 3);
        assert_eq!(field_7 / field_8, expected_div_field);

        // equality
        let field_1 = FiniteField::new(6, 3);
        let field_2 = FiniteField::new(6, 3);
        assert_eq!(field_1, field_2);

        // equality
        let field_3 = FiniteField::new(9, 3);
        let field_4 = FiniteField::new(5, 3);
        assert_ne!(field_3, field_4)
    }

    #[test]
    #[should_panic(expected = "Add Operation should be within the same FiniteField")]
    fn test_add_fails() {
        let field_1 = FiniteField::new(15, 10);
        let field_2 = FiniteField::new(12, 9);
        let expected_field = FiniteField::new(7, 10);
        assert_eq!(field_1 + field_2, expected_field);
    }

    #[test]
    #[should_panic(expected = "Mul Operation should be within the same FiniteField")]
    fn test_mul_fails() {
        let field_1 = FiniteField::new(15, 10);
        let field_2 = FiniteField::new(3, 7);
        let expected_field = FiniteField::new(5, 10);
        assert_eq!(field_1 * field_2, expected_field)
    }
}
