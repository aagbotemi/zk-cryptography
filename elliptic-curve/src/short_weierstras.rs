use field::field::{Field, FieldTrait};

use crate::utils::{bit, bits};

// An elliptic curve is denoted by the equation
// y^2 = x^3 + ax + b -> Short weierstras
// y^2 + a1xy + a3y = x^3 + a2x^2 + a4x + a6 -> General weierstras
// If 4a^3 + 27b^2 != 0 the curve is non singular
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct EllipticCurve {
    pub a: Field,
    pub b: Field,
}

#[derive(Debug)]
pub enum EllipticCurveError {
    InvalidPoint(ECPoint),
    InvalidScalar(usize),
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ECPoint {
    pub x: Field,
    pub y: Field,
    pub curve: EllipticCurve,
    pub is_infinity: bool,
}

pub trait EllipticCurveTrait {
    fn ec_point(&self, x: Field, y: Field) -> ECPoint;
    fn is_on_curve(&self, ec_point: &ECPoint) -> bool;
    fn add(&self, point_a: &ECPoint, point_b: &ECPoint) -> Result<ECPoint, EllipticCurveError>;
    fn double(&self, point_a: &ECPoint) -> Result<ECPoint, EllipticCurveError>;
    fn scalar_multiplication(
        &self,
        point: &ECPoint,
        scalar: usize,
    ) -> Result<ECPoint, EllipticCurveError>;
}

impl ECPoint {
    pub fn new(x: Field, y: Field, curve: EllipticCurve) -> Self {
        Self {
            x,
            y,
            curve,
            is_infinity: false,
        }
    }
}

impl EllipticCurve {
    pub fn new(a: Field, b: Field) -> EllipticCurve {
        EllipticCurve { a, b }
    }
    pub fn zero(curve: &EllipticCurve) -> ECPoint {
        ECPoint {
            x: Field::new(0, curve.a.modulus()),
            y: Field::new(0, curve.a.modulus()),
            curve: *curve,
            is_infinity: true,
        }
    }

    /// x3 = λ^2 - x1 - x2 mod p
    /// y3 = λ(x1 - x3) - y1 mod p
    ///
    fn compute_x3_y3_add(
        &self,
        slope: &Field,
        x1: &Field,
        y1: &Field,
        x2: &Field,
    ) -> (Field, Field) {
        // x3
        let λ2 = slope.pow(2);
        let λ2_minus_x1 = λ2 - *x1;
        let x3 = λ2_minus_x1 - *x2;

        // y3
        let x1_minus_x3 = *x1 - x3;
        let λ_multiplied_by_x1_minus_x3 = *slope * x1_minus_x3;
        let y3 = λ_multiplied_by_x1_minus_x3 - *y1;

        (x3, y3)
    }

    /// x3 = λ^2 - 2x1 mod p
    /// y3 = λ(x1 - x3) - y1 mod p
    ///
    fn compute_x3_y3_double(
        &self,
        slope: &Field,
        x1: &Field,
        y1: &Field,
        two: &Field,
    ) -> (Field, Field) {
        // x3
        let λ2 = slope.pow(2);
        let two_x1 = *two * *x1;
        let x3 = λ2 - two_x1;

        // y3
        let x1_minus_x3 = *x1 - x3;
        let λ_multiplied_by_x1_minus_x3 = *slope * x1_minus_x3;
        let y3 = λ_multiplied_by_x1_minus_x3 - *y1;

        (x3, y3)
    }
}

impl EllipticCurveTrait for EllipticCurve {
    fn ec_point(&self, x: Field, y: Field) -> ECPoint {
        ECPoint::new(x, y, *self)
    }

    /// Check if a point A = (x,y) belongs to the elliptic curve:
    ///
    /// if `y^2 = x^3 + ax + b mod p` then returns `true`, if not, returns `false`.
    ///
    fn is_on_curve(&self, ec_point: &ECPoint) -> bool {
        if ec_point.is_infinity {
            return true;
        }

        let y2 = ec_point.y * ec_point.y;
        let x3 = ec_point.x * ec_point.x * ec_point.x;
        let ax = self.a * ec_point.x;

        let rhs = x3 + ax + self.b;

        y2 == rhs
    }

    // C(x3,y3) = A(x1,y1) + B(x2,y2)
    //
    // B != A => λ = (y2 - y1) / (x2 - x1) mod p
    // B == A => λ = (3 * x1^2 + a) / (2 * y1) mod p // point doubling
    //
    //     (y2 - y1)
    // λ = --------- mod p
    //     (x2 - x1)
    //
    fn add(&self, point_a: &ECPoint, point_b: &ECPoint) -> Result<ECPoint, EllipticCurveError> {
        if !self.is_on_curve(&point_a) {
            return Err(EllipticCurveError::InvalidPoint(*point_a));
        }

        if !self.is_on_curve(&point_b) {
            return Err(EllipticCurveError::InvalidPoint(*point_b));
        }

        if point_a.is_infinity {
            return Ok(*point_b);
        }
        if point_b.is_infinity {
            return Ok(*point_a);
        }

        if point_a == point_b {
            return self.double(point_a);
        }

        // calculate the slope
        let numerator = point_b.y - point_a.y;
        let denominator = point_b.x - point_a.x;
        let slope = numerator / denominator;

        let (x3, y3) = self.compute_x3_y3_add(&slope, &point_a.x, &point_a.y, &point_b.x);

        let new_point = ECPoint::new(x3, y3, *self);
        assert!(self.is_on_curve(&new_point));

        Ok(new_point)
    }

    //     (3 * x1^2 + a)
    // λ = -------------- mod p
    //        (2 * y1)
    fn double(&self, point_a: &ECPoint) -> Result<ECPoint, EllipticCurveError> {
        if !self.is_on_curve(&point_a) {
            return Err(EllipticCurveError::InvalidPoint(*point_a));
        }

        if point_a.is_infinity {
            return Ok(*point_a);
        }

        let x1_2 = point_a.x.pow(2);
        let two_as_a_field = Field::new(2, x1_2.modulus());
        let three_as_a_field = Field::new(3, x1_2.modulus());

        // numerator
        let three_multiplied_by_x1_2 = three_as_a_field * x1_2;
        let numerator = three_multiplied_by_x1_2 + self.a;

        // denominator
        let denominator = two_as_a_field * point_a.y;

        // slope
        let slope = numerator / denominator;
        let (x3, y3) = self.compute_x3_y3_double(&slope, &point_a.x, &point_a.y, &two_as_a_field);

        let new_point = ECPoint::new(x3, y3, *self);
        assert!(self.is_on_curve(&new_point));

        Ok(new_point)
    }

    fn scalar_multiplication(
        &self,
        point: &ECPoint,
        scalar: usize,
    ) -> Result<ECPoint, EllipticCurveError> {
        if !self.is_on_curve(point) {
            return Err(EllipticCurveError::InvalidPoint(*point));
        }

        if point.is_infinity {
            return Ok(*point);
        }

        if scalar == 0 {
            return Err(EllipticCurveError::InvalidScalar(scalar));
        }

        let mut new_point = *point;

        for i in (0..bits(scalar) - 1).rev() {
            new_point = self.double(&new_point)?;
            if bit(scalar, i) {
                new_point = self.add(&new_point, point)?;
            }
        }

        assert!(self.is_on_curve(&new_point));

        Ok(new_point)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn setup_curve() -> EllipticCurve {
        let a = Field::new(2, 17);
        let b = Field::new(2, 17);
        EllipticCurve::new(a, b)
    }

    #[test]
    fn test_point_in_curve() {
        // y^2 = x^3 + 2x + 2 mod 17
        let ec_curve = setup_curve();

        // (6,3) + (5,1) = (10,6)
        let p1 = ec_curve.ec_point(Field::new(6, 17), Field::new(3, 17));
        let p2 = ec_curve.ec_point(Field::new(5, 17), Field::new(1, 17));
        let p3 = ec_curve.ec_point(Field::new(10, 17), Field::new(6, 17));

        assert!(ec_curve.is_on_curve(&p1));
        assert!(ec_curve.is_on_curve(&p2));
        assert!(ec_curve.is_on_curve(&p3));

        let p4 = ec_curve.ec_point(Field::new(4, 17), Field::new(1, 17));
        let p5 = ec_curve.ec_point(Field::new(1, 17), Field::new(1, 17));
        let p6 = ec_curve.ec_point(Field::new(0, 17), Field::new(1, 17));

        assert!(!ec_curve.is_on_curve(&p4));
        assert!(!ec_curve.is_on_curve(&p5));
        assert!(!ec_curve.is_on_curve(&p6));
    }

    #[test]
    fn test_point_addition() {
        // y^2 = x^3 + 2x + 2 mod 17
        let ec_curve = setup_curve();

        // (6,3) + (5,1) = (10,6)
        let point_1 = ec_curve.ec_point(Field::new(6, 17), Field::new(3, 17));
        let point_2 = ec_curve.ec_point(Field::new(5, 17), Field::new(1, 17));
        let point_result = ec_curve.ec_point(Field::new(10, 17), Field::new(6, 17));

        let result_1 = ec_curve.add(&point_1, &point_2);
        assert_eq!(result_1.unwrap(), point_result);

        // (16,13) + (5,1) = (0,6)
        let point_3 = ec_curve.ec_point(Field::new(16, 17), Field::new(13, 17));
        let point_4 = ec_curve.ec_point(Field::new(5, 17), Field::new(1, 17));
        let point_result_2 = ec_curve.ec_point(Field::new(0, 17), Field::new(6, 17));
        let point_result_3 = ec_curve.ec_point(Field::new(0, 17), Field::new(11, 17));

        let result_2 = ec_curve.add(&point_4, &point_3);
        let result_3 = ec_curve.add(&point_3, &point_3);

        assert_eq!(result_2.unwrap(), point_result_2);
        assert_eq!(result_3.unwrap(), point_result_3);
    }

    #[test]
    fn test_point_double() {
        // y^2 = x^3 + 2x + 2 mod 17
        let ec_curve = setup_curve();

        // (5,1) + (5,1) = 2 (5, 1) = (6,3)
        let point_1 = ec_curve.ec_point(Field::new(5, 17), Field::new(1, 17));
        let point_result = ec_curve.ec_point(Field::new(6, 17), Field::new(3, 17));

        let result = ec_curve.double(&point_1);
        assert_eq!(result.unwrap(), point_result);
    }

    #[test]
    fn test_scalar_multiplication() {
        let ec_curve = setup_curve();

        // 2 * (5,1) = (6,3)
        let point = ec_curve.ec_point(Field::new(5, 17), Field::new(1, 17));
        let point_result = ec_curve.ec_point(Field::new(6, 17), Field::new(3, 17));
        let point_result_2 = ec_curve.ec_point(Field::new(3, 17), Field::new(16, 17));

        let result = ec_curve.scalar_multiplication(&point, 2);
        let result_2 = ec_curve.scalar_multiplication(&point, 15);
        assert_eq!(result.unwrap(), point_result);
        assert_eq!(result_2.unwrap(), point_result_2);
    }

    #[test]
    fn test_point_at_infinity() {
        let modulus = 17;
        let curve = setup_curve();

        // Create a regular point on the curve
        let x = Field::new(6, modulus);
        let y = Field::new(3, modulus);
        let point = ECPoint::new(x, y, curve);

        // Create the point at infinity
        let infinity = EllipticCurve::zero(&curve);
        assert!(infinity.is_infinity);
        assert_eq!(infinity.x, Field::new(0, curve.a.modulus()));
        assert_eq!(infinity.y, Field::new(0, curve.a.modulus()));

        // Test that the point at infinity is on the curve
        assert!(curve.is_on_curve(&infinity));

        // Test addition with point at infinity
        let result = curve.add(&point, &infinity).expect("Addition failed");
        assert_eq!(result, point);

        let result = curve.add(&infinity, &point).expect("Addition failed");
        assert_eq!(result, point);

        // Test doubling the point at infinity
        let result = curve.double(&infinity).expect("Point Doubling failed");
        assert_eq!(result, infinity);

        // Test scalar multiplication with point at infinity
        let result = curve.scalar_multiplication(&infinity, 5).unwrap();
        assert_eq!(result, infinity);
    }
}
