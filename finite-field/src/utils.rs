pub fn mod_pow(base: usize, exponent: usize, modulus: usize) -> usize {
    if modulus == 1 {
        return 0;
    }

    let mut result = 1;
    let mut base = base % modulus;
    let mut exp = exponent;

    while exp > 0 {
        if exp % 2 == 1 {
            result = (result * base) % modulus;
        }
        exp = exp >> 1;
        base = (base * base) % modulus;
    }

    result
}

pub fn check_is_less_than(a: usize, b: usize) -> bool {
    if a < b {
        true
    } else {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mod_pow() {
        let result_1 = mod_pow(2, 4, 5);
        assert_eq!(result_1, 1);

        let result_2 = mod_pow(3, 2, 4);
        assert_eq!(result_2, 1);

        let result_3 = mod_pow(5, 4, 7);
        assert_eq!(result_3, 2);
    }

    #[test]
    fn test_check_is_less_than() {
        let result_1 = check_is_less_than(15, 10);
        assert!(!result_1);
        let result_2 = check_is_less_than(12, 18);
        assert!(result_2);
        let result_3 = check_is_less_than(1505, 1030);
        assert!(!result_3);
    }
}
