pub fn size_of_mle_n_var_at_each_layer(layer_index: usize) -> usize {
    if layer_index == 0 {
        return 1 << 3;
    }

    let layer_index_plus_one = layer_index + 1;
    let number_of_variable = layer_index + (2 * layer_index_plus_one);

    1 << number_of_variable
}

pub fn transform_label_to_binary_and_to_decimal(
    layer_index: usize,
    a: usize,
    b: usize,
    c: usize,
) -> usize {
    let a_binary_string: String = binary_string(a, layer_index);
    let b_binary_string: String = binary_string(b, layer_index + 1);
    let c_binary_string: String = binary_string(c, layer_index + 1);

    let combined_binary_string = a_binary_string + &b_binary_string + &c_binary_string;

    usize::from_str_radix(&combined_binary_string, 2).unwrap_or(0)
}

/// Convert a number to a binary string of a given size
pub fn binary_string(index: usize, mut bit_count: usize) -> String {
    if bit_count == 0 {
        bit_count = 1;
    }
    let binary = format!("{:b}", index);
    "0".repeat(bit_count.saturating_sub(binary.len())) + &binary
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_size_of_mle_n_var() {
        assert_eq!(size_of_mle_n_var_at_each_layer(0), 8);
        assert_eq!(size_of_mle_n_var_at_each_layer(1), 32);
        assert_eq!(size_of_mle_n_var_at_each_layer(2), 256);
        assert_ne!(size_of_mle_n_var_at_each_layer(2), 128);
        assert_ne!(size_of_mle_n_var_at_each_layer(3), 512);
        assert_eq!(size_of_mle_n_var_at_each_layer(3), 2048);
        assert_eq!(size_of_mle_n_var_at_each_layer(4), 16384);
    }

    #[test]
    fn test_transform_binary_and_to_decimal() {
        // a at layer 0, b & c at layer 1
        assert_eq!(transform_label_to_binary_and_to_decimal(1, 1, 2, 3), 27);
        assert_eq!(transform_label_to_binary_and_to_decimal(2, 1, 2, 3), 83);
    }

    #[test]
    fn test_binary_string() {
        assert_eq!(binary_string(0, 0), "0");
        assert_eq!(binary_string(0, 1), "0");
        assert_eq!(binary_string(0, 2), "00");
        assert_eq!(binary_string(5, 3), "101");
    }
}
