pub fn size_of_mle_n_var_at_each_layer(layer_index: usize) -> usize {
    if layer_index == 0 {
        return 1 << 3;
    }

    let layer_index_plus_one = layer_index + 1;
    let number_of_variable = layer_index + (2 * layer_index_plus_one);

    1 << number_of_variable
}

pub fn transform_label_to_binary_and_to_decimal(a: usize, b: usize, c: usize) -> usize {
    let a_binary: String = format!("{:b}", a);
    let b_binary = format!("{:02b}", b);
    let c_binary = format!("{:02b}", c);

    let combined_binary = format!("{}{}{}", a_binary, b_binary, c_binary);
    usize::from_str_radix(&combined_binary, 2).unwrap_or(0)
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
        assert_eq!(transform_label_to_binary_and_to_decimal(0, 0, 1), 1);
        // a at layer 1, b & c at layer 2
        assert_eq!(transform_label_to_binary_and_to_decimal(1, 2, 3), 27);
        // a at layer 2, b & c at layer 3
        assert_eq!(transform_label_to_binary_and_to_decimal(3, 6, 7), 247); 
    }
}
