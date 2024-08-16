pub fn bits(scalar: usize) -> usize {
    // (0..(format!("{:02b}", scalar).len() - 1)).rev()
    if scalar == 0 {
        return 0;
    }

    let bits_per_usize = std::mem::size_of::<usize>() * 8;
    let leading_zeros = scalar.leading_zeros() as usize;

    bits_per_usize - leading_zeros
}

/// Returns whether the bit in the given position is set
pub fn bit(scalar: usize, index: usize) -> bool {
    let index_per_usize = std::mem::size_of::<usize>() * 8;
    let digit_index = index / index_per_usize;
    let bit_offset = index % index_per_usize;

    if digit_index == 0 {
        (scalar & (1 << bit_offset)) != 0
    } else {
        false
    }
}
