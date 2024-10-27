use ark_ff::PrimeField;

#[derive(Clone, PartialEq, Eq, Default, Debug)]
pub struct Domain<F: PrimeField> {
    /// This is a const size of the domain
    pub(crate) size: u64,
    /// This is the generator of the domain, ofter regarded as the root of unity (omega)
    pub(crate) generator: F,
    /// This is the inverse of the group generator
    pub(crate) group_gen_inverse: F,
    /// This is the inverse of the group size
    pub(crate) group_size_inverse: F,
}

impl<F: PrimeField> Domain<F> {
    pub fn new(num_of_coeffs: usize) -> Domain<F> {
        let size = if num_of_coeffs.is_power_of_two() {
            num_of_coeffs
        } else {
            num_of_coeffs.checked_next_power_of_two().unwrap()
        } as u64;

        let generator = F::get_root_of_unity(size).unwrap();
        let group_gen_inverse = generator.inverse().unwrap();
        let group_size_inverse = F::from(size).inverse().unwrap();

        Domain {
            size,
            generator,
            group_gen_inverse,
            group_size_inverse,
        }
    }
}
