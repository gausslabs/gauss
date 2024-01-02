pub trait ShoupRepresentationFq {
    fn shoup_representation_fq(&self, q: Self) -> u64;
}

impl ShoupRepresentationFq for u64 {
    fn shoup_representation_fq(&self, q: Self) -> u64 {
        ((*self as u128 * (1u128 << 64)) / q as u128) as u64
    }
}
