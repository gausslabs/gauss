use super::{num::UnsignedInteger, ring::Matrix};

pub trait RandomKeyDist<T> {}

pub trait RandomGaussianDist<T: UnsignedInteger> {
    fn random_vec_in_modulus(modulus: u64) -> Vec<T>;
    fn random_ring_poly(moduli_chain: &[u64]) -> Vec<T>;
}

pub trait RandomUniformDist<T: UnsignedInteger, M: Matrix<T>> {
    fn random_vec_in_modulus(&self, modulus: u64, size: usize) -> Vec<T>;
    fn random_ring_poly(&self, moduli_chain: &[u64], ring_size: usize) -> M;
}
