use std::{marker::PhantomData, os::unix::thread};

use crate::core_crypto::{
    matrix::{Matrix, MatrixMut},
    num::UnsignedInteger,
    random::RandomUniformDist,
};
use itertools::Itertools;
use ndarray::{iter::Iter, Array2, ArrayBase, Axis, Dim, IndexLonger, ViewRepr};
use rand::{distributions::Uniform, thread_rng, Rng};

pub fn random_vec_in_fq<T: UnsignedInteger + rand::distributions::uniform::SampleUniform>(
    size: usize,
    q: T,
) -> Vec<T> {
    let rng = thread_rng();
    rng.sample_iter(Uniform::new(T::zero(), q))
        .take(size)
        .collect_vec()
}

pub struct TestRng {}

impl RandomUniformDist<u64, Vec<Vec<u64>>> for TestRng {
    fn random_ring_poly(&mut self, moduli_chain: &[u64], ring_size: usize) -> Vec<Vec<u64>> {
        let mut rng = thread_rng();
        moduli_chain
            .iter()
            .map(|qi| {
                (&mut rng)
                    .sample_iter(Uniform::new(0, *qi))
                    .take(ring_size)
                    .collect_vec()
            })
            .collect_vec()
    }

    fn random_vec_in_modulus(&mut self, modulus: u64, size: usize) -> Vec<u64> {
        let rng = thread_rng();
        rng.sample_iter(Uniform::new(0, modulus))
            .take(size)
            .collect_vec()
    }
}

#[cfg(test)]
mod tests {
    use itertools::izip;
    use rand::{CryptoRng, RngCore};

    use super::*;
}
