use std::{marker::PhantomData, os::unix::thread};

use crate::{
    core_crypto::{
        matrix::{Matrix, MatrixMut},
        num::UnsignedInteger,
        random::RandomUniformDist,
    },
    keys::SecretKey,
    schemes::ops::generate_ternery_secret_with_hamming_weight,
};
use itertools::Itertools;
use ndarray::{iter::Iter, Array2, ArrayBase, Axis, Dim, IndexLonger, ViewRepr};
use rand::{distributions::Uniform, thread_rng, CryptoRng, Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

pub fn random_vec_in_fq<T: UnsignedInteger + rand::distributions::uniform::SampleUniform>(
    size: usize,
    q: T,
) -> Vec<T> {
    let rng = ChaCha8Rng::from_seed([0u8; 32]);
    rng.sample_iter(Uniform::new(T::zero(), q))
        .take(size)
        .collect_vec()
}

pub struct TestTernarySecret {
    values: Vec<i32>,
}

impl TestTernarySecret {
    pub fn new<
        R: RandomUniformDist<[u8], Parameters = u8>
            + RandomUniformDist<usize, Parameters = usize>
            + CryptoRng,
    >(
        rng: &mut R,
        ring_size: usize,
    ) -> Self {
        let values = generate_ternery_secret_with_hamming_weight(rng, ring_size >> 1, ring_size);
        TestTernarySecret { values }
    }
}

impl SecretKey for TestTernarySecret {
    type Scalar = i32;
    fn values(&self) -> &[Self::Scalar] {
        self.values.as_slice()
    }
}
