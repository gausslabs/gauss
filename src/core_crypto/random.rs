use std::{borrow::Borrow, cell::RefCell};

use num_traits::Zero;
use rand::{
    distributions::uniform::SampleUniform, thread_rng, CryptoRng, Rng, RngCore, SeedableRng,
};
use rand_chacha::{ChaCha8Core, ChaCha8Rng};

use super::{
    matrix::{Matrix, MatrixMut, Row, RowMut},
    num::UnsignedInteger,
};

pub trait RandomSecretByteGenerator {
    fn random_bytes(&mut self, size: usize) -> Vec<u8>;
}

pub trait RandomSecretValueGenerator<T> {
    fn random_value_in_range(&mut self, range: T) -> T;
}

pub trait RandomGaussianDist<S, P> {
    fn random_vec_in_modulus(&self, modulus: S, size: usize) -> Vec<S>;
    fn random_ring_poly(&self, moduli_chain: &[S], ring_size: usize) -> P;
}

pub trait RandomUniformDist<S, P> {
    fn random_vec_in_modulus(&self, modulus: S, size: usize) -> Vec<S>;
    fn random_ring_poly(&self, moduli_chain: &[S], ring_size: usize) -> P;
}

pub trait WithLocal {
    fn with_local<F, R>(func: F) -> R
    where
        F: Fn(&Self) -> R;

    fn with_local_mut<F, R>(func: F) -> R
    where
        F: Fn(&mut Self) -> R;
}

pub(crate) struct DefaultU64SeededRandomGenerator {
    rng: ChaCha8Rng,
}

thread_local! {
    pub(crate) static DEFAULT_U64_SEEDED_RNG: RefCell<DefaultU64SeededRandomGenerator> = RefCell::new(DefaultU64SeededRandomGenerator::new());
}

impl WithLocal for DefaultU64SeededRandomGenerator {
    fn with_local<F, R>(func: F) -> R
    where
        F: Fn(&Self) -> R,
    {
        DEFAULT_U64_SEEDED_RNG.with_borrow(|value| func(value))
    }

    fn with_local_mut<F, R>(func: F) -> R
    where
        F: Fn(&mut Self) -> R,
    {
        DEFAULT_U64_SEEDED_RNG.with_borrow_mut(|value| func(value))
    }
}

impl DefaultU64SeededRandomGenerator {
    pub fn new() -> Self {
        DefaultU64SeededRandomGenerator {
            rng: ChaCha8Rng::from_rng(thread_rng()).unwrap(),
        }
    }

    /// Seed the rng of thread local generator with new seed
    pub fn new_with_seed(seed: <ChaCha8Rng as SeedableRng>::Seed) {
        DefaultU64SeededRandomGenerator::with_local_mut(|generator| {
            let rng = ChaCha8Rng::from_seed(seed);
            generator.rng = rng;
        });
    }
}

impl RandomSecretByteGenerator for DefaultU64SeededRandomGenerator {
    fn random_bytes(&mut self, size: usize) -> Vec<u8> {
        let mut bytes = vec![0u8; size];
        self.rng.fill_bytes(&mut bytes);
        bytes
    }
}

impl<T: SampleUniform + PartialOrd + Zero> RandomSecretValueGenerator<T>
    for DefaultU64SeededRandomGenerator
{
    fn random_value_in_range(&mut self, range: T) -> T {
        self.rng.gen_range(T::zero()..range)
    }
}

impl<P: Matrix> RandomGaussianDist<u64, P> for DefaultU64SeededRandomGenerator {
    fn random_vec_in_modulus(&self, modulus: u64, size: usize) -> Vec<u64> {
        vec![0u64; size]
    }

    fn random_ring_poly(&self, moduli_chain: &[u64], ring_size: usize) -> P {
        P::zeros(moduli_chain.len(), ring_size)
    }
}
impl<P: Matrix> RandomUniformDist<u64, P> for DefaultU64SeededRandomGenerator {
    fn random_ring_poly(&self, moduli_chain: &[u64], ring_size: usize) -> P {
        P::zeros(moduli_chain.len(), ring_size)
    }
    fn random_vec_in_modulus(&self, modulus: u64, size: usize) -> Vec<u64> {
        vec![0u64; size]
    }
}
impl CryptoRng for DefaultU64SeededRandomGenerator {}
