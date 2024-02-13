use std::{borrow::Borrow, cell::RefCell};

use itertools::{izip, Itertools};
use num_traits::Zero;
use rand::{
    distributions::{uniform::SampleUniform, Uniform},
    thread_rng, CryptoRng, Rng, RngCore, SeedableRng,
};
use rand_chacha::{ChaCha8Core, ChaCha8Rng};

use super::matrix::{Matrix, MatrixMut, RowMut};

pub trait RandomSecretByteGenerator {
    fn random_bytes(&mut self, size: usize) -> Vec<u8>;
}

pub trait RandomSecretValueGenerator<T> {
    fn random_value_in_range(&mut self, range: T) -> T;
}

pub trait RandomSeed {
    type Seed;
    fn random_seed(&mut self) -> Self::Seed;
}

pub trait RandomGaussianDist<M> {
    type Parameters: ?Sized;
    fn random_fill(&mut self, parameters: &Self::Parameters, container: &mut M);
}

pub trait RandomUniformDist<M> {
    type Parameters: ?Sized;
    fn random_fill(&mut self, parameters: &Self::Parameters, container: &mut M);
}

pub trait InitWithSeed {
    type Seed;
    fn init_with_seed(seed: Self::Seed) -> Self;
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

    /// Seed the rng of default thread local generator with new seed
    pub fn new_with_seed(seed: <ChaCha8Rng as SeedableRng>::Seed) {
        DefaultU64SeededRandomGenerator::with_local_mut(|generator| {
            let rng = ChaCha8Rng::from_seed(seed);
            generator.rng = rng;
        });
    }
}

impl InitWithSeed for DefaultU64SeededRandomGenerator {
    type Seed = <ChaCha8Rng as SeedableRng>::Seed;
    fn init_with_seed(seed: Self::Seed) -> Self {
        let rng = ChaCha8Rng::from_seed(seed);
        DefaultU64SeededRandomGenerator { rng }
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

impl<M: MatrixMut<MatElement = u64>> RandomUniformDist<M> for DefaultU64SeededRandomGenerator
where
    <M as Matrix>::R: RowMut,
{
    type Parameters = [u64];
    fn random_fill(&mut self, parameters: &Self::Parameters, container: &mut M) {
        debug_assert!(
            parameters.len() == container.dimension().0,
            "Matrix rows do not equal moduli chain: {} != {}",
            container.dimension().0,
            parameters.len()
        );

        izip!(container.iter_rows_mut(), parameters.iter()).for_each(|(r, qi)| {
            izip!(
                r.as_mut().iter_mut(),
                (&mut self.rng).sample_iter(Uniform::new(0, *qi))
            )
            .for_each(|(r_el, random_el)| *r_el = random_el);
        });
    }
}

impl<M: MatrixMut<MatElement = u64>> RandomGaussianDist<M> for DefaultU64SeededRandomGenerator
where
    <M as Matrix>::R: RowMut,
{
    type Parameters = [u64];
    fn random_fill(&mut self, parameters: &Self::Parameters, container: &mut M) {
        debug_assert!(
            parameters.len() == container.dimension().0,
            "Matrix rows do not equal moduli chain: {} != {}",
            container.dimension().0,
            parameters.len()
        );

        // izip!(container.iter_rows_mut(), parameters.iter()).for_each(|(r,
        // qi)| {     izip!(
        //         r.as_mut().iter_mut(),
        //         (&mut self.rng).sample_iter(Uniform::new(0, *qi))
        //     )
        //     .for_each(|(r_el, random_el)| *r_el = random_el);
        // });
    }
}

impl CryptoRng for DefaultU64SeededRandomGenerator {}
