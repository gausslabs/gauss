use crate::core_crypto::num::UnsignedInteger;
use itertools::Itertools;
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
