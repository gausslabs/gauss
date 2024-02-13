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
