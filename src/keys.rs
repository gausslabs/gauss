use std::usize;

use crate::{ciphertext::Ciphertext, core_crypto::num::UnsignedInteger};

pub trait SecretKey {
    type Scalar;
    fn values(&self) -> &[Self::Scalar];
}

pub trait SecretWithHammingWeight: SecretKey {
    fn new_with_hamming_weight(hw: usize, ring_size: usize) -> Self;
}

pub trait Encryptor<M, C> {
    fn encrypt(&self, message: M) -> C;
}

pub trait Decryptor<M, C> {
    fn decrypt(&self, c: C) -> M;
}
