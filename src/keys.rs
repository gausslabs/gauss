use std::usize;

use crate::{ciphertext::Ciphertext, core_crypto::num::UnsignedInteger};

pub trait SecretKey {
    type Scalar;
    fn values(&self) -> &[Self::Scalar];
}

pub trait Encryptor<M: ?Sized, C> {
    fn encrypt(&self, message: &M) -> C;
}

pub trait Decryptor<M, C> {
    fn decrypt(&self, c: &C) -> M;
}

/// Encodes self to plaintext P
pub trait LevelEncoder<P> {
    fn encode(&self, level: usize) -> P;
}

/// Decodes self to message M
pub trait LevelDecoder<M> {
    fn decode(&self, level: usize) -> M;
}
