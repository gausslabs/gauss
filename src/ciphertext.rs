use std::{cell::RefCell, sync::Arc};

use crate::core_crypto::{
    matrix::{Matrix, MatrixMut, Row, RowMut},
    modulus::NativeModulusBackend,
    ntt::NativeNTTBackend,
    num::UnsignedInteger,
};

#[derive(PartialEq, Clone, Copy)]
pub enum Representation {
    Coefficient,
    Evaluation,
}

pub trait Ciphertext {
    type Scalar: UnsignedInteger;
    type Row: Row<Element = Self::Scalar> + RowMut<Element = Self::Scalar>;
    type Poly: Matrix<MatElement = Self::Scalar>
        + MatrixMut<MatElement = Self::Scalar, R = Self::Row>;

    fn representation(&self) -> Representation;
    fn representation_mut(&mut self) -> &mut Representation;
}

pub trait SeededCiphertext {
    type Seed;
    fn seed(&self) -> Self::Seed;
    fn seed_mut(&mut self) -> &mut Self::Seed;
}

pub trait InitSeededLevelledCiphertext {
    type C;
    type Seed;
    fn new(c: Self::C, seed: Self::Seed, level: usize, representation: Representation) -> Self;
}

pub trait InitLevelledCiphertext {
    type C;
    fn new(c: Self::C, level: usize, representation: Representation) -> Self;
}

pub trait RlweCiphertext: Ciphertext {
    fn c_partq(&self) -> &[Self::Poly];
    fn c_partq_mut(&mut self) -> &mut [Self::Poly];
    fn level(&self) -> usize;
    fn level_mut(&mut self) -> &mut usize;
    fn is_lazy(&self) -> bool;
    fn is_lazy_mut(&mut self) -> &mut bool;
}

pub trait ExtendedRlweCiphertext: RlweCiphertext {
    fn c_partp(&self) -> &[Self::Poly];
}

pub struct BfvCiphertextU64 {
    c_partq: Vec<u64>,
    level: usize,
}
