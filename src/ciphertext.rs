use crate::core_crypto::ring::{Matrix, MatrixMut};

#[derive(PartialEq)]
pub enum Representation {
    Coefficient,
    Evaluation,
}

pub trait Ciphertext {
    type Poly: Matrix<u64>;

    fn representation(&self) -> Representation;
}

pub trait BfvCiphertext: Ciphertext {
    fn degree(&self) -> usize;
    fn c_basisq(&self) -> &[Self::Poly];
    fn level(&self) -> usize;
}

pub trait ExtendedBfvCiphertext: BfvCiphertext {
    fn c_basisp(&self) -> &[Self::Poly];
}
