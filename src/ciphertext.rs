use crate::core_crypto::{
    matrix::{Matrix, MatrixMut, Row, RowMut},
    num::UnsignedInteger,
};

#[derive(PartialEq)]
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
}

pub trait InitialiseLevelledCiphertext {
    type C;
    fn new(c: Self::C, level: usize, representation: Representation) -> Self;
}

pub trait BfvCiphertext: Ciphertext {
    fn degree(&self) -> usize;
    fn c_partq(&self) -> &[Self::Poly];
    fn c_partq_mut(&self) -> &mut [Self::Poly];
    fn level(&self) -> usize;
}

pub trait ExtendedBfvCiphertext: BfvCiphertext {
    fn c_partp(&self) -> &[Self::Poly];
}
