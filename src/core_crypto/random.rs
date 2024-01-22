use super::{
    matrix::{Matrix, MatrixMut, Row, RowMut},
    num::UnsignedInteger,
};

pub trait RandomKeyDist<T> {}

pub trait RandomGaussianDist {
    type Scalar: UnsignedInteger;
    type Row: Row + RowMut;
    type Poly: Matrix<MatElement = Self::Scalar, R = Self::Row>;

    fn random_vec_in_modulus(&self, modulus: Self::Scalar, size: usize) -> Vec<Self::Scalar>;
    fn random_ring_poly(&self, moduli_chain: &[Self::Scalar], ring_size: usize) -> Self::Poly;
}

pub trait RandomUniformDist {
    type Scalar: UnsignedInteger;
    type Row: Row + RowMut;
    type Poly: Matrix<MatElement = Self::Scalar> + MatrixMut<R = Self::Row>;

    fn random_vec_in_modulus(&self, modulus: Self::Scalar, size: usize) -> Vec<Self::Scalar>;
    fn random_ring_poly(&self, moduli_chain: &[Self::Scalar], ring_size: usize) -> Self::Poly;
}
