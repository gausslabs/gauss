use std::fmt::Debug;

use aligned_vec::{avec, AVec};
use itertools::Itertools;

use crate::core_crypto::num::UnsignedInteger;

pub trait Matrix: AsRef<[Self::R]> + Debug {
    type MatElement;
    type R: Row<Element = Self::MatElement>;

    fn dimension(&self) -> (usize, usize);

    fn zeros(row: usize, col: usize) -> Self;

    fn get_row(&self, row_idx: usize) -> &Self::R {
        &self.as_ref()[row_idx]
    }

    fn get_col_iter(&self, column_idx: usize) -> impl Iterator<Item = &Self::MatElement> {
        self.as_ref().iter().map(move |r| &r.as_ref()[column_idx])
    }

    fn iter_rows(&self) -> impl Iterator<Item = &Self::R> {
        self.as_ref().iter().map(move |r| r)
    }

    fn get_element(&self, row_idx: usize, column_idx: usize) -> &Self::MatElement {
        &self.as_ref()[row_idx].as_ref()[column_idx]
    }
}

pub trait MatrixMut: Matrix + AsMut<[<Self as Matrix>::R]>
where
    <Self as Matrix>::R: RowMut,
{
    fn get_row_mut(&mut self, row_idx: usize) -> &mut [<Self as Matrix>::MatElement] {
        self.as_mut()[row_idx].as_mut()
    }

    #[inline]
    fn set(&mut self, row_idx: usize, column_idx: usize, val: <Self as Matrix>::MatElement) {
        self.as_mut()[row_idx].as_mut()[column_idx] = val;
    }

    fn get_col_iter_mut(
        &mut self,
        column_idx: usize,
    ) -> impl Iterator<Item = &mut Self::MatElement> {
        let res = self
            .as_mut()
            .iter_mut()
            .map(move |r| r.as_mut().get_mut(column_idx).unwrap());
        res
    }

    fn iter_rows_mut(&mut self) -> impl Iterator<Item = &mut Self::R> {
        self.as_mut().iter_mut().map(move |r| r)
    }

    fn get_element_mut(
        &mut self,
        row_idx: usize,
        column_idx: usize,
    ) -> &mut <Self as Matrix>::MatElement {
        &mut self.as_mut()[row_idx].as_mut()[column_idx]
    }
}

pub trait Drop2Dimension: Matrix {
    /// Drop 2 dimensionality of matrix and reduces it self to a single row
    /// vector
    fn drop(self) -> Vec<Self::MatElement>;
}

impl<T: UnsignedInteger> Drop2Dimension for Vec<Vec<T>> {
    fn drop(self) -> Self::R {
        self.into_iter().flatten().collect_vec()
    }
}

///This part is borrowed from zama fhe's containter type without modification
///https://github.com/zama-ai/tfhe-rs/blob/main/tfhe/src/core_crypto/commons/traits/container.rs
pub trait Row: AsRef<[Self::Element]> {
    type Element;
}
pub trait RowMut: Row + AsMut<[<Self as Row>::Element]> {}

impl<T> Row for AVec<T> {
    type Element = T;
}
impl<T> Row for Vec<T> {
    type Element = T;
}

impl<T> RowMut for Vec<T> {}

impl<T> Row for &[T] {
    type Element = T;
}

impl<T> Row for &mut [T] {
    type Element = T;
}

impl<T> RowMut for &mut [T] {}

impl<T> Row for aligned_vec::ABox<[T]> {
    type Element = T;
}

impl<T> RowMut for aligned_vec::ABox<[T]> {}

impl<T> Row for Box<[T]> {
    type Element = T;
}

impl<T> RowMut for Box<[T]> {}
impl<T> RowMut for aligned_vec::AVec<T> {}

pub trait IntoRowOwned: Row + AsMut<[Self::Element]> {
    fn collect<I: Iterator<Item = Self::Element>>(iter: I) -> Self;
}

impl<T> IntoRowOwned for aligned_vec::ABox<[T]> {
    fn collect<I: Iterator<Item = Self::Element>>(iter: I) -> Self {
        aligned_vec::AVec::<T, _>::from_iter(0, iter).into_boxed_slice()
    }
}

impl<T> Matrix for AVec<AVec<T>>
where
    T: UnsignedInteger,
{
    type MatElement = T;
    type R = AVec<T>;

    fn zeros(row: usize, col: usize) -> Self {
        avec![avec![T::zero(); col]; row]
    }

    fn dimension(&self) -> (usize, usize) {
        (self.len(), self[0].len())
    }
}

impl<T> MatrixMut for AVec<AVec<T>> where T: UnsignedInteger {}

impl<T> Matrix for Vec<Vec<T>>
where
    T: UnsignedInteger,
{
    type MatElement = T;
    type R = Vec<T>;

    fn zeros(row: usize, col: usize) -> Self {
        vec![vec![T::zero(); col]; row]
    }

    fn dimension(&self) -> (usize, usize) {
        (self.len(), self[0].len())
    }
}
impl<T> MatrixMut for Vec<Vec<T>> where T: UnsignedInteger {}

#[cfg(test)]
mod test {

    use super::*;
    use crate::core_crypto::matrix::MatrixMut;
    use aligned_vec::avec;

    #[test]
    fn test_matrix_avec() {
        let v1 = avec![1_u64, 2_u64];
        let v2 = avec![3_u64, 4_u64];
        let v3 = avec![v1, v2];
        let r = v3.get_row(0);
        assert_eq!(r, &v3[0]);
        let c: Vec<_> = v3.get_col_iter(0usize).collect();
        for (idx, element_ref) in c.iter().enumerate() {
            assert_eq!(*element_ref, &v3[idx][0]);
        }
        let e = v3.get_element(1, 1);
        assert_eq!(4, *e);
    }
    #[test]
    fn test_matrix_mut_avec() {
        let v1 = avec![1_u64, 2_u64];
        let v2 = avec![3_u64, 4_u64];
        let mut v3 = avec![v1, v2];
        let r = v3.get_row_mut(0);
        r[0] = 2;
        assert_eq!(2, v3[0][0]);
        let c: Vec<_> = v3.get_col_iter_mut(1usize).collect();
        for i in c {
            *i = 0;
        }
        for r in v3.iter() {
            assert_eq!(r[1], 0);
        }

        v3.set(1, 1, 0);
        assert_eq!(0, v3[1][1]);
    }
    #[test]
    fn test_matrix_vec() {
        let v1 = vec![1_u64, 2_u64];
        let v2 = vec![3_u64, 4_u64];
        let v3 = vec![v1, v2];
        let r = v3.get_row(0);
        assert_eq!(r, &v3[0]);
        let c: Vec<_> = v3.get_col_iter(0usize).collect();
        for (idx, element_ref) in c.iter().enumerate() {
            assert_eq!(*element_ref, &v3[idx][0]);
        }
        let e = v3.get_element(1, 1);
        assert_eq!(4, *e);
    }
    #[test]
    fn test_matrix_mut_vec() {
        let v1 = vec![1_u64, 2_u64];
        let v2 = vec![3_u64, 4_u64];
        let mut v3 = vec![v1, v2];
        let r = v3.get_row_mut(0);
        r[0] = 2;
        assert_eq!(2, v3[0][0]);
        let c: Vec<_> = v3.get_col_iter_mut(1usize).collect();
        for i in c {
            *i = 0;
        }
        for r in v3.iter() {
            assert_eq!(r[1], 0);
        }

        v3.set(1, 1, 0);
        assert_eq!(0, v3[1][1]);
    }
}
