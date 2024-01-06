use std::{marker::PhantomData, os::unix::thread};

use crate::core_crypto::{
    num::UnsignedInteger,
    random::RandomUniformDist,
    ring::{Matrix, MatrixMut, MatrixRef},
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

pub struct TestMatrix {
    data: ndarray::Array2<u64>,
    rows: usize,
    cols: usize,
}

pub struct TestMatrixIteratorRef<'a> {
    data: ArrayBase<ViewRepr<&'a u64>, Dim<[usize; 1]>>,
    index: usize,
    length: usize,
}

impl<'a> Iterator for TestMatrixIteratorRef<'a> {
    type Item = &'a u64;
    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.length {
            return None;
        }

        let prev = self.index;
        self.index = prev + 1;
        Some(self.data.index(prev))
    }
}

impl<'a> From<ArrayBase<ViewRepr<&'a u64>, Dim<[usize; 1]>>> for TestMatrixIteratorRef<'a> {
    fn from(value: ArrayBase<ViewRepr<&'a u64>, Dim<[usize; 1]>>) -> Self {
        Self {
            data: value,
            index: 0,
            length: value.len(),
        }
    }
}

pub struct TestMatrixIteratorMutRef<'a> {
    data: ArrayBase<ViewRepr<&'a mut u64>, Dim<[usize; 1]>>,
    index: usize,
    length: usize,
}

impl<'a> Iterator for TestMatrixIteratorMutRef<'a> {
    type Item = &'a mut u64;
    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.length {
            return None;
        }

        let prev = self.index;
        self.index = prev + 1;
        let value = unsafe { &mut *self.data.get_mut_ptr(0).unwrap() };
        Some(value)
    }
}

impl<'a> From<ArrayBase<ViewRepr<&'a mut u64>, Dim<[usize; 1]>>> for TestMatrixIteratorMutRef<'a> {
    fn from(value: ArrayBase<ViewRepr<&'a mut u64>, Dim<[usize; 1]>>) -> Self {
        Self {
            length: value.len(),
            data: value,
            index: 0,
        }
    }
}

impl Matrix<u64> for TestMatrix {
    fn zeros(rows: usize, cols: usize) -> Self {
        TestMatrix {
            data: ndarray::Array2::<u64>::zeros((rows, cols)),
            rows,
            cols,
        }
    }

    fn from_values(rows: usize, cols: usize, values: Vec<u64>) -> Self {
        let data = ndarray::Array2::<u64>::from_shape_vec((rows, cols), values).unwrap();
        TestMatrix { data, rows, cols }
    }

    fn dimension(&self) -> (usize, usize) {
        (self.rows, self.cols)
    }
}

impl<'a> MatrixRef<'a, u64> for TestMatrix {
    type IteratorRef = TestMatrixIteratorRef<'a>;

    fn get_col(&'a self, index: usize) -> Self::IteratorRef {
        self.data.column(index).into()
    }

    fn get_row(&'a self, index: usize) -> Self::IteratorRef {
        self.data.row(index).into()
    }

    fn get(&'a self, row: usize, col: usize) -> &'a u64 {
        self.data.get((row, col)).unwrap()
    }
}

impl<'a> MatrixMut<'a, u64> for TestMatrix {
    type IteratorMutRef = TestMatrixIteratorMutRef<'a>;

    fn get_col_mut(&'a mut self, index: usize) -> Self::IteratorMutRef {
        self.data.column_mut(index).into()
    }

    fn get_row_mut(&'a mut self, index: usize) -> Self::IteratorMutRef {
        self.data.column_mut(index).into()
    }

    fn get_mut(&'a mut self, row: usize, col: usize) -> &'a u64 {
        self.data.get_mut((row, col)).unwrap()
    }

    fn set(&mut self, row: usize, col: usize, value: u64) {
        *self.data.get_mut((row, col)).unwrap() = value;
    }
}

pub struct TestRng {}

impl RandomUniformDist<u64, TestMatrix> for TestRng {
    fn random_ring_poly(&self, moduli_chain: &[u64], ring_size: usize) -> TestMatrix {
        let mut rng = thread_rng();
        let els = moduli_chain
            .iter()
            .flat_map(|qi| {
                (&mut rng)
                    .sample_iter(Uniform::new(0, *qi))
                    .take(ring_size)
                    .collect_vec()
            })
            .collect_vec();

        TestMatrix::from_values(moduli_chain.len(), ring_size, els)
    }

    fn random_vec_in_modulus(&self, modulus: u64, size: usize) -> Vec<u64> {
        let rng = thread_rng();
        rng.sample_iter(Uniform::new(0, modulus))
            .take(size)
            .collect_vec()
    }
}

#[cfg(test)]
mod tests {
    use itertools::izip;
    use rand::{CryptoRng, RngCore};

    use super::*;

    #[test]
    fn native_vector_works() {
        let rows = 1;
        let cols = 100;
        let a = thread_rng()
            .sample_iter(Uniform::new(0, u64::MAX))
            .take(rows * cols)
            .collect_vec();

        let mat = TestMatrix::from_values(rows, cols, a.clone());

        // Check iteration over all columns
        for i in 0..cols {
            mat.get_col(i).for_each(|v| {
                assert!(a[i] == *v);
            });
        }

        // Check iteration over the 1st row
        for i in 0..rows {
            izip!(a.iter(), mat.get_row(i)).for_each(|(a0, v0)| {
                assert_eq!(a0, v0);
            });
        }
    }
}
