use crate::core_crypto::{num::UnsignedInteger, ring::Matrix};
use itertools::Itertools;
use ndarray::{Array2, ArrayBase, Axis, Dim, IndexLonger, ViewRepr};
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

pub struct NativeVector {
    data: ndarray::Array2<u64>,
    rows: usize,
    cols: usize,
}

struct NativeVectorIterator<'a> {
    data: ArrayBase<ViewRepr<&'a u64>, Dim<[usize; 1]>>,
    index: usize,
    length: usize,
}

impl<'a> Iterator for NativeVectorIterator<'a> {
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

impl<'a> From<ArrayBase<ViewRepr<&'a u64>, Dim<[usize; 1]>>> for NativeVectorIterator<'a> {
    fn from(value: ArrayBase<ViewRepr<&'a u64>, Dim<[usize; 1]>>) -> Self {
        Self {
            data: value,
            index: 0,
            length: value.len(),
        }
    }
}

impl<'a> Matrix<'a, u64, NativeVectorIterator<'a>> for NativeVector {
    fn new(rows: usize, cols: usize) -> Self {
        NativeVector {
            data: ndarray::Array2::<u64>::zeros((rows, cols)),
            rows,
            cols,
        }
    }

    fn from_values(rows: usize, cols: usize, values: Vec<u64>) -> Self {
        let data = ndarray::Array2::<u64>::from_shape_vec((rows, cols), values).unwrap();
        NativeVector { data, rows, cols }
    }

    fn iter_cols(&self) -> NativeVectorIterator<'a> {
        todo!()
    }
    fn iter_cols_mut(&self) -> NativeVectorIterator<'a> {
        todo!()
    }

    fn iter_rows(&self) -> NativeVectorIterator<'a> {
        todo!()
    }
    fn iter_rows_mut(&self) -> NativeVectorIterator<'a> {
        todo!()
    }

    fn get_col(&self, index: usize) -> NativeVectorIterator<'_> {
        self.data.column(index).into()
    }
    fn get_row(&self, index: usize) -> NativeVectorIterator<'_> {
        self.data.row(index).into()
    }

    fn get_index(&self, x: usize, y: usize) -> &u64 {
        self.data.get((x, y)).unwrap()
    }
    fn get_index_mut(&mut self, x: usize, y: usize) -> &mut u64 {
        self.data.get_mut((x, y)).unwrap()
    }

    fn set_index(&mut self, x: usize, y: usize, value: u64) {
        *self.data.get_mut((x, y)).unwrap() = value;
    }

    fn dimension(&self) -> (usize, usize) {
        (self.rows, self.cols)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn native_vector_works() {
        let rows = 1;
        let cols = 100;
        let a = thread_rng()
            .sample_iter(Uniform::new(0, u64::MAX))
            .take(rows * cols)
            .collect_vec();

        let mat = NativeVector::from_values(rows, cols, a.clone());

        // Check iteration over a given column
        for i in 0..cols {
            mat.get_col(1).for_each(|v| {
                assert!(a[i] == *v);
            });
        }

        // Check iteration over the 1st row
        for i in 0..rows {
            
        }
    }
}
