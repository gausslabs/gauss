use std::{iter::Map, ops::Neg, path::Iter};

use aligned_vec::{AVec, CACHELINE_ALIGN};
use itertools::{izip, Itertools};
use num_bigint::{BigInt, BigUint, ToBigInt};
use num_traits::{ToPrimitive, Zero};

use crate::{
    core_crypto::{
        matrix::{Matrix, MatrixMut, RowMut},
        num::big_float::BigFloat,
    },
    parameters::{self, Parameters},
};

use super::{mod_inverse, moduli_chain_to_biguint};

pub trait TryConvertFrom<T: ?Sized> {
    type Parameters: ?Sized;

    fn try_convert_from(value: &T, parameters: &Self::Parameters) -> Self;
}

pub trait TryConvertFromIter<'a, T>
where
    T: 'a,
{
    type Parameters: ?Sized;

    fn try_convert_from(iter: impl Iterator<Item = &'a T>, parameters: &Self::Parameters) -> Self;
}

pub trait TryConvertFromParts<T>: TryConvertFrom<T> {
    fn try_convert_with_two_parts(
        value_part0: &T,
        value_part1: &T,
        parameters_part0: &Self::Parameters,
        parameters_part1: &Self::Parameters,
    ) -> Self;
}

impl TryConvertFrom<[i32]> for Vec<Vec<u64>> {
    type Parameters = [u64];
    fn try_convert_from(value: &[i32], parameters: &Self::Parameters) -> Self {
        parameters
            .iter()
            .map(|qi| {
                value
                    .iter()
                    .map(|signed_x| {
                        let v = signed_x.abs() as u64;
                        let v = v % qi;
                        if signed_x.is_negative() {
                            qi - v
                        } else {
                            v
                        }
                    })
                    .collect_vec()
            })
            .collect_vec()
    }
}

impl TryConvertFrom<[i32]> for AVec<AVec<u64>> {
    type Parameters = [u64];
    fn try_convert_from(value: &[i32], parameters: &Self::Parameters) -> Self {
        let res_iter = parameters.iter().map(|qi| {
            let row_iter = value.iter().map(|signed_x| {
                let v = signed_x.abs() as u64;
                let v = v % qi;
                if signed_x.is_negative() {
                    qi - v
                } else {
                    v
                }
            });
            AVec::from_iter(CACHELINE_ALIGN, row_iter)
        });
        AVec::from_iter(CACHELINE_ALIGN, res_iter)
    }
}

impl<'a> TryConvertFromIter<'a, f64> for Vec<Vec<u64>> {
    type Parameters = [u64];
    // TODO we have to transpose the matrix, becasue we can only use the iter once.
    fn try_convert_from(
        iter: impl Iterator<Item = &'a f64>,
        parameters: &Self::Parameters,
    ) -> Self {
        iter.map(|signed_x| {
            parameters
                .iter()
                .map(|qi| {
                    let v = signed_x.abs() as u64;
                    let v = v % qi;
                    if signed_x.is_negative() {
                        qi - v
                    } else {
                        v
                    }
                })
                .collect_vec()
        })
        .collect_vec()
    }
}

impl TryConvertFrom<[BigUint]> for Vec<Vec<u64>> {
    type Parameters = [u64];
    fn try_convert_from(value: &[BigUint], parameters: &Self::Parameters) -> Self {
        parameters
            .iter()
            .map(|qi| {
                value
                    .iter()
                    .map(|big_x| u64::try_from(big_x % *qi).unwrap())
                    .collect_vec()
            })
            .collect_vec()
    }
}

impl TryConvertFrom<[BigUint]> for AVec<AVec<u64>> {
    type Parameters = [u64];
    fn try_convert_from(value: &[BigUint], parameters: &Self::Parameters) -> Self {
        let res_iter = parameters.iter().map(|qi| {
            let row_iter = value
                .iter()
                .map(|big_x| u64::try_from(big_x % *qi).unwrap());
            AVec::from_iter(CACHELINE_ALIGN, row_iter)
        });
        AVec::from_iter(CACHELINE_ALIGN, res_iter)
    }
}

impl<M> TryConvertFrom<M> for Vec<BigUint>
where
    M: Matrix<MatElement = u64>,
{
    type Parameters = [u64];

    fn try_convert_from(value: &M, parameters: &Self::Parameters) -> Self {
        let big_q = moduli_chain_to_biguint(parameters);

        // q/q_i
        let mut q_over_qi_vec = vec![];
        // [[q/q_i]^{-1}]_q_i
        let mut q_over_qi_inv_modqi_vec = vec![];
        parameters.iter().for_each(|qi| {
            let q_over_qi = &big_q / qi;
            let q_over_qi_inv_modqi =
                BigUint::from(mod_inverse((&q_over_qi % qi).to_u64().unwrap(), *qi));
            q_over_qi_vec.push(q_over_qi);
            q_over_qi_inv_modqi_vec.push(q_over_qi_inv_modqi);
        });

        let (_, ring_size) = value.dimension();

        let mut out_biguint_coeffs = vec![];
        for ri in 0..ring_size {
            let mut x = BigUint::zero();
            value.get_col_iter(ri).enumerate().for_each(|(i, xi)| {
                x += xi * &q_over_qi_vec[i] * &q_over_qi_inv_modqi_vec[i];
            });
            out_biguint_coeffs.push(x % &big_q);
        }

        out_biguint_coeffs
    }
}

impl<M> TryConvertFrom<M> for Vec<BigInt>
where
    M: Matrix<MatElement = u64>,
{
    type Parameters = [u64];

    fn try_convert_from(value: &M, parameters: &Self::Parameters) -> Self {
        let big_q = moduli_chain_to_biguint(parameters);

        // q/q_i
        let mut q_over_qi_vec = vec![];
        // [[q/q_i]^{-1}]_q_i
        let mut q_over_qi_inv_modqi_vec = vec![];
        parameters.iter().for_each(|qi| {
            let q_over_qi = &big_q / qi;
            let q_over_qi_inv_modqi =
                BigUint::from(mod_inverse((&q_over_qi % qi).to_u64().unwrap(), *qi));
            q_over_qi_vec.push(q_over_qi);
            q_over_qi_inv_modqi_vec.push(q_over_qi_inv_modqi);
        });

        let (_, ring_size) = value.dimension();

        let mut out_bigint_coeffs = vec![];
        for ri in 0..ring_size {
            let mut x = BigUint::zero();
            value.get_col_iter(ri).enumerate().for_each(|(i, xi)| {
                x += xi * &q_over_qi_vec[i] * &q_over_qi_inv_modqi_vec[i];
            });
            x = x % &big_q;

            // convert x from unsigned representation to signed representation
            if x >= &big_q >> 1 {
                out_bigint_coeffs.push((&big_q - x).to_bigint().unwrap() * -1);
            } else {
                out_bigint_coeffs.push(x.to_bigint().unwrap());
            }
        }

        out_bigint_coeffs
    }
}

impl<M> TryConvertFrom<M> for Vec<BigFloat>
where
    M: Matrix<MatElement = u64>,
{
    type Parameters = [u64];

    fn try_convert_from(value: &M, parameters: &Self::Parameters) -> Self {
        let big_q = moduli_chain_to_biguint(parameters);

        // q/q_i
        let mut q_over_qi_vec = vec![];
        // [[q/q_i]^{-1}]_q_i
        let mut q_over_qi_inv_modqi_vec = vec![];
        parameters.iter().for_each(|qi| {
            let q_over_qi = &big_q / qi;
            let q_over_qi_inv_modqi =
                BigUint::from(mod_inverse((&q_over_qi % qi).to_u64().unwrap(), *qi));
            q_over_qi_vec.push(q_over_qi);
            q_over_qi_inv_modqi_vec.push(q_over_qi_inv_modqi);
        });

        let (_, ring_size) = value.dimension();

        let mut out_coeffs = vec![];
        for ri in 0..ring_size {
            let mut x = BigUint::zero();
            value.get_col_iter(ri).enumerate().for_each(|(i, xi)| {
                x += xi * &q_over_qi_vec[i] * &q_over_qi_inv_modqi_vec[i];
            });
            x = x % &big_q;

            // convert x from unsigned representation to signed representation
            if x >= &big_q >> 1 {
                out_coeffs.push(((&big_q - x).to_bigint().unwrap().neg()).into());
            } else {
                out_coeffs.push(x.to_bigint().unwrap().to_f64().unwrap().into());
            }
        }

        out_coeffs
    }
}

impl<M> TryConvertFromParts<M> for Vec<BigUint>
where
    M: Matrix<MatElement = u64>,
{
    fn try_convert_with_two_parts(
        value_part0: &M,
        value_part1: &M,
        parameters_part0: &Self::Parameters,
        parameters_part1: &Self::Parameters,
    ) -> Self {
        let big_q = moduli_chain_to_biguint(parameters_part0);
        let big_p = moduli_chain_to_biguint(parameters_part1);
        let big_qp = &big_p * &big_q;

        // qp/q_i
        let mut qp_over_qi_vec = vec![];
        // [[qp/q_i]^{-1}]_q_i
        let mut qp_over_qi_inv_modqi_vec = vec![];
        parameters_part0.iter().for_each(|qi| {
            let qp_over_qi = &big_qp / qi;
            let qp_over_qi_inv_modqi =
                BigUint::from(mod_inverse((&qp_over_qi % qi).to_u64().unwrap(), *qi));
            qp_over_qi_vec.push(qp_over_qi);
            qp_over_qi_inv_modqi_vec.push(qp_over_qi_inv_modqi);
        });

        // qp/p_j
        let mut qp_over_pj_vec = vec![];
        // [[qp/p_j]^{-1}]_p_j
        let mut qp_over_pj_inv_modpj_vec = vec![];
        parameters_part1.iter().for_each(|pj| {
            let qp_over_pj = &big_qp / pj;
            let qp_over_pj_inv_modpj =
                BigUint::from(mod_inverse((&qp_over_pj % pj).to_u64().unwrap(), *pj));
            qp_over_pj_vec.push(qp_over_pj);
            qp_over_pj_inv_modpj_vec.push(qp_over_pj_inv_modpj);
        });

        let (_, ring_size) = value_part0.dimension();

        let mut out_biguint_coeffs = vec![];
        for ri in 0..ring_size {
            let mut x = BigUint::zero();

            // q part
            value_part0
                .get_col_iter(ri)
                .enumerate()
                .for_each(|(i, xi)| {
                    x += xi * &qp_over_qi_vec[i] * &qp_over_qi_inv_modqi_vec[i];
                });

            // p part
            value_part1
                .get_col_iter(ri)
                .enumerate()
                .for_each(|(i, xi)| {
                    x += xi * &qp_over_pj_vec[i] * &qp_over_pj_inv_modpj_vec[i];
                });

            out_biguint_coeffs.push(x % &big_qp);
        }

        out_biguint_coeffs
    }
}

// TODO(Jay): I wanted to implement `try_convert_from` from [f64] to Matrix.
// Couldn't figure out the right way to do so. After implementing this we can
// `CastToZp<Uint> for F` trait bound (CastToZp trait itself) in `simd_encode`
// and instead use TryConvertFrom<[f64]> for Matrix

// pub trait TryConvertFromMut<T: ?Sized> {
//     type Parameters: ?Sized;

//     fn try_convert_from_mut(value: &T, parameters: &Self::Parameters, out:
// &mut Self); }
// impl<M> TryConvertFromMut<dyn Iterator<Item = f64>> for M
// where
//     M: MatrixMut,
//     <M as Matrix>::R: RowMut,
// {
//     type Parameters = [u64];
//     fn try_convert_from_mut(
//         value: &dyn Iterator<Item = f64>,
//         parameters: &Self::Parameters,
//         out: &mut Self,
//     ) {
//         // TODO(Jay): I don't think calculating mbig_q here is critical but
// not sure.         // Recheck
//         let big_q = moduli_chain_to_biguint(&parameters);
//         let dim = out.dimension();

//         debug_assert!(
//             dim.0 == parameters.len(),
//             "Expected Matrix to have {} rows but has {}",
//             parameters.len(),
//             dim.0
//         );

//         let ring_size = dim.1;

//         izip!((0..ring_size).into_iter(), *value).for_each(|(col_index, v)|
// {});     }
// }

#[cfg(test)]
mod tests {
    use aligned_vec::AVec;
    use itertools::Itertools;
    use num_bigint::{BigUint, RandBigInt};
    use num_traits::One;
    use rand::thread_rng;

    use super::{TryConvertFrom, TryConvertFromParts};

    use crate::core_crypto::prime::generate_primes_vec;

    #[test]
    fn convert_from_and_to_biguint_for_u64_moduli_chain_works() {
        let ring_size = 1 << 4;
        let q_chain_bits = vec![50, 50];
        let q_chain = generate_primes_vec(&q_chain_bits, ring_size, &[]);

        let mut big_q = BigUint::one();
        q_chain.iter().for_each(|x| {
            big_q *= *x;
        });

        let mut rng = thread_rng();
        let poly_big = (0..ring_size)
            .map(|_| rng.gen_biguint_below(&big_q))
            .collect_vec();

        // decompose biguint poly into q_i's
        let poly_decomposed =
            <Vec<Vec<u64>> as TryConvertFrom<[BigUint]>>::try_convert_from(&poly_big, &q_chain);

        // recompose decomposed poly into biguint
        let poly_big_back = Vec::<BigUint>::try_convert_from(&poly_decomposed, &q_chain);

        assert_eq!(poly_big, poly_big_back);
    }

    #[test]
    fn convert_from_and_to_biguint_for_u64_two_parts_moduli_chain_works_for_vec() {
        let ring_size = 1 << 4;
        let q_chain = generate_primes_vec(&[50, 50], ring_size, &[]);
        let p_chain = generate_primes_vec(&[50, 50], ring_size, &q_chain);

        let mut big_q = BigUint::one();
        q_chain.iter().for_each(|x| {
            big_q *= *x;
        });
        let mut big_p = BigUint::one();
        p_chain.iter().for_each(|x| {
            big_p *= *x;
        });
        let big_qp = &big_p * &big_q;

        let mut rng = thread_rng();
        let poly_big = (0..ring_size)
            .map(|_| rng.gen_biguint_below(&big_qp))
            .collect_vec();

        // decompose biguint poly into q_i's
        let poly_decomposed_part0 =
            <Vec<Vec<u64>> as TryConvertFrom<[BigUint]>>::try_convert_from(&poly_big, &q_chain);
        let poly_decomposed_part1 =
            <Vec<Vec<u64>> as TryConvertFrom<[BigUint]>>::try_convert_from(&poly_big, &p_chain);

        // recompose decomposed poly into biguint
        let poly_big_back = Vec::<BigUint>::try_convert_with_two_parts(
            &poly_decomposed_part0,
            &poly_decomposed_part1,
            &q_chain,
            &p_chain,
        );

        // assert_eq!(poly_big, poly_big_back);
    }

    #[test]
    fn convert_from_and_to_biguint_for_u64_two_parts_moduli_chain_works_for_avec() {
        let ring_size = 1 << 4;
        let q_chain = generate_primes_vec(&[50, 50], ring_size, &[]);
        let p_chain = generate_primes_vec(&[50, 50], ring_size, &q_chain);

        let mut big_q = BigUint::one();
        q_chain.iter().for_each(|x| {
            big_q *= *x;
        });
        let mut big_p = BigUint::one();
        p_chain.iter().for_each(|x| {
            big_p *= *x;
        });
        let big_qp = &big_p * &big_q;

        let mut rng = thread_rng();
        let poly_big = (0..ring_size)
            .map(|_| rng.gen_biguint_below(&big_qp))
            .collect_vec();

        // decompose biguint poly into q_i's
        let poly_decomposed_part0 =
            <AVec<AVec<u64>> as TryConvertFrom<[BigUint]>>::try_convert_from(&poly_big, &q_chain);
        let poly_decomposed_part1 =
            <AVec<AVec<u64>> as TryConvertFrom<[BigUint]>>::try_convert_from(&poly_big, &p_chain);

        // recompose decomposed poly into biguint
        let poly_big_back = Vec::<BigUint>::try_convert_with_two_parts(
            &poly_decomposed_part0,
            &poly_decomposed_part1,
            &q_chain,
            &p_chain,
        );

        // assert_eq!(poly_big, poly_big_back);
    }
}
