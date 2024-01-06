use itertools::{izip, Itertools};
use num_bigint::BigUint;
use num_traits::{ToPrimitive, Zero};

use crate::core_crypto::ring::{Matrix, MatrixRef};

use super::{mod_inverse, moduli_chain_to_biguint};

pub trait TryConvertFrom<T> {
    type Parameters;

    fn try_convert_from(value: T, parameters: Self::Parameters) -> Self;
}

pub trait TryConvertFromParts<T> {
    type Parameters;

    fn try_convert_with_one_part(value: T, parameters: Self::Parameters) -> Self;
    fn try_convert_with_two_parts(
        value_part0: T,
        value_part1: T,
        parameters_part0: Self::Parameters,
        parameters_part1: Self::Parameters,
    ) -> Self;
}

impl<'a, M> TryConvertFrom<&'a [BigUint]> for M
where
    M: Matrix<u64>,
{
    type Parameters = &'a Vec<u64>;
    fn try_convert_from(value: &'a [BigUint], parameters: Self::Parameters) -> Self {
        let values = parameters
            .iter()
            .flat_map(|qi| value.iter().map(|big_x| (big_x % *qi).to_u64().unwrap()))
            .collect_vec();

        Matrix::from_values(parameters.len(), value.len(), values)
    }
}

impl<'a, M> TryConvertFromParts<&'a M> for Vec<BigUint>
where
    M: MatrixRef<'a, u64>,
{
    type Parameters = &'a Vec<u64>;

    fn try_convert_with_one_part(value: &'a M, parameters: Self::Parameters) -> Self {
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
            value.get_col(ri).enumerate().for_each(|(i, xi)| {
                x += xi * &q_over_qi_vec[i] * &q_over_qi_inv_modqi_vec[i];
            });
            out_biguint_coeffs.push(x % &big_q);
        }

        out_biguint_coeffs
    }

    fn try_convert_with_two_parts(
        value_part0: &'a M,
        value_part1: &'a M,
        parameters_part0: Self::Parameters,
        parameters_part1: Self::Parameters,
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
            value_part0.get_col(ri).enumerate().for_each(|(i, xi)| {
                x += xi * &qp_over_qi_vec[i] * &qp_over_qi_inv_modqi_vec[i];
            });

            // p part
            value_part1.get_col(ri).enumerate().for_each(|(i, xi)| {
                x += xi * &qp_over_pj_vec[i] * &qp_over_pj_inv_modpj_vec[i];
            });

            out_biguint_coeffs.push(x % &big_qp);
        }

        out_biguint_coeffs
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;
    use num_bigint::{BigUint, RandBigInt};
    use num_traits::One;
    use rand::{thread_rng, Rng};

    use crate::{
        core_crypto::{
            prime::generate_primes_vec,
            ring::{self, Matrix},
        },
        utils::test_utils::TestMatrix,
    };

    use super::{TryConvertFrom, TryConvertFromParts};

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
        let poly_decomposed = TestMatrix::try_convert_from(&poly_big, &q_chain);

        // recompose decomposed poly into biguint
        let poly_big_back = Vec::<BigUint>::try_convert_with_one_part(&poly_decomposed, &q_chain);

        assert_eq!(poly_big, poly_big_back);
    }

    #[test]
    fn convert_from_and_to_biguint_for_u64_two_parts_moduli_chain_works() {
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
        let poly_decomposed_part0 = TestMatrix::try_convert_from(&poly_big, &q_chain);
        let poly_decomposed_part1 = TestMatrix::try_convert_from(&poly_big, &p_chain);

        // recompose decomposed poly into biguint
        let poly_big_back = Vec::<BigUint>::try_convert_with_two_parts(
            &poly_decomposed_part0,
            &poly_decomposed_part1,
            &q_chain,
            &p_chain,
        );

        assert_eq!(poly_big, poly_big_back);
    }
}
