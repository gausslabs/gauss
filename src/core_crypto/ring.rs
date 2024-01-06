use std::usize;

use super::modulus::{BarrettBackend, MontgomeryBackend, MontgomeryScalar};
use itertools::{izip, Itertools};

pub trait MatrixRef<'a, T: 'a>: Matrix<T> {
    type IteratorRef: Iterator<Item = &'a T>;

    fn get_col(&'a self, index: usize) -> Self::IteratorRef;
    fn get_row(&'a self, index: usize) -> Self::IteratorRef;

    fn get(&'a self, row: usize, col: usize) -> &T;
}

pub trait MatrixMut<'a, T: 'a> {
    type IteratorMutRef: Iterator<Item = &'a mut T>;

    fn get_col_mut(&'a mut self, index: usize) -> Self::IteratorMutRef;
    fn get_row_mut(&'a mut self, index: usize) -> Self::IteratorMutRef;

    fn get_mut(&'a mut self, row: usize, col: usize) -> &T;

    fn set(&mut self, row: usize, col: usize, value: T);
}

pub trait Matrix<T> {
    fn zeros(rows: usize, cols: usize) -> Self;
    fn from_values(rows: usize, cols: usize, values: Vec<T>) -> Self;

    fn dimension(&self) -> (usize, usize);
}

/// Given input $x \in R_{QP}$ calculates and returns $[\lceil \frac{t}{P} \cdot x \rfloor]_{Q}$
///
/// We implement "Complex Scaling in CRT representation" of section 2.4 in [HPS](https://eprint.iacr.org/2018/117.pdf)
///
/// TODO (Jay): describe params
pub fn scale_and_round<
    'a,
    MRef: MatrixRef<'a, u64>,
    MMut: MatrixMut<'a, u64>,
    ModOps: MontgomeryBackend<u64, u128> + BarrettBackend<u64, u128>,
>(
    q_out: &mut MMut,
    q_in: &'a MRef,
    p_in: &'a MRef,
    modq_operators: &[ModOps],
    modp_operators: &[ModOps],
    qp_over_pj_inv_modpj_times_tq_modqi_rational: &[Vec<MontgomeryScalar<u64>>],
    qp_over_pj_inv_modpj_times_tq_fractional: &[f64],
    qp_over_qi_inv_modqi_times_tq_over_qi_modqi: &[MontgomeryScalar<u64>],
    q_size: usize,
    p_size: usize,
    qp_size: usize,
    ring_size: usize,
) {
    for n in 0..ring_size {
        // summation for fractional can be done without modular reduction per `qi`
        let mut sum_fractional = 0f64;
        p_in.get_col(n).enumerate().for_each(|(j, px_j)| {
            // px_j * \theta_j
            // TODO (Jay): This will likely result in low precision. A better will be to split fractional into fractional high and fractional low
            sum_fractional += (*px_j as f64) * qp_over_pj_inv_modpj_times_tq_fractional[j];
        });
        let sum_fractional = sum_fractional as u64;

        for i in 0..q_size {
            let modqi_op = modq_operators.get(i).unwrap();

            // convert px_j to montgomery space of \mod{q_i}
            //
            // One thing to note here is \log{pxj} can be much greater than \log{q_i}. For example,
            // when p_j is 60 bits and q_i is 30 bits. Mapping to montgomery space in $q_j$ should still work
            // because the function accepts input in range [0, r) where r = 2^{64}.
            let pxjs_in_mont = p_in
                .get_col(n)
                .map(|x| modqi_op.normal_to_mont_space(*x))
                .collect_vec();

            // \sum px_j * \omega_j
            let mut sum_rational = modqi_op.mont_fma(
                &pxjs_in_mont,
                &qp_over_pj_inv_modpj_times_tq_modqi_rational[i],
            );

            // qx_i * [(qp/q_i)^{-1}]_{q_i} * (t * q)/q_i \mod qi
            let mut qx_i = modqi_op.normal_to_mont_space(*q_in.get(i, n));
            qx_i = modqi_op.mont_mul(qx_i, qp_over_qi_inv_modqi_times_tq_over_qi_modqi[i]);
            sum_rational = modqi_op.mont_add(sum_rational, qx_i);

            let mut sum_rational = modqi_op.mont_to_normal(sum_rational);

            // add fractional part
            sum_rational =
                modqi_op.add_mod_fast(sum_rational, modqi_op.barrett_reduce(sum_fractional));

            q_out.set(i, n, sum_rational);
        }
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;
    use num_traits::{FromPrimitive, ToPrimitive};

    use crate::{
        core_crypto::{
            modulus::{ModulusBackendConfig, NativeModulusBackend},
            prime::generate_primes_vec,
            random::RandomUniformDist,
        },
        utils::{
            convert::TryConvertFromParts,
            mod_inverse, moduli_chain_to_biguint,
            test_utils::{TestMatrix, TestRng},
        },
    };

    use super::*;

    #[test]
    fn scale_and_round_works() {
        let n = 1 << 4;
        let t = 65537u64;
        let q_chain = generate_primes_vec(&[50; 10], n, &[]);
        let p_chain = generate_primes_vec(&[50; 10], n, &q_chain);

        // Q and P mod operators
        let modq_operators = q_chain
            .iter()
            .map(|qi| NativeModulusBackend::initialise(*qi))
            .collect_vec();
        let modp_operators = p_chain
            .iter()
            .map(|pi| NativeModulusBackend::initialise(*pi))
            .collect_vec();

        let big_q = moduli_chain_to_biguint(&q_chain);
        let big_p = moduli_chain_to_biguint(&p_chain);
        let big_qp = &big_p * &big_q;

        // v = [[qp / p_j]^{-1}_{p_j} * tq]
        // \omega_j = floor(v/p_j)
        // \theta_j = (v%p_j)/p_j
        let v = p_chain
            .iter()
            .map(|p_j| {
                &big_q
                    * t
                    * BigUint::from_u64(mod_inverse(
                        ((&big_qp / *p_j) % p_j).to_u64().unwrap(),
                        *p_j,
                    ))
                    .unwrap()
            })
            .collect_vec();
        let omega = izip!(v.iter(), p_chain.iter())
            .map(|(a, p_j)| a / p_j)
            .collect_vec();
        let theta = izip!(v.iter(), p_chain.iter())
            .map(|(a, p_j)| (a % p_j).to_f64().unwrap() / (*p_j) as f64)
            .collect_vec();

        // [\omega_j]_q_i in montgomery space
        let omega_qis = modq_operators
            .iter()
            .map(|modqi_op| {
                omega
                    .iter()
                    .map(|o| {
                        let in_normal_space = (o % modqi_op.modulus).to_u64().unwrap();
                        modqi_op.normal_to_mont_space(in_normal_space)
                    })
                    .collect_vec()
            })
            .collect_vec();

        // [(qp / q_i)^{-1}]_q_i * (t*q)/q_i
        let q_values = modq_operators
            .iter()
            .map(|modqi| {
                let qi = modqi.modulus;

                let out_in_normal_space = ((((&big_q * t) / qi)
                    * BigUint::from(mod_inverse(((&big_qp / qi) % qi).to_u64().unwrap(), qi)))
                    % qi)
                    .to_u64()
                    .unwrap();

                modqi.normal_to_mont_space(out_in_normal_space)
            })
            .collect_vec();

        let test_rng = TestRng {};

        // Random polynomial in QP
        let poly0_q_part = test_rng.random_ring_poly(&q_chain, n);
        let poly0_p_part = test_rng.random_ring_poly(&p_chain, n);

        let mut poly_out = TestMatrix::zeros(q_chain.len(), n);

        scale_and_round(
            &mut poly_out,
            &poly0_q_part,
            &poly0_p_part,
            &modq_operators,
            &modp_operators,
            &omega_qis,
            &theta,
            &q_values,
            q_chain.len(),
            p_chain.len(),
            q_chain.len() + p_chain.len(),
            n,
        );

        let poly_in_big = Vec::<BigUint>::try_convert_with_two_parts(
            &poly0_q_part,
            &poly0_p_part,
            &q_chain,
            &p_chain,
        );
        let poly_out_big = Vec::<BigUint>::try_convert_with_one_part(&poly_out, &q_chain);

        // (t/P (x)) \mod{q} where x \in [0, qp)
        let poly_out_big_expected = poly_in_big.iter().map(|x| {
            if x >= &(&big_qp >> 1usize) {
                &big_q - (((((&big_qp - x) * t) + (&big_p >> 1)) / &big_p) % &big_q)
            } else {
                (((x * t) + (&big_p >> 1)) / &big_p) % &big_q
            }
        });

        izip!(poly_out_big.iter(), poly_out_big_expected).for_each(|(p0, p1)| {
            let bits = if p0 < &p1 {
                (p1 - p0).bits()
            } else {
                (p0 - p1).bits()
            };
            assert!(bits <= 1);
            // dbg!(bits);
        });
    }
}
