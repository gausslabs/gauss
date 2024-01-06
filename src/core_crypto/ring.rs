use std::usize;

use super::modulus::{BarrettBackend, MontgomeryBackend, MontgomeryScalar};
use itertools::{izip, Itertools};

pub trait Matrix<'a, T: 'a> {
    type IteratorRef: Iterator<Item = &'a T>;
    type IteratorMutRef: Iterator<Item = &'a mut T>;

    fn new(rows: usize, cols: usize) -> Self;
    fn from_values(rows: usize, cols: usize, values: Vec<T>) -> Self;

    fn get_col(&'a self, index: usize) -> Self::IteratorRef;
    fn get_row(&'a self, index: usize) -> Self::IteratorRef;

    fn get_col_mut(&'a mut self, index: usize) -> Self::IteratorMutRef;
    fn get_row_mut(&'a mut self, index: usize) -> Self::IteratorMutRef;

    fn get(&'a self, row: usize, col: usize) -> &'a T;
    fn get_mut(&'a mut self, row: usize, col: usize) -> &'a T;

    fn set(&mut self, row: usize, col: usize, value: T);

    fn dimension(&self) -> (usize, usize);
}

// pub fn random_ring_polynomial_in_uniform_dist<
//     R: CryptoRng + RngCore,
//     I: Iterator<Item = u64>,
//     ModOps: ModulusRandomVecInDistGenerator<u64, R, I>,
// >(
//     modq_operators: &[ModOps],
//     ring_size: usize,
// ) {
//     // let v =
// }

/// Given input $x \in R_{QP}$ calculates and returns $[\lceil \frac{t}{P} \cdot x \rfloor]_{Q}$
///
/// We implement "Complex Scaling in CRT representation" of section 2.4 in [HPS](https://eprint.iacr.org/2018/117.pdf)
///
/// TODO (Jay): describe params
pub fn scale_and_round<
    'a,
    M: Matrix<'a, u64>,
    ModOps: MontgomeryBackend<u64, u128> + BarrettBackend<u64, u128>,
>(
    q_out: &mut M,
    q_in: &'a M,
    p_in: &'a M,
    modq_operators: &[ModOps],
    modp_operators: &[ModOps],
    qp_over_pj_inv_modpj_mul_q_modqi_rational: Vec<&[MontgomeryScalar<u64>]>,
    qp_over_pj_inv_modpj_mul_q_fractional: &[f64],
    qp_over_qi_inv_modqi_mul_qp_over_qi_modqi: &[MontgomeryScalar<u64>],
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
            sum_fractional += ((*px_j as f64) * qp_over_pj_inv_modpj_mul_q_fractional[j]);
        });
        let sum_fractional = sum_fractional as u64;

        let xjs_in_mont = izip!(p_in.get_col(n), modp_operators.iter())
            .map(|(x, mod_op)| mod_op.normal_to_mont_space(*x))
            .collect_vec();

        for i in 0..q_size {
            let qi_mod_op = modq_operators.get(i).unwrap();

            // \sum px_j * \omega_j
            let mut sum_rational =
                qi_mod_op.mont_fma(&xjs_in_mont, &qp_over_pj_inv_modpj_mul_q_modqi_rational[i]);

            // qx_i * [(qp/q_i)^{-1}]_{q_i} * qp/q_i \mod qi
            let mut qx_i = qi_mod_op.normal_to_mont_space(*q_in.get(n, i));
            qx_i = qi_mod_op.mont_mul(qx_i, qp_over_qi_inv_modqi_mul_qp_over_qi_modqi[i]);
            sum_rational = qi_mod_op.mont_add(sum_rational, qx_i);

            let mut sum_rational = qi_mod_op.mont_to_normal(sum_rational);

            // add fractional part
            sum_rational =
                qi_mod_op.add_mod_fast(sum_rational, qi_mod_op.barrett_reduce(sum_fractional));

            q_out.set(n, i, sum_rational);
        }
    }
}

#[cfg(test)]
mod tests {

    fn dwdwa() {}
}
