use super::{
    modulus::{BarrettBackend, MontgomeryBackend, MontgomeryScalar},
    num::UnsignedInteger,
};
use itertools::{izip, Itertools};

pub trait Matrix<'a, T: UnsignedInteger> {
    fn iter_cols(&self) -> std::slice::Iter<&[T]>;
    fn iter_cols_mut(&self) -> std::slice::Iter<&mut [T]>;

    fn iter_rows(&self) -> std::slice::Iter<&[T]>;
    fn iter_rows_mut(&self) -> std::slice::Iter<&mut [T]>;

    fn get_col(&self, index: usize) -> &'a [T];
    fn get_row(&self, index: usize) -> &'a [T];

    fn get_index(&self, x: usize, y: usize) -> &'a T;
    fn get_index_mut(&mut self, x: usize, y: usize) -> &'a mut T;

    fn set_index(&mut self, x: usize, y: usize, value: T);

    fn dimension() -> (usize, usize);
}

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
    q_in: &M,
    p_in: &M,
    q_mod_operators: &[ModOps],
    p_mod_operators: &[ModOps],
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
        p_in.get_col(n).iter().enumerate().for_each(|(j, px_j)| {
            // px_j * \theta_j
            // TODO (Jay): This will likely result in low precision. A better will be to split fractional into fractional high and fractional low
            sum_fractional += ((*px_j as f64) * qp_over_pj_inv_modpj_mul_q_fractional[j]);
        });
        let sum_fractional = sum_fractional as u64;

        let xjs_in_mont = izip!(p_in.get_col(n).iter(), p_mod_operators.iter())
            .map(|(x, mod_op)| mod_op.normal_to_mont_space(*x))
            .collect_vec();

        for i in 0..q_size {
            let qi_mod_op = q_mod_operators.get(i).unwrap();

            // \sum px_j * \omega_j
            let mut sum_rational =
                qi_mod_op.mont_fma(&xjs_in_mont, &qp_over_pj_inv_modpj_mul_q_modqi_rational[i]);

            // qx_i * [(qp/q_i)^{-1}]_{q_i} * qp/q_i \mod qi
            let mut qx_i = qi_mod_op.normal_to_mont_space(*q_in.get_index(n, i));
            qx_i = qi_mod_op.mont_mul(qx_i, qp_over_qi_inv_modqi_mul_qp_over_qi_modqi[i]);
            sum_rational = qi_mod_op.mont_add(sum_rational, qx_i);

            let mut sum_rational = qi_mod_op.mont_to_normal(sum_rational);

            // add fractional part
            sum_rational =
                qi_mod_op.add_mod_fast(sum_rational, qi_mod_op.barrett_reduce(sum_fractional));

            q_out.set_index(n, i, sum_rational);
        }
    }
}
