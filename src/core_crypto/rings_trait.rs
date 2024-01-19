use std::os::unix::raw::mode_t;
use crate::core_crypto::traits::align::{Matrix,MatrixMut};
use super::modulus::{BarrettBackend, MontgomeryBackend, MontgomeryScalar};
use aligned_vec::AVec;
use itertools::{izip, Itertools};

type InnerMat<T> = AVec<AVec<T>>;

/// Given input polnyomial x \in Q outputs $[\lceil \frac{P \cdot x}{Q} \rfloor]_P$
///
/// We implement "Modulus Switching between Arbitrary RNS Bases" presented in
/// Appendix E of [2021/204](https://eprint.iacr.org/2021/204.pdf).
pub fn fast_convert_p_over_q<
    ModOps: MontgomeryBackend<u64, u128> + BarrettBackend<u64, u128>,
>(
    p_out: &mut InnerMat<u64>,
    q_in: &mut InnerMat<u64>,
    neg_p_times_q_over_qi_inv_modqi: &[u64],
    qi_inv_per_modpj: &[Vec<MontgomeryScalar<u64>>],
    one_over_qi: &[f64],
    modq_operators: &[ModOps],
    modp_operators: &[ModOps],
    q_size: usize,
    p_size: usize,
    ring_size: usize,
) {
    for ri in 0..ring_size {
        let mut mu = 0.5f64;
        // qxi_values = for each i: qxi * [-p(q/qi)^{-1}]_q_i \mod{q_i}
        let qxi_values = izip!(
            q_in.get_col_iter(ri),
            neg_p_times_q_over_qi_inv_modqi.iter(),
            one_over_qi.iter(),
            modq_operators.iter()
        )
        .map(|(qxi, op, one_over_qi_v, modqi)| {
            // qxi * [-p(q/qi)^{-1}]_q_i \mod{q_i}
            let value = modqi.mul_mod_fast(*qxi, *op);

            // To estimate \mu: value * \frac{1}{q_i}
            mu += (value as f64 * one_over_qi_v);

            value
        })
        .collect_vec();

        let mu = mu as u64;

        for j in 0..p_size {
            let modpj = &modp_operators[j];
            let qxi_values_mont = qxi_values
                .iter()
                .map(|x| modpj.normal_to_mont_space(*x))
                .collect_vec();

            // \sum qxi_values[i] * q_i^{-1} \mod{p_j}
            let out_in_pj = modpj.mont_fma(&qxi_values_mont, &qi_inv_per_modpj[j]);
            let mut out_in_pj = modpj.mont_to_normal(out_in_pj);

            // subtract oveflow \mu
            out_in_pj = modpj.sub_mod_fast(out_in_pj, mu);

            p_out.set(j, ri, out_in_pj);
        }
    }
}


/// Switches basis of a polynomial in Q basis to P basis
pub fn switch_crt_basis<
    ModOps: MontgomeryBackend<u64, u128> + BarrettBackend<u64, u128>,
>(
    p_out: &mut InnerMat<u64>,
    q_in: &InnerMat<u64>,
    q_over_qi_inv_modqi: &[u64],
    q_over_qi_per_modpj: &[Vec<MontgomeryScalar<u64>>],
    mu_times_q_per_modpj: &[Vec<MontgomeryScalar<u64>>],
    one_over_qi: &[f64],
    modq_operators: &[ModOps],
    modp_operators: &[ModOps],
    q_size: usize,
    p_size: usize,
    ring_size: usize,
) {
    for ri in 0..ring_size {
        let mut mu = 0.5f64;
        // q_values = for each i: qx_i * {q/q_i}^{-1}_q_i \mod{q_i}
        let q_values = izip!(
            q_in.get_col_iter(ri),
            q_over_qi_inv_modqi.iter(),
            modq_operators.iter(),
            one_over_qi.iter()
        )
        .map(|(qx_i, v, modqi, one_over_qi_value)| {
            // qx_i * {q/q_i}^{-1}_q_i \mod{q_i}
            let v = modqi.mul_mod_fast(*qx_i, *v);

            // To calculate \mu:
            // qx_i * {q/q_i}^{-1}_q_i \mod{q_i} * (1/q_i)
            mu += v as f64 * one_over_qi_value;

            return v;
        })
        .collect_vec();

        let mu = mu as usize;

        for j in 0..p_size {
            let modpj = &modp_operators[j];

            // map q_values to p_j mont space
            // [q_values_i]_p_j
            let qxi_times_q_hat_inv_modpj = q_values
                .iter()
                .map(|x| modpj.normal_to_mont_space(*x))
                .collect_vec();

            // [[input]_Q + \mu Q]_p_j =  \sum [q_valules_i]_p_j * q/q_i \mod{p_j}
            let mut out_in_pj = modpj.mont_fma(&qxi_times_q_hat_inv_modpj, &q_over_qi_per_modpj[j]);

            // subtract overflow: \mu Q
            out_in_pj = modpj.mont_sub(out_in_pj, mu_times_q_per_modpj[j][mu]);
            let input_modpj = modpj.mont_to_normal(out_in_pj);

            p_out.set(j, ri, input_modpj);
        }
    }
}

/// Simple scale and round procedure. Given input polynomial x \in R_Q it outputs scaled polynomial
/// $\frac{t}{Q}\cdot x \in R_t$
///
/// We implement simple scaling procedure as outlined in [HPS] along with digit decomposition technique
/// to reduce error accumulation due to precision loss when $max(log{q_i}) > 51$ bits as outlined in https://eprint.iacr.org/2021/204.pdf.
///
/// Usually $q_i < 2^{61}$. Thus it suffices to limit k = 2 and decomposition base \beta to $2^(max(log{q_i})/2)$.
pub fn simple_scale_and_round<

    ModOps: MontgomeryBackend<u64, u128> + BarrettBackend<u64, u128>,
>(
    t_out: &mut InnerMat<u64>,
    q_in: &InnerMat<u64>,
    q_over_qi_inv_modqi_times_t_over_qi_modt: &[MontgomeryScalar<u64>],
    beta_times_q_over_qi_inv_modqi_times_t_over_qi_modt: &[MontgomeryScalar<u64>],
    q_over_qi_inv_modqi_times_t_over_qi_fractional: &[f64],
    beta_times_q_over_qi_inv_modqi_times_t_over_qi_fractional: &[f64],
    log_beta: usize,
    modt_operator: ModOps,
    q_size: usize,
    ring_size: usize,
) {
    for ri in 0..ring_size {
        let mut fractional_lo = 0.5f64;
        let mut fractional_hi = 0.5f64;
        let mut qxis_lo_mont = Vec::with_capacity(q_size);
        let mut qxis_hi_mont = Vec::with_capacity(q_size);
        izip!(
            q_in.get_col_iter(ri),
            q_over_qi_inv_modqi_times_t_over_qi_fractional.iter(),
            beta_times_q_over_qi_inv_modqi_times_t_over_qi_fractional.iter()
        )
        .for_each(|(qxi, fr_lo, fr_hi)| {
            let qxi_hi = *qxi >> log_beta;
            let qxi_lo = *qxi - (qxi_hi << log_beta);

            fractional_lo += (qxi_lo as f64 * fr_lo);
            fractional_hi += (qxi_hi as f64 * fr_hi);

            qxis_hi_mont.push(modt_operator.normal_to_mont_space(qxi_hi));
            qxis_lo_mont.push(modt_operator.normal_to_mont_space(qxi_lo));
        });

        let out_lo_modt =
            modt_operator.mont_fma(&qxis_lo_mont, q_over_qi_inv_modqi_times_t_over_qi_modt);
        let out_hi_modt = modt_operator.mont_fma(
            &qxis_hi_mont,
            &beta_times_q_over_qi_inv_modqi_times_t_over_qi_modt,
        );

        let out_modt = modt_operator.mont_add(out_hi_modt, out_lo_modt);
        let mut out_modt = modt_operator.mont_to_normal(out_modt);

        // add fractional parts
        out_modt = modt_operator
            .add_mod_fast(out_modt, modt_operator.barrett_reduce(fractional_hi as u64));
        out_modt = modt_operator
            .add_mod_fast(out_modt, modt_operator.barrett_reduce(fractional_lo as u64));

        t_out.set(0, ri, out_modt);
    }
}

/// Given input $x \in R_{QP}$ calculates and returns $[\lceil \frac{t}{P} \cdot x \rfloor]_{Q}$
///
/// We implement "Complex Scaling in CRT representation" of section 2.4 in [HPS](https://eprint.iacr.org/2018/117.pdf)
///
/// TODO (Jay): describe params
pub fn scale_and_round<
    ModOps: MontgomeryBackend<u64, u128> + BarrettBackend<u64, u128>,
>(
    q_out: &mut InnerMat<u64>,
    q_in: &InnerMat<u64>,
    p_in: &InnerMat<u64>,
    modq_operators: &[ModOps],
    modp_operators: &[ModOps],
    qp_over_pj_inv_modpj_times_tq_per_modqi_rational: &[Vec<MontgomeryScalar<u64>>],
    qp_over_pj_inv_modpj_times_tq_fractional: &[f64],
    qp_over_qi_inv_modqi_times_tq_over_qi_modqi: &[MontgomeryScalar<u64>],
    q_size: usize,
    p_size: usize,
    qp_size: usize,
    ring_size: usize,
) {
    for ri in 0..ring_size {
        // summation for fractional can be done without modular reduction per `qi`
        let mut sum_fractional = 0.5f64;
        p_in.get_col_iter(ri).enumerate().for_each(|(j, px_j)| {
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
                .get_col_iter(ri)
                .map(|x| modqi_op.normal_to_mont_space(*x))
                .collect_vec();

            // \sum px_j * \omega_j
            let mut sum_rational = modqi_op.mont_fma(
                &pxjs_in_mont,
                &qp_over_pj_inv_modpj_times_tq_per_modqi_rational[i],
            );

            // qx_i * [(qp/q_i)^{-1}]_{q_i} * (t * q)/q_i \mod qi
            let mut qx_i = modqi_op.normal_to_mont_space(q_in[i][ri]);
            qx_i = modqi_op.mont_mul(qx_i, qp_over_qi_inv_modqi_times_tq_over_qi_modqi[i]);
            sum_rational = modqi_op.mont_add(sum_rational, qx_i);

            let mut sum_rational = modqi_op.mont_to_normal(sum_rational);

            // add fractional part
            sum_rational =
                modqi_op.add_mod_fast(sum_rational, modqi_op.barrett_reduce(sum_fractional));

            q_out.set(i, ri, sum_rational);
        }
    }
}
