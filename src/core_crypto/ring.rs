use super::{
    matrix::{Matrix, MatrixMut, Row, RowMut},
    modulus::{
        BarrettBackend, ModulusArithmeticBackend, ModulusVecBackend, MontgomeryBackend,
        MontgomeryScalar,
    },
    ntt::{Ntt, NttConfig},
    num::UnsignedInteger,
};
use itertools::{izip, Itertools};
use num_traits::AsPrimitive;

pub fn foward_lazy<
    Scalar: UnsignedInteger,
    M: MatrixMut<MatElement = Scalar>,
    N: Ntt<Scalar = Scalar>,
>(
    p: &mut M,
    ntt_ops: &[N],
) where
    <M as Matrix>::R: RowMut,
{
    izip!(p.iter_rows_mut(), ntt_ops.iter()).for_each(|(r, nttop)| nttop.forward_lazy(r.as_mut()));
}

pub fn foward<Scalar: UnsignedInteger, M: MatrixMut<MatElement = Scalar>, N: Ntt<Scalar = Scalar>>(
    p: &mut M,
    ntt_ops: &[N],
) where
    <M as Matrix>::R: RowMut,
{
    izip!(p.iter_rows_mut(), ntt_ops.iter()).for_each(|(r, nttop)| nttop.forward(r.as_mut()));
}

pub fn backward<
    Scalar: UnsignedInteger,
    M: MatrixMut<MatElement = Scalar>,
    N: Ntt<Scalar = Scalar>,
>(
    p: &mut M,
    ntt_ops: &[N],
) where
    <M as Matrix>::R: RowMut,
{
    izip!(p.iter_rows_mut(), ntt_ops.iter()).for_each(|(r, nttop)| nttop.backward(r.as_mut()));
}

pub fn backward_lazy<
    Scalar: UnsignedInteger,
    M: MatrixMut<MatElement = Scalar>,
    N: Ntt<Scalar = Scalar>,
>(
    p: &mut M,
    ntt_ops: &[N],
) where
    <M as Matrix>::R: RowMut,
{
    izip!(p.iter_rows_mut(), ntt_ops.iter()).for_each(|(r, nttop)| nttop.backward_lazy(r.as_mut()));
}

/// Input polynomials coefficents are in rnage [0, 2qi). Reduces coefficients to
/// [0,qi) in-place.
pub fn reduce_from_lazy_mut<
    Scalar: UnsignedInteger,
    MMut: MatrixMut<MatElement = Scalar>,
    ModOps: ModulusVecBackend<Scalar>,
>(
    q: &mut MMut,
    modq_ops: &[ModOps],
) where
    <MMut as Matrix>::R: RowMut,
{
    izip!(q.iter_rows_mut(), modq_ops.iter()).for_each(|(r, modqi)| {
        modqi.reduce_from_lazy_vec(r.as_mut());
    });
}

/// Inputs are assumed to be in Evaluation representation
pub fn mul_lazy_mut<
    T: UnsignedInteger,
    MRef: Matrix<MatElement = T>,
    MMut: MatrixMut<MatElement = T>,
    ModOps: ModulusVecBackend<T>,
>(
    q0: &mut MMut,
    q1: &MRef,
    modq_ops: &[ModOps],
) where
    <MMut as Matrix>::R: RowMut,
{
    debug_assert!(
        q0.dimension() == q1.dimension(),
        "Inputs matrices have unequal dimensions: {:#?}!={:#?}",
        q0.dimension(),
        q1.dimension()
    );

    izip!(q0.iter_rows_mut(), q1.iter_rows(), modq_ops.iter()).for_each(|(r0, r1, modqi)| {
        modqi.mul_lazy_mod_vec(r0.as_mut(), r1.as_ref());
    });
}

pub fn add_lazy_mut<
    Scalar: UnsignedInteger,
    MRef: Matrix<MatElement = Scalar>,
    MMut: MatrixMut<MatElement = Scalar>,
    ModOps: ModulusVecBackend<Scalar>,
>(
    q0: &mut MMut,
    q1: &MRef,
    modq_ops: &[ModOps],
) where
    <MMut as Matrix>::R: RowMut,
{
    debug_assert!(
        q0.dimension() == q1.dimension(),
        "Inputs matrices have unequal dimensions: {:#?}!={:#?}",
        q0.dimension(),
        q1.dimension()
    );

    izip!(q0.iter_rows_mut(), q1.iter_rows(), modq_ops.iter()).for_each(|(r0, r1, modqi)| {
        modqi.add_lazy_mod_vec(r0.as_mut(), r1.as_ref());
    });
}

pub fn add_mut<
    Scalar: UnsignedInteger,
    MRef: Matrix<MatElement = Scalar>,
    MMut: MatrixMut<MatElement = Scalar>,
    ModOps: ModulusVecBackend<Scalar>,
>(
    q0: &mut MMut,
    q1: &MRef,
    modq_ops: &[ModOps],
) where
    <MMut as Matrix>::R: RowMut,
{
    debug_assert!(
        q0.dimension() == q1.dimension(),
        "Input matrices have unequal dimensions: {:#?}!={:#?}",
        q0.dimension(),
        q1.dimension()
    );

    izip!(q0.iter_rows_mut(), q1.iter_rows(), modq_ops.iter()).for_each(|(r0, r1, modqi)| {
        modqi.add_mod_vec(r0.as_mut(), r1.as_ref());
    });
}

pub fn sub_mut<
    Scalar: UnsignedInteger,
    MRef: Matrix<MatElement = Scalar>,
    MMut: MatrixMut<MatElement = Scalar>,
    ModOps: ModulusVecBackend<Scalar>,
>(
    q0: &mut MMut,
    q1: &MRef,
    modq_ops: &[ModOps],
) where
    <MMut as Matrix>::R: RowMut,
{
    debug_assert!(
        q0.dimension() == q1.dimension(),
        "Input matrices have unequal dimensions: {:#?}!={:#?}",
        q0.dimension(),
        q1.dimension()
    );

    izip!(q0.iter_rows_mut(), q1.iter_rows(), modq_ops.iter()).for_each(|(r0, r1, modqi)| {
        modqi.sub_mod_vec(r0.as_mut(), r1.as_ref());
    });
}

pub fn sub_lazy_mut<
    Scalar: UnsignedInteger,
    MRef: Matrix<MatElement = Scalar>,
    MMut: MatrixMut<MatElement = Scalar>,
    ModOps: ModulusVecBackend<Scalar>,
>(
    q0: &mut MMut,
    q1: &MRef,
    modq_ops: &[ModOps],
) where
    <MMut as Matrix>::R: RowMut,
{
    debug_assert!(
        q0.dimension() == q1.dimension(),
        "Input matrices have unequal dimensions: {:#?}!={:#?}",
        q0.dimension(),
        q1.dimension()
    );

    izip!(q0.iter_rows_mut(), q1.iter_rows(), modq_ops.iter()).for_each(|(r0, r1, modqi)| {
        modqi.sub_lazy_mod_vec(r0.as_mut(), r1.as_ref());
    });
}

pub fn neg_mut<
    Scalar: UnsignedInteger,
    MMut: MatrixMut<MatElement = Scalar>,
    ModOps: ModulusVecBackend<Scalar>,
>(
    q: &mut MMut,
    modq_ops: &[ModOps],
) where
    <MMut as Matrix>::R: RowMut,
{
    izip!(q.iter_rows_mut(), modq_ops.iter()).for_each(|(r, modqi)| {
        modqi.neg_mod_vec(r.as_mut());
    });
}

pub fn neg_lazy_mut<
    Scalar: UnsignedInteger,
    MMut: MatrixMut<MatElement = Scalar>,
    ModOps: ModulusVecBackend<Scalar>,
>(
    q: &mut MMut,
    modq_ops: &[ModOps],
) where
    <MMut as Matrix>::R: RowMut,
{
    izip!(q.iter_rows_mut(), modq_ops.iter()).for_each(|(r, modqi)| {
        modqi.neg_lazy_mod_vec(r.as_mut());
    });
}

/// Given input polnyomial x \in Q outputs $[\lceil \frac{P \cdot x}{Q}
/// \rfloor]_P$
///
/// We implement "Modulus Switching between Arbitrary RNS Bases" presented in
/// Appendix E of [2021/204](https://eprint.iacr.org/2021/204.pdf).
pub fn fast_convert_p_over_q<
    MRef: Matrix<MatElement = u64>,
    MMut: MatrixMut<MatElement = u64>,
    ModOps: MontgomeryBackend<u64, u128> + BarrettBackend<u64, u128> + ModulusArithmeticBackend<u64>,
>(
    p_out: &mut MMut,
    q_in: &MRef,
    neg_p_times_q_over_qi_inv_modqi: &[u64],
    qi_inv_per_modpj: &[Vec<MontgomeryScalar<u64>>],
    one_over_qi: &[f64],
    modq_operators: &[ModOps],
    modp_operators: &[ModOps],
    q_size: usize,
    p_size: usize,
    ring_size: usize,
) where
    <MMut as Matrix>::R: RowMut,
{
    debug_assert!(q_in.dimension() == (q_size, ring_size));
    debug_assert!(p_out.dimension() == (p_size, ring_size));
    debug_assert!(modq_operators.len() == q_size);
    debug_assert!(modp_operators.len() == p_size);
    debug_assert!(neg_p_times_q_over_qi_inv_modqi.len() == q_size);
    debug_assert!(qi_inv_per_modpj.len() == p_size);
    debug_assert!(one_over_qi.len() == q_size);

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
    MRef: Matrix<MatElement = u64>,
    MMut: MatrixMut<MatElement = u64>,
    ModOps: MontgomeryBackend<u64, u128>
        + BarrettBackend<u64, u128>
        + ModulusArithmeticBackend<MontgomeryScalar<u64>>,
>(
    p_out: &mut MMut,
    q_in: &MRef,
    q_over_qi_inv_modqi: &[u64],
    q_over_qi_per_modpj: &[Vec<MontgomeryScalar<u64>>],
    mu_times_q_per_modpj: &[Vec<MontgomeryScalar<u64>>],
    one_over_qi: &[f64],
    modq_operators: &[ModOps],
    modp_operators: &[ModOps],
    q_size: usize,
    p_size: usize,
    ring_size: usize,
) where
    <MMut as Matrix>::R: RowMut,
{
    debug_assert!(p_out.dimension() == (p_size, ring_size));
    debug_assert!(q_in.dimension() == (q_size, ring_size));
    debug_assert!(q_over_qi_inv_modqi.len() == q_size);
    debug_assert!(q_over_qi_per_modpj.len() == p_size);
    debug_assert!(mu_times_q_per_modpj.len() == p_size);
    debug_assert!(one_over_qi.len() == q_size);
    debug_assert!(modq_operators.len() == q_size);
    debug_assert!(modp_operators.len() == p_size);

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
            out_in_pj = modpj.sub_mod_fast(out_in_pj, mu_times_q_per_modpj[j][mu]);
            let input_modpj = modpj.mont_to_normal(out_in_pj);

            p_out.set(j, ri, input_modpj);
        }
    }
}

/// Input q_in can have lazy coefficients. Output p_out_lazy has lazy
/// coefficients
pub fn approximate_switch_crt_basis<
    // `num_traits::AsPrimitive<u128> + num_traits::PrimInt,` are a consquenece of BarrettBackend
    // trait. Check issue #12
    S: UnsignedInteger + num_traits::AsPrimitive<u128> + num_traits::PrimInt,
    M: Matrix<MatElement = S>,
    MMut: MatrixMut<MatElement = S>,
    ModOp: ModulusVecBackend<S> + BarrettBackend<S, u128> + MontgomeryBackend<S, u128>,
>(
    p_out_lazy: &mut MMut,
    q_in: &M,
    q_over_qi_inv_modqi: &[S],
    q_over_qi_per_modpj: &[Vec<MontgomeryScalar<S>>],
    modq_operators: &[ModOp],
    modp_operators: &[ModOp],
    ring_size: usize,
    p_size: usize,
    q_size: usize,
) where
    <MMut as Matrix>::R: RowMut,
    u128: AsPrimitive<S>,
{
    debug_assert!(p_out_lazy.dimension() == (p_size, ring_size));
    debug_assert!(q_in.dimension() == (q_size, ring_size));
    debug_assert!(modq_operators.len() == q_size);
    debug_assert!(modp_operators.len() == p_size);
    debug_assert!(q_over_qi_inv_modqi.len() == q_size);
    debug_assert!(q_over_qi_per_modpj.len() == p_size);

    for ri in 0..ring_size {
        let q_values = izip!(
            modq_operators.iter(),
            q_in.get_col_iter(ri),
            q_over_qi_inv_modqi.iter()
        )
        .map(|(modqi, x_qi, q_over_qi_inv)| modqi.mul_mod_fast_lazy(*x_qi, *q_over_qi_inv))
        .collect_vec();

        for j in 0..p_size {
            let modpj = &modp_operators[j];

            // map q_values to mont space \mod pj
            let q_vals_in_mont = q_values
                .iter()
                .map(|v| modpj.normal_to_mont_space_lazy(*v))
                .collect_vec();

            // fma q_over_qi_modpj
            let x = modpj.mont_fma(&q_vals_in_mont, &q_over_qi_per_modpj[j]);

            // map from mont to normal
            let x = modpj.mont_to_normal_lazy(x);

            p_out_lazy.set(j, ri, x);
        }
    }
}

/// Input polynomials can have lazy coefficients. We assume input polynomial
/// part of subbasis Q is in evaluation form and part of subbasis P is in
/// coefficient form. Output poylnomial in basis Q will have lazy coefficients
/// and is in evaluation representation
///
/// Input polynmial x \in R_QP outputs its scaled representation ((1/P) x) \in
/// R_Q.
pub fn approximate_mod_down<
    // `num_traits::AsPrimitive<u128> + num_traits::PrimInt,` are a consquenece of BarrettBackend
    // trait. Check issue #12
    S: UnsignedInteger + num_traits::AsPrimitive<u128> + num_traits::PrimInt,
    M: Matrix<MatElement = S>,
    MMut: MatrixMut<MatElement = S>,
    ModOp: ModulusVecBackend<S> + BarrettBackend<S, u128> + MontgomeryBackend<S, u128>,
    NttOp: Ntt<Scalar = S>,
>(
    q_out_eval_lazy: &mut MMut,
    q_in_eval: &M,
    p_in: &M,
    p_over_pj_inv_modpj: &[S],
    p_over_pj_per_modqi: &[Vec<MontgomeryScalar<S>>],
    p_inv_modqi: &[S],
    modq_operators: &[ModOp],
    modp_operators: &[ModOp],
    q_nttops: &[NttOp],
    ring_size: usize,
    q_size: usize,
    p_size: usize,
) where
    <MMut as Matrix>::R: RowMut,
    u128: AsPrimitive<S>,
{
    debug_assert!(q_out_eval_lazy.dimension() == (q_size, ring_size));
    debug_assert!(q_in_eval.dimension() == (q_size, ring_size));
    debug_assert!(p_in.dimension() == (p_size, ring_size));
    debug_assert!(q_nttops.len() == q_size);
    debug_assert!(modq_operators.len() == q_size);
    debug_assert!(modp_operators.len() == p_size);
    debug_assert!(p_inv_modqi.len() == q_size);
    debug_assert!(p_over_pj_inv_modpj.len() == p_size);
    debug_assert!(p_over_pj_per_modqi.len() == q_size);

    // Switch P subbasis to Q
    approximate_switch_crt_basis(
        q_out_eval_lazy,
        p_in,
        p_over_pj_inv_modpj,
        p_over_pj_per_modqi,
        modp_operators,
        modq_operators,
        ring_size,
        q_size,
        p_size,
    );
    foward_lazy(q_out_eval_lazy, q_nttops);

    neg_lazy_mut(q_out_eval_lazy, modq_operators);
    add_lazy_mut(q_out_eval_lazy, q_in_eval, modq_operators);

    izip!(
        modq_operators.iter(),
        q_out_eval_lazy.iter_rows_mut(),
        p_inv_modqi.iter()
    )
    .for_each(|(modqi, x_qi, p_inv_qi)| modqi.scalar_mul_lazy_mod_vec(x_qi.as_mut(), *p_inv_qi));
}

/// Simple scale and round procedure. Given input polynomial x \in R_Q it
/// outputs scaled polynomial $\frac{t}{Q}\cdot x \in R_t$
///
/// We implement simple scaling procedure as outlined in [HPS] along with digit
/// decomposition technique to reduce error accumulation due to precision loss when $max(log{q_i}) > 51$ bits as outlined in https://eprint.iacr.org/2021/204.pdf.
///
/// Usually $q_i < 2^{61}$. Thus it suffices to limit k = 2 and decomposition
/// base \beta to $2^(max(log{q_i})/2)$.
pub fn simple_scale_and_round<
    MRef: Matrix<MatElement = u64>,
    MMut: MatrixMut<MatElement = u64>,
    ModOps: MontgomeryBackend<u64, u128>
        + BarrettBackend<u64, u128>
        + ModulusArithmeticBackend<u64>
        + ModulusArithmeticBackend<MontgomeryScalar<u64>>,
>(
    t_out: &mut MMut,
    q_in: &MRef,
    q_over_qi_inv_modqi_times_t_over_qi_modt: &[MontgomeryScalar<u64>],
    beta_times_q_over_qi_inv_modqi_times_t_over_qi_modt: &[MontgomeryScalar<u64>],
    q_over_qi_inv_modqi_times_t_over_qi_fractional: &[f64],
    beta_times_q_over_qi_inv_modqi_times_t_over_qi_fractional: &[f64],
    log_beta: usize,
    modt_operator: &ModOps,
    q_size: usize,
    ring_size: usize,
) where
    <MMut as Matrix>::R: RowMut,
{
    debug_assert!(t_out.dimension() == (1, ring_size));
    debug_assert!(q_in.dimension() == (q_size, ring_size));
    debug_assert!(q_over_qi_inv_modqi_times_t_over_qi_modt.len() == q_size);
    debug_assert!(beta_times_q_over_qi_inv_modqi_times_t_over_qi_modt.len() == q_size);
    debug_assert!(q_over_qi_inv_modqi_times_t_over_qi_fractional.len() == q_size);
    debug_assert!(beta_times_q_over_qi_inv_modqi_times_t_over_qi_fractional.len() == q_size);

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

        let out_modt = modt_operator.add_mod_fast(out_hi_modt, out_lo_modt);
        let mut out_modt = modt_operator.mont_to_normal(out_modt);

        // add fractional parts
        out_modt = modt_operator
            .add_mod_fast(out_modt, modt_operator.barrett_reduce(fractional_hi as u64));
        out_modt = modt_operator
            .add_mod_fast(out_modt, modt_operator.barrett_reduce(fractional_lo as u64));

        t_out.set(0, ri, out_modt);
    }
}

/// Given input $x \in R_{QP}$ calculates and returns $[\lceil \frac{t}{P} \cdot
/// x \rfloor]_{Q}$
///
/// We implement "Complex Scaling in CRT representation" of section 2.4 in [HPS](https://eprint.iacr.org/2018/117.pdf)
///
/// TODO (Jay): describe params
pub fn scale_and_round<
    MRef: Matrix<MatElement = u64>,
    MMut: MatrixMut<MatElement = u64>,
    ModOps: MontgomeryBackend<u64, u128>
        + BarrettBackend<u64, u128>
        + ModulusArithmeticBackend<u64>
        + ModulusArithmeticBackend<MontgomeryScalar<u64>>,
>(
    q_out: &mut MMut,
    q_in: &MRef,
    p_in: &MRef,
    modq_operators: &[ModOps],
    modp_operators: &[ModOps],
    qp_over_pj_inv_modpj_times_tq_per_modqi_rational: &[Vec<MontgomeryScalar<u64>>],
    qp_over_pj_inv_modpj_times_tq_fractional: &[f64],
    qp_over_qi_inv_modqi_times_tq_over_qi_modqi: &[MontgomeryScalar<u64>],
    q_size: usize,
    p_size: usize,
    qp_size: usize,
    ring_size: usize,
) where
    <MMut as Matrix>::R: RowMut,
{
    debug_assert!(q_out.dimension() == (q_size, ring_size));
    debug_assert!(q_in.dimension() == (q_size, ring_size));
    debug_assert!(p_in.dimension() == (p_size, ring_size));
    debug_assert!(modq_operators.len() == q_size);
    debug_assert!(modp_operators.len() == p_size);
    debug_assert!(qp_over_pj_inv_modpj_times_tq_per_modqi_rational.len() == q_size);
    debug_assert!(qp_over_pj_inv_modpj_times_tq_fractional.len() == p_size);
    debug_assert!(qp_over_qi_inv_modqi_times_tq_over_qi_modqi.len() == q_size);

    for ri in 0..ring_size {
        // summation for fractional can be done without modular reduction per `qi`
        let mut sum_fractional = 0.5f64;
        p_in.get_col_iter(ri).enumerate().for_each(|(j, px_j)| {
            // px_j * \theta_j
            // TODO (Jay): This will likely result in low precision. A better way is to
            // split fractional into fractional high and fractional low
            sum_fractional += (*px_j as f64) * qp_over_pj_inv_modpj_times_tq_fractional[j];
        });
        let sum_fractional = sum_fractional as u64;

        for i in 0..q_size {
            let modqi_op = modq_operators.get(i).unwrap();

            // convert px_j to montgomery space of \mod{q_i}
            //
            // One thing to note here is \log{pxj} can be much greater than \log{q_i}. For
            // example, when p_j is 60 bits and q_i is 30 bits. Mapping to
            // montgomery space in $q_j$ should still work because the function
            // accepts input in range [0, r) where r = 2^{64}.
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
            let mut qx_i = modqi_op.normal_to_mont_space(*q_in.get_element(i, ri));
            qx_i = modqi_op.mont_mul(qx_i, qp_over_qi_inv_modqi_times_tq_over_qi_modqi[i]);
            sum_rational = modqi_op.add_mod_fast(sum_rational, qx_i);

            let mut sum_rational = modqi_op.mont_to_normal(sum_rational);

            // add fractional part
            sum_rational =
                modqi_op.add_mod_fast(sum_rational, modqi_op.barrett_reduce(sum_fractional));

            q_out.set(i, ri, sum_rational);
        }
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;
    use num_traits::{FromPrimitive, ToPrimitive};

    use crate::{
        core_crypto::{
            matrix::MatrixEntity,
            modulus::{ModulusBackendConfig, NativeModulusBackend},
            ntt::NativeNTTBackend,
            prime::generate_primes_vec,
            random::{DefaultU64SeededRandomGenerator, RandomUniformDist, DEFAULT_U64_SEEDED_RNG},
        },
        utils::{
            convert::{TryConvertFrom, TryConvertFromParts},
            mod_inverse, mod_inverse_big_unit, moduli_chain_to_biguint,
        },
    };

    use super::*;

    #[test]
    fn fast_convert_p_over_q_works() {
        let n = 1 << 4;
        let q_chain = generate_primes_vec(&[60; 10], n, &[]);
        let p_chain = generate_primes_vec(&[60; 10], n, &q_chain);

        let big_q = moduli_chain_to_biguint(&q_chain);
        let big_p = moduli_chain_to_biguint(&p_chain);

        // we will sample radnom polynomial x \in R_Q and calculate [\frac{P}{Q} \cdot
        // x]_P
        let modq_operators = q_chain
            .iter()
            .map(|qi| NativeModulusBackend::initialise(*qi))
            .collect_vec();
        let modp_operators = p_chain
            .iter()
            .map(|pi| NativeModulusBackend::initialise(*pi))
            .collect_vec();

        // precomputes
        let neg_p_times_q_over_qi_inv_modqi = q_chain
            .iter()
            .map(|qi| {
                // (q/q_i)^{-1}_q_i
                let q_over_qi_inv_modqi = mod_inverse(((&big_q / qi) % qi).to_u64().unwrap(), *qi);
                let neg_q_over_qi_inv_modqi_times_p = (qi - ((&big_p * q_over_qi_inv_modqi) % qi))
                    .to_u64()
                    .unwrap();

                neg_q_over_qi_inv_modqi_times_p
            })
            .collect_vec();
        let qi_inv_per_modpj = modp_operators
            .iter()
            .map(|modpj| {
                q_chain
                    .iter()
                    .map(|qi| {
                        modpj.normal_to_mont_space(mod_inverse(*qi % modpj.modulus, modpj.modulus))
                    })
                    .collect_vec()
            })
            .collect_vec();
        let one_over_qi = q_chain.iter().map(|qi| 1f64 / *qi as f64).collect_vec();

        let mut test = DefaultU64SeededRandomGenerator::new();

        let mut poly_q_in = <Vec<Vec<u64>> as MatrixEntity>::zeros(q_chain.len(), n);
        RandomUniformDist::<Vec<Vec<u64>>>::random_fill(&mut test, &q_chain, &mut poly_q_in);
        let mut poly_p_out = Vec::<Vec<u64>>::zeros(p_chain.len(), n);

        fast_convert_p_over_q(
            &mut poly_p_out,
            &poly_q_in,
            &neg_p_times_q_over_qi_inv_modqi,
            &qi_inv_per_modpj,
            &one_over_qi,
            &modq_operators,
            &modp_operators,
            q_chain.len(),
            p_chain.len(),
            n,
        );

        let poly_q_biguint = Vec::<BigUint>::try_convert_from(&poly_q_in, &q_chain);
        let poly_p_biguint = Vec::<BigUint>::try_convert_from(&poly_p_out, &p_chain);

        let poly_p_biguint_expected = poly_q_biguint
            .iter()
            .map(|xi| {
                if xi >= &(&big_q >> 1) {
                    &big_p - ((((&big_p * (&big_q - xi)) + (&big_q >> 1)) / &big_q) % &big_p)
                } else {
                    (((&big_p * xi) + (&big_q >> 1)) / &big_q) % &big_p
                }
            })
            .collect_vec();

        izip!(poly_p_biguint.iter(), poly_p_biguint_expected).for_each(|(p0, p1)| {
            let bits = if p0 < &p1 {
                (p1 - p0).bits()
            } else {
                (p0 - p1).bits()
            };

            assert!(bits == 0);
        });
    }

    #[test]
    fn switch_crt_basis_works() {
        let n = 1 << 4;
        let q_chain = generate_primes_vec(&[60; 10], n, &[]);
        let p_chain = generate_primes_vec(&[60; 10], n, &q_chain);

        let big_q = moduli_chain_to_biguint(&q_chain);
        let big_p = moduli_chain_to_biguint(&p_chain);

        // we will switch a polynomial from basis Q to basis P

        let modq_operators = q_chain
            .iter()
            .map(|qi| NativeModulusBackend::initialise(*qi))
            .collect_vec();
        let modp_operators = p_chain
            .iter()
            .map(|pi| NativeModulusBackend::initialise(*pi))
            .collect_vec();

        // precomputes
        // q/qi and [{q/qi}^{-1}]_q_i
        let (q_over_qi, q_over_qi_inv_modqi): (Vec<BigUint>, Vec<u64>) = q_chain
            .iter()
            .map(|qi| {
                let q_over_qi = &big_q / qi;
                let q_over_qi_inv_modqi = mod_inverse((&q_over_qi % qi).to_u64().unwrap(), *qi);
                (q_over_qi, q_over_qi_inv_modqi)
            })
            .unzip();
        let q_over_qi_per_modpj = modp_operators
            .iter()
            .map(|modpj| {
                q_over_qi
                    .iter()
                    .map(|v| modpj.normal_to_mont_space((v % modpj.modulus).to_u64().unwrap()))
                    .collect_vec()
            })
            .collect_vec();
        let one_over_qi = q_chain.iter().map(|qi| 1f64 / *qi as f64).collect_vec();
        let mu_times_q = (0..q_chain.len()).map(|index| &big_q * index).collect_vec();
        let mu_times_q_per_modpj = modp_operators
            .iter()
            .map(|modpj| {
                mu_times_q
                    .iter()
                    .map(|mu_times_q| {
                        modpj.normal_to_mont_space((mu_times_q % modpj.modulus).to_u64().unwrap())
                    })
                    .collect_vec()
            })
            .collect_vec();

        let mut test_rng = DefaultU64SeededRandomGenerator::new();

        let mut poly_q_in = <Vec<Vec<u64>> as MatrixEntity>::zeros(q_chain.len(), n);
        RandomUniformDist::<Vec<Vec<u64>>>::random_fill(&mut test_rng, &q_chain, &mut poly_q_in);

        let mut poly_p_out = <Vec<Vec<u64>> as MatrixEntity>::zeros(p_chain.len(), n);

        switch_crt_basis(
            &mut poly_p_out,
            &poly_q_in,
            &q_over_qi_inv_modqi,
            &q_over_qi_per_modpj,
            &mu_times_q_per_modpj,
            &one_over_qi,
            &modq_operators,
            &modp_operators,
            q_chain.len(),
            p_chain.len(),
            n,
        );

        let poly_q_biguint = Vec::<BigUint>::try_convert_from(&poly_q_in, &q_chain);
        let poly_p_biguint = Vec::<BigUint>::try_convert_from(&poly_p_out, &p_chain);

        let poly_p_biguint_expected = poly_q_biguint
            .iter()
            .map(|qxi| {
                if qxi >= &(&big_q >> 1) {
                    &big_p - ((&big_q - qxi) % &big_p)
                } else {
                    qxi % &big_p
                }
            })
            .collect_vec();

        izip!(poly_p_biguint.iter(), poly_p_biguint_expected).for_each(|(p0, p1)| {
            let bits = if p0 < &p1 {
                (p1 - p0).bits()
            } else {
                (p0 - p1).bits()
            };
            assert!(bits == 0);
        });
    }

    #[test]
    fn approximate_mod_down_works() {
        let n = 1 << 4;
        let q_chain = generate_primes_vec(&[60; 3], n, &[]);
        let p_chain = generate_primes_vec(&[60; 6], n, &q_chain);

        let big_q = moduli_chain_to_biguint(&q_chain);
        let big_p = moduli_chain_to_biguint(&p_chain);

        // we will scale a polnoymial x \in R_PQ by 1/Q and get its representation 1/Q x
        // \in R_P
        let modq_operators = q_chain
            .iter()
            .map(|qi| NativeModulusBackend::initialise(*qi))
            .collect_vec();
        let modp_operators = p_chain
            .iter()
            .map(|pi| NativeModulusBackend::initialise(*pi))
            .collect_vec();
        let modq_nttops = q_chain
            .iter()
            .map(|qi| NativeNTTBackend::init(*qi, n))
            .collect_vec();
        let modp_nttops = p_chain
            .iter()
            .map(|pi| NativeNTTBackend::init(*pi, n))
            .collect_vec();

        // precomputes
        // approximate_switch_crt_basis
        // q/qi and [{q/qi}^{-1}]_q_i
        let (q_over_qi, q_over_qi_inv_modqi): (Vec<BigUint>, Vec<u64>) = q_chain
            .iter()
            .map(|qi| {
                let q_over_qi = &big_q / qi;
                let q_over_qi_inv_modqi = mod_inverse((&q_over_qi % qi).to_u64().unwrap(), *qi);
                (q_over_qi, q_over_qi_inv_modqi)
            })
            .unzip();
        let q_over_qi_per_modpj = modp_operators
            .iter()
            .map(|modpj| {
                q_over_qi
                    .iter()
                    .map(|v| modpj.normal_to_mont_space((v % modpj.modulus).to_u64().unwrap()))
                    .collect_vec()
            })
            .collect_vec();
        // approximate mod down
        let q_inv_modpj = p_chain
            .iter()
            .map(|pj| {
                mod_inverse_big_unit(&big_q, &BigUint::from_u64(*pj).unwrap())
                    .to_u64()
                    .unwrap()
            })
            .collect_vec();

        let mut test_rng = DefaultU64SeededRandomGenerator::new();

        // P U Q
        let mut poly_q_in = <Vec<Vec<u64>> as MatrixEntity>::zeros(q_chain.len(), n);
        let mut poly_p_in = <Vec<Vec<u64>> as MatrixEntity>::zeros(p_chain.len(), n);
        RandomUniformDist::<Vec<Vec<u64>>>::random_fill(&mut test_rng, &q_chain, &mut poly_q_in);
        RandomUniformDist::<Vec<Vec<u64>>>::random_fill(&mut test_rng, &p_chain, &mut poly_p_in);

        let mut poly_p_in_eval_lazy = poly_p_in.clone();
        foward_lazy(&mut poly_p_in_eval_lazy, &modp_nttops);

        let mut poly_p_out_eval_lazy = Vec::<Vec<u64>>::zeros(p_chain.len(), n);

        approximate_mod_down(
            &mut poly_p_out_eval_lazy,
            &poly_p_in_eval_lazy,
            &poly_q_in,
            &q_over_qi_inv_modqi,
            &q_over_qi_per_modpj,
            &q_inv_modpj,
            &modp_operators,
            &modq_operators,
            &modp_nttops,
            n,
            p_chain.len(),
            q_chain.len(),
        );
        // reduce_from_lazy_mut(&mut poly_p_out_eval_lazy, &modp_operators);
        backward(&mut poly_p_out_eval_lazy, &modp_nttops);
        let poly_p_out = Vec::<BigUint>::try_convert_from(&poly_p_out_eval_lazy, &p_chain);

        let poly_pq_biguint =
            Vec::<BigUint>::try_convert_with_two_parts(&poly_p_in, &poly_q_in, &p_chain, &q_chain);

        // [1/Q * x]_P
        let big_pq = &big_p * &big_q;
        let poly_p_biguint_expected = poly_pq_biguint
            .iter()
            .map(|qxi| {
                if qxi >= &(&big_pq >> 1) {
                    let v = ((&big_pq - qxi) + (&big_q >> 1)) / &big_q;
                    &big_p - (v % &big_p)
                } else {
                    ((qxi + (&big_q >> 1)) / &big_q) % &big_p
                }
            })
            .collect_vec();

        izip!(poly_p_out.iter(), poly_p_biguint_expected).for_each(|(p0, p1)| {
            let bits = if p0 < &p1 {
                (p1 - p0).bits()
            } else {
                (p0 - p1).bits()
            };
            assert!(bits <= 3);
        });
    }

    #[test]
    fn simple_scale_and_round_works() {
        let n = 1 << 4;
        let t = 65537u64;
        let q_chain = generate_primes_vec(&[60; 10], n, &[]);

        let modt_operator = NativeModulusBackend::initialise(t);
        let big_q = moduli_chain_to_biguint(&q_chain);

        // we will scale and round polynomial x \in R_Q to \frac{t \cdot x}{Q} \in R_t

        // precomputes
        let log_beta = ((64 - q_chain.iter().max().unwrap().leading_zeros()) / 2) as usize;
        let beta = 1u64 << log_beta;
        let mut q_over_qi_inv_mod_qi_times_t_over_qi_modt_vec = vec![];
        let mut beta_times_q_over_qi_inv_mod_qi_times_t_over_qi_modt_vec = vec![];
        let mut q_over_qi_inv_mod_qi_times_t_over_qi_fractional_vec = vec![];
        let mut beta_times_q_over_qi_inv_mod_qi_times_t_over_qi_fractional_vec = vec![];
        q_chain.iter().for_each(|qi| {
            let q_over_qi_inv_mod_qi = mod_inverse(((&big_q / qi) % qi).to_u64().unwrap(), *qi);

            let q_over_qi_inv_mod_qi_times_t = BigUint::from(t) * q_over_qi_inv_mod_qi;

            // v_i = ((q/q_i)^{-1}_q_i * t)
            // rational part: v_i / q_i \mod{t}

            q_over_qi_inv_mod_qi_times_t_over_qi_modt_vec.push(modt_operator.normal_to_mont_space(
                ((&q_over_qi_inv_mod_qi_times_t / qi) % t).to_u64().unwrap(),
            ));
            // fractional part: (v_i % q_i) / q_i
            q_over_qi_inv_mod_qi_times_t_over_qi_fractional_vec
                .push(((&q_over_qi_inv_mod_qi_times_t % qi).to_f64().unwrap()) / (*qi as f64));

            // \beta * v_i = (\beta * (q/q_i)^{-1}_q_i * t) / q_i \mod{t}
            // rational part: \beta * v_i / q_i \mod{t}
            beta_times_q_over_qi_inv_mod_qi_times_t_over_qi_modt_vec.push(
                modt_operator.normal_to_mont_space(
                    (((beta * &q_over_qi_inv_mod_qi_times_t) / qi) % t)
                        .to_u64()
                        .unwrap(),
                ),
            );
            // fractional part: ((\beta * v_i) % q_i) / q_i
            beta_times_q_over_qi_inv_mod_qi_times_t_over_qi_fractional_vec.push(
                ((beta * &q_over_qi_inv_mod_qi_times_t) % qi)
                    .to_f64()
                    .unwrap()
                    / (*qi as f64),
            );
        });

        let mut test_rng = DefaultU64SeededRandomGenerator::new();

        let mut poly_q_in = <Vec<Vec<u64>> as MatrixEntity>::zeros(q_chain.len(), n);
        RandomUniformDist::<Vec<Vec<u64>>>::random_fill(&mut test_rng, &q_chain, &mut poly_q_in);

        let mut poly_t_out = Vec::<Vec<u64>>::zeros(1, n);

        simple_scale_and_round(
            &mut poly_t_out,
            &poly_q_in,
            &q_over_qi_inv_mod_qi_times_t_over_qi_modt_vec,
            &beta_times_q_over_qi_inv_mod_qi_times_t_over_qi_modt_vec,
            &q_over_qi_inv_mod_qi_times_t_over_qi_fractional_vec,
            &beta_times_q_over_qi_inv_mod_qi_times_t_over_qi_fractional_vec,
            log_beta,
            &modt_operator,
            q_chain.len(),
            n,
        );

        let poly_q_biguint = Vec::<BigUint>::try_convert_from(&poly_q_in, &q_chain);
        let poly_t_biguint = Vec::<BigUint>::try_convert_from(&poly_t_out, &vec![t]);

        let poly_t_biguint_expected = poly_q_biguint.iter().map(|xi| {
            if xi >= &(&big_q >> 1) {
                t - ((((t * (&big_q - xi)) + (&big_q >> 1)) / &big_q) % t)
            } else {
                (((t * xi) + (&big_q >> 1)) / &big_q) % t
            }
        });

        izip!(poly_t_biguint.iter(), poly_t_biguint_expected).for_each(|(p0, p1)| {
            let bits = if p0 < &p1 {
                (p1 - p0).bits()
            } else {
                (p0 - p1).bits()
            };
            assert!(bits <= 1);
        });
    }

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

        let mut test_rng = DefaultU64SeededRandomGenerator::new();

        // Random polynomial in QP
        let mut poly0_q_part = <Vec<Vec<u64>> as MatrixEntity>::zeros(q_chain.len(), n);
        let mut poly0_p_part = <Vec<Vec<u64>> as MatrixEntity>::zeros(p_chain.len(), n);
        RandomUniformDist::<Vec<Vec<u64>>>::random_fill(&mut test_rng, &q_chain, &mut poly0_q_part);
        RandomUniformDist::<Vec<Vec<u64>>>::random_fill(&mut test_rng, &p_chain, &mut poly0_p_part);

        let mut poly_out = Vec::<Vec<u64>>::zeros(q_chain.len(), n);

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
        let poly_out_big = Vec::<BigUint>::try_convert_from(&poly_out, &q_chain);

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
