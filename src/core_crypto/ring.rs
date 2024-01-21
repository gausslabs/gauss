use super::{
    matrix::{Matrix, MatrixMut, RowMut},
    modulus::{BarrettBackend, MontgomeryBackend, MontgomeryScalar},
};
use itertools::{izip, Itertools};

// pub trait MatrixRef<'a, T: 'a>: Matrix<T> {
//     type IteratorRef: Iterator<Item = &'a T>;

//     fn get_col(&'a self, index: usize) -> Self::IteratorRef;
//     fn get_row(&'a self, index: usize) -> Self::IteratorRef;

//     fn get(&'a self, row: usize, col: usize) -> &T;
// }

// pub trait MatrixMut<'a, T: 'a> {
//     type IteratorMutRef: Iterator<Item = &'a mut T>;

//     fn get_col_mut(&'a mut self, index: usize) -> Self::IteratorMutRef;
//     fn get_row_mut(&'a mut self, index: usize) -> Self::IteratorMutRef;

//     fn get_mut(&'a mut self, row: usize, col: usize) -> &T;

//     fn set(&mut self, row: usize, col: usize, value: T);
// }

// pub trait Matrix<T> {
//     fn zeros(rows: usize, cols: usize) -> Self;
//     fn from_values(rows: usize, cols: usize, values: Vec<T>) -> Self;

//     fn dimension(&self) -> (usize, usize);
// }

// pub fn add_mut<
//     'a,
//     MRef: MatrixRef<'a, u64>,
//     MMut: MatrixMut<'a, u64>,
//     ModOps: MontgomeryBackend<u64, u128> + BarrettBackend<u64, u128>,
// >() {
// }

// pub fn add<
//     'a,
//     MRef: MatrixRef<'a, u64>,
//     ModOps: MontgomeryBackend<u64, u128> + BarrettBackend<u64, u128>,
// >() {
// }

/// Given input polnyomial x \in Q outputs $[\lceil \frac{P \cdot x}{Q}
/// \rfloor]_P$
///
/// We implement "Modulus Switching between Arbitrary RNS Bases" presented in
/// Appendix E of [2021/204](https://eprint.iacr.org/2021/204.pdf).
pub fn fast_convert_p_over_q<
    MRef: Matrix<MatElement = u64>,
    MMut: MatrixMut<MatElement = u64>,
    ModOps: MontgomeryBackend<u64, u128> + BarrettBackend<u64, u128>,
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
    ModOps: MontgomeryBackend<u64, u128> + BarrettBackend<u64, u128>,
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
    ModOps: MontgomeryBackend<u64, u128> + BarrettBackend<u64, u128>,
>(
    t_out: &mut MMut,
    q_in: &MRef,
    q_over_qi_inv_modqi_times_t_over_qi_modt: &[MontgomeryScalar<u64>],
    beta_times_q_over_qi_inv_modqi_times_t_over_qi_modt: &[MontgomeryScalar<u64>],
    q_over_qi_inv_modqi_times_t_over_qi_fractional: &[f64],
    beta_times_q_over_qi_inv_modqi_times_t_over_qi_fractional: &[f64],
    log_beta: usize,
    modt_operator: ModOps,
    q_size: usize,
    ring_size: usize,
) where
    <MMut as Matrix>::R: RowMut,
{
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

/// Given input $x \in R_{QP}$ calculates and returns $[\lceil \frac{t}{P} \cdot
/// x \rfloor]_{Q}$
///
/// We implement "Complex Scaling in CRT representation" of section 2.4 in [HPS](https://eprint.iacr.org/2018/117.pdf)
///
/// TODO (Jay): describe params
pub fn scale_and_round<
    MRef: Matrix<MatElement = u64>,
    MMut: MatrixMut<MatElement = u64>,
    ModOps: MontgomeryBackend<u64, u128> + BarrettBackend<u64, u128>,
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
    for ri in 0..ring_size {
        // summation for fractional can be done without modular reduction per `qi`
        let mut sum_fractional = 0.5f64;
        p_in.get_col_iter(ri).enumerate().for_each(|(j, px_j)| {
            // px_j * \theta_j
            // TODO (Jay): This will likely result in low precision. A better will be to
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
            sum_rational = modqi_op.mont_add(sum_rational, qx_i);

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
            modulus::{ModulusBackendConfig, NativeModulusBackend},
            prime::generate_primes_vec,
            random::RandomUniformDist,
        },
        utils::{
            convert::TryConvertFromParts, mod_inverse, moduli_chain_to_biguint, test_utils::TestRng,
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

        let test = TestRng {};

        let poly_q_in = test.random_ring_poly(&q_chain, n);
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

        let poly_q_biguint = Vec::<BigUint>::try_convert_with_one_part(&poly_q_in, &q_chain);
        let poly_p_biguint = Vec::<BigUint>::try_convert_with_one_part(&poly_p_out, &p_chain);

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

        let test_rng = TestRng {};

        let poly_q_in = test_rng.random_ring_poly(&q_chain, n);
        let mut poly_p_out = Vec::<Vec<u64>>::zeros(p_chain.len(), n);

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

        let poly_q_biguint = Vec::<BigUint>::try_convert_with_one_part(&poly_q_in, &q_chain);
        let poly_p_biguint = Vec::<BigUint>::try_convert_with_one_part(&poly_p_out, &p_chain);

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

        let test_rng = TestRng {};

        let poly_q_in = test_rng.random_ring_poly(&q_chain, n);
        let mut poly_t_out = Vec::<Vec<u64>>::zeros(1, n);

        simple_scale_and_round(
            &mut poly_t_out,
            &poly_q_in,
            &q_over_qi_inv_mod_qi_times_t_over_qi_modt_vec,
            &beta_times_q_over_qi_inv_mod_qi_times_t_over_qi_modt_vec,
            &q_over_qi_inv_mod_qi_times_t_over_qi_fractional_vec,
            &beta_times_q_over_qi_inv_mod_qi_times_t_over_qi_fractional_vec,
            log_beta,
            modt_operator,
            q_chain.len(),
            n,
        );

        let poly_q_biguint = Vec::<BigUint>::try_convert_with_one_part(&poly_q_in, &q_chain);
        let poly_t_biguint = Vec::<BigUint>::try_convert_with_one_part(&poly_t_out, &vec![t]);

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

        let test_rng = TestRng {};

        // Random polynomial in QP
        let poly0_q_part = test_rng.random_ring_poly(&q_chain, n);
        let poly0_p_part = test_rng.random_ring_poly(&p_chain, n);

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
