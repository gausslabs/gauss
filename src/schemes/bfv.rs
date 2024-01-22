use std::task::Poll;

use itertools::izip;
use rand::{CryptoRng, RngCore};

use crate::{
    ciphertext::{BfvCiphertext, Ciphertext, Representation},
    core_crypto::{
        matrix::{Matrix, MatrixMut, RowMut},
        modulus::ModulusVecBackend,
        ntt::Ntt,
        num::UnsignedInteger,
        random::{RandomGaussianDist, RandomUniformDist},
        ring::{
            self, add_mut, fast_convert_p_over_q, mul_lazy_mut, neg_mut, scale_and_round,
            switch_crt_basis,
        },
    },
    keys::{EncodedMessage, SecretKey},
    parameters::{
        BfvEncryptionParameters, BfvMultiplicationAlgorithm2Parameters, PolyModulusOpParameters,
        PolyNttOpParameters,
    },
    utils::convert::{TryConvertFrom, TryConvertFromParts},
};

fn foward_lazy<
    Scalar: UnsignedInteger,
    Poly: MatrixMut<MatElement = Scalar>,
    N: Ntt<Scalar = Scalar>,
>(
    p: &mut Poly,
    ntt_ops: &[N],
) where
    <Poly as Matrix>::R: RowMut,
{
    izip!(p.iter_rows_mut(), ntt_ops.iter()).for_each(|(r, nttop)| nttop.forward_lazy(r.as_mut()));
}

fn backward<
    Scalar: UnsignedInteger,
    Poly: MatrixMut<MatElement = Scalar>,
    N: Ntt<Scalar = Scalar>,
>(
    p: &mut Poly,
    ntt_ops: &[N],
) where
    <Poly as Matrix>::R: RowMut,
{
    izip!(p.iter_rows_mut(), ntt_ops.iter()).for_each(|(r, nttop)| nttop.backward(r.as_mut()));
}

pub fn secret_key_encryption<
    'a,
    Scalar: 'a + UnsignedInteger,
    Poly: MatrixMut<MatElement = Scalar>,
    S: SecretKey<Scalar = i32>,
    E: EncodedMessage<Scalar = Scalar>,
    P: BfvEncryptionParameters<Scalar = Scalar> + PolyNttOpParameters,
    C: BfvCiphertext<Poly = Poly>,
    R: RandomUniformDist<Scalar = Scalar, Poly = Poly>
        + RandomGaussianDist<Scalar = Scalar, Poly = Poly>
        + CryptoRng,
>(
    secret: &'a S,
    message: &E,
    parameters: &'a P,
    rng: &R,
    level: usize,
) -> C
where
    Poly: TryConvertFrom<&'a [i32], Parameters = &'a [Scalar]>,
    <Poly as Matrix>::R: RowMut,
{
    let q_size = level + 1;
    let ring_size = parameters.ring_size();

    // m(X)
    let mut encoded_m = Vec::<Scalar>::with_capacity(ring_size);
    encoded_m.copy_from_slice(message.value());

    // [Qm(X)]_t
    let modt = parameters.modt_operator();
    let q_modt = parameters.q_modt();
    modt.scalar_mul_mod_vec(&mut encoded_m, q_modt);

    // \Delta m = t^{-1}*[Qm(X)]_t \mod Q
    let modq_ops = parameters.modq_operators_at_level(level);
    let mut m = <C::Poly>::zeros(q_size, ring_size);
    izip!(
        m.iter_rows_mut(),
        parameters.t_inv_modqi_at_level(level).iter(),
        modq_ops.iter()
    )
    .for_each(|(row_i, t_inv_modqi, modqi)| {
        row_i.as_mut().copy_from_slice(&encoded_m);
        modqi.scalar_mul_mod_vec(row_i.as_mut(), *t_inv_modqi);
    });

    let q_moduli_chain = parameters.q_moduli_chain_at_level(level);
    let mut s = <C::Poly>::try_convert_from(secret.values(), q_moduli_chain);
    let mut a = RandomUniformDist::random_ring_poly(rng, q_moduli_chain, ring_size);

    // a*s
    let ntt_ops = parameters.basisq_ntt_ops_at_level(level);
    foward_lazy(&mut s, ntt_ops);
    foward_lazy(&mut a, ntt_ops);
    mul_lazy_mut(&mut s, &a, modq_ops);
    backward(&mut s, ntt_ops);

    // a*s + e
    let e = RandomGaussianDist::random_ring_poly(rng, q_moduli_chain, ring_size);
    add_mut(&mut s, &e, modq_ops);

    // a*s + e + \Delta m
    add_mut(&mut s, &m, modq_ops);

    // -a
    neg_mut(&mut a, modq_ops);

    let c = vec![s, a];
    BfvCiphertext::new(c, level)
}

pub fn ciphertext_add<
    Scalar: UnsignedInteger,
    C: BfvCiphertext<Scalar = Scalar>,
    P: PolyModulusOpParameters<Scalar = Scalar>,
>(
    c0: &mut C,
    c1: &C,
    parameters: &P,
) {
    debug_assert_eq!(
        c0.level(),
        c1.level(),
        "Ciphertext levels are unequal: {} != {}",
        c0.level(),
        c1.level()
    );
    debug_assert!(
        c0.representation() == c1.representation(),
        "Ciphertext c0 and c1 representations are not equal"
    );

    let modq_ops = parameters.modq_vec_ops_at_level(c0.level());
    izip!(c0.c_partq_mut().iter_mut(), c1.c_partq().iter(),).for_each(|(p0, p1)| {
        izip!(p0.iter_rows_mut(), p1.iter_rows(), modq_ops.iter()).for_each(|(r0, r1, modqi)| {
            modqi.add_mod_vec(r0.as_mut(), r1.as_ref());
        });
    });
}

pub fn ciphertext_sub<
    Scalar: UnsignedInteger,
    C: BfvCiphertext<Scalar = Scalar>,
    P: PolyModulusOpParameters<Scalar = Scalar>,
>(
    c0: &mut C,
    c1: &C,
    parameters: &P,
) {
    debug_assert_eq!(
        c0.level(),
        c1.level(),
        "Ciphertext levels are unequal: {} != {}",
        c0.level(),
        c1.level()
    );
    debug_assert!(
        c0.representation() == c1.representation(),
        "Ciphertext c0 and c1 representations are not equal"
    );

    let modq_ops = parameters.modq_vec_ops_at_level(c0.level());
    izip!(c0.c_partq_mut().iter_mut(), c1.c_partq().iter(),).for_each(|(p0, p1)| {
        izip!(p0.iter_rows_mut(), p1.iter_rows(), modq_ops.iter()).for_each(|(r0, r1, modqi)| {
            modqi.sub_mod_vec(r0.as_mut(), r1.as_ref());
        });
    });
}

pub fn ciphertext_mul<
    C: BfvCiphertext<Scalar = u64>,
    P: PolyModulusOpParameters<Scalar = u64>
        + BfvMultiplicationAlgorithm2Parameters<Scalar = u64>
        + PolyNttOpParameters<Scalar = u64>,
>(
    c0: &mut C,
    c1: &C,
    parameters: &P,
) where
    C::Poly: Clone,
{
    debug_assert!(c0.c_partq().len() == 2, "Ciphetext c0 has degree > 2");
    debug_assert!(c1.c_partq().len() == 2, "Ciphetext c1 has degree > 2");
    debug_assert_eq!(
        c0.level(),
        c1.level(),
        "Ciphertext levels are unequal: {} != {}",
        c0.level(),
        c1.level()
    );
    debug_assert!(
        c0.representation() == Representation::Coefficient,
        "Ciphertext c0 must be in Coefficient representation"
    );
    debug_assert!(
        c1.representation() == Representation::Coefficient,
        "Ciphertext c1 must be in Coefficient representation"
    );

    let level = c0.level();

    // In BFV multiplication algorithm 2. Size of Q basis is assumed to equal
    // P basis
    let q_size = level + 1;
    let p_size = level + 1;

    let ring_size = parameters.ring_size();

    // Scale and round c1 in basis Q by \frac{P}{Q} and output c0 in basis P
    let mut poverq_c10_partp = <C::Poly as Matrix>::zeros(p_size, ring_size);
    fast_convert_p_over_q(
        &mut poverq_c10_partp,
        &c1.c_partq()[0],
        parameters.neg_p_times_q_over_qi_inv_modqi_at_level(level),
        parameters.qi_inv_per_modpj_at_level(level),
        parameters.one_over_qi_at_level(level),
        parameters.modq_operators_at_level(level),
        parameters.modp_operators_at_level(level),
        q_size,
        p_size,
        ring_size,
    );
    let mut poverq_c11_partp = <C::Poly>::zeros(p_size, ring_size);
    fast_convert_p_over_q(
        &mut poverq_c11_partp,
        &c1.c_partq()[1],
        parameters.neg_p_times_q_over_qi_inv_modqi_at_level(level),
        parameters.qi_inv_per_modpj_at_level(level),
        parameters.one_over_qi_at_level(level),
        parameters.modq_operators_at_level(level),
        parameters.modp_operators_at_level(level),
        q_size,
        p_size,
        ring_size,
    );

    // Switch $P/Q * c1$ in basis P to basis Q
    let mut poverq_c10_partq = <C::Poly>::zeros(q_size, ring_size);
    switch_crt_basis(
        &mut poverq_c10_partq,
        &poverq_c10_partp,
        parameters.p_over_pj_inv_modpj_at_level(level),
        parameters.p_over_pj_per_modqi_at_level(level),
        parameters.mu_times_p_per_modqi_at_level(level),
        parameters.one_over_pj_at_level(level),
        parameters.modq_operators_at_level(level),
        parameters.modp_operators_at_level(level),
        q_size,
        p_size,
        ring_size,
    );
    let mut poverq_c11_partq = <C::Poly>::zeros(q_size, ring_size);
    switch_crt_basis(
        &mut poverq_c11_partq,
        &poverq_c11_partp,
        parameters.p_over_pj_inv_modpj_at_level(level),
        parameters.p_over_pj_per_modqi_at_level(level),
        parameters.mu_times_p_per_modqi_at_level(level),
        parameters.one_over_pj_at_level(level),
        parameters.modq_operators_at_level(level),
        parameters.modp_operators_at_level(level),
        q_size,
        p_size,
        ring_size,
    );

    // Switch c0 in basis Q to basis P
    let mut c00_partp = <C::Poly>::zeros(p_size, ring_size);
    switch_crt_basis(
        &mut c00_partp,
        &c1.c_partq()[0],
        parameters.q_over_qi_inv_modqi_at_level(level),
        parameters.q_over_qi_per_modpj_at_level(level),
        parameters.mu_times_q_per_modpj_at_level(level),
        parameters.one_over_qi_at_level(level),
        parameters.modq_operators_at_level(level),
        parameters.modp_operators_at_level(level),
        q_size,
        p_size,
        ring_size,
    );
    let mut c01_partp = <C::Poly>::zeros(p_size, ring_size);
    switch_crt_basis(
        &mut c01_partp,
        &c1.c_partq()[1],
        parameters.q_over_qi_inv_modqi_at_level(level),
        parameters.q_over_qi_per_modpj_at_level(level),
        parameters.mu_times_q_per_modpj_at_level(level),
        parameters.one_over_qi_at_level(level),
        parameters.modq_operators_at_level(level),
        parameters.modp_operators_at_level(level),
        q_size,
        p_size,
        ring_size,
    );

    // Convert polynomials in Coefficient representation to lazy Evaluation
    // representation
    {
        // Part Q
        // c0
        let basisq_nttops = parameters.basisq_ntt_ops_at_level(level);
        foward_lazy(&mut c0.c_partq_mut()[0], &basisq_nttops);
        foward_lazy(&mut c0.c_partq_mut()[1], &basisq_nttops);
        // P/Q c1
        foward_lazy(&mut poverq_c10_partq, &basisq_nttops);
        foward_lazy(&mut poverq_c11_partq, &basisq_nttops);

        // Part P
        let basisp_nttops = parameters.basisp_ntt_ops_at_level(level);
        // c0
        foward_lazy(&mut c00_partp, &basisp_nttops);
        foward_lazy(&mut c01_partp, &basisp_nttops);
        // P/Q c1
        foward_lazy(&mut poverq_c10_partp, &basisp_nttops);
        foward_lazy(&mut poverq_c11_partp, &basisp_nttops);
    }

    // part Q tensor
    let mut c_res_partq = {
        let modq_oprators = parameters.modq_vec_ops_at_level(level);

        let mut c00_c10_partq = poverq_c10_partq.clone();
        mul_lazy_mut(&mut c00_c10_partq, &c0.c_partq()[0], modq_oprators);

        // c00 * c11 + c01 * c10
        let mut c00c11_c01c10_partq = poverq_c10_partq;
        mul_lazy_mut(&mut c00c11_c01c10_partq, &c0.c_partq()[1], modq_oprators);
        mul_lazy_mut(&mut c0.c_partq_mut()[1], &poverq_c11_partq, modq_oprators);
        add_mut(&mut c00c11_c01c10_partq, &c0.c_partq()[1], modq_oprators);

        let mut c01_c11_partq = poverq_c11_partq.clone();
        mul_lazy_mut(&mut c01_c11_partq, &c0.c_partq()[1], modq_oprators);

        [c00_c10_partq, c00c11_c01c10_partq, c01_c11_partq]
    };

    // Basis P tensor
    let mut c_res_partp = {
        let modp_oprators = parameters.modp_vec_ops_at_level(level);

        let mut c00_c10_partp = poverq_c10_partp.clone();
        mul_lazy_mut(&mut c00_c10_partp, &c00_partp, modp_oprators);

        // c00 * c11 + c01 * c10
        let mut c00c11_c01c10_partp = poverq_c10_partp;
        mul_lazy_mut(&mut c00c11_c01c10_partp, &c01_partp, modp_oprators);
        mul_lazy_mut(&mut c00_partp, &poverq_c11_partp, modp_oprators);
        add_mut(&mut c00c11_c01c10_partp, &c00_partp, modp_oprators);

        let mut c01_c11_partp = poverq_c11_partp.clone();
        mul_lazy_mut(&mut c01_c11_partp, &c01_partp, modp_oprators);

        [c00_c10_partp, c00c11_c01c10_partp, c01_c11_partp]
    };

    // Convert polynomials in lazy Evaluation representation to Coefficient
    {
        let basisq_nttops = parameters.basisq_ntt_ops_at_level(level);
        c_res_partq
            .iter_mut()
            .for_each(|p| backward(p, &basisq_nttops));

        let basisp_nttops = parameters.basisp_ntt_ops_at_level(level);
        c_res_partp
            .iter_mut()
            .for_each(|p| backward(p, &basisp_nttops))
    }

    // Scale and round polynomials in basis QP by t/P  and output in basis Q
    let c_res_basisq = c0.c_partq_mut();
    izip!(
        c_res_basisq.iter_mut(),
        c_res_partp.iter(),
        c_res_partq.iter()
    )
    .for_each(|(basisq, partp, partq)| {
        scale_and_round(
            basisq,
            partq,
            partp,
            parameters.modq_operators_at_level(level),
            parameters.modp_operators_at_level(level),
            parameters.qp_over_pj_inv_modpj_times_tq_per_modqi_rational_at_level(level),
            parameters.qp_over_pj_inv_modpj_times_tq_fractional_at_level(level),
            parameters.qp_over_qi_inv_modqi_times_tq_over_qi_modqi_at_level(level),
            q_size,
            p_size,
            q_size + p_size,
            ring_size,
        );
    });
}
