use itertools::{izip, Itertools};
use num_traits::Signed;
use rand::{CryptoRng, RngCore};

use crate::{
    ciphertext::{Ciphertext, InitialiseLevelledCiphertext, Representation, RlweCiphertext},
    core_crypto::{
        matrix::{Matrix, MatrixMut, RowMut},
        modulus::ModulusVecBackend,
        ntt::Ntt,
        num::UnsignedInteger,
        random::{
            RandomGaussianDist, RandomSecretByteGenerator, RandomSecretValueGenerator,
            RandomUniformDist,
        },
        ring::{
            add_lazy_mut, add_mut, backward, fast_convert_p_over_q, foward_lazy, mul_lazy_mut,
            neg_mut, scale_and_round, simple_scale_and_round, switch_crt_basis,
        },
    },
    keys::SecretKey,
    parameters::{
        BfvDecryptionParameters, BfvEncodingDecodingParameters, BfvEncryptionParameters,
        BfvMultiplicationAlgorithm2Parameters, PolyModulusOpParameters, PolyNttOpParameters,
    },
    utils::convert::TryConvertFrom,
};

pub fn generate_ternery_secret_with_hamming_weight<
    Scalar: Signed + Clone,
    R: RandomSecretByteGenerator + RandomSecretValueGenerator<usize>,
>(
    rng: &mut R,
    hamming_weight: usize,
    ring_size: usize,
) -> Vec<Scalar> {
    let mut bytes = rng.random_bytes(hamming_weight.div_ceil(8));

    let mut secret = vec![Scalar::zero(); ring_size];
    let mut secret_indices = (0..ring_size).into_iter().collect_vec();
    let mut bit_index = 0;
    let mut byte_index = 0;
    for _ in 0..hamming_weight {
        let secret_index = rng.random_value_in_range(secret_indices.len());

        if bytes[byte_index] & 1u8 == 1 {
            secret[secret_indices[secret_index]] = Scalar::one();
        } else {
            secret[secret_indices[secret_index]] = Scalar::neg(Scalar::one());
        }

        bytes[byte_index] >>= 1;

        // remove secret_index from secret_indices
        secret_indices[secret_index] = *secret_indices.last().unwrap();
        secret_indices.truncate(secret_indices.len() - 1);

        if bit_index == 7 {
            bit_index = 0;
            byte_index += 1;
        } else {
            bit_index += 1;
        }
    }

    secret
}

pub fn simd_encode_message<
    Scalar: UnsignedInteger,
    P: BfvEncodingDecodingParameters<Scalar = Scalar>,
>(
    m: &[Scalar],
    parameters: &P,
) -> Vec<Scalar> {
    debug_assert!(
        m.len() <= parameters.ring_size(),
        "Message length {} > ring size {}",
        m.len(),
        parameters.ring_size()
    );

    let ring_size = parameters.ring_size();

    let matrix_representation_index = parameters.matrix_representation_index_map();
    let mut message = vec![Scalar::zero(); ring_size];
    m.iter().enumerate().for_each(|(index, v)| {
        message[matrix_representation_index[index]] = *v;
    });

    parameters.modt_op().reduce_vec(&mut message);
    parameters.t_ntt_op().backward(&mut message);

    message
}

pub fn simd_decode_message<
    Scalar: UnsignedInteger,
    P: BfvEncodingDecodingParameters<Scalar = Scalar>,
>(
    m: &[Scalar],
    parameters: &P,
) -> Vec<Scalar> {
    let mut m = m.to_vec();
    parameters.t_ntt_op().forward(&mut m);

    let matrix_representation_index = parameters.matrix_representation_index_map();
    let mut message = vec![Scalar::zero(); parameters.ring_size()];
    message
        .iter_mut()
        .enumerate()
        .for_each(|(index, v)| *v = m[matrix_representation_index[index]]);
    message
}

pub fn secret_key_encryption<
    Scalar: UnsignedInteger,
    Poly: MatrixMut<MatElement = Scalar>,
    S: SecretKey<Scalar = i32>,
    P: BfvEncryptionParameters<Scalar = Scalar>,
    C: RlweCiphertext<Poly = Poly> + InitialiseLevelledCiphertext<C = Vec<Poly>>,
    R: RandomUniformDist<Poly, Parameters = [Scalar]>
        + RandomGaussianDist<Poly, Parameters = [Scalar]>
        + CryptoRng,
>(
    secret: &S,
    message: &[Scalar],
    parameters: &P,
    rng: &mut R,
    level: usize,
) -> C
where
    Poly: TryConvertFrom<[i32], Parameters = [Scalar]>,
    <Poly as Matrix>::R: RowMut,
{
    let q_size = parameters.max_level() - level + 1;
    let ring_size = parameters.ring_size();

    // m(X)
    let mut encoded_m = message.to_vec();

    // [Qm(X)]_t
    let modt = parameters.modt_op();
    let q_modt = parameters.q_modt_at_level(level);
    modt.scalar_mul_mod_vec(&mut encoded_m, *q_modt);

    // \Delta m = (-t^{-1}*[Qm(X)]_t) \mod Q
    let modq_ops = parameters.modq_ops_at_level(level);
    let mut m = <C::Poly>::zeros(q_size, ring_size);
    izip!(
        m.iter_rows_mut(),
        parameters.neg_t_inv_modqi_at_level(level).iter(),
        modq_ops.iter()
    )
    .for_each(|(row_i, neg_t_inv_modqi, modqi)| {
        row_i.as_mut().copy_from_slice(&encoded_m);
        modqi.scalar_mul_mod_vec(row_i.as_mut(), *neg_t_inv_modqi);
    });

    let q_moduli_chain = parameters.q_moduli_chain_at_level(level);
    let mut s = Poly::try_convert_from(&secret.values(), &q_moduli_chain);
    let mut a_eval = Poly::zeros(q_moduli_chain.len(), ring_size);
    RandomUniformDist::random_fill(rng, q_moduli_chain, &mut a_eval);

    // a*s
    let ntt_ops = parameters.basisq_ntt_ops_at_level(level);
    foward_lazy(&mut s, ntt_ops);
    mul_lazy_mut(&mut s, &a_eval, modq_ops);

    backward(&mut s, ntt_ops);

    // a*s + e
    let mut e = Poly::zeros(q_moduli_chain.len(), ring_size);
    RandomGaussianDist::random_fill(rng, q_moduli_chain, &mut e);
    add_mut(&mut s, &e, modq_ops);

    // a*s + e + \Delta m
    add_mut(&mut s, &m, modq_ops);

    // -a
    backward(&mut a_eval, ntt_ops);
    let mut neg_a = a_eval;
    neg_mut(&mut neg_a, modq_ops);

    let c = vec![s, neg_a];
    InitialiseLevelledCiphertext::new(c, level, Representation::Coefficient)
}

pub fn secret_key_decryption<
    C: RlweCiphertext<Scalar = u64>,
    P: BfvDecryptionParameters<Scalar = u64>,
    S: SecretKey<Scalar = i32>,
>(
    c: &C,
    secret: &S,
    parameters: &P,
) -> C::Poly
where
    <C as Ciphertext>::Poly: Clone + TryConvertFrom<[i32], Parameters = [u64]>,
{
    let level = c.level();

    let basisq_ntt_ops = parameters.basisq_ntt_ops_at_level(level);

    let q_moduli_chain = parameters.q_moduli_chain_at_level(level);

    let mut s = C::Poly::try_convert_from(&secret.values(), &q_moduli_chain);
    foward_lazy(&mut s, basisq_ntt_ops);

    let mut c0 = c.c_partq()[0].clone();
    if c.representation() == Representation::Coefficient {
        foward_lazy(&mut c0, basisq_ntt_ops);
    }

    let s_clone = s.clone();
    let modq_ops = parameters.modq_ops_at_level(level);
    for i in 1..c.c_partq().len() {
        let mut ci = if c.representation() == Representation::Coefficient {
            let mut ci = c.c_partq()[i].clone();
            foward_lazy(&mut ci, basisq_ntt_ops);
            ci
        } else {
            c.c_partq()[i].clone()
        };

        mul_lazy_mut(&mut ci, &s, modq_ops);
        add_lazy_mut(&mut c0, &ci, modq_ops);
        mul_lazy_mut(&mut s, &s_clone, modq_ops);
    }

    backward(&mut c0, basisq_ntt_ops);

    let ring_size = parameters.ring_size();
    let q_size = parameters.max_level() - level + 1;
    let mut t_out = C::Poly::zeros(1, ring_size);
    simple_scale_and_round(
        &mut t_out,
        &c0,
        parameters.q_over_qi_inv_modqi_times_t_over_qi_modt_at_level(level),
        parameters.beta_times_q_over_qi_inv_modqi_times_t_over_qi_modt_at_level(level),
        parameters.q_over_qi_inv_modqi_times_t_over_qi_fractional_at_level(level),
        parameters.beta_times_q_over_qi_inv_modqi_times_t_over_qi_fractional_at_level(level),
        parameters.log_beta(),
        parameters.modt_op(),
        q_size,
        ring_size,
    );
    t_out
}

pub fn ciphertext_add<
    Scalar: UnsignedInteger,
    C: RlweCiphertext<Scalar = Scalar>,
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
    C: RlweCiphertext<Scalar = Scalar>,
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
    C: RlweCiphertext<Scalar = u64>,
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
    let q_size = parameters.max_level() - level + 1;
    let p_size = parameters.max_level() - level + 1;

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
