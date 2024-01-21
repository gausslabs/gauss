use itertools::izip;

use crate::{
    ciphertext::{BfvCiphertext, Ciphertext, Representation},
    core_crypto::{
        matrix::{Matrix, MatrixMut},
        modulus::ModulusVecBackend,
        ring::{fast_convert_p_over_q, switch_crt_basis},
    },
    parameters::{BfvMultiplicationAlgorithm2Parameters, PolyModulusOpParameters},
};

pub fn ciphertext_add<C: BfvCiphertext, P: PolyModulusOpParameters>(
    c0: &C,
    c1: &C,
    parameters: &P,
) -> C {
    debug_assert_eq!(
        c0.level(),
        c1.level(),
        "Ciphertext levels are unequal: {} != {}",
        c0.level(),
        c1.level()
    );

    let mod_ops = parameters.mod_vec_ops_at_level(c0.level());
    izip!(c0.c_basisq().iter(), c1.c_basisq().iter()).map(|(p0s, p1s)| {
        // TODO
    });

    todo!()
}

pub fn ciphertext_mul<
    C: BfvCiphertext,
    P: PolyModulusOpParameters + BfvMultiplicationAlgorithm2Parameters,
>(
    c0: &C,
    c1: &C,
    parameters: &P,
) {
    debug_assert!(c0.c_basisq().len() == 2, "Ciphetext c0 has degree > 2");
    debug_assert!(c1.c_basisq().len() == 2, "Ciphetext c1 has degree > 2");
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

    // // In BFV multiplication algorithm 2. Size of Q basis is assumed to equal
    // P // basis
    // let q_size = level + 1;
    // let p_size = level + 1;

    // let ring_size = parameters.ring_size();

    // // Scale and round c0 in basis Q by \frac{P}{Q} and output c0 in basis P
    // let (basisp_rows, basisp_cols) =
    // parameters.basisp_dimension_at_level(level); let mut c00_basisp =
    // <C::Poly as Matrix<u64>>::zeros(p_size, ring_size);
    // fast_convert_p_over_q(
    //     &mut c00_basisp,
    //     &c0.c_basisq()[0],
    //     parameters.neg_p_times_q_over_qi_inv_modqi_at_level(level),
    //     parameters.qi_inv_per_modpj_at_level(level),
    //     parameters.one_over_qi_at_level(level),
    //     parameters.modq_operators_at_level(level),
    //     parameters.modp_operators_at_level(level),
    //     q_size,
    //     p_size,
    //     ring_size,
    // );
    // let mut c01_basisp = <C::Poly as Matrix<u64>>::zeros(p_size, ring_size);
    // fast_convert_p_over_q(
    //     &mut c01_basisp,
    //     &c0.c_basisq()[0],
    //     parameters.neg_p_times_q_over_qi_inv_modqi_at_level(level),
    //     parameters.qi_inv_per_modpj_at_level(level),
    //     parameters.one_over_qi_at_level(level),
    //     parameters.modq_operators_at_level(level),
    //     parameters.modp_operators_at_level(level),
    //     q_size,
    //     p_size,
    //     ring_size,
    // );

    // Switch c0 in basis P to basis Q
    // let mut c00_basisq = <C::Poly as Matrix<u64>>::zeros(q_size, ring_size);
    // switch_crt_basis(
    //     &mut c00_basisq,
    //     &c00_basisp,
    //     parameters.p_over_pj_inv_modpj_at_level(level),
    //     parameters.p_over_pj_per_modqi_at_level(level),
    //     parameters.mu_times_p_per_modqi_at_level(level),
    //     parameters.one_over_pj_at_level(level),
    //     parameters.modq_operators_at_level(level),
    //     parameters.modp_operators_at_level(level),
    //     q_size,
    //     p_size,
    //     ring_size,
    // );
}
