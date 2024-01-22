use crate::core_crypto::{
    modulus::{
        BarrettBackend, ModulusArithmeticBackend, ModulusVecBackend, MontgomeryBackend,
        MontgomeryScalar,
    },
    ntt::Ntt,
    num::UnsignedInteger,
};

pub trait Parameters {
    type Scalar: UnsignedInteger;
}

pub trait PolyNttOpParameters: Parameters {
    type NttOp: Ntt<Scalar = Self::Scalar>;
    fn basisq_ntt_ops_at_level(&self, level: usize) -> &[Self::NttOp];
    fn basisp_ntt_ops_at_level(&self, level: usize) -> &[Self::NttOp];
}

pub trait PolyModulusOpParameters: Parameters {
    type ModVecOp: ModulusVecBackend<Self::Scalar>;

    fn modq_vec_ops_at_level(&self, level: usize) -> &[Self::ModVecOp];
    fn modp_vec_ops_at_level(&self, level: usize) -> &[Self::ModVecOp];

    fn ring_size(&self) -> usize;

    fn basisq_dimension_at_level(&self, level: usize) -> (usize, usize);
    fn basisp_dimension_at_level(&self, level: usize) -> (usize, usize);
}

// Encrytion parameters

// Levelling Down Parameters

// SIMD encoding parameters

pub trait BfvMultiplicationAlgorithm2Parameters: Parameters<Scalar = u64> {
    type ModOp: MontgomeryBackend<u64, u128>
        + BarrettBackend<u64, u128>
        + ModulusArithmeticBackend<u64>
        + ModulusArithmeticBackend<MontgomeryScalar<u64>>;

    fn modq_operators_at_level(&self, level: usize) -> &[Self::ModOp];
    fn modp_operators_at_level(&self, level: usize) -> &[Self::ModOp];

    // Switch polynomial basis from Q to P //
    fn q_over_qi_inv_modqi_at_level(&self, level: usize) -> &[Self::Scalar];
    fn q_over_qi_per_modpj_at_level(&self, level: usize) -> &[Vec<MontgomeryScalar<Self::Scalar>>];
    fn mu_times_q_per_modpj_at_level(&self, level: usize)
        -> &[Vec<MontgomeryScalar<Self::Scalar>>];
    fn one_over_qi_at_level(&self, level: usize) -> &[f64];

    // Switch polynomial basis from P to Q //
    fn p_over_pj_inv_modpj_at_level(&self, level: usize) -> &[Self::Scalar];
    fn p_over_pj_per_modqi_at_level(&self, level: usize) -> &[Vec<MontgomeryScalar<Self::Scalar>>];
    fn mu_times_p_per_modqi_at_level(&self, level: usize)
        -> &[Vec<MontgomeryScalar<Self::Scalar>>];
    fn one_over_pj_at_level(&self, level: usize) -> &[f64];

    // Fast convert polynomial in basis Q to basis in P after scaling by \frac{P}{Q}
    // //
    fn neg_p_times_q_over_qi_inv_modqi_at_level(&self, level: usize) -> &[Self::Scalar];
    fn qi_inv_per_modpj_at_level(&self, level: usize) -> &[Vec<MontgomeryScalar<Self::Scalar>>];
    // fn one_over_qi_at_level(&self, level: usize) -> &[f64];

    // Scale and round polynomial in basis QP by \frac{t}{P} and output in basis Q
    // //
    fn qp_over_pj_inv_modpj_times_tq_per_modqi_rational_at_level(
        &self,
        level: usize,
    ) -> &[Vec<MontgomeryScalar<Self::Scalar>>];
    fn qp_over_pj_inv_modpj_times_tq_fractional_at_level(&self, level: usize) -> &[f64];
    fn qp_over_qi_inv_modqi_times_tq_over_qi_modqi_at_level(
        &self,
        level: usize,
    ) -> &[MontgomeryScalar<u64>];
}
