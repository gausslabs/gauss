use itertools::Itertools;
use num_bigint::BigUint;
use num_traits::{FromPrimitive, One, ToPrimitive};

use crate::{
    core_crypto::{
        modulus::{
            BarrettBackend, ModulusArithmeticBackend, ModulusBackendConfig, ModulusVecBackend,
            MontgomeryBackend, MontgomeryScalar,
        },
        ntt::{Ntt, NttConfig},
        num::{BFloat, ComplexNumber, UnsignedInteger},
    },
    utils::{mod_inverse, mod_inverse_big_unit},
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

pub trait BfvMultiplicationAlgorithm2Parameters: Parameters<Scalar = u64> {
    type ModOp: MontgomeryBackend<u64, u128>
        + BarrettBackend<u64, u128>
        + ModulusArithmeticBackend<u64>
        + ModulusArithmeticBackend<MontgomeryScalar<u64>>;

    fn max_level(&self) -> usize;

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

// Encrytion parameters
pub trait BfvEncryptionParameters: Parameters {
    type ModOp: ModulusVecBackend<Self::Scalar>;
    type NttOp: Ntt<Scalar = Self::Scalar>;

    fn max_level(&self) -> usize;
    fn modq_ops_at_level(&self, level: usize) -> &[Self::ModOp];
    fn basisq_ntt_ops_at_level(&self, level: usize) -> &[Self::NttOp];
    fn modt_op(&self) -> &Self::ModOp;
    fn ring_size(&self) -> usize;
    fn q_moduli_chain_at_level(&self, level: usize) -> &[Self::Scalar];

    fn q_modt_at_level(&self, level: usize) -> &Self::Scalar;
    fn neg_t_inv_modqi_at_level(&self, level: usize) -> &[Self::Scalar];
}

pub trait BfvDecryptionParameters: Parameters {
    type NttOp: Ntt<Scalar = Self::Scalar>;
    type ModOps: ModulusVecBackend<Self::Scalar>
        + MontgomeryBackend<u64, u128>
        + BarrettBackend<u64, u128>
        + ModulusArithmeticBackend<Self::Scalar>
        + ModulusArithmeticBackend<MontgomeryScalar<Self::Scalar>>;

    fn basisq_ntt_ops_at_level(&self, level: usize) -> &[Self::NttOp];
    fn q_moduli_chain_at_level(&self, level: usize) -> &[Self::Scalar];
    fn modq_ops_at_level(&self, level: usize) -> &[Self::ModOps];
    fn modt_op(&self) -> &Self::ModOps;
    fn ring_size(&self) -> usize;
    fn max_level(&self) -> usize;

    fn q_over_qi_inv_modqi_times_t_over_qi_modt_at_level(
        &self,
        level: usize,
    ) -> &[MontgomeryScalar<Self::Scalar>];
    fn beta_times_q_over_qi_inv_modqi_times_t_over_qi_modt_at_level(
        &self,
        level: usize,
    ) -> &[MontgomeryScalar<Self::Scalar>];
    fn q_over_qi_inv_modqi_times_t_over_qi_fractional_at_level(&self, level: usize) -> &[f64];
    fn beta_times_q_over_qi_inv_modqi_times_t_over_qi_fractional_at_level(
        &self,
        level: usize,
    ) -> &[f64];
    fn log_beta(&self) -> usize;
}

// Levelling Down Parameters

// SIMD encoding parameters

pub trait BfvEncodingDecodingParameters: Parameters {
    type ModOp: ModulusVecBackend<Self::Scalar>;
    type NttOp: Ntt<Scalar = Self::Scalar>;

    fn ring_size(&self) -> usize;
    fn matrix_representation_index_map(&self) -> &[usize];
    fn t_ntt_op(&self) -> &Self::NttOp;
    fn modt_op(&self) -> &Self::ModOp;
}

// CKKS encoding decoding parameters
pub trait CkksEncDecParameters: Parameters {
    type F: BFloat;
    type Complex: ComplexNumber<Self::F>;
    type BU;
    type ModOp: ModulusVecBackend<Self::Scalar>;
    type NttOp: Ntt<Scalar = Self::Scalar>;

    fn delta(&self) -> Self::F;
    fn psi_powers(&self) -> &[Self::Complex];
    fn rot_group(&self) -> &[usize];
    fn ring_size(&self) -> usize;
    fn bigq_at_level(&self, level: usize) -> &Self::BU;
    fn q_moduli_chain_at_level(&self, level: usize) -> &[Self::Scalar];

    fn q_modops_at_level(&self, level: usize) -> &[Self::ModOp];
    fn q_nttops_at_level(&self, level: usize) -> &[Self::NttOp];
}
