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

// Encrytion parameters
pub trait BfvEncryptionParameters: Parameters {
    type ModOp: ModulusVecBackend<Self::Scalar>;

    fn modq_operators_at_level(&self, level: usize) -> &[Self::ModOp];
    fn modt_operator(&self) -> Self::ModOp;
    fn ring_size(&self) -> usize;
    fn q_moduli_chain_at_level(&self, level: usize) -> &[Self::Scalar];

    fn q_modt(&self) -> Self::Scalar;
    fn t_inv_modqi_at_level(&self, level: usize) -> &[Self::Scalar];
}

// Levelling Down Parameters

// SIMD encoding parameters

pub trait BfvEncodingDecodingParameters: Parameters {
    type ModOp: ModulusVecBackend<Self::Scalar>;
    type NttOp: Ntt<Scalar = Self::Scalar>;

    fn ring_size(&self) -> usize;
    fn matrix_representation_index_map(&self) -> Vec<usize>;
    fn t_ntt_op(&self) -> &Self::NttOp;
    fn modt_op(&self) -> &Self::ModOp;
}

// TODO (Jay)
// // Taken from [SEAL](https://github.com/microsoft/SEAL/blob/82b07db635132e297282649e2ab5908999089ad2/native/src/seal/batchencoder.cpp)
//         let row = degree >> 1;
//         let m = degree << 1;
//         let gen = 3;
//         let mut pos = 1;
//         let mut matrix_reps_index_map = vec![0usize; degree];
//         for i in 0..row {
//             let index1 = (pos - 1) >> 1;
//             let index2 = (m - pos - 1) >> 1;
//             matrix_reps_index_map[i] = index1.reverse_bits() >>
// (degree.leading_zeros() + 1);             matrix_reps_index_map[i | row] =
// index2.reverse_bits() >> (degree.leading_zeros() + 1);             pos *=
// gen;             pos &= m - 1;
//         }
