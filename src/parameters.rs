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
        num::UnsignedInteger,
        ring,
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

    fn modq_ops_at_level(&self, level: usize) -> &[Self::ModOp];
    fn modt_op(&self) -> &Self::ModOp;
    fn ring_size(&self) -> usize;
    fn q_moduli_chain_at_level(&self, level: usize) -> &[Self::Scalar];

    fn q_modt_at_level(&self, level: usize) -> &Self::Scalar;
    fn t_inv_modqi_at_level(&self, level: usize) -> &[Self::Scalar];
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

pub struct BfvClientParametersForScalarU64<M, N> {
    ring_size: usize,
    q_moduli_chain: Vec<u64>,
    levels_plus_1: usize,
    modq_ops: Vec<M>,
    modt_op: M,
    t_ntt_op: N,
    q_modt_at_level: Vec<u64>,
    t_inv_modqi_at_level: Vec<Vec<u64>>,
    matrix_representation_index_map: Vec<usize>,
}

impl<M, N> Parameters for BfvClientParametersForScalarU64<M, N> {
    type Scalar = u64;
}

impl<M: ModulusVecBackend<u64> + ModulusBackendConfig<u64>, N: NttConfig<Scalar = u64> + Ntt>
    BfvClientParametersForScalarU64<M, N>
{
    fn new(
        q_moduli_chain: &[u64],
        p_moduli_chain: &[u64],
        specialp_moduli_chain: &[u64],
        t: u64,
        ring_size: usize,
    ) -> Self {
        let levels = q_moduli_chain.len() - 1;

        let mut q = BigUint::one();
        for qi in q_moduli_chain.iter() {
            q *= *qi;
        }

        let mut q_clone = q.clone();
        let mut q_modt_at_level = vec![];
        let mut t_inv_modqi_at_level = vec![];
        let big_t = BigUint::from_u64(t).unwrap();
        for i in 0..levels + 1 {
            q_modt_at_level.push((&q_clone % t).to_u64().unwrap());

            let t_inv_modq = mod_inverse_big_unit(&big_t, &q_clone);
            let mut t_inv_modqi = vec![];
            for j in 0..i + 1 {
                t_inv_modqi.push((&t_inv_modq % q_moduli_chain[j]).to_u64().unwrap());
            }
            t_inv_modqi_at_level.push(t_inv_modqi);

            q_clone /= q_moduli_chain[levels - i];
        }

        // Taken from [SEAL](https://github.com/microsoft/SEAL/blob/82b07db635132e297282649e2ab5908999089ad2/native/src/seal/batchencoder.cpp)
        let row = ring_size >> 1;
        let m = ring_size << 1;
        let gen = 3;
        let mut pos = 1;
        let mut matrix_representation_index_map = vec![0usize; ring_size];
        for i in 0..row {
            let index1 = (pos - 1) >> 1;
            let index2 = (m - pos - 1) >> 1;
            matrix_representation_index_map[i] =
                index1.reverse_bits() >> (ring_size.leading_zeros() + 1);
            matrix_representation_index_map[i | row] =
                index2.reverse_bits() >> (ring_size.leading_zeros() + 1);
            pos *= gen;
            pos &= m - 1;
        }

        let modt_op = M::initialise(t);
        let t_ntt_op = N::init(t, ring_size);

        let modq_ops = q_moduli_chain
            .iter()
            .map(|qi| M::initialise(*qi))
            .collect_vec();

        Self {
            ring_size,
            q_moduli_chain: q_moduli_chain.to_vec(),
            levels_plus_1: q_moduli_chain.len(),
            modq_ops,
            modt_op,
            t_ntt_op,
            q_modt_at_level,
            t_inv_modqi_at_level,
            matrix_representation_index_map,
        }
    }
}

impl<M: ModulusVecBackend<u64> + ModulusBackendConfig<u64>, N: NttConfig<Scalar = u64> + Ntt>
    BfvEncodingDecodingParameters for BfvClientParametersForScalarU64<M, N>
{
    type ModOp = M;
    type NttOp = N;

    fn matrix_representation_index_map(&self) -> &[usize] {
        self.matrix_representation_index_map.as_slice()
    }

    fn modt_op(&self) -> &Self::ModOp {
        &self.modt_op
    }

    fn ring_size(&self) -> usize {
        self.ring_size
    }

    fn t_ntt_op(&self) -> &Self::NttOp {
        &self.t_ntt_op
    }
}

impl<M: ModulusVecBackend<u64> + ModulusBackendConfig<u64>, N: NttConfig<Scalar = u64> + Ntt>
    BfvEncryptionParameters for BfvClientParametersForScalarU64<M, N>
{
    type ModOp = M;

    fn modq_ops_at_level(&self, level: usize) -> &[Self::ModOp] {
        &self.modq_ops[..self.levels_plus_1 - level]
    }

    fn modt_op(&self) -> &Self::ModOp {
        &self.modt_op
    }

    fn q_modt_at_level(&self, level: usize) -> &Self::Scalar {
        &self.q_modt_at_level[level]
    }

    fn q_moduli_chain_at_level(&self, level: usize) -> &[Self::Scalar] {
        &self.q_moduli_chain[..self.levels_plus_1 - level]
    }

    fn ring_size(&self) -> usize {
        self.ring_size
    }

    fn t_inv_modqi_at_level(&self, level: usize) -> &[Self::Scalar] {
        &self.t_inv_modqi_at_level[level]
    }
}
