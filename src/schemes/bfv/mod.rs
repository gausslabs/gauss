use crate::{
    ciphertext::{BfvCiphertext, Ciphertext, InitialiseLevelledCiphertext, Representation},
    core_crypto::{
        matrix::{Matrix, MatrixMut, Row, RowMut},
        modulus::{
            BarrettBackend, ModulusArithmeticBackend, ModulusBackendConfig, ModulusVecBackend,
            MontgomeryBackend, MontgomeryScalar, NativeModulusBackend,
        },
        ntt::{NativeNTTBackend, Ntt, NttConfig},
        random::{DefaultU64SeededRandomGenerator, WithLocal},
    },
    keys::{Encryptor, SecretKey, SecretWithHammingWeight},
    parameters::{
        BfvDecryptionParameters, BfvEncodingDecodingParameters, BfvEncryptionParameters, Parameters,
    },
    utils::{convert::TryConvertFrom, mod_inverse, mod_inverse_big_unit},
};
use itertools::Itertools;
use num_bigint::BigUint;
use num_traits::{FromPrimitive, One, ToPrimitive};
use std::sync::OnceLock;

use self::ops::{secret_key_encryption, simd_encode_message};

pub(crate) mod ops;

static BFV_CLIENT_PARAMETERS_U64: OnceLock<
    BfvClientParametersForScalarU64<NativeModulusBackend, NativeNTTBackend>,
> = OnceLock::new();

trait WithGlobal {
    fn with_global<F, R>(func: F) -> R
    where
        F: Fn(&Self) -> R;
}

pub struct BfvSecretKey {
    values: Vec<i32>,
}

impl SecretKey for BfvSecretKey {
    type Scalar = i32;
    fn values(&self) -> &[Self::Scalar] {
        &self.values
    }
}

impl SecretWithHammingWeight for BfvSecretKey {
    fn new_with_hamming_weight(hw: usize, ring_size: usize) -> Self {
        todo!()
    }
}

impl<'a, P> Encryptor<&'a [u64], BfvCiphertextScalarU64<P>> for BfvSecretKey
where
    P: TryConvertFrom<&'a [i32], Parameters = &'a [u64]> + MatrixMut<MatElement = u64>,
    <P as Matrix>::R: RowMut,
{
    fn encrypt(&self, message: &[u64]) -> BfvCiphertextScalarU64<P> {
        DefaultU64SeededRandomGenerator::with_local(|random_generator| {
            BfvClientParametersForScalarU64::with_global(|parameters| {
                let encoded_message = simd_encode_message(message, parameters);
                secret_key_encryption(self, &encoded_message, parameters, random_generator, 0)
            })
        })
    }
}

pub struct BfvCiphertextScalarU64<P> {
    c_partq: Vec<P>,
    level: usize,
    representation: Representation,
}

impl<P, R> Ciphertext for BfvCiphertextScalarU64<P>
where
    R: Row<Element = u64> + RowMut,
    P: Matrix<MatElement = u64> + MatrixMut<MatElement = u64, R = R>,
{
    type Scalar = u64;
    type Poly = P;
    type Row = R;

    fn representation(&self) -> Representation {
        self.representation
    }
}

impl<P> InitialiseLevelledCiphertext for BfvCiphertextScalarU64<P> {
    type C = Vec<P>;

    fn new(c: Self::C, level: usize, representation: Representation) -> Self {
        Self {
            c_partq: c,
            level,
            representation,
        }
    }
}

impl<P, R> BfvCiphertext for BfvCiphertextScalarU64<P>
where
    R: Row<Element = u64> + RowMut,
    P: Matrix<MatElement = u64> + MatrixMut<MatElement = u64, R = R>,
{
    fn c_partq(&self) -> &[Self::Poly] {
        &self.c_partq
    }

    fn c_partq_mut(&mut self) -> &mut [Self::Poly] {
        &mut self.c_partq
    }

    fn level(&self) -> usize {
        self.level
    }
}

/// Parameters that can only be used for encryption and decryption with Scalar
/// set to u64. Intended to be used only on the client side.
pub struct BfvClientParametersForScalarU64<M, N> {
    ring_size: usize,
    q_moduli_chain: Vec<u64>,
    levels_plus_1: usize,
    modq_ops: Vec<M>,
    modt_op: M,
    t_ntt_op: N,
    basisq_ntt_ops: Vec<N>,

    // Encryption precomputes
    q_modt_at_level: Vec<u64>,
    t_inv_modqi_at_level: Vec<Vec<u64>>,
    matrix_representation_index_map: Vec<usize>,

    // Decryption precomputes
    q_over_qi_inv_modqi_times_t_over_qi_modt_at_level: Vec<Vec<MontgomeryScalar<u64>>>,
    beta_times_q_over_qi_inv_modqi_times_t_over_qi_modt_at_level: Vec<Vec<MontgomeryScalar<u64>>>,
    q_over_qi_inv_modqi_times_t_over_qi_fractional_at_level: Vec<Vec<f64>>,
    beta_times_q_over_qi_inv_modqi_times_t_over_qi_fractional_at_level: Vec<Vec<f64>>,
    log_beta_at_level: Vec<usize>,
}

impl<M: Default, N: Default> Default for BfvClientParametersForScalarU64<M, N> {
    fn default() -> Self {
        Self {
            ring_size: 0,
            q_moduli_chain: vec![],
            levels_plus_1: 0,
            modq_ops: vec![],
            modt_op: M::default(),
            t_ntt_op: N::default(),
            basisq_ntt_ops: vec![],
            q_modt_at_level: vec![],
            t_inv_modqi_at_level: vec![],
            matrix_representation_index_map: vec![],
            q_over_qi_inv_modqi_times_t_over_qi_modt_at_level: vec![],
            beta_times_q_over_qi_inv_modqi_times_t_over_qi_modt_at_level: vec![],
            q_over_qi_inv_modqi_times_t_over_qi_fractional_at_level: vec![],
            beta_times_q_over_qi_inv_modqi_times_t_over_qi_fractional_at_level: vec![],
            log_beta_at_level: vec![],
        }
    }
}

impl<M, N> Parameters for BfvClientParametersForScalarU64<M, N> {
    type Scalar = u64;
}

impl<
        M: ModulusVecBackend<u64> + ModulusBackendConfig<u64> + MontgomeryBackend<u64, u128>,
        N: NttConfig<Scalar = u64> + Ntt + NttConfig,
    > BfvClientParametersForScalarU64<M, N>
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

        // Precomputes for encryption, decrption, encoding, decoding
        let mut q_clone = q.clone();
        let mut q_modt_at_level = vec![];
        let mut t_inv_modqi_at_level = vec![];
        let big_t = BigUint::from_u64(t).unwrap();
        for i in 0..levels + 1 {
            q_modt_at_level.push((&q_clone % t).to_u64().unwrap());

            // t^{-1} \mod Q
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

        // Precomputes for Simple Scale and Round procedure
        let mut big_q = q.clone();
        let mut q_over_qi_inv_modqi_times_t_over_qi_modt_at_level = vec![];
        let mut beta_times_q_over_qi_inv_modqi_times_t_over_qi_modt_at_level = vec![];
        let mut q_over_qi_inv_modqi_times_t_over_qi_fractional_at_level = vec![];
        let mut beta_times_q_over_qi_inv_modqi_times_t_over_qi_fractional_at_level = vec![];
        let mut log_beta_at_level = vec![];
        for i in 0..levels + 1 {
            let q_chain_at_level_i = &q_moduli_chain[..(levels + 1 - i)];

            let log_beta =
                ((64 - q_chain_at_level_i.iter().max().unwrap().leading_zeros()) / 2) as usize;
            let beta = 1u64 << log_beta;
            let mut q_over_qi_inv_modqi_times_t_over_qi_modt_vec = vec![];
            let mut beta_times_q_over_qi_inv_modqi_times_t_over_qi_modt_vec = vec![];
            let mut q_over_qi_inv_modqi_times_t_over_qi_fractional_vec = vec![];
            let mut beta_times_q_over_qi_inv_modqi_times_t_over_qi_fractional_vec = vec![];
            q_chain_at_level_i.iter().for_each(|qi| {
                let q_over_qi_inv_mod_qi = mod_inverse(((&big_q / qi) % qi).to_u64().unwrap(), *qi);

                let q_over_qi_inv_mod_qi_times_t = BigUint::from(t) * q_over_qi_inv_mod_qi;

                // v_i = ((q/q_i)^{-1}_q_i * t)
                // rational part: v_i / q_i \mod{t}
                q_over_qi_inv_modqi_times_t_over_qi_modt_vec.push(modt_op.normal_to_mont_space(
                    ((&q_over_qi_inv_mod_qi_times_t / qi) % t).to_u64().unwrap(),
                ));
                // fractional part: (v_i % q_i) / q_i
                q_over_qi_inv_modqi_times_t_over_qi_fractional_vec
                    .push(((&q_over_qi_inv_mod_qi_times_t % qi).to_f64().unwrap()) / (*qi as f64));

                // \beta * v_i = (\beta * (q/q_i)^{-1}_q_i * t) / q_i \mod{t}
                // rational part: \beta * v_i / q_i \mod{t}
                beta_times_q_over_qi_inv_modqi_times_t_over_qi_modt_vec.push(
                    modt_op.normal_to_mont_space(
                        (((beta * &q_over_qi_inv_mod_qi_times_t) / qi) % t)
                            .to_u64()
                            .unwrap(),
                    ),
                );
                // fractional part: ((\beta * v_i) % q_i) / q_i
                beta_times_q_over_qi_inv_modqi_times_t_over_qi_fractional_vec.push(
                    ((beta * &q_over_qi_inv_mod_qi_times_t) % qi)
                        .to_f64()
                        .unwrap()
                        / (*qi as f64),
                );
            });

            // insert
            q_over_qi_inv_modqi_times_t_over_qi_modt_at_level
                .push(q_over_qi_inv_modqi_times_t_over_qi_modt_vec);
            beta_times_q_over_qi_inv_modqi_times_t_over_qi_modt_at_level
                .push(beta_times_q_over_qi_inv_modqi_times_t_over_qi_modt_vec);
            q_over_qi_inv_modqi_times_t_over_qi_fractional_at_level
                .push(q_over_qi_inv_modqi_times_t_over_qi_fractional_vec);
            beta_times_q_over_qi_inv_modqi_times_t_over_qi_fractional_at_level
                .push(beta_times_q_over_qi_inv_modqi_times_t_over_qi_fractional_vec);
            log_beta_at_level.push(log_beta);

            big_q /= q_moduli_chain[levels - i];
        }

        let basisq_ntt_ops = q_moduli_chain
            .iter()
            .map(|qi| N::init(*qi, ring_size))
            .collect_vec();

        Self {
            ring_size,
            q_moduli_chain: q_moduli_chain.to_vec(),
            levels_plus_1: q_moduli_chain.len(),
            basisq_ntt_ops,
            modq_ops,
            modt_op,
            t_ntt_op,
            q_modt_at_level,
            t_inv_modqi_at_level,
            matrix_representation_index_map,
            q_over_qi_inv_modqi_times_t_over_qi_fractional_at_level,
            beta_times_q_over_qi_inv_modqi_times_t_over_qi_fractional_at_level,
            q_over_qi_inv_modqi_times_t_over_qi_modt_at_level,
            beta_times_q_over_qi_inv_modqi_times_t_over_qi_modt_at_level,
            log_beta_at_level,
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
    type NttOp = N;

    fn modq_ops_at_level(&self, level: usize) -> &[Self::ModOp] {
        &self.modq_ops[..self.levels_plus_1 - level]
    }

    fn basisq_ntt_ops_at_level(&self, level: usize) -> &[Self::NttOp] {
        &self.basisq_ntt_ops[..self.levels_plus_1 - level]
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

impl<
        M: ModulusVecBackend<u64>
            + BarrettBackend<u64, u128>
            + MontgomeryBackend<u64, u128>
            + ModulusArithmeticBackend<MontgomeryScalar<u64>>
            + ModulusArithmeticBackend<u64>,
        N: NttConfig<Scalar = u64> + Ntt,
    > BfvDecryptionParameters for BfvClientParametersForScalarU64<M, N>
{
    type ModOps = M;
    type NttOp = N;

    fn basisq_ntt_ops_at_level(&self, level: usize) -> &[Self::NttOp] {
        &self.basisq_ntt_ops[..self.levels_plus_1 - level]
    }

    fn modq_ops_at_level(&self, level: usize) -> &[Self::ModOps] {
        &self.modq_ops[..self.levels_plus_1 - level]
    }

    fn modt_op(&self) -> &Self::ModOps {
        &self.modt_op
    }

    fn q_moduli_chain_at_level(&self, level: usize) -> &[Self::Scalar] {
        &self.q_moduli_chain[..self.levels_plus_1 - level]
    }

    fn ring_size(&self) -> usize {
        self.ring_size
    }

    fn beta_times_q_over_qi_inv_modqi_times_t_over_qi_fractional_at_level(
        &self,
        level: usize,
    ) -> &[f64] {
        &self.beta_times_q_over_qi_inv_modqi_times_t_over_qi_fractional_at_level[level]
    }

    fn q_over_qi_inv_modqi_times_t_over_qi_fractional_at_level(&self, level: usize) -> &[f64] {
        &self.q_over_qi_inv_modqi_times_t_over_qi_fractional_at_level[level]
    }

    fn beta_times_q_over_qi_inv_modqi_times_t_over_qi_modt_at_level(
        &self,
        level: usize,
    ) -> &[MontgomeryScalar<Self::Scalar>] {
        &self.beta_times_q_over_qi_inv_modqi_times_t_over_qi_modt_at_level[level]
    }

    fn q_over_qi_inv_modqi_times_t_over_qi_modt_at_level(
        &self,
        level: usize,
    ) -> &[MontgomeryScalar<Self::Scalar>] {
        &self.q_over_qi_inv_modqi_times_t_over_qi_modt_at_level[level]
    }

    fn log_beta(&self) -> usize {
        self.log_beta_at_level[0]
    }
}

impl WithGlobal for BfvClientParametersForScalarU64<NativeModulusBackend, NativeNTTBackend> {
    fn with_global<F, R>(func: F) -> R
    where
        Self: Sized,
        F: FnOnce(&Self) -> R,
    {
        func(BFV_CLIENT_PARAMETERS_U64.get().unwrap())
    }
}
