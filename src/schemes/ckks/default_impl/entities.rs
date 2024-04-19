use std::{clone, iter::Once, sync::OnceLock};

use itertools::Itertools;
use num_bigint::BigUint;
use num_complex::Complex;
use num_traits::{zero, Zero};
use rand::SeedableRng;
use rand_chacha::{ChaCha8Core, ChaCha8Rng};

use crate::{
    ciphertext::{Ciphertext, Representation, RlweCiphertext, SeededCiphertext},
    core_crypto::{
        matrix::{Matrix, MatrixEntity, MatrixMut, RowMut},
        modulus::{ModulusBackendConfig, ModulusVecBackend, NativeModulusBackend},
        ntt::{NativeNTTBackend, Ntt, NttConfig},
        num::{big_float::BigFloat, BFloat, ComplexNumber},
        prime::generate_primes_vec,
        random::{DefaultU64SeededRandomGenerator, WithLocal},
        ring::{backward, foward_lazy},
    },
    keys::{Decryptor, Encryptor, LevelDecoder, LevelEncoder, SecretKey},
    parameters::{CkksEncDecParameters, Parameters},
    schemes::{
        ckks::ops::{
            secret_key_decryption, secret_key_encryption, simd_decode, simd_encode,
            ScaledCkksCiphertext,
        },
        ops::generate_ternery_secret_with_hamming_weight,
        WithGlobal,
    },
    utils::{convert::TryConvertFrom, moduli_chain_to_biguint, psi_powers},
};

use super::CkksCiphertext;

type DefaultBigFloat = crate::core_crypto::num::big_float::BigFloat;
type DefaultComplex = Complex<DefaultBigFloat>;

pub static NATIVE_CKKS_CLIENT_PARAMETERS_U64: OnceLock<
    CkksClientParametersU64<
        NativeModulusBackend,
        NativeNTTBackend,
        DefaultComplex,
        DefaultBigFloat,
    >,
> = OnceLock::new();

impl<M: MatrixMut<MatElement = u64> + MatrixEntity> LevelEncoder<M> for Vec<DefaultComplex>
where
    <M as Matrix>::R: RowMut,
{
    fn encode(&self, level: usize) -> M {
        CkksClientParametersU64::with_global(|params| {
            // encode
            let mut m_poly = M::zeros(params.q_moduli_chain_len, params.ring_size());
            simd_encode(&mut m_poly, self, params, level, params.delta());
            m_poly
        })
    }
}

impl<M: MatrixMut<MatElement = u64>> LevelDecoder<Vec<DefaultComplex>> for M
where
    <M as Matrix>::R: RowMut,
{
    fn decode(&self, level: usize) -> Vec<DefaultComplex> {
        CkksClientParametersU64::with_global(|params| {
            // decode
            let mut m_out = vec![DefaultComplex::zero(); params.ring_size() >> 1];
            simd_decode(self, params, level, params.delta(), &mut m_out);
            m_out
        })
    }
}

pub struct CkksSecretKey {
    values: Vec<i32>,
}

impl SecretKey for CkksSecretKey {
    type Scalar = i32;
    fn values(&self) -> &[Self::Scalar] {
        &self.values
    }
}

impl CkksSecretKey {
    pub fn new() -> CkksSecretKey {
        CkksClientParametersU64::with_global(|parameters| {
            DefaultU64SeededRandomGenerator::with_local_mut(|rng| {
                let values = generate_ternery_secret_with_hamming_weight(
                    rng,
                    parameters.secret_hw,
                    parameters.ring_size,
                );
                CkksSecretKey { values }
            })
        })
    }
}

impl<
        M: MatrixMut<MatElement = u64>
            + MatrixEntity
            + Clone
            + TryConvertFrom<[i32], Parameters = [u64]>,
    > Encryptor<[DefaultComplex], CkksCiphertextGenericStorage<M>> for CkksSecretKey
where
    <M as Matrix>::R: RowMut,
{
    fn encrypt(&self, message: &[DefaultComplex]) -> CkksCiphertextGenericStorage<M> {
        CkksClientParametersU64::with_global(|params| {
            DefaultU64SeededRandomGenerator::with_local_mut(|rng| {
                let delta = params.delta();

                // encode
                let mut m_poly = M::zeros(params.q_moduli_chain_len, params.ring_size());
                simd_encode(&mut m_poly, message, params, 0, delta);

                // `encrypt` function wants m_poly in Evaluation representation.
                foward_lazy(&mut m_poly, params.q_nttops_at_level(0));

                // encrypt
                let mut c_out = CkksCiphertextGenericStorage {
                    c: vec![M::zeros(params.q_moduli_chain_len, params.ring_size()); 2],
                    level: 0,
                    is_lazy: false,
                    seed: <ChaCha8Rng as SeedableRng>::Seed::default(),
                    representation: Representation::Evaluation,
                    scale: delta.clone(),
                };
                secret_key_encryption(&mut c_out, &m_poly, self, params, rng, 0);

                c_out
            })
        })
    }
}

impl<
        M: MatrixMut<MatElement = u64>
            + Clone
            + MatrixEntity
            + TryConvertFrom<[i32], Parameters = [u64]>,
    > Decryptor<Vec<DefaultComplex>, CkksCiphertextGenericStorage<M>> for CkksSecretKey
where
    <M as Matrix>::R: RowMut,
{
    fn decrypt(&self, c: &CkksCiphertextGenericStorage<M>) -> Vec<DefaultComplex> {
        CkksClientParametersU64::with_global(|params| {
            // decrypt
            let mut m_poly = M::zeros(
                params.q_moduli_chain_at_level(c.level()).len(),
                params.ring_size(),
            );
            secret_key_decryption(c, &mut m_poly, self, params);

            // decrypt outputs m_poly in Evaluation representation. Map it back to
            // Coefficient representation.
            backward(&mut m_poly, params.q_nttops_at_level(c.level()));

            // decode
            let mut m_out = vec![DefaultComplex::zero(); params.ring_size() >> 1];
            simd_decode(&m_poly, params, c.level(), c.scale(), &mut m_out);
            m_out
        })
    }
}

pub struct CkksClientParametersU64<ModOp, NttOp, Comp, F> {
    q_moduli_chain: Vec<u64>,
    ring_size: usize,
    delta: F,
    bigq_at_level: Vec<BigUint>,
    q_moduli_chain_len: usize,
    secret_hw: usize,

    psi_powers: Vec<Comp>,
    rot_group: Vec<usize>,

    q_modops: Vec<ModOp>,
    q_nttops: Vec<NttOp>,
}

impl<
        ModOp: ModulusBackendConfig<u64>,
        NttOp: NttConfig<Scalar = u64>,
        Comp: ComplexNumber<F>,
        F: BFloat,
    > CkksClientParametersU64<ModOp, NttOp, Comp, F>
{
    fn new(ring_size: usize, q_moduli_chain: Vec<u64>, delta: F) -> Self {
        let m = ring_size << 1;
        let n = ring_size;
        let l = n >> 1;

        let mut a = 1usize;
        let mut rot_group = vec![];
        for _ in 0..l {
            rot_group.push(a);
            a = (a * 5) % m;
        }

        let psi_powers = psi_powers(m as u32);
        let mut bigq_at_level = vec![];
        for l in 0..q_moduli_chain.len() {
            bigq_at_level.push(moduli_chain_to_biguint(
                &q_moduli_chain[..q_moduli_chain.len() - l],
            ));
        }

        let q_modops = q_moduli_chain
            .iter()
            .map(|qi| ModOp::initialise(*qi))
            .collect_vec();
        let q_nttops = q_moduli_chain
            .iter()
            .map(|qi| NttOp::init(*qi, ring_size))
            .collect_vec();

        CkksClientParametersU64 {
            delta,
            psi_powers,
            rot_group,
            ring_size,
            q_moduli_chain_len: q_moduli_chain.len(),
            q_moduli_chain,
            bigq_at_level,
            q_modops,
            q_nttops,
            secret_hw: ring_size >> 1,
        }
    }
}

impl<M, N, C, F> Parameters for CkksClientParametersU64<M, N, C, F> {
    type Scalar = u64;
}

impl<
        ModOp: ModulusVecBackend<u64>,
        NttOp: Ntt<Scalar = u64>,
        Comp: ComplexNumber<F>,
        F: BFloat,
    > CkksEncDecParameters for CkksClientParametersU64<ModOp, NttOp, Comp, F>
{
    type BU = BigUint;
    type Complex = Comp;
    type F = F;
    type ModOp = ModOp;
    type NttOp = NttOp;

    fn bigq_at_level(&self, level: usize) -> &Self::BU {
        &self.bigq_at_level[level]
    }
    fn delta(&self) -> &Self::F {
        &self.delta
    }
    fn psi_powers(&self) -> &[Self::Complex] {
        &self.psi_powers
    }
    fn ring_size(&self) -> usize {
        self.ring_size
    }
    fn rot_group(&self) -> &[usize] {
        &self.rot_group
    }

    fn q_moduli_chain_at_level(&self, level: usize) -> &[Self::Scalar] {
        &self.q_moduli_chain[..self.q_moduli_chain_len - level]
    }
    fn q_modops_at_level(&self, level: usize) -> &[Self::ModOp] {
        &self.q_modops[..self.q_moduli_chain_len - level]
    }
    fn q_nttops_at_level(&self, level: usize) -> &[Self::NttOp] {
        &self.q_nttops[..self.q_moduli_chain_len - level]
    }
}

impl WithGlobal
    for CkksClientParametersU64<NativeModulusBackend, NativeNTTBackend, Complex<BigFloat>, BigFloat>
{
    fn with_global<F, R>(func: F) -> R
    where
        F: Fn(&Self) -> R,
    {
        func(NATIVE_CKKS_CLIENT_PARAMETERS_U64.get().unwrap())
    }
}

pub struct CkksCiphertextGenericStorage<M> {
    c: Vec<M>,
    level: usize,
    is_lazy: bool,
    seed: <ChaCha8Rng as SeedableRng>::Seed,
    representation: Representation,
    scale: DefaultBigFloat,
}

impl<M> Ciphertext for CkksCiphertextGenericStorage<M>
where
    M: MatrixMut<MatElement = u64> + MatrixEntity,
    <M as Matrix>::R: RowMut,
{
    type Poly = M;
    type Scalar = u64;
    type Row = M::R;

    fn representation(&self) -> crate::ciphertext::Representation {
        self.representation
    }

    fn representation_mut(&mut self) -> &mut crate::ciphertext::Representation {
        &mut self.representation
    }
}

impl<M> SeededCiphertext for CkksCiphertextGenericStorage<M> {
    type Seed = <ChaCha8Rng as SeedableRng>::Seed;

    fn seed(&self) -> Self::Seed {
        self.seed
    }
    fn seed_mut(&mut self) -> &mut Self::Seed {
        &mut self.seed
    }
}

impl<M> RlweCiphertext for CkksCiphertextGenericStorage<M>
where
    M: MatrixMut<MatElement = u64> + MatrixEntity,
    <M as Matrix>::R: RowMut,
{
    fn c_partq(&self) -> &[Self::Poly] {
        &self.c
    }

    fn c_partq_mut(&mut self) -> &mut [Self::Poly] {
        &mut self.c
    }

    fn level(&self) -> usize {
        self.level
    }

    fn level_mut(&mut self) -> &mut usize {
        &mut self.level
    }

    fn is_lazy(&self) -> bool {
        self.is_lazy
    }

    fn is_lazy_mut(&mut self) -> &mut bool {
        &mut self.is_lazy
    }
}

impl<M: MatrixMut<MatElement = u64>> ScaledCkksCiphertext for CkksCiphertextGenericStorage<M>
where
    <M as Matrix>::R: RowMut,
{
    type F = DefaultBigFloat;
    fn scale(&self) -> &Self::F {
        &self.scale
    }

    fn scale_mut(&mut self) -> &mut Self::F {
        &mut self.scale
    }
}

pub fn build_parameters(q_moduli_chain_sizes: &[usize], delta: DefaultBigFloat, ring_size: usize) {
    let q_moduli_chain = generate_primes_vec(q_moduli_chain_sizes, ring_size, &[]);

    // Client side parameters
    let parameters = CkksClientParametersU64::<
        NativeModulusBackend,
        NativeNTTBackend,
        DefaultComplex,
        DefaultBigFloat,
    >::new(ring_size, q_moduli_chain, delta);

    NATIVE_CKKS_CLIENT_PARAMETERS_U64.get_or_init(|| parameters);
}
