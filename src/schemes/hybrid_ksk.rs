// key gen
// key switch

use std::thread::panicking;

use itertools::{izip, Itertools};
use num_bigint::BigUint;

use crate::{
    core_crypto::{
        matrix::{Matrix, MatrixMut, RowMut},
        modulus::ModulusVecBackend,
        random::{InitWithSeed, RandomGaussianDist, RandomSeed, RandomUniformDist},
    },
    parameters::Parameters,
};

pub trait HybridKskKeyGenParameters: Parameters {
    type ModOp: ModulusVecBackend<Self::Scalar>;

    fn dnum(&self) -> usize;
    fn ring_size(&self) -> usize;
    fn specialp_moduli_chain(&self) -> &[Self::Scalar];
    fn primes_at_level(&self) -> usize;

    fn modq_ops_at_level(&self, level: usize) -> &[Self::ModOp];
    fn modspecialp_ops_at_level(&self, level: usize) -> &[Self::ModOp];
    fn q_moduli_chain_at_level(&self, level: usize) -> &[Self::Scalar];
    fn specialp_moduli_chain_at_level(&self, level: usize) -> &[Self::Scalar];
    fn big_specialp(&self, level: usize) -> BigUint;
    fn big_qjs_at_level(&self, level: usize) -> &[BigUint];
    fn gammak_modqi_at_level(&self, level: usize) -> &[Vec<Self::Scalar>];
    fn gammak_modspecialpj_at_level(&self, level: usize) -> &[Vec<Self::Scalar>];
    fn alpha_at_level(&self, level: usize) -> usize;
}

fn generate_key<
    P: HybridKskKeyGenParameters,
    M: Matrix<MatElement = P::Scalar> + MatrixMut + Clone,
    R: RandomSeed + RandomGaussianDist<M>,
    NR: InitWithSeed<Seed = R::Seed> + RandomUniformDist<M>,
>(
    params: P,
    level: usize,
    p: M,
    rng: &mut R,
) where
    <M as Matrix>::R: RowMut,
{
    debug_assert!(p.dimension() == (params.primes_at_level(), params.ring_size()));

    let dnum = params.dnum();

    let alpha = (dnum as f64 / (level as f64 + 1.0)).ceil() as usize;
    let q_moduli_chain = params.q_moduli_chain_at_level(level);
    let specialp_moduli_chain = params.specialp_moduli_chain();

    let big_specialp = params.big_specialp(level);
    let ring_size = params.ring_size();

    // q_moduli_chain.iter().chunks(dnum)
    let seed = rng.random_seed();

    let gammak_modqi_at_level = params.gammak_modqi_at_level(level);
    let gammak_modspecialpj_at_level = params.gammak_modspecialpj_at_level(level);
    let modq_ops = params.modq_ops_at_level(level);
    let modspecialp_ops = params.modspecialp_ops_at_level(level);
    let mut prng = NR::init_with_seed(seed);
    for k in 0..alpha {
        let mut p_clone = p.clone();

        // part Q
        // let gammak = &gammak_modqi_at_level[k];
        // let mut ak_partq =
        //     RandomUniformDist::random_ring_poly(&mut prng, q_moduli_chain,
        // ring_size); let mut ek_partq =
        //     RandomGaussianDist::random_ring_poly(&mut rng, q_moduli_chain,
        // ring_size);

        // izip!(modq_ops.iter(), gammak.iter(),
        // p_clone.iter_rows_mut()).for_each(     |(modqi, gki, ak_modqi)| {
        //         modqi.scalar_mul_mod_vec(ak_modqi.as_mut(), *gki);
        //     },
        // );

        // part specialP
    }
}
fn keyswitch() {}
