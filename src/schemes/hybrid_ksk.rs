// key gen
// key switch

use std::thread::panicking;

use itertools::{izip, Itertools};
use num_bigint::BigUint;

use crate::{
    core_crypto::{
        matrix::{Matrix, MatrixMut, RowMut},
        modulus::{ModulusArithmeticBackend, ModulusVecBackend},
        ntt::Ntt,
        random::{InitWithSeed, RandomGaussianDist, RandomSeed, RandomUniformDist},
        ring::{add_lazy_mut, foward_lazy, mul_lazy_mut},
    },
    keys::SecretKey,
    parameters::Parameters,
    utils::convert::TryConvertFrom,
};

pub trait HybridKskKeyGenParameters: Parameters {
    type ModOp: ModulusVecBackend<Self::Scalar>;
    type NttOp: Ntt<Scalar = Self::Scalar>;

    fn dnum(&self) -> usize;
    fn ring_size(&self) -> usize;
    fn specialp_moduli_chain(&self) -> &[Self::Scalar];
    fn primes_at_level(&self) -> usize;

    fn q_modops_at_level(&self, level: usize) -> &[Self::ModOp];
    fn specialp_modops_at_level(&self, level: usize) -> &[Self::ModOp];

    fn q_nttops_at_level(&self, level: usize) -> &[Self::NttOp];
    fn specialp_nttops_at_level(&self, level: usize) -> &[Self::NttOp];

    fn q_moduli_chain_at_level(&self, level: usize) -> &[Self::Scalar];
    fn specialp_moduli_chain_at_level(&self, level: usize) -> &[Self::Scalar];

    fn big_specialp(&self, level: usize) -> BigUint;
    fn big_qjs_at_level(&self, level: usize) -> &[BigUint];

    fn gammak_modqi_at_level(&self, level: usize) -> &[Vec<Self::Scalar>];
    fn alpha_at_level(&self, level: usize) -> usize;
}

fn generate_key<
    P: HybridKskKeyGenParameters,
    M: Matrix<MatElement = P::Scalar> + MatrixMut,
    R: RandomSeed + RandomGaussianDist<[P::Scalar], Parameters = P::Scalar>,
    NR: InitWithSeed<Seed = R::Seed> + RandomUniformDist<[P::Scalar], Parameters = P::Scalar>,
    S: SecretKey<Scalar = i32>,
>(
    params: P,
    level: usize,
    p: M,
    secret: &S,
    ksk_polys: &mut [M],
    rng: &mut R,
) where
    <M as Matrix>::R: RowMut + Clone,
    M: TryConvertFrom<[i32], Parameters = [P::Scalar]>,
{
    let dnum = params.dnum();
    let alpha = (dnum as f64 / (level as f64 + 1.0)).ceil() as usize;

    debug_assert!(p.dimension() == (params.primes_at_level(), params.ring_size()));
    debug_assert!(
        ksk_polys.len() == alpha * 2,
        "KSK polynomials supplied {} != expected polynomials alpha*2 {}",
        ksk_polys.len(),
        alpha * 2
    );

    let q_moduli_chain = params.q_moduli_chain_at_level(level);
    let specialp_moduli_chain = params.specialp_moduli_chain();

    let ring_size = params.ring_size();

    // q_moduli_chain.iter().chunks(dnum)
    let seed = rng.random_seed();

    let gammak_modqi_at_level = params.gammak_modqi_at_level(level);
    let q_modops = params.q_modops_at_level(level);
    let specialp_modops = params.specialp_modops_at_level(level);
    let mut prng = NR::init_with_seed(seed);

    let q_nttops = params.q_nttops_at_level(level);
    let specialp_nttops = params.specialp_nttops_at_level(level);

    let mut s_partq = M::try_convert_from(secret.values(), q_moduli_chain);
    let mut s_partspecialp = M::try_convert_from(secret.values(), specialp_moduli_chain);
    foward_lazy(&mut s_partq, q_nttops);
    foward_lazy(&mut s_partspecialp, specialp_nttops);

    let q_size = q_moduli_chain.len();
    let specialp_size = specialp_moduli_chain.len();

    for (k, ck) in izip!(0..alpha, ksk_polys.chunks_exact_mut(2)) {
        let gammak = &gammak_modqi_at_level[k];

        let (c0_k, c1_k) = ck.split_at_mut(1);
        let c0_k = &mut c0_k[0];
        let c1_k = &mut c1_k[0];

        debug_assert!(c0_k.dimension() == (q_size + specialp_size, ring_size));
        debug_assert!(c1_k.dimension() == (q_size + specialp_size, ring_size));

        // part Q
        izip!(
            gammak.iter(),
            q_moduli_chain.iter(),
            q_modops.iter(),
            q_nttops.iter(),
            p.iter_rows(),
            c0_k.iter_rows_mut().take(q_size),
            c1_k.iter_rows_mut().take(q_size),
            s_partq.iter_rows()
        )
        .for_each(
            |(gammak_modqi, qi, qi_modop, qi_nttop, p_modqi, c0k_modqi, c1k_modqi, s_modqi)| {
                // Assume c_{1,k} is sampled in evaluation form
                RandomUniformDist::random_fill(&mut prng, qi, c1k_modqi.as_mut());

                // c1_k * s \mod qi
                let mut c1k_s = c1k_modqi.clone();
                qi_modop.mul_lazy_mod_vec(c1k_s.as_mut(), s_modqi.as_ref());

                // e \mod qi
                RandomGaussianDist::random_fill(rng, qi, c0k_modqi.as_mut());
                qi_nttop.forward_lazy(c0k_modqi.as_mut());

                // e + c1_k * s \mod qi
                qi_modop.add_lazy_mod_vec(c0k_modqi.as_mut(), c1k_s.as_ref());

                // gamma_k * p \mod qi
                let mut p_modqi_clone = p_modqi.clone();
                qi_modop.scalar_mul_mod_vec(p_modqi_clone.as_mut(), *gammak_modqi);

                // e + c1_k * s + gamma_k * p \mod qi
                qi_modop.add_lazy_mod_vec(c0k_modqi.as_mut(), p_modqi_clone.as_ref());
            },
        );

        // part specialP
        izip!(
            specialp_moduli_chain.iter(),
            specialp_modops.iter(),
            specialp_nttops.iter(),
            c0_k.iter_rows_mut().skip(q_size),
            c1_k.iter_rows_mut().skip(q_size),
            s_partspecialp.iter_rows()
        )
        .for_each(|(pj, pj_modop, pj_nttop, c0k_modpj, c1k_modpj, s_modpj)| {
            // Assume c_{1, k} is sampled in evaluation form
            RandomUniformDist::random_fill(&mut prng, pj, c1k_modpj.as_mut());

            // c1_k * s \mod pj
            let mut c1k_s = c1k_modpj.clone();
            pj_modop.mul_lazy_mod_vec(c1k_s.as_mut(), s_modpj.as_ref());

            // e \mod pj
            RandomGaussianDist::random_fill(rng, pj, c0k_modpj.as_mut());
            pj_nttop.forward_lazy(c0k_modpj.as_mut());

            // e + c1_k * s \mod pj
            pj_modop.add_lazy_mod_vec(c0k_modpj.as_mut(), c1k_s.as_ref());
        });
    }
}
fn keyswitch() {}
