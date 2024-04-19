use itertools::{izip, Itertools};
use num_traits::{AsPrimitive, PrimInt, Unsigned};

use crate::{
    core_crypto::{
        matrix::{Matrix, MatrixEntity, MatrixMut, RowMut},
        modulus::{BarrettBackend, ModulusVecBackend, MontgomeryBackend, MontgomeryScalar},
        ntt::Ntt,
        num::UnsignedInteger,
        random::{InitWithSeed, RandomGaussianDist, RandomUniformDist},
        ring::{approximate_mod_down, backward, backward_lazy, foward_lazy, reduce_from_lazy_mut},
    },
    keys::SecretKey,
    parameters::Parameters,
    utils::convert::TryConvertFrom,
};

pub trait HybridKskKeyGenParameters: Parameters {
    type ModOp: ModulusVecBackend<Self::Scalar>;
    type NttOp: Ntt<Scalar = Self::Scalar>;

    fn dnum(&self) -> usize;
    fn alpha_at_level(&self, level: usize) -> usize;
    fn ring_size(&self) -> usize;

    fn q_modops_at_level(&self, level: usize) -> &[Self::ModOp];
    fn specialp_modops_at_level(&self, level: usize) -> &[Self::ModOp];

    fn q_nttops_at_level(&self, level: usize) -> &[Self::NttOp];
    fn specialp_nttops_at_level(&self, level: usize) -> &[Self::NttOp];

    fn q_moduli_chain_at_level(&self, level: usize) -> &[Self::Scalar];
    fn specialp_moduli_chain_at_level(&self) -> &[Self::Scalar];

    fn gammak_modqi_at_level(&self, level: usize) -> &[Vec<Self::Scalar>];
}

pub trait HybridKskRuntimeParameters: Parameters
where
    // TODO(Jay): [Remove]`num_traits::AsPrimitive<u128>
    // + num_traits::PrimInt,` are a consquenece of
    // BarrettBackend trait. Check issue #12
    Self::Scalar: num_traits::AsPrimitive<u128> + num_traits::PrimInt,
    u128: AsPrimitive<Self::Scalar>,
{
    type ModOp: ModulusVecBackend<Self::Scalar>
        + BarrettBackend<Self::Scalar, u128>
        + MontgomeryBackend<Self::Scalar, u128>;
    type NttOp: Ntt<Scalar = Self::Scalar>;

    fn dnum(&self) -> usize;
    fn ring_size(&self) -> usize;

    fn q_modops_at_level(&self, level: usize) -> &[Self::ModOp];
    fn specialp_modops_at_level(&self, level: usize) -> &[Self::ModOp];

    fn q_nttops_at_level(&self, level: usize) -> &[Self::NttOp];
    fn specialp_nttops_at_level(&self, level: usize) -> &[Self::NttOp];

    fn q_moduli_chain_at_level(&self, level: usize) -> &[Self::Scalar];
    fn specialp_moduli_chain_at_level(&self, level: usize) -> &[Self::Scalar];

    fn specialp_over_specialpj_inv_modspecialpj_at_level(&self, level: usize) -> &[Self::Scalar];
    fn specialp_over_specialpj_per_modqi_at_level(
        &self,
        level: usize,
    ) -> &[Vec<MontgomeryScalar<Self::Scalar>>];
    fn specialp_inv_modqi_at_level(&self, level: usize) -> &[Self::Scalar];

    fn qj_over_qji_inv_modqji_level(&self, level: usize) -> &[Vec<Self::Scalar>];
    fn qj_over_qji_per_modqi_at_level(
        &self,
        level: usize,
    ) -> &[Vec<Vec<MontgomeryScalar<Self::Scalar>>>];
    fn qj_over_qji_per_modspecialpj_at_level(
        &self,
        level: usize,
    ) -> &[Vec<Vec<MontgomeryScalar<Self::Scalar>>>];
}

/// Input polynomial p must be in Evaluation form
pub fn generate_key<
    P: HybridKskKeyGenParameters,
    M: Matrix<MatElement = P::Scalar> + MatrixMut,
    Se: Clone,
    R: RandomUniformDist<Se, Parameters = u8>
        + RandomUniformDist<M, Parameters = [P::Scalar]>
        + RandomGaussianDist<M, Parameters = [P::Scalar]>
        + InitWithSeed<Seed = Se>,
    S: SecretKey<Scalar = i32>,
>(
    params: &P,
    level: usize,
    p_eval: &M,
    secret: &S,
    ksk_polys_out: &mut [M],
    seed_out: &mut Se,
    rng: &mut R,
) where
    <M as Matrix>::R: RowMut + Clone,
    M: TryConvertFrom<[i32], Parameters = [P::Scalar]>,
{
    let alpha = params.alpha_at_level(level);

    debug_assert!(
        p_eval.dimension()
            == (
                params.q_moduli_chain_at_level(level).len(),
                params.ring_size()
            )
    );
    debug_assert!(
        ksk_polys_out.len() == alpha * 2,
        "KSK polynomials supplied: {} != expected polynomials alpha*2: {}",
        ksk_polys_out.len(),
        alpha * 2
    );

    let q_moduli_chain = params.q_moduli_chain_at_level(level);
    let specialp_moduli_chain = params.specialp_moduli_chain_at_level();

    // TODO(Jay): q and spepcialp as single vector is only required to sample
    // pseudorandom part & error of ksks in QP_s. We can remove once we have a way
    // to pass a submatrix as mutable matrix
    let mut q_and_specuialp_moduli_chain = q_moduli_chain.to_vec();
    q_and_specuialp_moduli_chain.extend(specialp_moduli_chain.iter());

    let ring_size = params.ring_size();

    // let seed = rng.
    RandomUniformDist::random_fill(rng, &0u8, seed_out);
    let mut prng = R::init_with_seed(seed_out.clone());

    let gammak_modqi_at_level = params.gammak_modqi_at_level(level);

    let q_modops = params.q_modops_at_level(level);
    let specialp_modops = params.specialp_modops_at_level(level);

    let q_nttops = params.q_nttops_at_level(level);
    let specialp_nttops = params.specialp_nttops_at_level(level);

    let mut s_partq = M::try_convert_from(secret.values(), q_moduli_chain);
    let mut s_partspecialp = M::try_convert_from(secret.values(), specialp_moduli_chain);
    foward_lazy(&mut s_partq, q_nttops);
    foward_lazy(&mut s_partspecialp, specialp_nttops);

    let q_size = q_moduli_chain.len();
    let specialp_size = specialp_moduli_chain.len();

    for (k, ck) in izip!(0..alpha, ksk_polys_out.chunks_exact_mut(2)) {
        let gammak = &gammak_modqi_at_level[k];

        let (c0_k, c1_k) = ck.split_at_mut(1);
        let c0_k = &mut c0_k[0];
        let c1_k = &mut c1_k[0];

        debug_assert!(c0_k.dimension() == (q_size + specialp_size, ring_size));
        debug_assert!(c1_k.dimension() == (q_size + specialp_size, ring_size));

        // c_{1,k} is sampled in coefficient form. Since evaluaton form is Ntt
        // dependent, one cannot sample a seeded polynomial directly in Evaluation form.
        // We store ciphertext tuple as (c0, c1) where c1 = -a (a being the pseudo
        // random part of the ciphertext).
        // TODO(Jay): If we had a way to pass submatrix as mutable matrix then we can
        // separate sampling for moduli chain q and p_s
        RandomUniformDist::random_fill(&mut prng, q_and_specuialp_moduli_chain.as_slice(), c1_k);

        // Sample error polynomial in R_{QP_s} and store in c_{0,k}
        RandomGaussianDist::random_fill(&mut prng, q_and_specuialp_moduli_chain.as_slice(), c0_k);

        {
            // TODO(Jay): I really wish there was a way to pass submatrix as matrix because
            // that will allow us to use forward_lazy directly.
            izip!(
                q_nttops.iter(),
                c0_k.iter_rows_mut().take(q_size),
                c1_k.iter_rows_mut().take(q_size)
            )
            .for_each(|(nttqi, c0_modqi, c1_modqi)| {
                nttqi.forward_lazy(c0_modqi.as_mut());
                nttqi.forward_lazy(c1_modqi.as_mut());
            });
            izip!(
                specialp_nttops.iter(),
                c0_k.iter_rows_mut().skip(q_size),
                c1_k.iter_rows_mut().skip(q_size)
            )
            .for_each(|(nttqi, c0_modqi, c1_modqi)| {
                nttqi.forward_lazy(c0_modqi.as_mut());
                nttqi.forward_lazy(c1_modqi.as_mut());
            });
        }

        // part Q
        izip!(
            gammak.iter(),
            q_moduli_chain.iter(),
            q_modops.iter(),
            q_nttops.iter(),
            p_eval.iter_rows(),
            c0_k.iter_rows_mut().take(q_size),
            c1_k.iter_rows().take(q_size),
            s_partq.iter_rows()
        )
        .for_each(
            |(gammak_modqi, qi, qi_modop, qi_nttop, p_modqi, c0k_modqi, c1k_modqi, s_modqi)| {
                // Note that c_{1,k} \mod qi equals evaluation representation of -a \mod qi
                // a * s \mod qi
                let mut a_s = c1k_modqi.clone();
                qi_modop.neg_lazy_mod_vec(a_s.as_mut());
                qi_modop.mul_lazy_mod_vec(a_s.as_mut(), s_modqi.as_ref());

                // e + a * s \mod qi
                qi_modop.add_lazy_mod_vec(c0k_modqi.as_mut(), a_s.as_ref());

                // gamma_k * p \mod qi
                let mut p_modqi_clone = p_modqi.clone();
                qi_modop.scalar_mul_lazy_mod_vec(p_modqi_clone.as_mut(), *gammak_modqi);

                // e + a * s + gamma_k * p \mod qi
                qi_modop.add_lazy_mod_vec(c0k_modqi.as_mut(), p_modqi_clone.as_ref());
            },
        );

        // part specialP
        izip!(
            specialp_moduli_chain.iter(),
            specialp_modops.iter(),
            specialp_nttops.iter(),
            c0_k.iter_rows_mut().skip(q_size),
            c1_k.iter_rows().skip(q_size),
            s_partspecialp.iter_rows()
        )
        .for_each(|(pj, pj_modop, pj_nttop, c0k_modpj, c1k_modpj, s_modpj)| {
            // a * s \mod pj
            let mut a_s = c1k_modpj.clone();
            pj_modop.neg_lazy_mod_vec(a_s.as_mut());
            pj_modop.mul_lazy_mod_vec(a_s.as_mut(), s_modpj.as_ref());

            // e + a * s \mod pj
            pj_modop.add_lazy_mod_vec(c0k_modpj.as_mut(), a_s.as_ref());
        });
    }
}

/// Input polynomial x must be in coefficient form. Outout ciphertext
/// polynomials c0_out and c1_out are in Evaluation representation
pub fn keyswitch<
    M: Matrix,
    MMut: MatrixMut<MatElement = M::MatElement> + MatrixEntity + Clone,
    P: HybridKskRuntimeParameters<Scalar = M::MatElement>,
>(
    x: &M,
    c0_out: &mut MMut,
    c1_out: &mut MMut,
    ksk_polys: &[M],
    params: &P,
    level: usize,
) where
    <MMut as Matrix>::R: RowMut,
    M::MatElement: UnsignedInteger,
    // TODO(Jay): [Remove] these additional bounds are a consquenece of
    // BarrettBackend & Montgomery trait. Check issue #12
    <P as Parameters>::Scalar: AsPrimitive<u128> + PrimInt,
    u128: AsPrimitive<<P as Parameters>::Scalar>,
{
    let dnum = params.dnum();

    let q_moduli_chain = params.q_moduli_chain_at_level(level);
    let specialp_moduli_chain = params.specialp_moduli_chain_at_level(level);

    let q_nttops = params.q_nttops_at_level(level);
    let specialp_nttops = params.specialp_nttops_at_level(level);

    let q_modops = params.q_modops_at_level(level);
    let specialp_modops = params.specialp_modops_at_level(level);

    let q_size = q_moduli_chain.len();
    let specialp_size = specialp_moduli_chain.len();
    let ring_size = params.ring_size();

    debug_assert!(
        x.dimension() == (q_size, ring_size),
        "Dimensions of input polynomial to key switch are incorrect: Expect {:?} got {:?}",
        x.dimension(),
        (q_size, ring_size),
    );

    let mut k = 0;
    let mut c0_sum_subbasisq = MMut::zeros(q_size, ring_size);
    let mut c1_sum_subbasisq = MMut::zeros(q_size, ring_size);
    let mut c0_sum_subbasisspecialp = MMut::zeros(specialp_size, ring_size);
    let mut c1_sum_subbasisspecialp = MMut::zeros(specialp_size, ring_size);

    let qj_over_qji_inv_modqji_all_k = params.qj_over_qji_inv_modqji_level(level);
    let qj_over_qji_per_modqi_all_k = params.qj_over_qji_per_modqi_at_level(level);
    let qj_over_qji_per_modspecialpj_all_k = params.qj_over_qji_per_modspecialpj_at_level(level);
    for i in (0..q_size).step_by(dnum) {
        let start = i;
        let end = std::cmp::min(i + dnum, q_size); // end is excluded from range

        // let qj_moduli_chain = &q_moduli_chain[start..end];
        // let qj_modops = &q_modops[start..end];

        let qj_over_qji_inv_modqji = &qj_over_qji_inv_modqji_all_k[k];
        debug_assert!(
            qj_over_qji_inv_modqji.len() == end - start,
            "Q_j/q_[{k},i]^[-1] switch precomputes are incorrect: Expected length {} but got {}",
            end - start,
            qj_over_qji_inv_modqji.len()
        );

        let qj_over_qji_per_modqi = &qj_over_qji_per_modqi_all_k[k];
        debug_assert!(
            qj_over_qji_per_modqi.len() == q_size - (end - start),
            "Q_j/q_[{k},i] mod q_i switch precomputes are incorrect: Expected length {} but got {}",
            q_size - (end - start),
            qj_over_qji_per_modqi.len()
        );

        let qj_over_qji_per_modspecialpj = &qj_over_qji_per_modspecialpj_all_k[k];
        debug_assert!(
            qj_over_qji_per_modspecialpj.len() == specialp_size,
            "Q_j/q_[{k},i] mod specialpj switch precomputes are incorrect: Expected length {} but got {}",
            specialp_size,
            qj_over_qji_per_modspecialpj.len()
        );

        // Switch basis x \in Qj from Qj to QP_s
        // FIXME(Jay): Consolidate all 3 into 1
        let mut q_till_start = MMut::zeros(start, ring_size);
        let mut q_from_end = MMut::zeros(q_size - end, ring_size);
        let mut specialp_all = MMut::zeros(specialp_size, ring_size);
        for ri in 0..ring_size {
            let mut q_values = Vec::with_capacity(end - start);
            for j in start..end {
                q_values.push(
                    q_modops[j].mul_mod_fast_lazy(
                        *x.get_element(j, ri),
                        qj_over_qji_inv_modqji[j - start],
                    ),
                );
            }

            // q_0..q_{start-1}
            for (index, (qi, modqi, qj_over_qji_modqi)) in izip!(
                &q_moduli_chain[..start],
                &q_modops[..start],
                &qj_over_qji_per_modqi[..start]
            )
            .enumerate()
            {
                // map q_values to mont space
                let q_in_mont_space = q_values
                    .iter()
                    .map(|v| modqi.normal_to_mont_space_lazy(*v))
                    .collect_vec();

                let tmp = modqi.mont_fma(&q_in_mont_space, &qj_over_qji_modqi);
                let tmp = modqi.mont_to_normal_lazy(tmp);
                q_till_start.set(index, ri, tmp)
            }

            // q_{end}..q_{size-1}
            for (index, (qi, modqi, qj_over_qji_modqi)) in izip!(
                &q_moduli_chain[end..],
                &q_modops[end..],
                &qj_over_qji_per_modqi[start..]
            )
            .enumerate()
            {
                // map q_values to mont space
                let q_in_mont_space = q_values
                    .iter()
                    .map(|v| modqi.normal_to_mont_space_lazy(*v))
                    .collect_vec();

                let tmp = modqi.mont_fma(&q_in_mont_space, &qj_over_qji_modqi);
                let tmp = modqi.mont_to_normal_lazy(tmp);
                q_from_end.set(index, ri, tmp)
            }

            // p_0..p_l
            for (index, (pj, modpj, qj_over_qji_modpj)) in izip!(
                specialp_moduli_chain,
                specialp_modops,
                qj_over_qji_per_modspecialpj
            )
            .enumerate()
            {
                // map q_values to mont space
                let q_in_mont_space = q_values
                    .iter()
                    .map(|v| modpj.normal_to_mont_space_lazy(*v))
                    .collect_vec();

                let tmp = modpj.mont_fma(&q_in_mont_space, &qj_over_qji_modpj);
                let tmp = modpj.mont_to_normal_lazy(tmp);
                specialp_all.set(index, ri, tmp);
            }
        }

        // Multiply x \in QP_s with corresponding ksk ciphertext

        // q_0..q_{start-1}
        foward_lazy(&mut q_till_start, &q_nttops[..start]);
        let mut q_till_start_c0 = q_till_start.clone();
        let mut q_till_start_c1 = q_till_start;
        for (modqi, x_qi_c0, x_qi_c1, c0_k_qi, c1_k_qi, c0_sum_qi, c1_sum_qi) in izip!(
            q_modops[..start].iter(),
            q_till_start_c0.iter_rows_mut(),
            q_till_start_c1.iter_rows_mut(),
            ksk_polys[2 * k].iter_rows().take(start),
            ksk_polys[2 * k + 1].iter_rows().take(start),
            c0_sum_subbasisq.iter_rows_mut().take(start),
            c1_sum_subbasisq.iter_rows_mut().take(start)
        ) {
            modqi.mul_lazy_mod_vec(x_qi_c0.as_mut(), c0_k_qi.as_ref());
            modqi.mul_lazy_mod_vec(x_qi_c1.as_mut(), c1_k_qi.as_ref());

            modqi.add_lazy_mod_vec(c0_sum_qi.as_mut(), x_qi_c0.as_ref());
            modqi.add_lazy_mod_vec(c1_sum_qi.as_mut(), x_qi_c1.as_ref());
        }

        // q_{start}..q_{end-1}
        let mut x_partqj = {
            // TODO (Jay): Require a better way to copy elements of sub-matrix
            let mut tmp = MMut::zeros(end - start, ring_size);
            izip!(
                tmp.iter_rows_mut(),
                x.iter_rows().skip(start).take(end - start)
            )
            .for_each(|(to, from)| to.as_mut().copy_from_slice(from.as_ref()));
            tmp
        };
        foward_lazy(&mut x_partqj, &q_nttops[start..end]);
        let mut x_partqj_c0 = x_partqj.clone();
        let mut x_partqj_c1 = x_partqj;
        for (modqi, x_qi_c0, x_qi_c1, c0_k_qi, c1_k_qi, c0_sum_qi, c1_sum_qi) in izip!(
            q_modops[start..end].iter(),
            x_partqj_c0.iter_rows_mut(),
            x_partqj_c1.iter_rows_mut(),
            ksk_polys[2 * k].iter_rows().skip(start),
            ksk_polys[2 * k + 1].iter_rows().skip(start),
            c0_sum_subbasisq.iter_rows_mut().skip(start),
            c1_sum_subbasisq.iter_rows_mut().skip(start)
        ) {
            modqi.mul_lazy_mod_vec(x_qi_c0.as_mut(), c0_k_qi.as_ref());
            modqi.mul_lazy_mod_vec(x_qi_c1.as_mut(), c1_k_qi.as_ref());

            modqi.add_lazy_mod_vec(c0_sum_qi.as_mut(), x_qi_c0.as_ref());
            modqi.add_lazy_mod_vec(c1_sum_qi.as_mut(), x_qi_c1.as_ref());
        }

        // q_end..q_{size-1}
        foward_lazy(&mut q_from_end, &q_nttops[end..]);
        let mut q_from_end_c0 = q_from_end.clone();
        let mut q_from_end_c1 = q_from_end;
        for (modqi, x_qi_c0, x_qi_c1, c0_k_qi, c1_k_qi, c0_sum_qi, c1_sum_qi) in izip!(
            q_modops[end..].iter(),
            q_from_end_c0.iter_rows_mut(),
            q_from_end_c1.iter_rows_mut(),
            ksk_polys[2 * k].iter_rows().skip(end),
            ksk_polys[2 * k + 1].iter_rows().skip(end),
            c0_sum_subbasisq.iter_rows_mut().skip(end),
            c1_sum_subbasisq.iter_rows_mut().skip(end)
        ) {
            modqi.mul_lazy_mod_vec(x_qi_c0.as_mut(), c0_k_qi.as_ref());
            modqi.mul_lazy_mod_vec(x_qi_c1.as_mut(), c1_k_qi.as_ref());

            modqi.add_lazy_mod_vec(c0_sum_qi.as_mut(), x_qi_c0.as_ref());
            modqi.add_lazy_mod_vec(c1_sum_qi.as_mut(), x_qi_c1.as_ref());
        }

        // p_0..p{size-1}
        foward_lazy(&mut specialp_all, specialp_nttops);
        let mut specialp_all_c0 = specialp_all.clone();
        let mut specialp_all_c1 = specialp_all;
        for (modpj, x_pj_c0, x_pj_c1, c0_k_pj, c1_k_pj, c0_sum_pj, c1_sum_pj) in izip!(
            specialp_modops.iter(),
            specialp_all_c0.iter_rows_mut(),
            specialp_all_c1.iter_rows_mut(),
            ksk_polys[2 * k].iter_rows().skip(q_size),
            ksk_polys[2 * k + 1].iter_rows().skip(q_size),
            c0_sum_subbasisspecialp.iter_rows_mut(),
            c1_sum_subbasisspecialp.iter_rows_mut()
        ) {
            modpj.mul_lazy_mod_vec(x_pj_c0.as_mut(), c0_k_pj.as_ref());
            modpj.mul_lazy_mod_vec(x_pj_c1.as_mut(), c1_k_pj.as_ref());

            modpj.add_lazy_mod_vec(c0_sum_pj.as_mut(), x_pj_c0.as_ref());
            modpj.add_lazy_mod_vec(c1_sum_pj.as_mut(), x_pj_c1.as_ref());
        }

        k += 1;
    }

    // Change representation of only subbasis P_s to Coefficient for approximate mod
    // down
    backward_lazy(&mut c0_sum_subbasisspecialp, specialp_nttops);
    backward_lazy(&mut c1_sum_subbasisspecialp, specialp_nttops);

    // Approximate Mod Down: v \in QP_s => (1/P)v \in Q
    let specialp_over_specialpj_inv_modspecialpj =
        params.specialp_over_specialpj_inv_modspecialpj_at_level(level);
    let specialp_over_specialpj_per_modqi =
        params.specialp_over_specialpj_per_modqi_at_level(level);
    let specialp_inv_modqi = params.specialp_inv_modqi_at_level(level);
    approximate_mod_down(
        c0_out,
        &c0_sum_subbasisq,
        &c0_sum_subbasisspecialp,
        specialp_over_specialpj_inv_modspecialpj,
        specialp_over_specialpj_per_modqi,
        specialp_inv_modqi,
        q_modops,
        specialp_modops,
        q_nttops,
        ring_size,
        q_size,
        specialp_size,
    );
    approximate_mod_down(
        c1_out,
        &c1_sum_subbasisq,
        &c1_sum_subbasisspecialp,
        specialp_over_specialpj_inv_modspecialpj,
        specialp_over_specialpj_per_modqi,
        specialp_inv_modqi,
        q_modops,
        specialp_modops,
        q_nttops,
        ring_size,
        q_size,
        specialp_size,
    );
    reduce_from_lazy_mut(c0_out, q_modops);
    reduce_from_lazy_mut(c1_out, q_modops);
}
