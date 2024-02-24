use itertools::Itertools;
use num_bigint::BigUint;
use num_traits::{FromPrimitive, One, ToPrimitive};

use crate::{
    core_crypto::{
        modulus::{
            ModulusBackendConfig, MontgomeryBackend, MontgomeryScalar, NativeModulusBackend,
        },
        ntt::{NativeNTTBackend, NttConfig},
    },
    parameters::Parameters,
    utils::mod_inverse_big_unit,
};

use self::ops::{HybridKskKeyGenParameters, HybridKskRuntimeParameters};

pub(crate) mod ops;

struct HybridKskParametersScalarU64 {
    q_modops: Vec<NativeModulusBackend>,
    specialp_modops: Vec<NativeModulusBackend>,

    q_nttops: Vec<NativeNTTBackend>,
    specialp_nttops: Vec<NativeNTTBackend>,

    q_moduli_chain: Vec<u64>,
    specialp_moduli_chain: Vec<u64>,

    // Key generation parametrs //
    gammak_modqi_at_level: Vec<Vec<Vec<u64>>>,

    // Ksk runtime parameters //
    // To switch x \in Qj to x \in QP_s //
    qj_over_qji_inv_modqji_all_k_at_level: Vec<Vec<Vec<u64>>>,
    qj_over_qji_per_modqi_all_k_at_level: Vec<Vec<Vec<Vec<MontgomeryScalar<u64>>>>>,
    qj_over_qji_per_modspecialpj_all_k_at_level: Vec<Vec<Vec<Vec<MontgomeryScalar<u64>>>>>,
    // To scale by 1/P_s and switch v \in QP_s to (1/P)v \in Q //
    specialp_over_specialpj_inv_modspecialpj: Vec<u64>,
    specialp_over_specialpj_per_modqi: Vec<Vec<MontgomeryScalar<u64>>>,
    specialp_inv_modqi: Vec<u64>,

    // Other //
    dnum: usize,
    q_moduli_chain_len: usize,
    ring_size: usize,
    alpha_at_level: Vec<usize>,
}

impl HybridKskParametersScalarU64 {
    pub fn new(
        q_moduli_chain: &[u64],
        specialp_moduli_chain: &[u64],
        dnum: usize,
        ring_size: usize,
    ) -> Self {
        // In general if q moduli chain for the second last level will be smaller than
        // bit size of special p. And it is possible to improve efficiency for hybrid
        // keyswitching by dropping 1-2 primes from special p moduli chain. But we don't
        // handle for this case. It should instead be handled at parameter generation
        // phase

        let mut big_specialp = BigUint::one();
        for pj in specialp_moduli_chain.iter() {
            big_specialp *= *pj;
        }

        let mut gammak_modqi_at_level = vec![];
        let mut alpha_at_level = vec![];
        let mut qj_over_qji_inv_modqji_all_k_at_level = vec![];
        let mut qj_over_qji_per_modqi_all_k_at_level = vec![];
        let mut qj_over_qji_per_modspecialpj_all_k_at_level = vec![];
        for lvl in 0..(q_moduli_chain.len()) {
            let q_moduli_chain_at_level = &q_moduli_chain[..q_moduli_chain.len() - lvl];

            let mut big_q = BigUint::one();
            for qi in q_moduli_chain_at_level.iter() {
                big_q *= *qi;
            }
            let primes_count_at_lvl = q_moduli_chain_at_level.len();

            let mut gammak_modqi = vec![];
            let mut qj_over_qji_inv_modqji_all_k = vec![];
            let mut qj_over_qji_per_modqi_all_k = vec![];
            let mut qj_over_qji_per_modspecialpj_all_k = vec![];
            for i in (0..primes_count_at_lvl).step_by(dnum) {
                let start = i;
                let end = std::cmp::min(i + dnum, primes_count_at_lvl);

                let mut big_qj = BigUint::one();
                for j in (start..end) {
                    big_qj *= q_moduli_chain_at_level[j];
                }

                // gamma = P * Q/Qj^{-1} * Q/Qj
                let q_over_qj = &big_q / &big_qj;
                let gamma = &big_specialp * mod_inverse_big_unit(&q_over_qj, &big_qj) * &q_over_qj;

                let mut tmp = vec![];
                for qi in q_moduli_chain_at_level.iter() {
                    tmp.push((&gamma % *qi).to_u64().unwrap());
                }
                gammak_modqi.push(tmp);

                // Ksk runtime parameters //

                // Parameters to switch x \in Qj to QP_s//
                let mut qj_over_qji_inv_modqji = vec![];
                let mut qj_over_qji_vec = vec![];
                for j in (start..end) {
                    let qji = q_moduli_chain_at_level[j];
                    let qj_over_qji = &big_qj / qji;

                    qj_over_qji_inv_modqji.push(
                        mod_inverse_big_unit(&qj_over_qji, &BigUint::from_u64(qji).unwrap())
                            .to_u64()
                            .unwrap(),
                    );

                    qj_over_qji_vec.push(qj_over_qji);
                }

                let mut qj_over_qji_per_modqi = vec![];
                // q_0..q_{start-1}
                for qi in q_moduli_chain_at_level[..start].iter() {
                    // TODO (Jay): Move montgomery arithemtic out of the trait as an independent
                    // function. Then remove this and call the function directly
                    let modqi = NativeModulusBackend::initialise(*qi);
                    qj_over_qji_per_modqi.push(
                        qj_over_qji_vec
                            .iter()
                            .map(|v| modqi.normal_to_mont_space((v % qi).to_u64().unwrap()))
                            .collect_vec(),
                    );
                }
                // q_end..q_{size-1}
                for qi in q_moduli_chain_at_level[end..].iter() {
                    // TODO (Jay): Move montgomery arithemtic out of the trait as an independent
                    // function. Then remove this and call the function directly
                    let modqi = NativeModulusBackend::initialise(*qi);
                    qj_over_qji_per_modqi.push(
                        qj_over_qji_vec
                            .iter()
                            .map(|v| modqi.normal_to_mont_space((v % qi).to_u64().unwrap()))
                            .collect_vec(),
                    );
                }
                let mut qj_over_qji_per_modspecialpj = vec![];
                for pj in specialp_moduli_chain.iter() {
                    let modpj = NativeModulusBackend::initialise(*pj);
                    qj_over_qji_per_modspecialpj.push(
                        qj_over_qji_vec
                            .iter()
                            .map(|v| modpj.normal_to_mont_space((v % pj).to_u64().unwrap()))
                            .collect_vec(),
                    );
                }

                qj_over_qji_inv_modqji_all_k.push(qj_over_qji_inv_modqji);
                qj_over_qji_per_modqi_all_k.push(qj_over_qji_per_modqi);
                qj_over_qji_per_modspecialpj_all_k.push(qj_over_qji_per_modspecialpj);
            }

            gammak_modqi_at_level.push(gammak_modqi);
            qj_over_qji_inv_modqji_all_k_at_level.push(qj_over_qji_inv_modqji_all_k);
            qj_over_qji_per_modqi_all_k_at_level.push(qj_over_qji_per_modqi_all_k);
            qj_over_qji_per_modspecialpj_all_k_at_level.push(qj_over_qji_per_modspecialpj_all_k);

            let alpha = (q_moduli_chain_at_level.len() as f64 / dnum as f64).ceil() as usize;
            alpha_at_level.push(alpha);
        }

        // Parameters to switch v \in QP_s to (1/P)v \in Q //
        // P_s is independent of level so the precomputes to switch and scale by 1/P_s a
        // polynomial in QP_s (where Q is modulus at level) to Q remain unchanged
        let mut specialp_over_specialpj = vec![];
        let mut specialp_over_specialpj_inv_modspecialpj = vec![];
        for pj in specialp_moduli_chain.iter() {
            let tmp = &big_specialp / *pj;
            specialp_over_specialpj_inv_modspecialpj.push(
                mod_inverse_big_unit(&tmp, &BigUint::from_u64(*pj).unwrap())
                    .to_u64()
                    .unwrap(),
            );
            specialp_over_specialpj.push(tmp);
        }
        let mut specialp_over_specialpj_per_modqi = vec![];
        let mut specialp_inv_modqi = vec![];
        for qi in q_moduli_chain.iter() {
            let modqi = NativeModulusBackend::initialise(*qi);
            let tmp = specialp_over_specialpj
                .iter()
                .map(|v| modqi.normal_to_mont_space((v % *qi).to_u64().unwrap()))
                .collect_vec();
            specialp_over_specialpj_per_modqi.push(tmp);

            specialp_inv_modqi.push(
                mod_inverse_big_unit(&big_specialp, &BigUint::from_u64(*qi).unwrap())
                    .to_u64()
                    .unwrap(),
            );
        }

        let q_modops = q_moduli_chain
            .iter()
            .map(|qi| NativeModulusBackend::initialise(*qi))
            .collect_vec();
        let specialp_modops = specialp_moduli_chain
            .iter()
            .map(|pj| NativeModulusBackend::initialise(*pj))
            .collect_vec();
        let q_nttops = q_moduli_chain
            .iter()
            .map(|qi| NativeNTTBackend::init(*qi, ring_size))
            .collect_vec();
        let specialp_nttops = specialp_moduli_chain
            .iter()
            .map(|pj| NativeNTTBackend::init(*pj, ring_size))
            .collect_vec();

        HybridKskParametersScalarU64 {
            q_modops,
            specialp_modops,
            q_nttops,
            specialp_nttops,

            // KSK key gen //
            gammak_modqi_at_level,

            // KSK runtime //
            qj_over_qji_inv_modqji_all_k_at_level,
            qj_over_qji_per_modqi_all_k_at_level,
            qj_over_qji_per_modspecialpj_all_k_at_level,
            specialp_over_specialpj_inv_modspecialpj,
            specialp_over_specialpj_per_modqi,
            specialp_inv_modqi,

            q_moduli_chain: q_moduli_chain.to_vec(),
            specialp_moduli_chain: specialp_moduli_chain.to_vec(),
            dnum,
            q_moduli_chain_len: q_moduli_chain.len(),
            ring_size,
            alpha_at_level,
        }
    }

    fn q_union_speciap_size_at_level(&self, level: usize) -> usize {
        self.q_moduli_chain_len - level + self.specialp_moduli_chain.len()
    }
}

impl Parameters for HybridKskParametersScalarU64 {
    type Scalar = u64;
}

impl HybridKskKeyGenParameters for HybridKskParametersScalarU64 {
    type ModOp = NativeModulusBackend;
    type NttOp = NativeNTTBackend;

    fn dnum(&self) -> usize {
        self.dnum
    }

    fn gammak_modqi_at_level(&self, level: usize) -> &[Vec<Self::Scalar>] {
        self.gammak_modqi_at_level[level].as_slice()
    }

    fn q_modops_at_level(&self, level: usize) -> &[Self::ModOp] {
        &self.q_modops[..self.q_moduli_chain_len - level]
    }

    fn q_moduli_chain_at_level(&self, level: usize) -> &[Self::Scalar] {
        &self.q_moduli_chain[..self.q_moduli_chain_len - level]
    }

    fn q_nttops_at_level(&self, level: usize) -> &[Self::NttOp] {
        &self.q_nttops[..self.q_moduli_chain_len - level]
    }

    fn specialp_modops_at_level(&self, level: usize) -> &[Self::ModOp] {
        &self.specialp_modops
    }

    fn specialp_nttops_at_level(&self, level: usize) -> &[Self::NttOp] {
        &self.specialp_nttops
    }

    fn ring_size(&self) -> usize {
        self.ring_size
    }

    fn specialp_moduli_chain_at_level(&self) -> &[Self::Scalar] {
        &self.specialp_moduli_chain
    }

    fn alpha_at_level(&self, level: usize) -> usize {
        self.alpha_at_level[level]
    }
}

impl HybridKskRuntimeParameters for HybridKskParametersScalarU64 {
    type ModOp = NativeModulusBackend;
    type NttOp = NativeNTTBackend;

    fn dnum(&self) -> usize {
        self.dnum
    }

    fn ring_size(&self) -> usize {
        self.ring_size
    }

    fn q_modops_at_level(&self, level: usize) -> &[Self::ModOp] {
        &self.q_modops[..self.q_moduli_chain_len - level]
    }

    fn q_moduli_chain_at_level(&self, level: usize) -> &[Self::Scalar] {
        &self.q_moduli_chain[..self.q_moduli_chain_len - level]
    }

    fn q_nttops_at_level(&self, level: usize) -> &[Self::NttOp] {
        &self.q_nttops[..self.q_moduli_chain_len - level]
    }

    fn specialp_modops_at_level(&self, level: usize) -> &[Self::ModOp] {
        &self.specialp_modops
    }

    fn specialp_nttops_at_level(&self, level: usize) -> &[Self::NttOp] {
        &self.specialp_nttops
    }

    fn qj_over_qji_inv_modqji_level(&self, level: usize) -> &[Vec<Self::Scalar>] {
        &self.qj_over_qji_inv_modqji_all_k_at_level[level]
    }

    fn qj_over_qji_per_modqi_at_level(
        &self,
        level: usize,
    ) -> &[Vec<Vec<MontgomeryScalar<Self::Scalar>>>] {
        &self.qj_over_qji_per_modqi_all_k_at_level[level]
    }

    fn qj_over_qji_per_modspecialpj_at_level(
        &self,
        level: usize,
    ) -> &[Vec<Vec<MontgomeryScalar<Self::Scalar>>>] {
        &self.qj_over_qji_per_modspecialpj_all_k_at_level[level]
    }

    fn specialp_inv_modqi_at_level(&self, level: usize) -> &[Self::Scalar] {
        &self.specialp_inv_modqi[..self.q_moduli_chain_len - level]
    }

    fn specialp_moduli_chain_at_level(&self, level: usize) -> &[Self::Scalar] {
        &self.specialp_moduli_chain
    }

    fn specialp_over_specialpj_inv_modspecialpj_at_level(&self, level: usize) -> &[Self::Scalar] {
        &self.specialp_over_specialpj_inv_modspecialpj
    }

    fn specialp_over_specialpj_per_modqi_at_level(
        &self,
        level: usize,
    ) -> &[Vec<MontgomeryScalar<Self::Scalar>>] {
        &self.specialp_over_specialpj_per_modqi[..self.q_moduli_chain_len - level]
    }
}

#[cfg(test)]
mod tests {
    use itertools::izip;
    use num_bigint::BigInt;
    use rand::SeedableRng;
    use rand_chacha::{rand_core::le, ChaCha8Rng};

    use crate::{
        core_crypto::{
            matrix::{Matrix, MatrixMut},
            modulus::ModulusVecBackend,
            prime::generate_primes_vec,
            random::{DefaultU64SeededRandomGenerator, InitWithSeed, RandomUniformDist},
            ring::{
                add_lazy_mut, backward, backward_lazy, foward, foward_lazy, mul_lazy_mut, neg_mut,
                sub_mut,
            },
        },
        keys::SecretKey,
        utils::{convert::TryConvertFrom, test_utils::TestTernarySecret},
    };

    use self::ops::{generate_key, keyswitch};

    use super::*;

    #[test]
    fn hybird_key_switching_works() {
        let hybrid_ksk_params = {
            // Setup //
            let ring_size = 1 << 3;
            let q_moduli_chain_sizes = vec![60, 60, 60, 60, 50, 40, 50];
            let dnum = 2;
            let specialp_moduli_chain_sizes = vec![60, 60, 50];

            let q_moduli_chain = &generate_primes_vec(&q_moduli_chain_sizes, ring_size, &[]);
            let specialp_moduli_chain =
                &generate_primes_vec(&specialp_moduli_chain_sizes, ring_size, &q_moduli_chain);

            HybridKskParametersScalarU64::new(
                q_moduli_chain,
                &specialp_moduli_chain,
                dnum,
                ring_size,
            )
        };

        let mut rng = DefaultU64SeededRandomGenerator::new();

        // Test key switching at all levels from 0 till max_level
        for level in (0..hybrid_ksk_params.q_moduli_chain.len()) {
            // Test KSK //
            {
                let q_moduli_chain_at_level =
                    HybridKskKeyGenParameters::q_moduli_chain_at_level(&hybrid_ksk_params, level);
                let specialp_moduli_chain_at_level =
                    HybridKskKeyGenParameters::specialp_moduli_chain_at_level(&hybrid_ksk_params);
                let ring_size = HybridKskKeyGenParameters::ring_size(&hybrid_ksk_params);

                let q_nttops_at_level =
                    HybridKskKeyGenParameters::q_nttops_at_level(&hybrid_ksk_params, level);
                let q_modops_at_level =
                    HybridKskKeyGenParameters::q_modops_at_level(&hybrid_ksk_params, level);

                let mut p1 =
                    <Vec<Vec<u64>> as Matrix>::zeros(q_moduli_chain_at_level.len(), ring_size);
                RandomUniformDist::random_fill(&mut rng, q_moduli_chain_at_level, &mut p1);

                let secret = TestTernarySecret::new(&mut rng, ring_size);

                // Create key switching keys
                let mut ksk_polys = (0..(hybrid_ksk_params.alpha_at_level(level) * 2))
                    .into_iter()
                    .map(|_| {
                        <Vec<Vec<u64>> as Matrix>::zeros(
                            hybrid_ksk_params.q_union_speciap_size_at_level(level),
                            ring_size,
                        )
                    })
                    .collect_vec();
                let mut seed = <ChaCha8Rng as SeedableRng>::Seed::default();
                let mut p1_eval = p1.clone();
                foward_lazy(&mut p1_eval, &q_nttops_at_level);
                generate_key(
                    &hybrid_ksk_params,
                    level,
                    &p1_eval,
                    &secret,
                    &mut ksk_polys,
                    &mut seed,
                    &mut rng,
                );

                // Check key switching keys correctly enrypt gamma values
                {
                    assert!(
                        ksk_polys.len() == hybrid_ksk_params.alpha_at_level(level) * 2,
                        "Ciphertext in Ksk:{} but expect 2*alpha:{} at level:{level}",
                        ksk_polys.len(),
                        2 * hybrid_ksk_params.alpha_at_level(level)
                    );

                    let specialp_nttops_at_level =
                        HybridKskKeyGenParameters::specialp_nttops_at_level(
                            &hybrid_ksk_params,
                            level,
                        );
                    let specialp_modops_at_level =
                        HybridKskKeyGenParameters::specialp_modops_at_level(
                            &hybrid_ksk_params,
                            level,
                        );

                    let qp_chain = q_moduli_chain_at_level
                        .to_vec()
                        .into_iter()
                        .chain(specialp_moduli_chain_at_level.to_vec().into_iter())
                        .collect_vec();
                    let qp_modops = q_modops_at_level
                        .to_vec()
                        .into_iter()
                        .chain(specialp_modops_at_level.to_vec())
                        .collect_vec();
                    // Note(Jay): Due to the way NativeNTTBackend handles generating psi, two
                    // separately intialised instances of Ntt with same modulus and ring are not
                    // compatible. Given ksk polynomials are in evaluation, we `Clone` the Ntt
                    // instances which were used in hybrid key gen phase and avoid creating new
                    // instances, because doing so will be incorrect.
                    let qp_nttops = q_nttops_at_level
                        .to_vec()
                        .into_iter()
                        .chain(specialp_nttops_at_level.to_vec())
                        .collect_vec();

                    for (gamma_k, c_k) in izip!(
                        hybrid_ksk_params.gammak_modqi_at_level(level).iter(),
                        ksk_polys.chunks_exact(2)
                    ) {
                        let mut c0 = c_k[0].clone();
                        let mut c1 = c_k[1].clone();

                        let mut s = <Vec<Vec<u64>>>::try_convert_from(secret.values(), &qp_chain);
                        foward_lazy(&mut s, &qp_nttops);
                        mul_lazy_mut(&mut c1, &s, &qp_modops);
                        add_lazy_mut(&mut c0, &c1, &qp_modops);
                        backward(&mut c0, &qp_nttops);

                        let mut gamma_p1 = p1.clone();
                        izip!(
                            q_modops_at_level.iter(),
                            gamma_p1.iter_rows_mut(),
                            gamma_k.iter()
                        )
                        .for_each(|(modqi, r, g)| {
                            modqi.scalar_mul_mod_vec(r.as_mut(), *g);
                        });

                        let mut c0_partq = c0[..q_moduli_chain_at_level.len()].to_vec();
                        sub_mut(&mut c0_partq, &gamma_p1, q_modops_at_level);
                        Vec::<BigInt>::try_convert_from(&c0_partq, q_moduli_chain_at_level)
                            .iter()
                            .for_each(|v| {
                                // Difference bits depend on error in key switching
                                // key. We yet don't have a way to tighly bound it, hence no way to
                                // say for sure what the maximum difference will be. For now, we
                                // satisfy ourselves with 5 bits since \sigma is usually 3.2 and
                                // probability of any error sample
                                // outside range -6 \sigma to 6 \sigma is very very less.
                                assert!(v.bits() < 5);
                            });

                        // part specialp of c0 is 0 since \gamma \mod pj = 0
                        let c0_partspecialp = c0[q_moduli_chain_at_level.len()..].to_vec();
                        Vec::<BigInt>::try_convert_from(
                            &c0_partspecialp,
                            specialp_moduli_chain_at_level,
                        )
                        .iter()
                        .for_each(|v| {
                            assert!(v.bits() < 5);
                        });
                    }
                }

                let mut p2 =
                    <Vec<Vec<u64>> as Matrix>::zeros(q_moduli_chain_at_level.len(), ring_size);
                RandomUniformDist::random_fill(&mut rng, q_moduli_chain_at_level, &mut p2);

                // Create key switch ciphertext polynomials
                let mut c0_out =
                    <Vec<Vec<u64>> as Matrix>::zeros(q_moduli_chain_at_level.len(), ring_size);
                let mut c1_out = c0_out.clone();
                // Key switch
                keyswitch(
                    &p2,
                    &mut c0_out,
                    &mut c1_out,
                    ksk_polys.as_slice(),
                    &hybrid_ksk_params,
                    level,
                );

                // Check key switch output ciphertext polynomials c0, c1 encrypt p1p2
                let mut s =
                    <Vec<Vec<u64>>>::try_convert_from(secret.values(), q_moduli_chain_at_level);
                foward_lazy(&mut s, q_nttops_at_level); // c0 + c1*s
                mul_lazy_mut(&mut c1_out, &s, q_modops_at_level);
                add_lazy_mut(&mut c0_out, &c1_out, q_modops_at_level);
                backward(&mut c0_out, q_nttops_at_level);
                let p1p2_out_big =
                    Vec::<BigUint>::try_convert_from(&c0_out, q_moduli_chain_at_level);

                // Expected p1*p2
                let mut p1p2_expected = p1.clone();
                foward_lazy(&mut p1p2_expected, q_nttops_at_level);
                foward_lazy(&mut p2, q_nttops_at_level);
                mul_lazy_mut(&mut p1p2_expected, &p2, q_modops_at_level);
                backward(&mut p1p2_expected, q_nttops_at_level);
                let p1p2_expected_big =
                    Vec::<BigUint>::try_convert_from(&p1p2_expected, q_moduli_chain_at_level);

                izip!(p1p2_out_big.iter(), p1p2_expected_big.iter()).for_each(|(p0, p1)| {
                    let bits = if p0 < &p1 {
                        (p1 - p0).bits()
                    } else {
                        (p0 - p1).bits()
                    };
                    // The difference depends on error in key switching keys + rounding error from
                    // $x \in R_QP -> (1/P x) \in R_Q$. Ideally it should not exceed more than a few
                    // bits
                    assert!(bits <= 4);
                });
            }
        }
    }
}
