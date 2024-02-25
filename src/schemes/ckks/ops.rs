use itertools::{izip, Itertools};
use std::{
    fmt::Debug,
    ops::{Rem, Sub},
};

use crate::{
    ciphertext::{Representation, RlweCiphertext, SeededCiphertext},
    core_crypto::{
        matrix::{Matrix, MatrixMut, RowMut},
        num::{BFloat, BInt, BUint, CastToZp, ComplexNumber, UnsignedInteger},
        random::{InitWithSeed, RandomGaussianDist, RandomUniformDist},
        ring::{
            add_lazy_mut, foward, foward_lazy, mul_lazy_mut, neg_lazy_mut, neg_mut,
            reduce_from_lazy_mut,
        },
    },
    keys::SecretKey,
    parameters::CkksEncDecParameters,
    utils::{bit_reverse_map, convert::TryConvertFrom},
};

// specialIFFT
pub fn special_inv_fft<F: BFloat, C: ComplexNumber<F> + Clone + Copy>(
    v: &mut [C],
    psi_powers: &[C],
    rot_group: &[usize],
) {
    debug_assert!(
        v.len().is_power_of_two(),
        "Expected length of input to be power of 2 but is {}",
        v.len()
    );
    debug_assert!(
        v.len() * 4 == psi_powers.len(),
        "psi_powers must have powers of psi for 0 <= j < 4l, but its length is {}",
        psi_powers.len()
    );

    let v_len = v.len();
    // l = M * 4, where M denotes M of M^{th} cyrclotomic polynomial
    let which_cyclo = 4 * v_len;
    let mut m = v_len;
    while m >= 1 {
        for i in (0..v_len).step_by(m) {
            let lenh = m >> 1;
            let lenq = m << 2;
            let gap = which_cyclo / lenq;

            for j in 0..lenh {
                let idx = (lenq - (rot_group[j] % lenq)) * gap;

                // X + Y
                let u = v[i + j] + &v[i + j + lenh];
                // (X - Y) \cdot \psi_{idx}
                let k = (v[i + j] - v[i + j + lenh]) * psi_powers[idx];

                v[i + j] = u;
                v[i + j + lenh] = k;
            }
        }

        m >>= 1;
    }

    bit_reverse_map(v);
    v.iter_mut()
        .for_each(|a| *a = *a / &F::from(v_len).unwrap());
}

pub fn special_fft<F: BFloat, C: ComplexNumber<F> + Copy>(
    v: &mut [C],
    psi_powers: &[C],
    rot_group: &[usize],
) {
    debug_assert!(
        v.len().is_power_of_two(),
        "Expected length of input to be power of 2 but is {}",
        v.len()
    );
    debug_assert!(
        v.len() * 4 == psi_powers.len(),
        "psi_powers must have powers of psi for 0 <= j < 4l, but its length is {}",
        psi_powers.len()
    );

    bit_reverse_map(v);

    let v_len = v.len();
    // l = M * 4, where M denotes M of M^{th} cyrclotomic polynomial
    let which_cyclo = 4 * v_len;
    let mut m = 2;
    while m <= v_len {
        for i in (0..v_len).step_by(m) {
            let lenh = m >> 1;
            let lenq = m << 2;
            let gap = which_cyclo / lenq;

            for j in 0..lenh {
                let idx = (rot_group[j] % lenq) * gap;

                // TODO(Jay): Remove bound of Copy
                let u = v[i + j];
                let k = v[i + j + lenh] * psi_powers[idx];

                // X + \psi Y
                v[i + j] = u + k;
                // X - \psi Y
                v[i + j + lenh] = u - k;
            }
        }
        m <<= 1;
    }
}

pub fn simd_encode<
    Scalar: UnsignedInteger + TryFrom<Uint>,
    Uint: BUint,
    F: BFloat + CastToZp<Uint>,
    C: ComplexNumber<F> + Clone + Copy,
    MMut: MatrixMut<MatElement = Scalar>,
    P: CkksEncDecParameters<F = F, Scalar = Scalar, BU = Uint, Complex = C>,
>(
    p: &mut MMut,
    m: &[C],
    params: &P,
    level: usize,
    delta: F,
) where
    <Scalar as TryFrom<Uint>>::Error: Debug,
    for<'a> &'a Uint: Rem<Scalar, Output = Uint> + Sub<&'a Uint, Output = Uint>,
    <MMut as Matrix>::R: RowMut,
{
    let ring_size = params.ring_size();
    let slots = ring_size >> 1;

    debug_assert!(
        m.len() == slots,
        "Expected m vector length {slots} but is {}",
        m.len()
    );

    let psi_powers = params.psi_powers();
    let rot_group = params.rot_group();
    debug_assert!(
        psi_powers.len() == ring_size * 2,
        "Expected psi^i for 0 <= i < {}(M) but psi_powers has length {}",
        ring_size * 2,
        psi_powers.len()
    );
    debug_assert!(
        rot_group.len() == slots,
        "Expected (5^j mod M) for 0 <= j < {}(l) but rot_group has length {}",
        slots,
        rot_group.len()
    );

    let mut m = m.to_vec();
    special_inv_fft(&mut m, &psi_powers, &rot_group);

    // scale by delta
    izip!(m.iter_mut()).for_each(|v| {
        *v = *v * &delta;
    });

    let q_moduli_chain = params.q_moduli_chain_at_level(level);
    let big_q = params.bigq_at_level(level);

    for ri in 0..ring_size >> 1 {
        let delta_m_ri = CastToZp::cast(&m[ri].re(), big_q);
        izip!(p.get_col_iter_mut(ri), q_moduli_chain.iter()).for_each(|(x_qi, qi)| {
            *x_qi = (&delta_m_ri % *qi).try_into().unwrap();
        });
    }

    for ri in ring_size >> 1..ring_size {
        let delta_m_ri = CastToZp::cast(&m[ri - (ring_size >> 1)].img(), big_q);
        izip!(p.get_col_iter_mut(ri), q_moduli_chain.iter()).for_each(|(x_qi, qi)| {
            *x_qi = (&delta_m_ri % *qi).try_into().unwrap();
        });
    }
}

pub fn simd_decode<
    Scalar: UnsignedInteger,
    F: BFloat,
    // TODO(Jay): Remove Copy bound
    C: ComplexNumber<F> + Copy,
    M: Matrix<MatElement = Scalar>,
    P: CkksEncDecParameters<F = F, Scalar = Scalar, Complex = C>,
>(
    p: &M,
    params: &P,
    level: usize,
    delta: F,
    m_out: &mut [C],
) where
    Vec<F>: TryConvertFrom<M, Parameters = [Scalar]>,
{
    let ring_size = params.ring_size();
    let slots = ring_size >> 1;
    debug_assert!(
        m_out.len() == ring_size >> 1,
        "Expected m_out to have {} slots but has {}",
        slots,
        m_out.len()
    );

    let q_moduli_chain = params.q_moduli_chain_at_level(level);
    // TODO(Jay): `try_convert_from` first computes recomposition factors. Hence is
    // quie expensive.
    let mut p_reals = Vec::<F>::try_convert_from(p, &q_moduli_chain);
    p_reals
        .iter_mut()
        .map(|v| {
            // scale by 1/delta
            *v = *v / delta;
        })
        .collect_vec();

    for k in 0..slots {
        m_out[k] = C::new(p_reals[k], p_reals[slots + k]);
    }

    let psi_powers = params.psi_powers();
    let rot_group = params.rot_group();
    special_fft(m_out, &psi_powers, &rot_group);
}

/// Encrypt message polynomial m_eval_lazy using secret key.
///
/// m_eval_lazy must be in Evaluation representation and can have lazy
/// coefficients.
pub fn secret_key_encryption<
    Scalar: UnsignedInteger,
    Seed,
    C: RlweCiphertext<Scalar = Scalar> + SeededCiphertext<Seed = Seed>,
    S: SecretKey,
    P: CkksEncDecParameters<Scalar = Scalar>,
    R: RandomGaussianDist<C::Poly, Parameters = [Scalar]>
        + RandomUniformDist<C::Poly, Parameters = [Scalar]>
        + RandomUniformDist<Seed, Parameters = u8>
        + InitWithSeed<Seed = Seed>,
>(
    c_out: &mut C,
    m_eval_lazy: &C::Poly,
    secret: &S,
    params: &P,
    rng: &mut R,
    level: usize,
) where
    C::Poly: TryConvertFrom<[S::Scalar], Parameters = [Scalar]> + Clone,
{
    let q_moduli_chain = params.q_moduli_chain_at_level(level);
    let ring_size = params.ring_size();
    debug_assert!(
        c_out.c_partq().len() == 2,
        "Expedted ciphertext to have 2 polynomials but has {}",
        c_out.c_partq().len()
    );
    debug_assert!(
        c_out.c_partq()[0].dimension() == (q_moduli_chain.len(), ring_size),
        "c0 poly should have {:?} dimension but has {:?}",
        (q_moduli_chain.len(), ring_size),
        c_out.c_partq()[0].dimension()
    );
    debug_assert!(
        c_out.c_partq()[1].dimension() == (q_moduli_chain.len(), ring_size),
        "c1 poly should have {:?} dimension but has {:?}",
        (q_moduli_chain.len(), ring_size),
        c_out.c_partq()[1].dimension()
    );

    let q_modops = params.q_modops_at_level(level);
    let q_nttops = params.q_nttops_at_level(level);

    {
        // Sample seeded -a in Coefficient form
        RandomUniformDist::<Seed>::random_fill(rng, &0u8, c_out.seed_mut());
        let mut prng = R::init_with_seed(c_out.seed());
        RandomUniformDist::<C::Poly>::random_fill(
            &mut prng,
            &q_moduli_chain,
            &mut c_out.c_partq_mut()[1],
        );
        foward_lazy(&mut c_out.c_partq_mut()[1], &q_nttops);
    }

    // a * s
    let mut a = c_out.c_partq_mut()[1].clone();
    neg_lazy_mut(&mut a, q_modops);
    let mut sa = C::Poly::try_convert_from(secret.values(), q_moduli_chain);
    foward_lazy(&mut sa, q_nttops);
    mul_lazy_mut(&mut sa, &a, q_modops);

    // sample e
    // RandomGaussianDist::random_fill(rng, &q_moduli_chain, &mut
    // c_out.c_partq_mut()[0]); foward_lazy(&mut c_out.c_partq_mut()[0],
    // q_nttops);

    // a*s + e + m
    add_lazy_mut(&mut c_out.c_partq_mut()[0], &sa, q_modops);
    add_lazy_mut(&mut c_out.c_partq_mut()[0], m_eval_lazy, q_modops);

    *c_out.representation_mut() = Representation::Evaluation;
    *c_out.level_mut() = level;
    *c_out.is_lazy_mut() = true;
}

/// Decrypt RLWE ciphertext using secret
///
/// Output message polynomial `m_eval_lazty` is in Evaluation representation
/// and has lazy coefficients
pub fn secret_key_decryption<
    Scalar: UnsignedInteger,
    C: RlweCiphertext<Scalar = Scalar>,
    S: SecretKey,
    P: CkksEncDecParameters<Scalar = Scalar>,
>(
    c: &C,
    m_eval_lazy: &mut C::Poly,
    secret: &S,
    params: &P,
) where
    C::Poly: Clone + TryConvertFrom<[S::Scalar], Parameters = [Scalar]>,
{
    let level = c.level();
    let q_moduli_chain = params.q_moduli_chain_at_level(level);

    debug_assert!(
        c.representation() == Representation::Evaluation,
        "Expected ciphertext to have Evaluation representation"
    );
    debug_assert!(
        c.c_partq().len() >= 2,
        "Ciphertext must of atleast degree 2"
    );
    debug_assert!(
        m_eval_lazy.dimension() == c.c_partq()[0].dimension(),
        "m_eval has incorrect dimensions. Expected {:?} but has {:?}",
        c.c_partq()[0].dimension(),
        m_eval_lazy.dimension(),
    );

    let q_modops = params.q_modops_at_level(level);
    let q_nttops = params.q_nttops_at_level(level);

    *m_eval_lazy = c.c_partq()[0].clone();

    let mut s = C::Poly::try_convert_from(secret.values(), q_moduli_chain);
    foward_lazy(&mut s, q_nttops);
    let s_clone = s.clone();

    for i in 1..c.c_partq().len() {
        let mut p = c.c_partq()[i].clone();
        mul_lazy_mut(&mut p, &s, q_modops);
        add_lazy_mut(m_eval_lazy, &p, q_modops);

        mul_lazy_mut(&mut s, &s_clone, q_modops);
    }
}

// specialFFT
// encoding
// decoeding
// how do you implement random round ?
// check encoding & decoding error

#[cfg(test)]
mod tests {
    use crate::{
        core_crypto::{
            modulus::NativeModulusBackend, ntt::NativeNTTBackend, prime::generate_primes_vec, ring,
        },
        parameters::Parameters,
        utils::{moduli_chain_to_biguint, psi_powers},
    };

    use super::*;
    use itertools::Itertools;
    use num_bigint::BigUint;
    use num_complex::{Complex, Complex64, ComplexDistribution};
    use num_traits::{zero, Zero};
    use rand::{distributions::Uniform, thread_rng, Rng};

    #[test]
    fn special_fft_round_trip() {
        // generate random complex values
        // create rot_group
        let m = 32;
        let n = m >> 1;
        let l = n >> 1;

        let mut a = 1usize;
        let mut rot_group = vec![];
        for _ in 0..l {
            rot_group.push(a);
            a = (a * 5) % m;
        }

        // vec of length l with random complex values
        let reals = Uniform::new(0.0, 100.0);
        let imags = Uniform::new(0.0, 100.0);
        let complex_distr = ComplexDistribution::new(reals, imags);
        let mut values = thread_rng()
            .sample_iter(complex_distr)
            .take(l)
            .collect_vec();

        let psi_powers = psi_powers(m as u32);

        let values_clone = values.clone();

        special_inv_fft(&mut values, &psi_powers, &rot_group);
        special_fft(&mut values, &psi_powers, &rot_group);

        // values after round trip will not equal values before due to errors
        // accumulated from floating operations. But the difference should be
        // negligible. I don't know yet how to quantify the difference, hence
        // can't hardcode a value. Therefore, we pretty print here
        // difference and check ourselves that the error is negligible
        // TODO(Jay): Pretty print the difference
        // dbg!(values_clone);
        // dbg!(values);
    }

    struct CkksParameters {
        delta: f64,
        psi_powers: Vec<Complex<f64>>,
        rot_group: Vec<usize>,
        ring_size: usize,
        q_moduli_chain: Vec<u64>,
        bigq: BigUint,
    }
    impl CkksParameters {
        fn new(ring_size: usize, q_moduli_chain: Vec<u64>, delta: f64) -> Self {
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
            let bigq = moduli_chain_to_biguint(&q_moduli_chain);

            CkksParameters {
                delta,
                psi_powers,
                rot_group,
                ring_size,
                q_moduli_chain,
                bigq,
            }
        }
    }

    impl Parameters for CkksParameters {
        type Scalar = u64;
    }

    impl CkksEncDecParameters for CkksParameters {
        type BU = BigUint;
        type Complex = Complex<f64>;
        type F = f64;
        type ModOp = NativeModulusBackend;
        type NttOp = NativeNTTBackend;

        fn bigq_at_level(&self, level: usize) -> &Self::BU {
            &self.bigq
        }
        fn delta(&self) -> Self::F {
            self.delta
        }
        fn psi_powers(&self) -> &[Self::Complex] {
            &self.psi_powers
        }
        fn q_moduli_chain_at_level(&self, level: usize) -> &[Self::Scalar] {
            &self.q_moduli_chain
        }
        fn ring_size(&self) -> usize {
            self.ring_size
        }
        fn rot_group(&self) -> &[usize] {
            &self.rot_group
        }

        fn q_modops_at_level(&self, level: usize) -> &[Self::ModOp] {
            todo!()
        }
        fn q_nttops_at_level(&self, level: usize) -> &[Self::NttOp] {
            todo!()
        }
    }

    #[test]
    fn encoding_decoding_works() {
        let ring_size = 1 << 4;
        let q_moduli_chain = generate_primes_vec(&[50, 50, 50], ring_size, &[]);
        let params = CkksParameters::new(ring_size, q_moduli_chain, 2.0_f64.powi(40i32));

        // vec of length l with random complex values
        let reals = Uniform::new(0.0, 100.0);
        let imags = Uniform::new(0.0, 100.0);
        let complex_distr = ComplexDistribution::new(reals, imags);
        let values = thread_rng()
            .sample_iter(complex_distr)
            .take(params.ring_size() >> 1)
            .collect_vec();

        let level = 0;
        let mut p = <Vec<Vec<u64>> as Matrix>::zeros(
            params.q_moduli_chain_at_level(0).len(),
            params.ring_size(),
        );

        simd_encode(&mut p, &values, &params, level, params.delta());

        let mut m_out = vec![Complex::<f64>::zero(); ring_size >> 1];
        simd_decode(&p, &params, level, params.delta(), &mut m_out);

        // dbg!(values);
        // dbg!(m_out);
    }
}
