use std::{
    clone,
    fmt::Debug,
    ops::{Div, Rem, Sub},
    process::Output,
};

use itertools::{izip, Itertools};
use num_traits::{Num, One};

use crate::{
    core_crypto::{
        matrix::{Matrix, MatrixMut, RowMut},
        num::{ComplexNumber, Float, UnsignedInteger},
        ring,
    },
    utils::bit_reverse_map,
};

// specialIFFT
fn special_inv_fft<F: Float, C: ComplexNumber<F> + Clone + Copy>(
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
        v.len() * 4 + 1 == psi_powers.len(),
        "psi_powers must have powers of psi for 0 <= j <= 4l, but its length is {}",
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
                let u = v[i + j] + v[i + j + lenh];
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
        .for_each(|a| *a = *a / <F as From<u32>>::from(v_len as u32));
}

fn special_fft<F: Float, C: ComplexNumber<F> + Copy + Clone>(
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
        v.len() * 4 + 1 == psi_powers.len(),
        "psi_powers must have powers of psi for 0 <= j <= 4l, but its length is {}",
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

/// Scale input float a by delta and randomly rounds delta_a. Maps delta_a from
/// signed interval to unsigned interval [0, big_q) where big_q is product of
/// smaller primes qi. Maps delta_a in unsigned interval to big_q's moduli chain
/// and returns
fn map_float_to_zq_moduli_chain<Scalar: UnsignedInteger, Uint, F: Float>(
    a: F,
    delta: F,
    big_q: &Uint,
    q_moduli_chain: &[Uint],
) -> Vec<Scalar>
where
    <Uint as TryFrom<F>>::Error: Debug,
    <Uint as TryFrom<Scalar>>::Error: Debug,
    <Scalar as TryFrom<Uint>>::Error: Debug,
    Uint: TryFrom<F> + TryFrom<Scalar> + Num + PartialOrd,
    for<'a> &'a Uint: Rem<&'a Uint, Output = Uint> + Sub<&'a Uint, Output = Uint>,
    Scalar: TryFrom<Uint>,
{
    let delta_a = a * delta;
    let delta_a = delta_a.round(); // TODO random round

    // convert signed to unsigned representation
    let delta_a_modq = if delta_a < F::zero() {
        big_q - &Uint::try_from(-delta).unwrap()
    } else {
        Uint::try_from(delta).unwrap()
    };

    q_moduli_chain
        .iter()
        .map(|qi| <Scalar as TryFrom<Uint>>::try_from(&delta_a_modq % qi).unwrap())
        .collect_vec()
}

pub fn simd_encode<
    Scalar: UnsignedInteger,
    Uint,
    F: Float,
    C: ComplexNumber<F> + Clone + Copy,
    MMut: MatrixMut<MatElement = Scalar>,
>(
    p: &mut MMut,
    m: &[C],
) where
    <Uint as TryFrom<F>>::Error: Debug,
    <Uint as TryFrom<Scalar>>::Error: Debug,
    <Scalar as TryFrom<Uint>>::Error: Debug,
    Uint: TryFrom<F> + TryFrom<Scalar> + Num + PartialOrd,
    Scalar: TryFrom<Uint>,
    for<'a> &'a Uint: Rem<&'a Uint, Output = Uint> + Sub<&'a Uint, Output = Uint>,
    <MMut as Matrix>::R: RowMut,
{
    let slots = 0usize;
    let ring_size = 0;
    debug_assert!(
        m.len() == slots,
        "Expected m vector length {slots} but is {}",
        m.len()
    );
    debug_assert!(
        ring_size == 2 * slots,
        "Expected ring_size = {} (2*slots), but is {}",
        2 * slots,
        ring_size
    );

    let psi_powers = vec![];
    let rot_group = vec![];

    let mut m = m.to_vec();
    special_inv_fft(&mut m, &psi_powers, &rot_group);

    let delta = F::zero();

    let q_moduli_chain: Vec<Scalar> = vec![];
    let q_moduli_chain_big: Vec<Uint> = q_moduli_chain
        .iter()
        .map(|v| (*v).try_into().unwrap())
        .collect_vec();
    let big_q = Uint::one();

    for ri in 0..ring_size >> 1 {
        izip!(
            p.get_col_iter_mut(ri),
            map_float_to_zq_moduli_chain::<Scalar, _, _>(
                m[ri].re(),
                delta,
                &big_q,
                &q_moduli_chain_big
            )
            .iter()
        )
        .for_each(|(to, from)| *to = *from);
    }

    for ri in ring_size >> 1..ring_size {
        izip!(
            p.get_col_iter_mut(ri),
            map_float_to_zq_moduli_chain::<Scalar, _, _>(
                m[ri].img(),
                delta,
                &big_q,
                &q_moduli_chain_big
            )
            .iter()
        )
        .for_each(|(to, from)| *to = *from);
    }
}

// specialFFT
// encoding
// decoeding
// how do you implement random round ?
// check encoding & decoding error

#[cfg(test)]
mod tests {
    use crate::utils::psi_powers;

    use super::*;
    use itertools::Itertools;
    use num_complex::{Complex64, ComplexDistribution};
    use rand::{distributions::Uniform, thread_rng, Rng};

    #[test]
    fn special_inv_fft_round_trip() {
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
    }
}
