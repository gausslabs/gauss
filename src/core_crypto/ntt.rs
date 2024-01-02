use itertools::Itertools;
use rand::thread_rng;

use crate::utils::mod_inverse;

use super::{
    modulus::{BarrettBackend, ModulusBackendConfig, NativeModulusBackend},
    prime::find_primitive_root,
};

pub fn forward_butterly(x: &u64, y: &u64, w: &u64, q: &u64) -> (u64, u64) {
    let v = ((*y as u128 * *w as u128) % *q as u128) as u64;
    ((x + v) % q, (x + (q - v)) % q)
}

pub unsafe fn inverse_butterfly(x: *mut u64, y: *mut u64, w_inv: &u64, q: &u64) {
    let u = (*x + *y) % q;
    *y = (((*x + (q - *y)) as u128 * *w_inv as u128) % *q as u128) as u64;
    *x = u;
}

pub fn ntt(a: &mut [u64], psi: &[u64], q: u64) {
    debug_assert!(a.len() == psi.len());

    let n = a.len();
    let mut t = n;

    let mut m = 1;
    while m < n {
        t >>= 1;

        for i in 0..m {
            let j_1 = 2 * i * t;
            let j_2 = j_1 + t;
            let w = unsafe { psi.get_unchecked(m + i) };
            for j in j_1..j_2 {
                unsafe {
                    let (x, y) =
                        forward_butterly(a.get_unchecked(j), a.get_unchecked(j + t), w, &q);
                    *a.get_unchecked_mut(j) = x;
                    *a.get_unchecked_mut(j + t) = y;
                }
            }
        }

        m <<= 1;
    }
}

pub fn ntt_inv(a: &mut [u64], psi_inv: &[u64], n_inv: u64, q: u64) {
    debug_assert!(a.len() == psi_inv.len());

    let mut m = a.len();
    let mut t = 1;
    while m > 1 {
        let mut j_1: usize = 0;
        let h = m >> 1;
        for i in 0..h {
            let j_2 = j_1 + t;
            unsafe {
                let w_inv = psi_inv.get_unchecked(h + i);

                for j in j_1..j_2 {
                    let x = a.get_unchecked_mut(j) as *mut u64;
                    let y = a.get_unchecked_mut(j + t) as *mut u64;
                    inverse_butterfly(x, y, w_inv, &q)
                }
            }
            j_1 = j_1 + 2 * t;
        }
        t *= 2;
        m >>= 1;
    }

    a.iter_mut()
        .for_each(|a0| *a0 = ((*a0 as u128 * n_inv as u128) % q as u128) as u64);
}

pub struct NativeNTTBackend {
    q: u64,
    n: u64,
    n_inv: u64,
    psi_powers_bo: Box<[u64]>,
    psi_inv_powers_bo: Box<[u64]>,
}

impl NativeNTTBackend {
    pub fn new(q: u64, n: u64) -> NativeNTTBackend {
        // \psi = 2n^{th} primitive root of unity in F_q
        let mut rng = thread_rng();
        let psi =
            find_primitive_root(q, n * 2, &mut rng).expect("Unable to find 2n^th root of unity");
        let psi_inv = mod_inverse(psi, q);
        dbg!(psi, psi_inv);
        let modulus = NativeModulusBackend::initialise(q);

        let mut psi_powers = Vec::with_capacity(n as usize);
        let mut psi_inv_powers = Vec::with_capacity(n as usize);
        let mut running_psi = 1;
        let mut running_psi_inv = 1;
        for _ in 0..n {
            psi_powers.push(running_psi);
            psi_inv_powers.push(running_psi_inv);

            running_psi = modulus.mul_mod_fast(running_psi, psi);
            running_psi_inv = modulus.mul_mod_fast(running_psi_inv, psi_inv);
        }

        // powers stored in bit reversed order
        let mut psi_powers_bo = vec![0u64; n as usize];
        let mut psi_inv_powers_bo = vec![0u64; n as usize];
        let shift_by = n.leading_zeros() + 1;
        for i in 0..n as usize {
            // i in bit reversed order
            let bo_index = i.reverse_bits() >> shift_by;

            psi_powers_bo[bo_index] = psi_powers[i];
            psi_inv_powers_bo[bo_index] = psi_inv_powers[i];
        }

        // n^{-1} \mod{q}
        let n_inv = mod_inverse(n, q);

        NativeNTTBackend {
            q,
            n,
            n_inv,
            psi_powers_bo: psi_powers_bo.into_boxed_slice(),
            psi_inv_powers_bo: psi_inv_powers_bo.into_boxed_slice(),
        }
    }

    pub fn ntt(&self, a: &mut [u64]) {
        ntt(a, &self.psi_powers_bo, self.q);
    }

    pub fn ntt_inv(&self, a: &mut [u64]) {
        ntt_inv(a, &self.psi_inv_powers_bo, self.n_inv, self.q);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core_crypto::num::UnsignedInteger;
    use rand::{distributions::Uniform, Rng};

    const Q_60_BITS: u64 = 1152921504606748673;
    const N: u64 = 1 << 4;

    const K: usize = 128;

    fn random_vec_in_fq<T: UnsignedInteger + rand::distributions::uniform::SampleUniform>(
        size: usize,
        q: T,
    ) -> Vec<T> {
        let rng = thread_rng();
        rng.sample_iter(Uniform::new(T::zero(), q))
            .take(size)
            .collect_vec()
    }

    #[test]
    fn ntt_works() {
        let ntt_backend = NativeNTTBackend::new(Q_60_BITS, N);
        for _ in 0..K {
            let mut a = random_vec_in_fq(N as usize, Q_60_BITS);
            let a_clone = a.clone();

            ntt_backend.ntt(&mut a);
            assert_ne!(a, a_clone);
            ntt_backend.ntt_inv(&mut a);
            assert_eq!(a, a_clone);
        }
    }
}
