use super::{barrett::BarrettBackend, ModulusBackendConfig, ModulusVecBackend};
use itertools::izip;

pub struct NativeModulusBackend {
    modulus: u64,
    barrett_constant: u64,
    barrett_alpha: usize,
    modulus_bits: usize,
}

impl ModulusBackendConfig<u64> for NativeModulusBackend {
    fn initialise(modulus: u64) -> NativeModulusBackend {
        let (alpha, mu) = <NativeModulusBackend>::precompute_alpha_and_barrett_constant(modulus);

        NativeModulusBackend {
            modulus,
            barrett_alpha: alpha,
            barrett_constant: mu,
            modulus_bits: 64 - modulus.leading_zeros() as usize,
        }
    }
}

impl BarrettBackend<u64, u128> for NativeModulusBackend {
    #[inline]
    fn barrett_alpha(&self) -> usize {
        self.barrett_alpha
    }
    #[inline]
    fn modulus(&self) -> u64 {
        self.modulus
    }
    #[inline]
    fn modulus_bits(&self) -> usize {
        self.modulus_bits
    }
    #[inline]
    fn barrett_constant(&self) -> u64 {
        self.barrett_constant
    }
}

impl ModulusVecBackend<u64> for NativeModulusBackend {
    fn add_mod_vec(&self, a: &mut [u64], b: &[u64]) {
        izip!(a.iter_mut(), b.iter()).for_each(|(a0, b0)| {
            *a0 = self.add_mod_fast(*a0, *b0);
        })
    }

    fn sub_mod_vec(&self, a: &mut [u64], b: &[u64]) {
        izip!(a.iter_mut(), b.iter()).for_each(|(a0, b0)| {
            *a0 = self.sub_mod_fast(*a0, *b0);
        })
    }

    fn mul_mod_vec(&self, a: &mut [u64], b: &[u64]) {
        izip!(a.iter_mut(), b.iter()).for_each(|(a0, b0)| {
            *a0 = self.mul_mod_fast(*a0, *b0);
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{thread_rng, Rng};

    const PRIME_60_BITS: u64 = 1152921504606748673;

    const N: usize = 1000;

    #[test]
    fn native_modulus_backend_works() {
        let p = PRIME_60_BITS;
        let mut rng = thread_rng();
        let modulus_backend = NativeModulusBackend::initialise(p);
        for _ in 0..N {
            // Case when a,b < p
            let a = rng.gen::<u64>() % p;
            let b = rng.gen::<u64>() % p;

            let c = modulus_backend.mul_mod_fast(a, b);
            let c_expected = ((a as u128 * b as u128) % p as u128) as u64;
            assert_eq!(c, c_expected);

            let c = modulus_backend.add_mod_fast(a, b);
            let c_expected = (a + b) % p;
            assert_eq!(c, c_expected);

            let c = modulus_backend.sub_mod_fast(a, b);
            let c_expected = if a > b { a - b } else { p - (b - a) };
            assert_eq!(c, c_expected);
        }
    }
}
