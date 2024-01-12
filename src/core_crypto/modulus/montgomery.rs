use crate::{
    core_crypto::modulus::BarrettBackend, core_crypto::num::UnsignedInteger,
    utils::FastModularInverse,
};
use itertools::izip;
use num_traits::AsPrimitive;

/// Wrapper around `MontgomeryScalar`s.
///
/// This is just to prevent cross-arithmatic between Scalars in Montgomery space and Scalars in Normal space.
#[derive(Debug, Clone, Copy)]
pub struct MontgomeryScalar<Scalar: UnsignedInteger>(pub(crate) Scalar);

impl<Scalar: UnsignedInteger> MontgomeryScalar<Scalar> {
    pub fn zero() -> Self {
        MontgomeryScalar(Scalar::zero())
    }
}

impl<Scalar: UnsignedInteger> std::fmt::Display for MontgomeryScalar<Scalar> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "MontgomeryScalar({})", self.0)?;
        Ok(())
    }
}

pub trait MontgomeryBackendConfig<Scalar, ScalarDoubled>
where
    Scalar: UnsignedInteger + AsPrimitive<ScalarDoubled> + FastModularInverse + 'static,
    ScalarDoubled: UnsignedInteger + AsPrimitive<Scalar> + 'static,
{
    fn initialise(modulus: Scalar) -> (Scalar, Scalar) {
        let n_inv_modr = Scalar::fast_inverse(modulus);

        // r^2 (mod n) where `r` is either 2^{32/64/128}
        let r_modn = (ScalarDoubled::one() << (Scalar::BITS as usize)) % modulus.as_();
        let r_square_modn = (r_modn * r_modn) % modulus.as_();
        (n_inv_modr, r_square_modn.as_())
    }
}

pub trait MontgomeryBackend<Scalar, ScalarDoubled>
where
    Scalar: UnsignedInteger + AsPrimitive<ScalarDoubled> + 'static,
    ScalarDoubled: UnsignedInteger + AsPrimitive<Scalar> + 'static,
{
    fn modulus(&self) -> Scalar;

    fn modulus_twice(&self) -> Scalar;

    /// n^{-1} (mod r) where n is the original modulus
    fn n_inverse_modr(&self) -> Scalar;

    // TODO (Bhargav)
    // Debug mod_inv value

    // montgomery reduction
    fn redc(&self, x: MontgomeryScalar<Scalar>) -> Scalar {
        let N = self.modulus().as_();
        let _N = self.n_inverse_modr();

        let x_hi = x.0.as_() >> Scalar::BITS;
        let m = (x.0.wrapping_mul(&_N)).as_();
        let mN = m * N;
        let mN_hi = mN >> Scalar::BITS;
        let mut tmp = x_hi + N;
        tmp = tmp - mN_hi;
        let result = if x_hi < mN_hi { tmp } else { x_hi - mN_hi };
        assert!(result < N);

        result.as_()
    }

    /// r^2 (mod n)
    fn r_square_modn(&self) -> Scalar;

    fn mont_sub(
        &self,
        a: MontgomeryScalar<Scalar>,
        b: MontgomeryScalar<Scalar>,
    ) -> MontgomeryScalar<Scalar> {
        debug_assert!(a.0 < self.modulus(), "Input {a} >= {}", self.modulus());
        debug_assert!(b.0 < self.modulus(), "Input {b} >= {}", self.modulus());

        let o = if a.0 >= b.0 {
            a.0 - b.0
        } else {
            (a.0 + self.modulus()) - b.0
        };

        MontgomeryScalar(o)
    }

    fn mont_add(
        &self,
        a: MontgomeryScalar<Scalar>,
        b: MontgomeryScalar<Scalar>,
    ) -> MontgomeryScalar<Scalar> {
        debug_assert!(a.0 < self.modulus(), "Input {a} >= {}", self.modulus());
        debug_assert!(b.0 < self.modulus(), "Input {b} >= {}", self.modulus());

        let mut c = a.0 + b.0;
        if c >= self.modulus() {
            c -= self.modulus();
        }
        MontgomeryScalar(c)
    }

    /// For a & b in Montgomery space, outputs `o = ab*r^{-1} (mod n)`,
    /// where `n` is the original modulus and `o \in [0,n)`
    ///
    /// In general, a*b should in range [0, qr). We assume q is smaller 2^{60}.
    /// If r = 2^64, \log{qr} = 64 + 60 = 124. Therefore a and b can be range [0, 2q)
    ///
    /// - [Reference 1](https://en.algorithmica.org/hpc/number-theory/montgomery/)
    /// - [Reference 2](https://jeffhurchalla.com/2022/04/29/optimized-montgomery-multiplication-with-smaller-modulus-sizes/)
    fn mont_mul(
        &self,
        a: MontgomeryScalar<Scalar>,
        b: MontgomeryScalar<Scalar>,
    ) -> MontgomeryScalar<Scalar> {
        let r = Scalar::BITS as usize;
        let n_inv = self.n_inverse_modr();
        let n = self.modulus();

        let ab = a.0.as_() * b.0.as_();

        // ab * n' (mod r)
        let q = ab.as_().wrapping_mul(&n_inv);

        // (q * n) / r
        let m = ((q.as_() * n.as_()) >> r).as_();

        // ab / r
        let ab_hi = (ab >> r).as_();

        let out = if ab_hi < m { ab_hi + n - m } else { ab_hi - m };
        MontgomeryScalar(out)
    }

    /// For a & b in Montgomery space, outputs `o = ab*r^{-1} (mod n)`,
    /// where `n` is the original modulus and `o \in [0,2n-2)`
    fn mont_mul_lazy(
        &self,
        a: MontgomeryScalar<Scalar>,
        b: MontgomeryScalar<Scalar>,
    ) -> MontgomeryScalar<Scalar> {
        let r = Scalar::BITS as usize;
        let n_inv = self.n_inverse_modr();
        let n = self.modulus();

        let ab = a.0.as_() * b.0.as_();

        // ab * n' (mod r)
        let q = ab.as_().wrapping_mul(&n_inv);

        // (q * n) / r
        let m = ((q.as_() * n.as_()) >> r).as_();

        // ab / r
        let ab_hi = (ab >> r).as_();

        MontgomeryScalar(ab_hi + n - m)
    }

    /// Caculates $\sum{a_i * b_i} \mod{q}$ in Montgomery space.
    ///
    /// Inputs should be in range [0, 2q) and output is in range [0, q)
    ///
    /// TODO (Jay):
    ///     - Apparently FMA can be optmised further. Check https://jeffhurchalla.com/2022/05/01/the-montgomery-multiply-accumulate/
    fn mont_fma(
        &self,
        a: &[MontgomeryScalar<Scalar>],
        b: &[MontgomeryScalar<Scalar>],
    ) -> MontgomeryScalar<Scalar> {
        debug_assert!(a.len() == b.len(), "Length of a and b are not equal");

        let q = self.modulus();
        let R = Scalar::from(Scalar::BITS).unwrap();
        let N = self.modulus();

        let mut sum = MontgomeryScalar::zero();
        izip!(a.iter(), b.iter()).for_each(|(&a0, &b0)| {
            let mul = ScalarDoubled::from(a0.0).unwrap() * ScalarDoubled::from(b0.0).unwrap();
            let u = Scalar::from(mul >> Scalar::BITS).unwrap();
            let v = Scalar::from(mul).unwrap();
            let c = sum.0;
            let w = if u < N - c { u + c } else { u + c - N };

            let s = w * R + v;
            sum = MontgomeryScalar(self.redc(MontgomeryScalar(s)));
        });
        sum
    }

    /// Transforms input scalar from normal space to montogmery space.
    /// Calculates and returns a * r (mod n), where n is original modulus.
    ///
    /// We precompute r^2 (mod n) and calculate
    /// (a * r^2)r^{-1} = a *r (mod n) instead
    ///
    /// Note that `mont_mul` accepts input in range [0, nr). Since (r^2 % n) < n,
    /// input `a` must be < r (i.e. 2^{64} if Scalae is u64)
    fn normal_to_mont_space(&self, a: Scalar) -> MontgomeryScalar<Scalar> {
        self.mont_mul(MontgomeryScalar(a), MontgomeryScalar(self.r_square_modn()))
    }

    /// Transforms input from montogmery space to normal space.
    /// Calculates a * r^{-1} (mod n).
    fn mont_to_normal(&self, a: MontgomeryScalar<Scalar>) -> Scalar {
        // since (a * 1) r^{-1} (mod n) = a * r^{-1} (mod n),
        // one can transform `a` from montgomery space to normal
        // space via calling `mont_mul` routine with `1` as second
        // operand.
        self.mont_mul(a, MontgomeryScalar(Scalar::one())).0
    }
}

#[cfg(test)]
mod tests {
    use super::{MontgomeryBackend, MontgomeryBackendConfig, MontgomeryScalar};
    use crate::{
        core_crypto::{
            modulus::NativeModulusBackend,
            num::{NumericConstants, UnsignedInteger},
        },
        utils::{mod_exponent, mod_inverse},
    };

    struct MontgomeryBackendTest {
        modulus: u64,
        n_inv_modr: u64,
        r_squared_modn: u64,
    }

    impl MontgomeryBackend<u64, u128> for MontgomeryBackendTest {
        fn modulus(&self) -> u64 {
            self.modulus
        }

        fn modulus_twice(&self) -> u64 {
            self.modulus.wrapping_mul(2)
        }

        fn n_inverse_modr(&self) -> u64 {
            self.n_inv_modr
        }

        fn r_square_modn(&self) -> u64 {
            self.r_squared_modn
        }
    }

    impl PartialEq for MontgomeryScalar<u64> {
        fn eq(&self, other: &Self) -> bool {
            self.0 == other.0
        }
    }

    #[test]
    fn test_mont_add() {
        let MODULUS = 11;

        let (n_inv_modr, r_squared_modn) =
            <NativeModulusBackend as MontgomeryBackendConfig<u64, u128>>::initialise(MODULUS);

        let mont_test = MontgomeryBackendTest {
            modulus: MODULUS,
            n_inv_modr,
            r_squared_modn,
        };

        (0..10).for_each(|x| {
            (0..10).for_each(|y| {
                let a = MontgomeryScalar(x);
                let b = MontgomeryScalar(y);
                assert_eq!(
                    mont_test.mont_add(a, b),
                    MontgomeryScalar((x + y) % MODULUS)
                );
            });
        });
    }

    #[test]
    fn test_mont_sub() {
        let MODULUS = 11;

        let (n_inv_modr, r_squared_modn) =
            <NativeModulusBackend as MontgomeryBackendConfig<u64, u128>>::initialise(MODULUS);

        let mont_test = MontgomeryBackendTest {
            modulus: MODULUS,
            n_inv_modr,
            r_squared_modn,
        };

        (0..10).for_each(|x| {
            (0..10).for_each(|y| {
                let a = MontgomeryScalar(x);
                let b = MontgomeryScalar(y);
                assert_eq!(
                    mont_test.mont_sub(a, b),
                    MontgomeryScalar((x + 10 * MODULUS).wrapping_sub(y) % MODULUS)
                );
            });
        });
    }

    #[test]
    fn test_mont_mul() {
        let MODULUS = 11;
        let r_inv = 9; // r = 2^64

        let (n_inv_modr, r_squared_modn) =
            <NativeModulusBackend as MontgomeryBackendConfig<u64, u128>>::initialise(MODULUS);

        let mont_test = MontgomeryBackendTest {
            modulus: MODULUS,
            n_inv_modr,
            r_squared_modn,
        };

        (0..10).for_each(|x| {
            (0..10).for_each(|y| {
                let a = MontgomeryScalar(x);
                let b = MontgomeryScalar(y);
                println!(
                    "{} * {} = {}, {}",
                    a,
                    b,
                    mont_test.mont_mul(a, b),
                    MontgomeryScalar((x * y * r_inv) % MODULUS)
                );
            });
        });
    }

    #[test]
    fn test_mont_redc() {
        let MODULUS = 11;

        let (n_inv_modr, r_squared_modn) =
            <NativeModulusBackend as MontgomeryBackendConfig<u64, u128>>::initialise(MODULUS);

        let mont_test = MontgomeryBackendTest {
            modulus: MODULUS,
            n_inv_modr,
            r_squared_modn,
        };

        let x = MontgomeryScalar(50_u64);
        let res = mont_test.redc(x);

        println!("{}", res);
    }

    #[test]
    fn test_inv_modr() {
        let MODULUS = 11;

        let (n_inv_modr, r_squared_modn) =
            <NativeModulusBackend as MontgomeryBackendConfig<u64, u128>>::initialise(MODULUS);

        println!("{}", MODULUS.wrapping_mul(n_inv_modr));
    }
}
