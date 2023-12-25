use num_traits::AsPrimitive;

use crate::{core_crypto::num::UnsignedInteger, utils::FastModularInverse};

struct MontgomeryScalar<Scalar: UnsignedInteger>(Scalar);

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

    /// n^{-1} (mod r) where n is the original modulus
    fn n_inverse_modr(&self) -> Scalar;

    fn mont_sub(&self, a: Scalar, b: Scalar) -> Scalar {
        debug_assert!(a < self.modulus(), "Input {a} >= {}", self.modulus());
        debug_assert!(b < self.modulus(), "Input {b} >= {}", self.modulus());

        if a >= b {
            a - b
        } else {
            (a + self.modulus()) - b
        }
    }

    fn mont_add(&self, a: Scalar, b: Scalar) -> Scalar {
        debug_assert!(a < self.modulus(), "Input {a} >= {}", self.modulus());
        debug_assert!(b < self.modulus(), "Input {b} >= {}", self.modulus());

        let mut c = a + b;
        if c >= self.modulus() {
            c -= self.modulus();
        }
        c
    }

    /// r^2 (mod n)
    fn r_square_modn(&self) -> Scalar;

    /// For a & b in Montgomery space, outputs `o = ab*r^{-1} (mod n)`,
    /// where `n` is the original modulus and `o \in [0,n)`
    ///
    /// - [Reference](https://en.algorithmica.org/hpc/number-theory/montgomery/)
    fn mont_mul(&self, a: Scalar, b: Scalar) -> Scalar {
        let r = Scalar::BITS as usize;
        let n_inv = self.n_inverse_modr();
        let n = self.modulus();

        let ab = a.as_() * b.as_();

        // ab * n' (mod r)
        let q = ab.as_().wrapping_mul(&n_inv);

        // (q * n) / r
        let m = ((q.as_() * n.as_()) >> r).as_();

        // ab / r
        let ab_hi = (ab >> r).as_();

        if ab_hi < m {
            ab_hi + n - m
        } else {
            ab_hi - m
        }
    }

    /// For a & b in Montgomery space, outputs `o = ab*r^{-1} (mod n)`,
    /// where `n` is the original modulus and `o \in [0,2n-2)`
    fn mont_mul_lazy(&self, a: Scalar, b: Scalar) -> Scalar {
        let r = Scalar::BITS as usize;
        let n_inv = self.n_inverse_modr();
        let n = self.modulus();

        let ab = a.as_() * b.as_();

        // ab * n' (mod r)
        let q = ab.as_().wrapping_mul(&n_inv);

        // (q * n) / r
        let m = ((q.as_() * n.as_()) >> r).as_();

        // ab / r
        let ab_hi = (ab >> r).as_();

        ab_hi + n - m
    }

    /// Transforms input scalar from normal space to montogmery space.
    /// Calculates and returns a * r (mod n), where n is original modulus.
    ///
    /// We precompute r^2 (mod n) and calculate
    /// (a * r^2)r^{-1} = a *r (mod n) instead
    fn normal_to_mont_space(&self, a: Scalar) -> Scalar {
        self.mont_mul(a, self.r_square_modn())
    }

    /// Transforms input from montogmery space to normal space.
    /// Calculates a * r^{-1} (mod n).
    fn mont_to_normal(&self, a: Scalar) -> Scalar {
        // since (a * 1) r^{-1} (mod n) = a * r^{-1} (mod n),
        // one can transform `a` from montgomery space to normal
        // space via calling `mont_mul` routine with `1` as second
        // operand.
        self.mont_mul(a, Scalar::one())
    }
}
