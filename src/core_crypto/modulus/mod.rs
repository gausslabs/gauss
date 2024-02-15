use super::num::UnsignedInteger;

pub mod barrett;
pub mod montgomery;
pub mod native_backend;
mod shoup;

pub use barrett::BarrettBackend;
pub use montgomery::{MontgomeryBackend, MontgomeryBackendConfig, MontgomeryScalar};
pub use native_backend::NativeModulusBackend;
use rand::{CryptoRng, Rng, RngCore};
pub use shoup::ShoupRepresentationFq;

pub trait ModulusBackendConfig<Scalar> {
    fn initialise(modulus: Scalar) -> Self;
}

pub trait ModulusArithmeticBackend<Scalar: UnsignedInteger> {
    fn modulus(&self) -> Scalar;

    fn twice_modulus(&self) -> Scalar;

    fn neg_mod_fast(&self, a: Scalar) -> Scalar {
        debug_assert!(
            a <= self.modulus(),
            "Input {a} > (modulus){}",
            self.modulus()
        );

        self.modulus() - a
    }

    fn neg_lazy_mod_fast(&self, a: Scalar) -> Scalar {
        debug_assert!(
            a < self.twice_modulus(),
            "Input {a} > (2*modulus){}",
            self.twice_modulus()
        );

        self.twice_modulus() - a
    }

    fn add_mod_fast(&self, a: Scalar, b: Scalar) -> Scalar {
        debug_assert!(
            a <= self.modulus(),
            "Input a:{a} > (modulus){}",
            self.modulus()
        );
        debug_assert!(
            b <= self.modulus(),
            "Input b:{b} > (modulus){}",
            self.modulus()
        );

        let mut c = a + b;
        if c >= self.modulus() {
            c -= self.modulus();
        }
        c
    }

    /// Lazy modular addition of a,b \in [0, 2q).
    ///
    /// Output is in range [0, 2q)
    fn add_lazy_mod_fast(&self, a: Scalar, b: Scalar) -> Scalar {
        debug_assert!(
            a <= self.twice_modulus(),
            "Input a:{a} > (2*modulus){}",
            self.twice_modulus()
        );
        debug_assert!(
            b <= self.twice_modulus(),
            "Input b:{b} > (2*modulus){}",
            self.twice_modulus()
        );

        let twice_modulus = self.twice_modulus();

        let mut c = a + b;
        if c >= twice_modulus {
            c -= twice_modulus;
        }
        c
    }

    fn sub_mod_fast(&self, a: Scalar, b: Scalar) -> Scalar {
        debug_assert!(
            a <= self.modulus(),
            "Input a:{a} > (modulus){}",
            self.modulus()
        );
        debug_assert!(
            b <= self.modulus(),
            "Input b:{b} > (modulus){}",
            self.modulus()
        );

        if a >= b {
            a - b
        } else {
            (a + self.modulus()) - b
        }
    }

    fn sub_lazy_mod_fast(&self, a: Scalar, b: Scalar) -> Scalar {
        debug_assert!(
            a <= self.twice_modulus(),
            "Input a:{a} > (2*modulus){}",
            self.twice_modulus()
        );
        debug_assert!(
            b <= self.twice_modulus(),
            "Input b:{b} > (2*modulus){}",
            self.twice_modulus()
        );

        if a >= b {
            a - b
        } else {
            (a + self.twice_modulus()) - b
        }
    }
}

pub trait ModulusVecBackend<Scalar>
where
    Scalar: UnsignedInteger,
{
    /// Elements in input vector `a` are in range [0, Scalar::MAX].
    /// Elements are reduce inplace to rnage [0, q)
    fn reduce_vec(&self, a: &mut [Scalar]);
    fn neg_mod_vec(&self, a: &mut [Scalar]);
    fn add_mod_vec(&self, a: &mut [Scalar], b: &[Scalar]);
    fn sub_mod_vec(&self, a: &mut [Scalar], b: &[Scalar]);
    fn mul_mod_vec(&self, a: &mut [Scalar], b: &[Scalar]);
    fn scalar_mul_mod_vec(&self, a: &mut [Scalar], b: Scalar);

    /// Inplace reduce elemnts of vector `a` in range [0, 2q) to
    /// [0, 2q)
    fn reduce_from_lazy_vec(&self, a: &mut [Scalar]);

    /// Inplace modular multiplication a=a*b. Input a nad b are in range [0, 2q)
    /// and output a is in range [0, 2q]
    fn mul_lazy_mod_vec(&self, a: &mut [Scalar], b: &[Scalar]);
    /// Inplace modular addition a=a+b. Input a and b are in range [0, 2q) and
    /// output a is in range [0, 2q]
    fn add_lazy_mod_vec(&self, a: &mut [Scalar], b: &[Scalar]);
    /// Inplace modular subtraction a=a+b. Input a and b are in range [0, 2q)
    /// and output a is in range [0, 2q]
    fn sub_lazy_mod_vec(&self, a: &mut [Scalar], b: &[Scalar]);
    /// Inplace scalar modular multiplication a=a*b. Input a and b are in range
    /// [0, 2q) and output a is in range [0, 2q]
    fn scalar_mul_lazy_mod_vec(&self, a: &mut [Scalar], b: Scalar);
    fn neg_lazy_mod_vec(&self, a: &mut [Scalar]);
}
