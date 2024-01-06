use super::num::UnsignedInteger;

mod barrett;
mod montgomery;
mod native_backend;
mod shoup;

pub use barrett::BarrettBackend;
pub use montgomery::{MontgomeryBackend, MontgomeryBackendConfig, MontgomeryScalar};
pub use native_backend::NativeModulusBackend;
use rand::{CryptoRng, Rng, RngCore};
pub use shoup::ShoupRepresentationFq;

pub trait ModulusBackendConfig<Scalar> {
    fn initialise(modulus: Scalar) -> Self;
}

pub trait ModulusVecBackend<Scalar>
where
    Scalar: UnsignedInteger,
{
    fn add_mod_vec(&self, a: &mut [Scalar], b: &[Scalar]);
    fn sub_mod_vec(&self, a: &mut [Scalar], b: &[Scalar]);
    fn mul_mod_vec(&self, a: &mut [Scalar], b: &[Scalar]);
}

pub trait ModulusRandomVecInDistGenerator<'a, Scalar, R>
where
    Scalar: UnsignedInteger,
    R: Rng + 'a,
{
    type IteratorUniform: Iterator<Item = Scalar>;
    type IteratorGaussian: Iterator<Item = Scalar>;

    fn random_vec_unifrom_dist_modulus(&self, size: usize, rng: &'a mut R)
        -> Self::IteratorUniform;
    fn random_vec_gaussian_dist_modulus(
        &self,
        std_dev: usize,
        size: usize,
        rng: &mut R,
    ) -> Self::IteratorGaussian;
}
