use super::num::UnsignedInteger;

mod barrett;
mod montgomery;
mod native_backend;
mod shoup;

pub use barrett::BarrettBackend;
pub use montgomery::{MontgomeryBackend, MontgomeryBackendConfig, MontgomeryScalar};
pub use native_backend::NativeModulusBackend;
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
