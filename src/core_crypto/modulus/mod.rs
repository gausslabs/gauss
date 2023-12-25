use super::num::UnsignedInteger;

mod barrett;
mod montgomery;
mod native_backend;

pub use native_backend::NativeModulusBackend;

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
