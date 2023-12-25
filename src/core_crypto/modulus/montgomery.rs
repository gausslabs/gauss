use crate::core_crypto::num::UnsignedInteger;

pub trait MontgomeryBackend<Scalar, ScalarDoubled>
where
    Scalar: UnsignedInteger,
    ScalarDoubled: UnsignedInteger,
{

    

    fn transform(a: Scalar) {

    }
}
