use crate::core_crypto::num::UnsignedInteger;

pub trait SecretKey {
    type Scalar;

    fn values(&self) -> &[Self::Scalar];
}

pub trait EncodedMessage {
    type Scalar: UnsignedInteger;
    fn encode() -> Self;
    fn value(&self) -> &[Self::Scalar];
}
