use num_traits::{
    AsPrimitive, CheckedShl, CheckedShr, Num, NumAssign, NumOps, One, PrimInt, ToPrimitive,
    WrappingAdd, WrappingMul, WrappingShl, WrappingShr, WrappingSub, Zero,
};
use std::{
    fmt::{Debug, Display},
    ops::{
        Add, AddAssign, Div, DivAssign, Mul, MulAssign, Rem, RemAssign, Shl, Shr, Sub, SubAssign,
    },
};

use crate::core_crypto::modulus::MontgomeryScalar;

use super::{NumericConstants, UnsignedInteger};

// TODO(Jay): Switch to macro that implements UnsignedInteger for any wrapper
// type

macro_rules! impl_num_ops_with_generic_scalar {
    ($trait:ident, $scalar: ident, $typ:ident, $op_name:ident, $op:tt) => {
        impl<$scalar: UnsignedInteger> $trait for $typ<$scalar> {
            type Output = Self;

            fn $op_name(self, other: Self) -> Self {
                Self(self.0 $op other.0)
            }
        }
    }
}

macro_rules! wrapping_impl_with_generic_scalar {
    ($trait_name:ident, $scalar: ident,$method:ident, $t:ident) => {
        impl<$scalar: UnsignedInteger> $trait_name for $t<$scalar> {
            #[inline]
            fn $method(&self, v: &Self) -> Self {
                Self(self.0.$method(&v.0))
            }
        }
    };

    ($trait_name:ident, $method:ident, $t:ty, $rhs:ty) => {
        impl<$scalar: UnsignedInteger> $trait_name<$rhs> for $t<$scalar> {
            #[inline]
            fn $method(&self, v: &$rhs) -> Self {
                <$t>::$method(*self, *v)
            }
        }
    };
}

impl<Scalar: UnsignedInteger> Shl<u32> for MontgomeryScalar<Scalar> {
    type Output = Self;
    fn shl(self, rhs: u32) -> Self::Output {
        Self(self.0 << rhs)
    }
}

impl<Scalar: UnsignedInteger> Shl<usize> for MontgomeryScalar<Scalar> {
    type Output = Self;
    fn shl(self, rhs: usize) -> Self::Output {
        Self(self.0 << rhs)
    }
}

impl<Scalar: UnsignedInteger> Shr<u32> for MontgomeryScalar<Scalar> {
    type Output = Self;
    fn shr(self, rhs: u32) -> Self::Output {
        Self(self.0 >> rhs)
    }
}

impl<Scalar: UnsignedInteger> Shr<usize> for MontgomeryScalar<Scalar> {
    type Output = Self;
    fn shr(self, rhs: usize) -> Self::Output {
        Self(self.0 >> rhs)
    }
}

impl<Scalar: UnsignedInteger> CheckedShr for MontgomeryScalar<Scalar> {
    fn checked_shr(&self, rhs: u32) -> Option<Self> {
        Some(*self >> rhs)
    }
}

impl<Scalar: UnsignedInteger> CheckedShl for MontgomeryScalar<Scalar> {
    fn checked_shl(&self, rhs: u32) -> Option<Self> {
        Some(*self << rhs)
    }
}

impl<Scalar: UnsignedInteger> WrappingShr for MontgomeryScalar<Scalar> {
    fn wrapping_shr(&self, rhs: u32) -> Self {
        Self(self.0 >> rhs)
    }
}

impl<Scalar: UnsignedInteger> WrappingShl for MontgomeryScalar<Scalar> {
    fn wrapping_shl(&self, rhs: u32) -> Self {
        Self(self.0 << rhs)
    }
}

impl<Scalar: UnsignedInteger> AddAssign for MontgomeryScalar<Scalar> {
    fn add_assign(&mut self, rhs: Self) {
        self.0 = self.0 + rhs.0
    }
}

impl<Scalar: UnsignedInteger> SubAssign for MontgomeryScalar<Scalar> {
    fn sub_assign(&mut self, rhs: Self) {
        self.0 = self.0 - rhs.0
    }
}

impl<Scalar: UnsignedInteger> MulAssign for MontgomeryScalar<Scalar> {
    fn mul_assign(&mut self, rhs: Self) {
        self.0 = self.0 - rhs.0
    }
}

impl<Scalar: UnsignedInteger> Rem for MontgomeryScalar<Scalar> {
    type Output = Self;
    fn rem(self, rhs: Self) -> Self::Output {
        Self(self.0 % rhs.0)
    }
}

impl<Scalar: UnsignedInteger> RemAssign for MontgomeryScalar<Scalar> {
    fn rem_assign(&mut self, rhs: Self) {
        self.0 = self.0 % rhs.0
    }
}

impl<Scalar: UnsignedInteger> Div for MontgomeryScalar<Scalar> {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        MontgomeryScalar(self.0 / rhs.0)
    }
}

impl<Scalar: UnsignedInteger> DivAssign for MontgomeryScalar<Scalar> {
    fn div_assign(&mut self, rhs: Self) {
        self.0 = self.0 / rhs.0
    }
}

impl<Scalar: UnsignedInteger> PartialEq for MontgomeryScalar<Scalar> {
    fn eq(&self, other: &Self) -> bool {
        self.0.eq(&other.0)
    }
}

impl<Scalar: UnsignedInteger> PartialOrd for MontgomeryScalar<Scalar> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl<Scalar: UnsignedInteger> One for MontgomeryScalar<Scalar> {
    fn one() -> Self {
        Self(Scalar::one())
    }
}

impl<Scalar: UnsignedInteger> Zero for MontgomeryScalar<Scalar> {
    fn zero() -> Self {
        Self(Scalar::zero())
    }

    fn is_zero(&self) -> bool {
        self.0 == Scalar::zero()
    }
}

impl<Scalar: UnsignedInteger> Num for MontgomeryScalar<Scalar> {
    type FromStrRadixErr = <Scalar as num_traits::Num>::FromStrRadixErr;
    fn from_str_radix(str: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        Scalar::from_str_radix(str, radix).map(|v| MontgomeryScalar(v))
    }
}

impl_num_ops_with_generic_scalar!(Add, Scalar, MontgomeryScalar, add, +);
impl_num_ops_with_generic_scalar!(Sub, Scalar, MontgomeryScalar, sub, -);
impl_num_ops_with_generic_scalar!(Mul, Scalar, MontgomeryScalar, mul, *);
wrapping_impl_with_generic_scalar!(WrappingAdd, Scalar, wrapping_add, MontgomeryScalar);
wrapping_impl_with_generic_scalar!(WrappingSub, Scalar, wrapping_sub, MontgomeryScalar);
wrapping_impl_with_generic_scalar!(WrappingMul, Scalar, wrapping_mul, MontgomeryScalar);

impl<Scalar: UnsignedInteger> UnsignedInteger for MontgomeryScalar<Scalar> {}

impl<Scalar: UnsignedInteger> NumericConstants for MontgomeryScalar<Scalar> {
    const BITS: u32 = Scalar::BITS;
    const MAX: Self = MontgomeryScalar(Scalar::MAX);
}
