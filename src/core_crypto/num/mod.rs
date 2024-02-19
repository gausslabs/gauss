use num_complex::Complex;
use num_traits::{
    AsPrimitive, CheckedShl, CheckedShr, MulAddAssign, Num, NumAssign, NumOps, One, Pow, PrimInt,
    ToPrimitive, WrappingAdd, WrappingMul, WrappingShl, WrappingShr, WrappingSub, Zero,
};
use std::{
    fmt::{Debug, Display},
    ops::{Add, AddAssign, Div, Mul, Sub, SubAssign},
};

pub(crate) mod impl_num;

pub trait UnsignedIntegerDoubled {}

pub trait NumericConstants {
    const BITS: u32;
    const MAX: Self;
}

pub trait UnsignedInteger:
    Num
    + NumAssign
    + WrappingMul
    + WrappingAdd
    + WrappingSub
    + WrappingShl
    + CheckedShl
    + WrappingShr
    + CheckedShr
    + PartialOrd
    + Display
    + NumericConstants
    + Debug
    + Copy
    + Clone
{
}

impl NumericConstants for u64 {
    const BITS: u32 = u64::BITS;
    const MAX: u64 = u64::MAX;
}

impl NumericConstants for u128 {
    const BITS: u32 = u128::BITS;
    const MAX: u128 = u128::MAX;
}

impl UnsignedInteger for u64 {}
impl UnsignedInteger for u128 {}

pub trait Float: num_traits::Float + From<u32> {
    const PI: Self;
}

impl Float for f64 {
    const PI: Self = std::f64::consts::PI * 2.0;
}

pub trait ComplexNumber<T>: Num + Div<T, Output = Self>
where
    T: Float,
{
    fn nth_root(n: u32) -> Self;
    fn powu(&self, exp: u32) -> Self;
    fn re(&self) -> T;
    fn img(&self) -> T;
}

impl<T: Debug + Float + Clone> ComplexNumber<T> for Complex<T> {
    fn nth_root(n: u32) -> Self {
        Complex::<T>::from_polar(T::one(), T::PI / <T as From<u32>>::from(n))
    }
    fn powu(&self, exp: u32) -> Self {
        self.powu(exp)
    }
    fn re(&self) -> T {
        self.re
    }
    fn img(&self) -> T {
        self.im
    }
}
