use astro_float::BigFloat;
use num_bigint::{BigInt, BigUint};
use num_complex::{Complex, ComplexFloat};
use num_traits::{
    AsPrimitive, CheckedShl, CheckedShr, Float, MulAddAssign, Num, NumAssign, NumCast, NumOps,
    NumRef, One, Pow, PrimInt, RefNum, ToPrimitive, WrappingAdd, WrappingMul, WrappingShl,
    WrappingShr, WrappingSub, Zero,
};
use std::{
    fmt::{Debug, Display},
    ops::{Add, AddAssign, Div, Mul, Sub, SubAssign},
};

pub(crate) mod big_float;
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

// pub trait RefOps: for<'r> Div<&'r T, Output = Self> {}

// TODO(Jay): How to add trait bound for impl Add<&T, Output = T> for {}. Doing
// this remove unecessary copy bounds in ckks/ops.rs
pub trait ComplexNumber<T>: NumRef
where
    T: BFloat,
{
    fn new(re: T, im: T) -> Self;
    fn nth_root(n: u32) -> Self;
    fn powu(&self, exp: u32) -> Self;
    fn re(&self) -> &T;
    fn img(&self) -> &T;
}

pub trait BUInt: Num + PartialOrd {}

pub trait BInt: Num + PartialOrd {}

pub trait BFloat: Num + PartialOrd {
    fn abs(&self) -> Self;

    /// BFloat log2 may be quite expensive
    fn log2(&self) -> Self;

    fn min_value() -> Self;

    fn max_value() -> Self;
}

// TODO(Jay): CastToZp<To> is unecessary since it is only used on simd_encode
// where it can easily be replaced with TryConvertFrom<[f64]> to M
pub trait CastToZp<To> {
    fn cast(&self, q: &To) -> To;
}

pub trait ModInverse {
    fn mod_inverse(&self, q: Self) -> Self;
}

impl ComplexNumber<f64> for Complex<f64> {
    fn new(re: f64, im: f64) -> Self {
        Complex::new(re, im)
    }
    fn nth_root(n: u32) -> Self {
        Complex::<f64>::from_polar(f64::one(), (2.0 * std::f64::consts::PI) / (n as f64))
    }
    fn powu(&self, exp: u32) -> Self {
        self.powu(exp)
    }
    fn re(&self) -> &f64 {
        &self.re
    }
    fn img(&self) -> &f64 {
        &self.im
    }
}

impl BFloat for f64 {
    fn abs(&self) -> Self {
        f64::abs(*self)
    }

    fn log2(&self) -> Self {
        f64::log2(*self)
    }

    fn max_value() -> Self {
        f64::MAX
    }

    fn min_value() -> Self {
        f64::MIN
    }
}

// impl BFloat for BigFloat {
//     const TWICE_PI: Self = todo!();
// }

impl BUInt for BigUint {}

impl CastToZp<BigUint> for f64 {
    fn cast(&self, q: &BigUint) -> BigUint {
        let v = self.round();
        if v < 0.0 {
            q - (v.abs().to_u64().unwrap() % q)
        } else {
            v.to_u64().unwrap() % q
        }
    }
}
