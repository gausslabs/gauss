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

// TODO(Jay): How to add trait bound for impl Add<&T, Output = T> for {}. Doing
// this remove unecessary copy bounds in ckks/ops.rs
pub trait ComplexNumber<T>:
    Num + NumRef + for<'r> Div<&'r T, Output = Self> + for<'r> Mul<&'r T, Output = Self>
where
    T: BFloat,
{
    fn new(re: T, im: T) -> Self;
    fn nth_root(n: u32) -> Self;
    fn powu(&self, exp: u32) -> Self;
    fn re(&self) -> T;
    fn img(&self) -> T;
}

pub trait BUint: Num + NumRef + PartialOrd {}

pub trait BInt: Num + NumRef + PartialOrd {}

pub trait BFloat: Num + NumRef + Float + NumCast {
    const TWICE_PI: Self;
}

// TODO(Jay): CastToZp<To> is unecessary since it is only used on simd_encode
// where it can easily be replaced with TryConvertFrom<[f64]> to M
pub trait CastToZp<To> {
    fn cast(&self, q: &To) -> To;
}

pub trait ModInverse {
    fn mod_inverse(&self, q: Self) -> Self;
}

impl<T: BFloat + Clone> ComplexNumber<T> for Complex<T> {
    fn new(re: T, im: T) -> Self {
        Complex::new(re, im)
    }
    fn nth_root(n: u32) -> Self {
        Complex::<T>::from_polar(T::one(), T::TWICE_PI / T::from(n).unwrap())
    }
    fn powu(&self, exp: u32) -> Self {
        self.powu(exp)
    }
    fn re(&self) -> T {
        self.re.clone()
    }
    fn img(&self) -> T {
        self.im.clone()
    }
}

impl BFloat for f64 {
    const TWICE_PI: Self = std::f64::consts::PI * 2.0;
}

impl BUint for BigUint {}

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
