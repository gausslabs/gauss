use core::fmt;
use std::{
    ops::{Add, Div, Mul, Neg, Rem, Sub},
    sync::OnceLock,
};

use astro_float::{BigFloat as AstroBFloat, Consts};
use num_bigint::{BigInt, BigUint};
use num_complex::Complex;
use num_traits::{Float, Num, NumOps, NumRef, One, Zero};

use super::{BFloat, CastToZp, ComplexNumber};

#[derive(Clone, Debug)]
pub struct BigFloat(AstroBFloat);

const PRECISION: usize = 256usize;
const ROUNDING_MODE: astro_float::RoundingMode = astro_float::RoundingMode::None;

fn astro_one() -> &'static AstroBFloat {
    static V: OnceLock<AstroBFloat> = OnceLock::new();
    V.get_or_init(|| AstroBFloat::from_u8(1, PRECISION))
}

fn astro_zero() -> &'static AstroBFloat {
    static V: OnceLock<AstroBFloat> = OnceLock::new();
    V.get_or_init(|| AstroBFloat::from_u8(0, PRECISION))
}

fn astro_twice_pi() -> &'static AstroBFloat {
    let mut consts = Consts::new().unwrap();
    let pi = consts.pi(PRECISION, ROUNDING_MODE);
    let twice_pi = pi.mul(
        &AstroBFloat::from_u8(2, PRECISION),
        PRECISION,
        ROUNDING_MODE,
    );
    static V: OnceLock<AstroBFloat> = OnceLock::new();
    V.get_or_init(|| twice_pi)
}

impl PartialEq<BigFloat> for BigFloat {
    fn eq(&self, other: &BigFloat) -> bool {
        self.0.eq(&other.0)
    }

    fn ne(&self, other: &BigFloat) -> bool {
        self.0.ne(&other.0)
    }
}

macro_rules! impl_ops {
    ($(impl $Impl:ident<$Other:ty> for $Self:ty, $method:ident;)*) => {$(
        impl $Impl<$Other> for $Self {
            type Output = BigFloat;

            #[inline]
            fn $method(self, other: $Other) -> BigFloat {
                BigFloat(self.0.$method(&other.0, PRECISION, ROUNDING_MODE))
            }
        }
    )*}
}

impl_ops! {
    impl Mul<BigFloat> for BigFloat, mul;
    impl Mul<BigFloat> for &BigFloat, mul;
    impl Mul<&BigFloat> for BigFloat, mul;
    impl Mul<&BigFloat> for &BigFloat, mul;

    impl Add<BigFloat> for BigFloat, add;
    impl Add<BigFloat> for &BigFloat, add;
    impl Add<&BigFloat> for BigFloat, add;
    impl Add<&BigFloat> for &BigFloat, add;

    impl Sub<BigFloat> for BigFloat, sub;
    impl Sub<BigFloat> for &BigFloat, sub;
    impl Sub<&BigFloat> for BigFloat, sub;
    impl Sub<&BigFloat> for &BigFloat, sub;

    impl Div<BigFloat> for BigFloat, div;
    impl Div<BigFloat> for &BigFloat, div;
    impl Div<&BigFloat> for BigFloat, div;
    impl Div<&BigFloat> for &BigFloat, div;
}

macro_rules! impl_rem {
    ($(impl Rem<$Other:ty> for $Self:ty;)*) => {$(
        impl Rem<$Other> for $Self {
            type Output = BigFloat;

            #[inline]
            fn rem(self, other: $Other) -> BigFloat {
                BigFloat(self.0.rem(&other.0))
            }
        }
    )*}
}

impl_rem! {
    impl Rem<BigFloat> for BigFloat;
    impl Rem<BigFloat> for &BigFloat;
    impl Rem<&BigFloat> for BigFloat;
    impl Rem<&BigFloat> for &BigFloat;
}

impl Num for BigFloat {
    type FromStrRadixErr = astro_float::Error;

    fn from_str_radix(str: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        todo!()
    }
}

impl One for BigFloat {
    fn is_one(&self) -> bool
    where
        Self: PartialEq,
    {
        self.0.eq(astro_one())
    }

    fn one() -> Self {
        BigFloat(astro_one().clone())
    }
}

impl Zero for BigFloat {
    fn is_zero(&self) -> bool {
        self.0.eq(astro_zero())
    }

    fn zero() -> Self {
        BigFloat(astro_zero().clone())
    }
}

impl Neg for &BigFloat {
    type Output = BigFloat;
    fn neg(self) -> Self::Output {
        BigFloat((&self.0).neg())
    }
}

impl Neg for BigFloat {
    type Output = BigFloat;
    fn neg(self) -> Self::Output {
        BigFloat((self.0).neg())
    }
}

impl PartialOrd for BigFloat {
    fn ge(&self, other: &Self) -> bool {
        self.0.ge(&other.0)
    }

    fn gt(&self, other: &Self) -> bool {
        self.0.gt(&other.0)
    }

    fn le(&self, other: &Self) -> bool {
        self.0.le(&other.0)
    }

    fn lt(&self, other: &Self) -> bool {
        self.0.lt(&other.0)
    }

    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl BFloat for BigFloat {
    fn abs(&self) -> Self {
        BigFloat(self.0.abs())
    }

    fn log2(&self) -> Self {
        let mut consts = Consts::new().unwrap();
        BigFloat(self.0.log2(PRECISION, ROUNDING_MODE, &mut consts))
    }

    fn max_value() -> Self {
        BigFloat(AstroBFloat::max_value(PRECISION))
    }

    fn min_value() -> Self {
        BigFloat(AstroBFloat::min_value(PRECISION))
    }
}

impl ComplexNumber<BigFloat> for Complex<BigFloat> {
    fn img(&self) -> &BigFloat {
        &self.im
    }

    fn re(&self) -> &BigFloat {
        &self.re
    }

    fn new(re: BigFloat, im: BigFloat) -> Self {
        Complex::new(re, im)
    }

    fn nth_root(n: u32) -> Self {
        let twice_pi = astro_twice_pi();
        let theta = twice_pi.div(
            &AstroBFloat::from_u32(n, PRECISION),
            PRECISION,
            ROUNDING_MODE,
        );

        let mut consts = Consts::new().unwrap();
        let re = BigFloat(theta.cos(PRECISION, ROUNDING_MODE, &mut consts));
        let im = BigFloat(theta.sin(PRECISION, ROUNDING_MODE, &mut consts));

        Complex::<BigFloat>::new(re, im)
    }

    fn powu(&self, exp: u32) -> Self {
        self.powu(exp)
    }
}

impl From<f64> for BigFloat {
    fn from(value: f64) -> Self {
        BigFloat(AstroBFloat::from_f64(value, PRECISION))
    }
}

impl From<u32> for BigFloat {
    fn from(value: u32) -> Self {
        BigFloat(AstroBFloat::from_u32(value, PRECISION))
    }
}

impl From<BigInt> for BigFloat {
    fn from(value: BigInt) -> Self {
        let mut consts = Consts::new().unwrap();
        let (sign, digits) = value.to_radix_be(16);
        let sign = if sign == num_bigint::Sign::Minus {
            astro_float::Sign::Neg
        } else {
            astro_float::Sign::Pos
        };
        BigFloat(AstroBFloat::convert_from_radix(
            sign,
            &digits,
            1,
            astro_float::Radix::Hex,
            PRECISION,
            ROUNDING_MODE,
            &mut consts,
        ))
    }
}

impl CastToZp<BigUint> for BigFloat {
    fn cast(&self, q: &BigUint) -> BigUint {
        let mut consts = Consts::new().unwrap();

        let (sign, radix_repr, ex) = self
            .0
            .round(PRECISION, ROUNDING_MODE)
            .convert_to_radix(astro_float::Radix::Hex, ROUNDING_MODE, &mut consts)
            .unwrap();

        assert!(sign.is_positive());

        BigUint::from_radix_le(&radix_repr, 16).unwrap() % q
    }
}

impl std::fmt::Display for BigFloat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        <AstroBFloat as std::fmt::Display>::fmt(&self.0, f)
    }
}
