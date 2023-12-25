use std::fmt::Display;

use num_traits::{
    AsPrimitive, CheckedShl, CheckedShr, Num, NumAssign, PrimInt, ToPrimitive, WrappingShl,
    WrappingShr,
};

pub trait UnsignedInteger:
    PrimInt
    + Num
    + NumAssign
    + WrappingShl
    + CheckedShl
    + WrappingShr
    + CheckedShr
    + PartialOrd
    + Display
{
}

impl UnsignedInteger for u64 {}
impl UnsignedInteger for u128 {}
