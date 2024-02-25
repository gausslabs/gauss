use itertools::{izip, Itertools};
use num_bigint::{BigInt, BigUint, ToBigUint};
use num_complex::Complex;
use num_traits::{FromPrimitive, One, ToBytes, ToPrimitive, Zero};

use crate::core_crypto::{
    modulus::{BarrettBackend, ModulusBackendConfig, NativeModulusBackend},
    num::{BFloat, ComplexNumber, UnsignedInteger},
};
use std::{
    mem,
    ops::{Add, Mul, Sub},
    sync::Arc,
};

use std::cell::RefCell;

pub(crate) mod convert;
#[cfg(test)]
pub(crate) mod test_utils;

struct Mdo<S> {
    f: S,
}

thread_local! {
    static X: RefCell<Arc<Mdo<u64>>> = panic!("!");
}

pub trait FastModularInverse {
    /// Calculates modular inverse of `a` in `Self`
    ///
    /// `a` must be odd
    ///
    /// - [Reference](https://jeffhurchalla.com/2022/04/25/a-faster-multiplicative-inverse-mod-a-power-of-2/)
    fn fast_inverse(a: Self) -> Self;
}

impl FastModularInverse for u64 {
    fn fast_inverse(a: Self) -> Self {
        assert!(a & 1 == 1, "Modulus inverse of {a} does not exit");
        let three_a = a.wrapping_mul(3);
        let x0 = three_a.wrapping_mul(three_a);
        let mut y = 1u64.wrapping_sub(x0.wrapping_mul(a));

        let x1 = x0.wrapping_mul(y.wrapping_add(1));
        y = y.wrapping_mul(y);

        let x2 = x1.wrapping_mul(y.wrapping_add(1));
        y = y.wrapping_mul(y);

        let x3 = x2.wrapping_mul(y.wrapping_add(1));
        y = y.wrapping_mul(y);

        let x4 = x3.wrapping_mul(y.wrapping_add(1));
        y = y.wrapping_mul(y);

        let x5 = x4.wrapping_mul(y.wrapping_add(1));
        y = y.wrapping_mul(y);

        let x6 = x5.wrapping_mul(y.wrapping_add(1));
        x6
    }
}

impl FastModularInverse for u32 {
    fn fast_inverse(a: Self) -> Self {
        assert!(a & 1 == 1, "Modulus inverse of {a} does not exit");

        let three_a = a.wrapping_mul(3);
        let x0 = three_a.wrapping_mul(three_a);
        let mut y = 1u32.wrapping_sub(x0.wrapping_mul(a));

        let x1 = x0.wrapping_mul(y.wrapping_add(1));
        y = y.wrapping_mul(y);

        let x2 = x1.wrapping_mul(y.wrapping_add(1));
        y = y.wrapping_mul(y);

        let x3 = x2.wrapping_mul(y.wrapping_add(1));
        y = y.wrapping_mul(y);

        let x4 = x3.wrapping_mul(y.wrapping_add(1));
        y = y.wrapping_mul(y);

        let x5 = x4.wrapping_mul(y.wrapping_add(1));
        x5
    }
}

/// Calculates a^n \mod{q} using binary exponentation
/// TODO (Jay): Add tests for modular expoents
pub fn mod_exponent(a: u64, mut n: u64, q: u64) -> u64 {
    // let mut a_prod = a;
    // let mut a_n = 1;

    // let modulus = NativeModulusBackend::initialise(q);

    // while n > 0 {
    //     if n & 1 == 1 {
    //         a_n = modulus.mul_mod_fast_lazy(a_prod, a_n);
    //     }
    //     a_prod = modulus.mul_mod_fast_lazy(a_prod, a_prod);

    //     n = n >> 1u32;
    // }

    //TODO(Jay): Implmentation above is buggy. Fix the implementation and add tests

    let expected_an = num_bigint_dig::BigUint::modpow(
        &num_bigint_dig::BigUint::from_u64(a).unwrap(),
        &num_bigint_dig::BigUint::from_u64(n).unwrap(),
        &num_bigint_dig::BigUint::from_u64(q).unwrap(),
    );

    // assert_eq!(expected_an.to_u64().unwrap(), a_n);

    expected_an.to_u64().unwrap()
}

/// Calculates modular inverse `a^{-1}` of `a` s.t. a * a^{-1} = 1 \mod{q}
///
/// TODO (Jay): Add tests and assert whether `q` is prime
pub fn mod_inverse(a: u64, q: u64) -> u64 {
    mod_exponent(a, q - 2, q)
}

pub fn mod_inverse_big_unit(a: &BigUint, b: &BigUint) -> BigUint {
    use num_bigint_dig::BigUint as BigUintDig;
    use num_bigint_dig::ModInverse;

    let a_dig = BigUintDig::from_bytes_le(&a.to_le_bytes());
    let b_dig = BigUintDig::from_bytes_le(&b.to_le_bytes());
    let inv = a_dig
        .mod_inverse(&b_dig)
        .expect("Modular inverse a^{-1} s.t. a * a^{-1} = 1 mod b does not exist")
        .to_biguint()
        .unwrap();
    BigUint::from_bytes_le(&inv.to_bytes_le())
}

/// Extended GCD algorithm. The funciton calculates the GCD of a & b
/// and two new variables x & y that satisy ax + by == gcd (i.e. Bezout's
/// identity)
///
/// Refer to attached docs for implementation details
pub fn extended_gcd(mut a: i64, mut b: i64) -> (i64, i64, i64) {
    let mut swapped = false;
    if a < b {
        swapped = true;
        mem::swap(&mut a, &mut b);
    }

    let mut r1 = a;
    let mut r2 = b;

    let mut old_x = 1;
    let mut old_y = 0;
    let mut curr_x = 0;
    let mut curr_y = 1;

    while r2 > 0 {
        let q = r1 / r2;
        let rem = r1 % r2;

        let new_x = old_x - q * curr_x;
        let new_y = old_y - q * curr_y;

        r1 = r2;
        r2 = rem;
        old_x = curr_x;
        old_y = curr_y;
        curr_x = new_x;
        curr_y = new_y;
    }

    if swapped {
        mem::swap(&mut old_x, &mut old_y)
    }

    (r1, old_x, old_y)
}

pub fn moduli_chain_to_biguint<T: UnsignedInteger>(moduli_chain: &[T]) -> BigUint
where
    BigUint: From<T>,
{
    let mut big_q = BigUint::one();
    moduli_chain
        .iter()
        .for_each(|qi| big_q *= BigUint::from(*qi));
    big_q
}

pub fn negacyclic_mul<T: UnsignedInteger, F: Fn(&T, &T) -> T>(
    a: &[T],
    b: &[T],
    mul: F,
    modulus: T,
) -> Vec<T> {
    let mut r = vec![T::zero(); a.len()];
    for i in 0..a.len() {
        for j in 0..i + 1 {
            // println!("i: {j} {}", i - j);
            r[i] = (r[i] + mul(&a[j], &b[i - j])) % modulus;
        }

        for j in i + 1..a.len() {
            // println!("i: {j} {}", a.len() - j + i);
            r[i] = (r[i] + modulus - mul(&a[j], &b[a.len() - j + i])) % modulus;
        }
        // println!("")
    }

    return r;
}

/// Maps elements of a to their bit reversed indices. Length of a must be power
/// of 2
pub fn bit_reverse_map<T: Clone>(a: &mut [T]) {
    debug_assert!(
        a.len().is_power_of_two(),
        "Expected length of input to be some power of 2: but is {}",
        a.len()
    );

    let l = a.len().leading_zeros() + 1;
    let a_clone = a.to_vec();
    for (index, el) in izip!(a_clone.iter()).enumerate() {
        let r_index = index.reverse_bits() >> l;
        a[r_index] = el.clone();
    }
}

/// Calculates M^th root of unity and returns the subgroup \psi^{0}, \psi^{1},
/// ..., \psi^{M}
pub fn psi_powers<F: BFloat, C: ComplexNumber<F>>(m: u32) -> Vec<C> {
    let m_root_unity = C::nth_root(m);
    (0..m)
        .into_iter()
        .map(|i| m_root_unity.powu(i as u32))
        .collect_vec()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{thread_rng, Rng};

    #[test]
    fn extended_gcd_works() {
        let mut rng = thread_rng();

        for _ in 0..1000 {
            let a = rng.gen_range(0..100000);
            let b = rng.gen_range(0..100000);
            let (gcd, x, y) = extended_gcd(a, b);

            // assert bezout's identity
            assert_eq!(
                (a.wrapping_mul(x).wrapping_add(b.wrapping_mul(y))),
                gcd,
                "
                Expected {a}x{x} + {b}x{y} == {gcd}
                "
            );
        }
    }

    #[test]
    fn fast_modular_inverse_word_size_works() {
        let mut rng = thread_rng();

        for _ in 0..1 {
            // multiply by 2 and add 1 to get make sure inputs are odd.
            let a_u64 = rng.gen::<u64>().wrapping_mul(2).wrapping_add(1);
            let a_u64_inv = u64::fast_inverse(a_u64);
            assert_eq!(
                a_u64.wrapping_mul(a_u64_inv),
                1,
                "{a_u64} x {a_u64_inv} (mod 2^64) != 1"
            );

            let a_u32 = rng.gen::<u32>().wrapping_mul(2).wrapping_add(1);
            let a_u32_inv = u32::fast_inverse(a_u32);
            assert_eq!(
                a_u32.wrapping_mul(a_u32_inv),
                1,
                "{a_u32} x {a_u32_inv} (mod 2^32) != 1"
            );
        }
    }
}
