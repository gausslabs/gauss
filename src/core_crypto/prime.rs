use rand::{Rng, RngCore};

use crate::utils::mod_exponent;

// trait PrimalityTest<T> {
//     fn is_prime(candidate: T) -> bool;
// }

// TODO (Jay): this is only a workaround. Add a propoer way to perform primality
// tests.
fn is_probably_prime(candidate: u64) -> bool {
    num_bigint_dig::prime::probably_prime(&num_bigint_dig::BigUint::from(candidate), 0)
}

pub fn generate_primes_vec(
    sizes: &[usize],
    polynomial_degree: usize,
    skip_list: &[u64],
) -> Vec<u64> {
    let mut primes = vec![];
    sizes.iter().for_each(|s| {
        let mut upper_bound = 1u64 << s;
        loop {
            if let Some(p) = generate_prime(*s, (2 * polynomial_degree) as u64, upper_bound) {
                if !primes.contains(&p) && !skip_list.contains(&p) {
                    primes.push(p);
                    break;
                } else {
                    upper_bound = p;
                }
            } else {
                panic!("Not enough primes");
            }
        }
    });
    primes
}

/// Finds prime that satisfy
/// - $prime \lt upper_bound$
/// - $\log{prime} = num_bits$
/// - `prime % modulo == 1`
pub fn generate_prime(num_bits: usize, modulo: u64, upper_bound: u64) -> Option<u64> {
    let leading_zeros = (64 - num_bits) as u32;

    let mut tentative_prime = upper_bound - 1;
    while tentative_prime % modulo != 1 && tentative_prime.leading_zeros() == leading_zeros {
        tentative_prime -= 1;
    }

    while !is_probably_prime(tentative_prime)
        && tentative_prime.leading_zeros() == leading_zeros
        && tentative_prime >= modulo
    {
        tentative_prime -= modulo;
    }

    if is_probably_prime(tentative_prime) && tentative_prime.leading_zeros() == leading_zeros {
        Some(tentative_prime)
    } else {
        None
    }
}

/// Find n^{th} root of unity in field F_q, if one exists
///
/// Note: n^{th} root of unity exists if and only if $q = 1 \mod{n}$
pub(crate) fn find_primitive_root<R: RngCore>(q: u64, n: u64, rng: &mut R) -> Option<u64> {
    assert!(n.is_power_of_two(), "{n} is not power of two");

    // n^th root of unity only exists if n|(q-1)
    assert!(q % n == 1, "{n}^th root of unity in F_{q} does not exists");

    let t = (q - 1) / n;

    for _ in 0..100 {
        let mut omega = rng.gen::<u64>() % q;

        // \omega = \omega^t. \omega is now n^th root of unity
        omega = mod_exponent(omega, t, q);

        #[cfg(debug)]
        {
            let check = mod_exponent(omega, n, q);
            debug_assert!(
                check == 1,
                "omega({omega}) is not n^th root of unity: Expected 1 but is {check}",
            );
        }

        // We restrict n to be power of 2. Thus checking whether \omega is primitive
        // n^th root of unity is as simple as checking: \omega^{n/2} != 1
        if mod_exponent(omega, n >> 1, q) == 1 {
            continue;
        } else {
            return Some(omega);
        }
    }

    None
}

#[cfg(test)]
mod test {
    use rand::thread_rng;

    use crate::core_crypto::modulus::{BarrettBackend, ModulusBackendConfig, NativeModulusBackend};

    use super::*;

    const Q_60_BITS: u64 = 1152921504606748673;
    const N: u64 = 1 << 15;

    #[test]
    fn find_primitive_root_works() {
        let mut rng = thread_rng();
        let root = find_primitive_root(Q_60_BITS, N, &mut rng).unwrap();

        // check root^n = 1
        // TODO (Jay): Hardcode tests here, probably with some python script
        let modulus = NativeModulusBackend::initialise(Q_60_BITS);
        let mut root_n = 1;
        for i in 0..N {
            root_n = modulus.mul_mod_fast(root_n, root);
        }

        assert!(
            root_n == 1,
            "Incorrect {N}^th root of unity: {root}^{N} != 1"
        );
    }
}
