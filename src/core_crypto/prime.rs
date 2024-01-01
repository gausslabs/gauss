use rand::{Rng, RngCore};

use crate::utils::mod_exponent;

/// Find n^{th} root of unity in field F_q, if one exists
///
/// Note: n^{th} root of unity exists if and only if $q = 1 \mod{n}$
fn find_primitive_root<R: RngCore>(q: u64, n: u64, rng: &mut R) -> Option<u64> {
    assert!(n.is_power_of_two(), "{n} is not power of two");

    // n^th root of unity only exists if n|(q-1)
    assert!(q % n == 1, "{n}^th root of unity in F_{q} does not exists");

    let t = (q - 1) / n;

    for _ in 0..100 {
        let mut omega = rng.gen::<u64>() % q;

        // \omega = \omega^t. \omega is now n^th root of unity
        omega = mod_exponent(omega, t, q);

        // We restrict n to be power of 2. Thus checking whether \omega is primitive n^th root of unity is a simple check: \omega^{n/2} != 1
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
