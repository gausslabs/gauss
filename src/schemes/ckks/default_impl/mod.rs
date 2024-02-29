mod entities;

pub type CkksSecretKey = entities::CkksSecretKey;
pub type CkksCiphertext =
    self::entities::CkksCiphertextGenericStorage<aligned_vec::AVec<aligned_vec::AVec<u64>>>;

#[cfg(test)]
mod tests {
    use itertools::Itertools;
    use num_complex::{Complex, ComplexDistribution};
    use num_traits::{zero, Zero};
    use rand::{thread_rng, Rng};
    use rand_distr::Uniform;

    use crate::{
        core_crypto::{
            num::big_float::BigFloat,
            random::{DefaultU64SeededRandomGenerator, RandomUniformDist},
        },
        keys::{Decryptor, Encryptor, LevelDecoder, LevelEncoder},
        utils::print_precision_stats,
    };

    use self::entities::build_parameters;

    use super::*;

    #[test]
    fn complex_vec_encoding_decoding_works() {
        let ring_size = 1 << 4;
        let delta = 2.0f64.powi(40i32);
        build_parameters(&[59], delta.into(), ring_size);

        let mut test_rng = DefaultU64SeededRandomGenerator::new();
        let mut m = vec![Complex::<BigFloat>::zero(); ring_size >> 1];
        <DefaultU64SeededRandomGenerator as RandomUniformDist<[Complex<BigFloat>]>>::random_fill(
            &mut test_rng,
            &(-1.0f64, 1.0f64),
            m.as_mut_slice(),
        );
        let encoded_m: Vec<Vec<u64>> = m.encode(0);
        let m_back: Vec<Complex<BigFloat>> = encoded_m.decode(0);

        print_precision_stats(&m_back, &m);
    }

    #[test]
    fn ckks_encryption_decryption_works() {
        let ring_size = 1 << 4;
        let delta = 2.0f64.powi(45i32);
        build_parameters(&[59, 59], delta.into(), ring_size);

        let mut test_rng = DefaultU64SeededRandomGenerator::new();
        let mut m = vec![Complex::<BigFloat>::zero(); ring_size >> 1];
        <DefaultU64SeededRandomGenerator as RandomUniformDist<[Complex<BigFloat>]>>::random_fill(
            &mut test_rng,
            &(-1.0f64, 1.0f64),
            m.as_mut_slice(),
        );
        let secret = CkksSecretKey::new();

        let ct: CkksCiphertext = secret.encrypt(&m);
        let m_back = secret.decrypt(&ct);

        print_precision_stats(&m_back, &m);
    }
}
