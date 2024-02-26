mod entities;

pub type CkksSecretKey = entities::CkksSecretKey;
pub type CkksCiphertext =
    self::entities::CkksCiphertextGenericStorage<aligned_vec::AVec<aligned_vec::AVec<u64>>>;

#[cfg(test)]
mod tests {
    use itertools::Itertools;
    use num_complex::{Complex, ComplexDistribution};
    use rand::{thread_rng, Rng};
    use rand_distr::Uniform;

    use crate::{
        keys::{Decryptor, Encryptor, LevelDecoder, LevelEncoder},
        utils::print_precision_stats,
    };

    use self::entities::build_parameters;

    use super::*;

    fn sample_random_complex(low: f64, high: f64, size: usize) -> Vec<Complex<f64>> {
        let reals = Uniform::new(low, high);
        let imags = Uniform::new(low, high);
        let complex_distr = ComplexDistribution::new(reals, imags);
        thread_rng()
            .sample_iter(complex_distr)
            .take(size)
            .collect_vec()
    }

    #[test]
    fn complex_vec_encoding_decoding_works() {
        let ring_size = 1 << 4;
        let delta = 2.0f64.powi(40i32);
        build_parameters(&[50], delta, ring_size);

        let m = sample_random_complex(-1.0, 1.0, ring_size >> 1);
        let encoded_m: Vec<Vec<u64>> = m.encode(0);
        let m_back: Vec<Complex<f64>> = encoded_m.decode(0);

        print_precision_stats(&m_back, &m);
    }

    #[test]
    fn ckks_encryption_decryption_works() {
        let ring_size = 1 << 4;
        let delta = 2.0f64.powi(40i32);
        build_parameters(&[50, 50], delta, ring_size);

        let m = sample_random_complex(-1.0, 1.0, ring_size >> 1);
        let secret = CkksSecretKey::new();

        let ct: CkksCiphertext = secret.encrypt(&m);
        let m_back = secret.decrypt(&ct);

        print_precision_stats(&m_back, &m);
    }
}
