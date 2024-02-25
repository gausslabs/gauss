pub(crate) mod bfv;
pub(crate) mod ckks;
pub(crate) mod hybrid_ksk;

pub trait WithGlobal {
    fn with_global<F, R>(func: F) -> R
    where
        F: Fn(&Self) -> R;
}

// scheme invariant ops
pub(crate) mod ops {
    use itertools::Itertools;
    use num_traits::Signed;
    use rand::CryptoRng;

    use crate::core_crypto::random::RandomUniformDist;

    pub fn generate_ternery_secret_with_hamming_weight<
        Scalar: Signed + Clone,
        R: RandomUniformDist<[u8], Parameters = u8>
            + RandomUniformDist<usize, Parameters = usize>
            + CryptoRng,
    >(
        rng: &mut R,
        hamming_weight: usize,
        ring_size: usize,
    ) -> Vec<Scalar> {
        let mut bytes = vec![0u8; (hamming_weight.div_ceil(8)) as usize];
        RandomUniformDist::random_fill(rng, &0u8, bytes.as_mut_slice());

        let mut secret = vec![Scalar::zero(); ring_size];
        let mut secret_indices = (0..ring_size).into_iter().collect_vec();
        let mut bit_index = 0;
        let mut byte_index = 0;
        for _ in 0..hamming_weight {
            let mut secret_index = 0usize;
            rng.random_fill(&secret_indices.len(), &mut secret_index);

            if bytes[byte_index] & 1u8 == 1 {
                secret[secret_indices[secret_index]] = Scalar::one();
            } else {
                secret[secret_indices[secret_index]] = Scalar::neg(Scalar::one());
            }

            bytes[byte_index] >>= 1;

            // remove secret_index from secret_indices
            secret_indices[secret_index] = *secret_indices.last().unwrap();
            secret_indices.truncate(secret_indices.len() - 1);

            if bit_index == 7 {
                bit_index = 0;
                byte_index += 1;
            } else {
                bit_index += 1;
            }
        }

        secret
    }
}
