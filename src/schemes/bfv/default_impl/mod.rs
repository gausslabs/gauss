use entities::BfvSecretKey;

mod entities;

pub type BfvCiphertext = entities::BfvCiphertextScalarU64GenericStorage<Vec<Vec<u64>>>;

#[cfg(test)]
mod tests {
    use crate::keys::{Decryptor, Encryptor, SecretKey};

    use self::entities::{build_parameters, BfvClientParametersForScalarU64, WithGlobal};

    use super::*;

    #[test]
    fn secret_key_has_correct_hamming_weight() {
        for n in [4, 8, 15] {
            let ring_size = 1 << n;
            build_parameters(&[50, 50], 65537, ring_size);
            let secret = BfvSecretKey::new();

            let mut hamming_weight = 0;
            secret.values().iter().for_each(|v| {
                if *v != 0 {
                    hamming_weight += 1;
                }
            });

            BfvClientParametersForScalarU64::with_global(|parameters| {
                assert_eq!(
                    parameters.secret_hw(),
                    hamming_weight,
                    "Expected Hamming weight {} but is {}",
                    parameters.secret_hw(),
                    hamming_weight
                );
            });
        }
    }

    #[test]
    fn encryption_decryption_works() {
        build_parameters(&[50, 50], 65537, 1 << 4);
        let secret = BfvSecretKey::new();

        let ct: BfvCiphertext = secret.encrypt(&[21]);
        let m = secret.decrypt(ct);
        dbg!(m);
    }
}
