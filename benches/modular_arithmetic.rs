use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use itertools::Itertools;
use rand::{thread_rng, Rng};

use gauss::core_crypto::modulus::{
    barrett::BarrettBackend, ModulusArithmeticBackend, ModulusBackendConfig, MontgomeryBackend,
    MontgomeryScalar, NativeModulusBackend,
};

const PRIMES: [u64; 6] = [
    1021,                // 10-bit
    1048361,             // 20-bit
    1073741101,          // 30-bit
    1099511626321,       // 40-bit
    1125899906842201,    // 50-bit
    1152921504606845161, // 60-bit
];

const K: usize = 16;

// Barrett Arithmetic
fn bench_barrett_add(modulus_backend: &NativeModulusBackend, a: &u64, b: &u64) {
    modulus_backend.add_mod_fast(*a, *b);
}

fn bench_barrett_sub(modulus_backend: &NativeModulusBackend, a: &u64, b: &u64) {
    modulus_backend.sub_mod_fast(*a, *b);
}

fn bench_barrett_mul(modulus_backend: &NativeModulusBackend, a: &u64, b: &u64) {
    modulus_backend.mul_mod_fast(*a, *b);
}

// Montgomery Arithmetic

fn bench_mont_add(
    modulus_backend: &NativeModulusBackend,
    a: &MontgomeryScalar<u64>,
    b: &MontgomeryScalar<u64>,
) {
    modulus_backend.mont_add(*a, *b);
}

fn bench_mont_sub(
    modulus_backend: &NativeModulusBackend,
    a: &MontgomeryScalar<u64>,
    b: &MontgomeryScalar<u64>,
) {
    modulus_backend.mont_sub(*a, *b);
}

fn bench_mont_mul(
    modulus_backend: &NativeModulusBackend,
    a: &MontgomeryScalar<u64>,
    b: &MontgomeryScalar<u64>,
) {
    modulus_backend.mont_mul(*a, *b);
}

// Montgomery Vector Arithmetic

fn bench_mont_fma(
    modulus_backend: &NativeModulusBackend,
    a: &Vec<MontgomeryScalar<u64>>,
    b: &Vec<MontgomeryScalar<u64>>,
) {
    modulus_backend.mont_fma(a, b);
}

fn criterion_modular_arithmetic(c: &mut Criterion) {
    let mut group = c.benchmark_group("Montgomery Modular Arithmetic Benchmarks");
    let mut rng = thread_rng();

    PRIMES.iter().for_each(|p| {
        let modulus_backend = <NativeModulusBackend as ModulusBackendConfig<u64>>::initialise(*p);

        let a = rng.gen::<u64>() % p;
        let b = rng.gen::<u64>() % p;
        let a_mont = modulus_backend.normal_to_mont_space(a);
        let b_mont = modulus_backend.normal_to_mont_space(b);

        // Addition
        let id = BenchmarkId::new("Barrett Addition", p);
        group.bench_with_input(id, &(a, b), |bx, (a, b)| {
            bx.iter(|| bench_barrett_add(black_box(&modulus_backend), black_box(a), black_box(b)))
        });

        let id = BenchmarkId::new("Montgomery Addition", p);
        group.bench_with_input(id, &(a_mont, b_mont), |bx, (a_mont, b_mont)| {
            bx.iter(|| {
                bench_mont_add(
                    black_box(&modulus_backend),
                    black_box(&a_mont),
                    black_box(&b_mont),
                )
            })
        });

        // Subtraction
        let id = BenchmarkId::new("Barrett Subtraction", p);
        group.bench_with_input(id, &(a, b), |bx, (a, b)| {
            bx.iter(|| bench_barrett_sub(black_box(&modulus_backend), black_box(a), black_box(b)))
        });

        let id = BenchmarkId::new("Montgomery Subtraction", p);
        group.bench_with_input(id, &(a_mont, b_mont), |bx, (a_mont, b_mont)| {
            bx.iter(|| {
                bench_mont_sub(
                    black_box(&modulus_backend),
                    black_box(&a_mont),
                    black_box(&b_mont),
                )
            })
        });

        // Multiplication
        let id = BenchmarkId::new("Barrett Multiplication", p);
        group.bench_with_input(id, &(a, b), |bx, (a, b)| {
            bx.iter(|| bench_barrett_mul(black_box(&modulus_backend), black_box(a), black_box(b)))
        });

        let id = BenchmarkId::new("Montgomery Multiplication", p);
        group.bench_with_input(id, &(a_mont, b_mont), |bx, (a_mont, b_mont)| {
            bx.iter(|| {
                bench_mont_mul(
                    black_box(&modulus_backend),
                    black_box(&a_mont),
                    black_box(&b_mont),
                )
            })
        });
    })
}

fn criterion_vector_arithmetic(c: &mut Criterion) {
    let mut rng = thread_rng();

    PRIMES.iter().for_each(|p| {
        let mut group = c.benchmark_group(format!("Montgomery Vector Arithmetic mod {}", p));
        let modulus_backend = <NativeModulusBackend as ModulusBackendConfig<u64>>::initialise(*p);

        (3..K).for_each(|n| {
            let size = 1 << n;
            let v = (0..size).map(|_| rng.gen::<u64>()).collect_vec();
            let a_vec_mont = v
                .iter()
                .map(|v| modulus_backend.normal_to_mont_space(*v))
                .collect_vec();
            let b_vec_mont = v
                .iter()
                .map(|v| modulus_backend.normal_to_mont_space(*v))
                .collect_vec();

            let id = BenchmarkId::new("FMA", size);
            group.bench_with_input(id, &(a_vec_mont, b_vec_mont), |b, (a_mont, b_mont)| {
                b.iter(|| {
                    bench_mont_fma(
                        black_box(&modulus_backend),
                        black_box(a_mont),
                        black_box(b_mont),
                    )
                })
            });
        });
    });
}

criterion_group!(
    benches,
    criterion_modular_arithmetic,
    criterion_vector_arithmetic
);
criterion_main!(benches);
