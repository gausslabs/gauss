use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use gauss::core_crypto::ntt::NativeNTTBackend;
use itertools::Itertools;
use rand::{thread_rng, Rng};

const K: usize = 8;

const PRIMES: [u64; 10] = [
    32257,               // 15-bit
    1047041,             // 20-bit
    33551873,            // 25-bit
    1073738753,          // 30-bit
    34359724033,         // 35-bit
    1099511603713,       // 40-bit
    35184372088321,      // 45-bit
    1125899906826241,    // 50-bit
    36028797018963457,   // 55-bit
    1152921504606844417, // 60-bit
];

fn bench_ntt(ntt_backend: &NativeNTTBackend, a: &mut Vec<u64>) {
    ntt_backend.ntt(a)
}

fn bench_intt(ntt_backend: &NativeNTTBackend, a: &mut Vec<u64>) {
    ntt_backend.ntt(a)
}

fn criterion_benchmark(c: &mut Criterion) {
    let mut rng = thread_rng();

    PRIMES.iter().for_each(|p| {
        let mut group = c.benchmark_group(format!("Montgomery Vector Arithmetic mod {}", p));

        (1..K).for_each(|n| {
            let size = 1 << n;
            let ntt_backend = NativeNTTBackend::new(*p, size);
            let v = (0..size).map(|_| rng.gen::<u64>()).collect_vec();

            let id = BenchmarkId::new("NTT", size);
            group.bench_with_input(id, &v, |b, v| {
                let mut _v = v.clone();
                b.iter(|| bench_ntt(black_box(&ntt_backend), black_box(&mut _v)))
            });

            let id = BenchmarkId::new("INTT", size);
            group.bench_with_input(id, &v, |b, v| {
                let mut _v = v.clone();
                b.iter(|| bench_intt(black_box(&ntt_backend), black_box(&mut _v)))
            });
        });
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
