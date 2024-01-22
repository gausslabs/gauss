use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use gauss::core_crypto::ntt::{NativeNTTBackend, Ntt, NttConfig};
use itertools::Itertools;
use rand::{distributions::Uniform, thread_rng, Rng};

const K: usize = 15;

const PRIMES: [u64; 5] = [
    1047041,             // 20-bit
    1073738753,          // 30-bit
    1099511603713,       // 40-bit
    1125899906826241,    // 50-bit
    1152921504606844417, // 60-bit
];

fn bench_forward(ntt_backend: &NativeNTTBackend, a: &mut Vec<u64>) {
    ntt_backend.forward(a)
}

fn bench_forward_lazy(ntt_backend: &NativeNTTBackend, a: &mut Vec<u64>) {
    ntt_backend.forward_lazy(a)
}

fn bench_backward(ntt_backend: &NativeNTTBackend, a: &mut Vec<u64>) {
    ntt_backend.backward(a)
}

fn bench_backward_lazy(ntt_backend: &NativeNTTBackend, a: &mut Vec<u64>) {
    ntt_backend.backward_lazy(a)
}

fn criterion_benchmark(c: &mut Criterion) {
    PRIMES.iter().for_each(|p| {
        let mut group = c.benchmark_group(format!("Number Theoretic Transforms mod {}", p));

        (4..=K).for_each(|n| {
            let rng = thread_rng();
            let size = 1 << n;
            let ntt_backend = NativeNTTBackend::init(*p, size);
            let v = rng.sample_iter(Uniform::new(0, p)).take(size).collect_vec();

            let id = BenchmarkId::new("Forward", size);
            group.bench_with_input(id, &v, |b, v| {
                let mut _v = v.clone();
                b.iter(|| bench_forward(black_box(&ntt_backend), black_box(&mut _v)))
            });

            let id = BenchmarkId::new("Forward Lazy", size);
            group.bench_with_input(id, &v, |b, v| {
                let mut _v = v.clone();
                b.iter(|| bench_forward_lazy(black_box(&ntt_backend), black_box(&mut _v)))
            });

            let id = BenchmarkId::new("Backward", size);
            group.bench_with_input(id, &v, |b, v| {
                let mut _v = v.clone();
                b.iter(|| bench_backward(black_box(&ntt_backend), black_box(&mut _v)))
            });

            let id = BenchmarkId::new("Backward Lazy", size);
            group.bench_with_input(id, &v, |b, v| {
                let mut _v = v.clone();
                b.iter(|| bench_backward_lazy(black_box(&ntt_backend), black_box(&mut _v)))
            });
        });
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
