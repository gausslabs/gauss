use std::{char::MAX, cmp::Ordering, ops::Index, thread::Thread};

use itertools::Itertools;
use num_traits::{
    float::{self, FloatCore},
    Float, Inv, Num, ToPrimitive,
};
use rand::{
    distributions::uniform::SampleUniform, distributions::Uniform, thread_rng, CryptoRng, Rng,
    RngCore, SeedableRng,
};

use crate::utils::get_bit_at;

use std::{borrow::Borrow, cell::RefCell};

use itertools::izip;
use num_traits::Zero;
use rand_chacha::{ChaCha8Core, ChaCha8Rng};

use super::matrix::{Matrix, MatrixMut, RowMut};

pub trait RandomSecretByteGenerator {
    fn random_bytes(&mut self, size: usize) -> Vec<u8>;
}

pub trait RandomSecretValueGenerator<T> {
    fn random_value_in_range(&mut self, range: T) -> T;
}

pub trait RandomGaussianDist<S, P> {
    fn random_vec_in_modulus(&mut self, modulus: S, size: usize) -> Vec<S>;
    fn random_ring_poly(&mut self, moduli_chain: &[S], ring_size: usize) -> P;
}

pub trait RandomUniformDist<S, P> {
    fn random_vec_in_modulus(&mut self, modulus: S, size: usize) -> Vec<S>;
    fn random_ring_poly(&mut self, moduli_chain: &[S], ring_size: usize) -> P;
}

pub trait WithLocal {
    fn with_local<F, R>(func: F) -> R
    where
        F: Fn(&Self) -> R;

    fn with_local_mut<F, R>(func: F) -> R
    where
        F: Fn(&mut Self) -> R;
}

pub(crate) struct DefaultU64SeededRandomGenerator {
    rng: ChaCha8Rng,
}

thread_local! {
    pub(crate) static DEFAULT_U64_SEEDED_RNG: RefCell<DefaultU64SeededRandomGenerator> = RefCell::new(DefaultU64SeededRandomGenerator::new());
}

impl WithLocal for DefaultU64SeededRandomGenerator {
    fn with_local<F, R>(func: F) -> R
    where
        F: Fn(&Self) -> R,
    {
        DEFAULT_U64_SEEDED_RNG.with_borrow(|value| func(value))
    }

    fn with_local_mut<F, R>(func: F) -> R
    where
        F: Fn(&mut Self) -> R,
    {
        DEFAULT_U64_SEEDED_RNG.with_borrow_mut(|value| func(value))
    }
}

impl DefaultU64SeededRandomGenerator {
    pub fn new() -> Self {
        DefaultU64SeededRandomGenerator {
            rng: ChaCha8Rng::from_rng(thread_rng()).unwrap(),
        }
    }

    /// Seed the rng of thread local generator with new seed
    pub fn new_with_seed(seed: <ChaCha8Rng as SeedableRng>::Seed) {
        DefaultU64SeededRandomGenerator::with_local_mut(|generator| {
            let rng = ChaCha8Rng::from_seed(seed);
            generator.rng = rng;
        });
    }
}

impl RandomSecretByteGenerator for DefaultU64SeededRandomGenerator {
    fn random_bytes(&mut self, size: usize) -> Vec<u8> {
        let mut bytes = vec![0u8; size];
        self.rng.fill_bytes(&mut bytes);
        bytes
    }
}

impl<T: SampleUniform + PartialOrd + Zero> RandomSecretValueGenerator<T>
    for DefaultU64SeededRandomGenerator
{
    fn random_value_in_range(&mut self, range: T) -> T {
        self.rng.gen_range(T::zero()..range)
    }
}

impl<P: Matrix> RandomGaussianDist<u64, P> for DefaultU64SeededRandomGenerator {
    fn random_vec_in_modulus(&mut self, modulus: u64, size: usize) -> Vec<u64> {
        vec![0u64; size]
    }

    fn random_ring_poly(&mut self, moduli_chain: &[u64], ring_size: usize) -> P {
        P::zeros(moduli_chain.len(), ring_size)
    }
}

impl<P: Matrix<MatElement = u64> + MatrixMut> RandomUniformDist<u64, P>
    for DefaultU64SeededRandomGenerator
where
    <P as Matrix>::R: RowMut,
{
    fn random_ring_poly(&mut self, moduli_chain: &[u64], ring_size: usize) -> P {
        let rows = moduli_chain.len();
        let mut poly = P::zeros(rows, ring_size);

        izip!(poly.iter_rows_mut(), moduli_chain.iter()).for_each(|(r, qi)| {
            izip!(
                r.as_mut().iter_mut(),
                (&mut self.rng).sample_iter(Uniform::new(0, *qi))
            )
            .for_each(|(r_el, random_el)| *r_el = random_el);
        });

        poly
    }
    fn random_vec_in_modulus(&mut self, modulus: u64, size: usize) -> Vec<u64> {
        vec![0u64; size]
    }
}

const MAX_TREE_DEPTH: usize = 64;
const MAX_LEVELS: usize = 4;
const PRECISION: u64 = 53;
const BERNOULLI_FLIPS: u64 = 23;
const CENTER_COUNT: u64 = 1024;
const SMOOTHING_PARAMTER: f64 = 6.;

#[derive(Clone)]
struct PeikertSampler {
    fin: i64,
    m_vals: Vec<f64>,
    b_mean: f64,
}

impl PeikertSampler {
    pub fn new(mean: f64, std: f64) -> Self {
        let acc = 1e-17;
        let fin = ((acc.log10() * -2.).sqrt() * std).ceil() as i64;

        let mean = if mean >= 0. {
            mean.floor()
        } else {
            mean.ceil()
        };

        let mean = mean - mean;
        let variance = std * std;

        let cumsum = ((-1 * fin)..=fin)
            .map(|_x| {
                let x = _x as f64;
                (-((x - mean) * (x - mean) / (variance * 2.))).exp()
            })
            .sum::<f64>();

        let mut m_vals = ((-1 * fin)..=fin)
            .map(|i| {
                let inner = (i as f64 - mean) * (i as f64 - mean) / (2. * variance) as f64;
                ((-inner).exp()) / cumsum
            })
            .collect_vec();

        (1..m_vals.len()).for_each(|i| {
            m_vals[i] += m_vals[i - 1];
        });

        PeikertSampler {
            b_mean: mean,
            fin,
            m_vals,
        }
    }

    pub fn random_integer(self) -> i64 {
        let mut rng = thread_rng();

        let seed = rng.gen_range(0. ..1.);
        let val = self.m_vals.binary_search_by(|&x| {
            if x < seed {
                Ordering::Less
            } else if x == seed {
                Ordering::Equal
            } else {
                Ordering::Greater
            }
        });

        let val = match val {
            Ok(x) => x,
            Err(x) => x,
        } as i64;

        val - self.fin + self.b_mean as i64
    }
}

#[derive(Clone)]
enum Sampler {
    BaseSampler(PeikertSampler),
    CombinedSampler(CombinedSampler),
}

impl Sampler {
    pub fn generate_integer(&mut self) -> i64 {
        match self {
            Self::BaseSampler(s) => <PeikertSampler as Clone>::clone(&(*s)).random_integer(),
            Self::CombinedSampler(s) => s.random_integer(),
        }
    }

    pub fn random_bit(&mut self) -> bool {
        self.generate_integer() & 1 > 0
    }
}

#[derive(Clone)]
struct CombinedSampler {
    x1: i64,
    x2: i64,
    rng1: Box<Sampler>,
    rng2: Box<Sampler>,
}

impl CombinedSampler {
    // SampleI in UCSD paper
    pub fn random_integer(&mut self) -> i64 {
        match (self.rng1.as_mut(), self.rng2.as_mut()) {
            (Sampler::BaseSampler(sampler1), Sampler::BaseSampler(sampler2)) => self
                .x1
                .wrapping_mul(<PeikertSampler as Clone>::clone(&sampler1).random_integer())
                .wrapping_add(
                    self.x2
                        .wrapping_mul(<PeikertSampler as Clone>::clone(&sampler2).random_integer()),
                ),
            (Sampler::CombinedSampler(sampler1), Sampler::CombinedSampler(sampler2)) => self
                .x1
                .wrapping_mul(sampler1.random_integer())
                .wrapping_add(self.x2.wrapping_mul(sampler2.random_integer())),
            (Sampler::BaseSampler(sampler1), Sampler::CombinedSampler(sampler2)) => self
                .x1
                .wrapping_mul(<PeikertSampler as Clone>::clone(&sampler1).random_integer())
                .wrapping_add(self.x2.wrapping_mul(sampler2.random_integer())),
            (Sampler::CombinedSampler(sampler1), Sampler::BaseSampler(sampler2)) => self
                .x1
                .wrapping_mul(sampler1.random_integer())
                .wrapping_add(
                    self.x2
                        .wrapping_mul(<PeikertSampler as Clone>::clone(&sampler2).random_integer()),
                ),
        }
    }
}

pub struct DiscreteGuassianGenerator {
    mean: f64,
    std: f64,

    base_samplers: Vec<Sampler>,
    wide_sampler: Sampler,
    mask: i64,
    sampler_variance: f64,
    wide_variance: f64,
    k: u64,
    log_base: u64,
}

impl DiscreteGuassianGenerator {
    pub fn new(mean: f64, std: f64, log_base: u64) -> Self {
        let mut x1 = 0;
        let mut x2 = 0;
        let n = SMOOTHING_PARAMTER;

        let k = ((PRECISION - BERNOULLI_FLIPS) as f64 / log_base as f64).ceil() as u64;
        let mask = (1 << log_base) - 1;

        let rngs = (0..CENTER_COUNT)
            .map(|_| Sampler::BaseSampler(PeikertSampler::new(mean, std)))
            .collect_vec();

        let base_variance = std * std;
        let mut wide_variance = base_variance;
        let mut wide_sampler = rngs[0].to_owned();

        for i in 0..MAX_LEVELS {
            x1 = (wide_variance / (2. * n * n)).sqrt().floor() as i64;
            x2 = (x1 - 1).max(1);

            // SampleI
            if i == 1 {
                let rng = Box::new(rngs[0].to_owned());

                wide_sampler = Sampler::CombinedSampler(CombinedSampler {
                    x1,
                    x2,
                    rng1: rng.clone(),
                    rng2: rng,
                });
            } else {
                wide_sampler = Sampler::CombinedSampler(CombinedSampler {
                    x1,
                    x2,
                    rng1: Box::new(wide_sampler.clone()),
                    rng2: Box::new(wide_sampler.clone()),
                })
            }
            wide_variance = (x1 * x1 + x2 * x2).to_f64().unwrap() * wide_variance;
        }

        let t = (1 / (1 << (2 * log_base))) as f64;
        let s = t.powi(k as i32);
        let sampler_variance = base_variance * (1. + s * (k as f64));

        Self {
            mean,
            std,
            base_samplers: rngs,
            wide_sampler,
            mask,
            sampler_variance,
            wide_variance,
            k,
            log_base,
        }
    }

    // Sample Z from UCSD paper
    pub fn generate_integer(&mut self) -> i64 {
        let mean = self.mean;
        let std = self.std;

        let variance = std * std;

        let x = self.wide_sampler.generate_integer() as f64; // x = SampleI(max)
        let rhs = ((variance - self.sampler_variance) / self.wide_variance).sqrt();
        let mut c = mean + x * rhs;
        let c_floor = c.floor();

        c -= c_floor;

        (c_floor as i64) + self.flip_and_round(c)
    }

    fn flip_and_round(&mut self, center: f64) -> i64 {
        let c = center as i64 * (1 << PRECISION);
        let base_c = c >> BERNOULLI_FLIPS;

        let mut random_bit = false;

        for i in (0..BERNOULLI_FLIPS).rev() {
            random_bit = self.base_samplers[0].random_bit();
            let c_i = get_bit_at(c, i as u8);
            if random_bit && !c_i {
                return self.sample_c(base_c);
            }
            if !random_bit && c_i {
                return self.sample_c(base_c + 1);
            }
        }
        self.sample_c(base_c + 1)
    }

    fn sample_c(&mut self, center: i64) -> i64 {
        let mut c = center;
        let mut sample: i64 = 0;

        (0..self.k).for_each(|_| {
            sample = self.base_samplers[(self.mask & c) as usize].generate_integer();
            if (self.mask & c) > 0 && c < 0 {
                sample -= 1;
            }
            (0..self.log_base).for_each(|_| {
                c /= 2;
            });
            c += sample;
        });

        c
    }
}

#[cfg(test)]
mod tests {
    use super::{DiscreteGuassianGenerator, PeikertSampler, CENTER_COUNT};
    const STD: f64 = 1.0;
    const MEAN: f64 = 0.0;

    #[test]
    fn test_discrete_gaussian() {
        let log_base = ((CENTER_COUNT as f64).log10() / (2_f64).log10()) as u64;
        let mut gen = DiscreteGuassianGenerator::new(MEAN, STD, log_base);

        for _ in 0..100 {
            let res = gen.generate_integer();
            assert!(res.abs() < 10);
        }
    }

    #[test]
    fn test_peikert_sampler() {
        for _ in 0..100 {
            let peikert = PeikertSampler::new(MEAN, STD);
            let x = peikert.random_integer();
            assert!(x.abs() < 10);
        }
    }
}

impl CryptoRng for DefaultU64SeededRandomGenerator {}
