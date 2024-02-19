use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use rand_distr::{Distribution, Normal};

/// `TruncatedDiscreteGaussian` represents the parameters of a truncated discrete Gaussian distribution with:
/// * standard deviation `sigma`
/// * bounds `[-bound, bound]`, where bound is equal to `(6*sigma).round()`
/// * mean 0
#[derive(Clone, Debug)]
struct TruncatedDiscreteGaussian {
    sigma: f64,
    bound: i64,
}

#[derive(Clone, Debug)]
struct GaussianSampler {
    params: TruncatedDiscreteGaussian,
    rng: ChaCha20Rng,
}

impl GaussianSampler {
    /// `new` returns a new `GaussianSampler` with the given parameters starting from the standard deviation `sigma`.
    pub fn new(sigma: f64) -> Self {
        assert!(sigma > 0.0, "sigma must be positive");
        let bound = (6.0 * sigma).round() as i64;
        let rng = ChaCha20Rng::from_entropy();
        Self {
            params: TruncatedDiscreteGaussian { sigma, bound },
            rng,
        }
    }

    /// `sample` returns a sample from the truncated discrete Gaussian distribution.
    /// The technique used is patterned after [lattigo](https://github.com/tuneinsight/lattigo/blob/c031b14be1fb3697945709d7afbed264fa845442/ring/sampler_gaussian.go#L71).
    /// In particular, `sampled_val` is sampled from a normal distribution with mean 0 and standard deviation `sigma`.
    /// If `sampled_val` is within the bounds, it is rounded and returned.
    pub fn sample(&mut self) -> i64 {
        let discrete_normal = Normal::new(0.0, self.params.sigma).unwrap();
        let sampled_val = discrete_normal.sample(&mut self.rng);
        if sampled_val.abs() < self.params.bound as f64 {
            sampled_val.round() as i64
        } else {
            self.sample()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;
    #[test]
    fn gaussian_sampler_truncation_limit() {
        // Assert that the samples are within the bounds of the truncated Gaussian distribution.
        let sigma = rand::thread_rng().gen_range(1.0..=100.0);
        let mut sampler = GaussianSampler::new(sigma);

        for _ in 0..5000000 {
            let sample = sampler.sample();
            assert!((-sampler.params.bound..=sampler.params.bound).contains(&sample));
        }
    }

    #[test]
    fn gaussian_sampler_mean() {
        // Assert that the mean of the samples is close to 0.
        let sigma = rand::thread_rng().gen_range(1.0..=100.0);
        let mut sampler = GaussianSampler::new(sigma);

        let mut sum = 0;
        let n = 5000000;
        for _ in 0..n {
            sum += sampler.sample();
        }
        let mean = sum as f64 / n as f64;
        assert!((mean - 0.1..=mean + 0.1).contains(&0.0));
    }

    #[test]
    fn gaussian_sampler_pmf_consistency() {
        // Assert that the probability mass function for a random `x` observed by performing sample on the `sampler`
        // is consistent with the expected distribution.
        let sigma = rand::thread_rng().gen_range(1.0..=100.0);
        let mut sampler = GaussianSampler::new(sigma);
        let val = rand::thread_rng().gen_range(-sampler.params.bound..=sampler.params.bound);

        let pmf_expected = (1.0 / ((2.0 * std::f64::consts::PI).sqrt() * sampler.params.sigma))
            * (-0.5 * (val as f64 / sampler.params.sigma).powi(2)).exp();

        let mut count = 0;
        let n = 5000000;
        for _ in 0..n {
            if sampler.sample() == val {
                count += 1;
            }
        }
        let pmf_observed = count as f64 / n as f64;
        assert!((pmf_observed - 0.001..=pmf_observed + 0.001).contains(&pmf_expected));
    }
}
