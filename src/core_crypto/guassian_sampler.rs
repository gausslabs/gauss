use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use rand_distr::{Distribution, Normal};

#[derive(Clone, Debug)]
struct TruncatedDiscreteGaussian {
    sigma: f64,
    bound: i64,
}

/// `TruncatedDiscreteGaussianSampler` represents an instance to sample values
/// from a trucated discrete Gaussian distribution with:
/// * standard deviation `sigma`
/// * bounds `[-bound, bound]`, where bound is equal to `(6*sigma).round()`
/// * mean 0
#[derive(Clone, Debug)]
struct TruncatedDiscreteGaussianSampler {
    sigma: f64,
    bound: i64,
    normal: Normal<f64>,
    rng: ChaCha20Rng,
}

impl TruncatedDiscreteGaussianSampler {
    /// `new` returns a new `TruncatedDiscreteGaussianSampler` with the given
    /// parameters starting from the standard deviation `sigma`.
    pub fn new(sigma: f64) -> Self {
        assert!(sigma > 0.0, "sigma must be positive");
        let bound = (6.0 * sigma).round() as i64;
        let rng = ChaCha20Rng::from_entropy();
        let normal = Normal::new(0.0, sigma).unwrap();
        Self {
            sigma,
            bound,
            normal,
            rng,
        }
    }

    /// `sample` returns a sample from the truncated discrete Gaussian
    /// distribution. The technique used is patterned after [lattigo](https://github.com/tuneinsight/lattigo/blob/c031b14be1fb3697945709d7afbed264fa845442/ring/sampler_gaussian.go#L71).
    /// In particular, `sampled_val` is sampled from a normal distribution with
    /// mean 0 and standard deviation `sigma`. If `sampled_val` is within
    /// the bounds, it is rounded and returned.
    pub fn sample(&mut self) -> i64 {
        let sampled_val = self.normal.sample(&mut self.rng);
        if sampled_val.abs() < self.bound as f64 {
            sampled_val.round() as i64
        } else {
            self.sample()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::TruncatedDiscreteGaussianSampler;
    use rand::Rng;
    #[test]
    fn gaussian_sampler_truncation_limit() {
        // Assert that the samples are integers within the bounds of the truncated
        // Gaussian distribution.
        let sigma = rand::thread_rng().gen_range(1.0..=100.0);
        let mut sampler = TruncatedDiscreteGaussianSampler::new(sigma);

        for _ in 0..5000000 {
            let sample = sampler.sample();
            assert!((-sampler.bound..=sampler.bound).contains(&sample));
        }
    }
}
