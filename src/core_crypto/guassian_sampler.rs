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

/// TODO: this struct is likely to contain further fields in the future
#[derive(Clone, Debug)]
struct GaussianSampler {
    params: TruncatedDiscreteGaussian,
}

impl GaussianSampler {
    /// `new` returns a new `GaussianSampler` with the given parameters starting from the standard deviation `sigma`.
    pub fn new(sigma: f64) -> Self {
        assert!(sigma > 0.0, "sigma must be positive");
        let bound = (6.0 * sigma).round() as i64;
        Self {
            params: TruncatedDiscreteGaussian { sigma, bound },
        }
    }
    /// `norm_f64` returns a normally distributed f64 in
    /// the range [-f64::MAX, f64::MAX], bounds included,
    /// with standard normal distribution (mean = 0, std_dev = 1).
    fn norm_f64(&self) -> f64 {
        let normal = Normal::new(0.0, 1.0).unwrap();
        normal.sample(&mut rand::thread_rng())
    }

    /// `sample` returns a sample from the truncated discrete Gaussian distribution.
    /// The technique used is patterned after [lattigo](https://github.com/tuneinsight/lattigo/blob/c031b14be1fb3697945709d7afbed264fa845442/ring/sampler_gaussian.go#L71).
    /// In particular, `norm` is sampled from the standard normal distribution of mean 0 and standard deviation 1.
    /// The `norm` is then scaled by `sigma` and, if the result is within the bounds, it is rounded and returned.
    pub fn sample(&self) -> i64 {
        let norm = self.norm_f64();
        let scaled = norm * self.params.sigma;
        if scaled.abs() < self.params.bound as f64 {
            scaled.round() as i64
        } else {
            self.sample()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;
    const SIGMA: f64 = 3.2;

    #[test]
    fn gaussian_sampler_truncation_limit() {
        let sampler = GaussianSampler::new(SIGMA);

        for _ in 0..(1 << 25) {
            let norm = sampler.norm_f64();
            assert!((f64::MIN..=f64::MAX).contains(&norm));
            let sample = sampler.sample();
            assert!((-sampler.params.bound..=sampler.params.bound).contains(&sample));
        }
    }

    #[test]
    fn gaussian_sampler_mean() {
        let sampler = GaussianSampler::new(SIGMA);

        let mut sum = 0;
        let n = 1 << 25;
        for _ in 0..n {
            sum += sampler.sample();
        }
        let mean = sum as f64 / n as f64;
        assert!((mean - 0.001..=mean + 0.001).contains(&0.0));
    }

    #[test]
    fn gaussian_pmf_consistency() {
        let sampler = GaussianSampler::new(SIGMA);
        let val = rand::thread_rng().gen_range(-sampler.params.bound..=sampler.params.bound);

        let pmf_expected = (1.0 / ((2.0 * std::f64::consts::PI).sqrt() * sampler.params.sigma))
            * (-0.5 * (val as f64 / sampler.params.sigma).powi(2)).exp();

        let mut count = 0;
        let n = 1 << 25;
        for _ in 0..n {
            if sampler.sample() == val {
                count += 1;
            }
        }
        let pmf_observed = count as f64 / n as f64;
        assert!((pmf_observed - 0.001..=pmf_observed + 0.001).contains(&pmf_expected));
    }
}
