data {
  int<lower=0> N;      // Number of data points
  int<lower=0> y[N];   // Observed counts
  real<lower=0> a;     // Shape parameter for gamma prior
  real<lower=0> b;     // Rate parameter for gamma prior
}

parameters {
  real<lower=0> phi; // Rate parameter of Poisson distribution
}

model {
  // Prior distribution
  phi ~ gamma(a, b);

  // Likelihood
  for (i in 1:N) {
    y[i] ~ poisson(phi);
  }
}
