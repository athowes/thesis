// epil.stan

data {
  int<lower=0> N;        // Number of patients
  int<lower=0> J;        // Number of clinic visits
  int<lower=0> K;        // Number of predictors (inc. intercept)
  matrix[N * J, K] X;    // Design matrix
  int<lower=0> y[N * J]; // Outcome variable
  matrix[N * J, N] E;    // Epsilon matrix
}

parameters {
  vector[K] beta;            // Vector of coefficients
  vector[N] epsilon;         // Patient specific errors
  vector[N * J] nu;          // Patient-visit errors
  real<lower=0> tau_epsilon; // Precision of epsilon
  real<lower=0> tau_nu;      // Precision of nu
}

transformed parameters {
  vector[N * J] eta = X * beta + nu + E * epsilon;
}

model {
  beta ~ normal(0, 100);
  tau_epsilon ~ gamma(0.001, 0.001);
  tau_nu ~ gamma(0.001, 0.001);
  epsilon ~ normal(0, sqrt(1 / tau_epsilon));
  nu ~ normal(0, sqrt(1 / tau_nu));
  y ~ poisson_log(eta);
}