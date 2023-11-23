// loalaozip.cpp

#include <TMB.hpp>
// #include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(y);             // Binomial response
  DATA_VECTOR(N);             // Binomial sample size
  DATA_SPARSE_MATRIX(design); // eta = design * W
  DATA_SCALAR(nu);            // Matern shape
  DATA_SCALAR(rho_u);
  DATA_SCALAR(rho_alpha);
  DATA_SCALAR(sigma_u);
  DATA_SCALAR(sigma_alpha);
  DATA_SCALAR(betaprec);
  
  DATA_MATRIX(D); // Distance matrix for calculating Matern
  
  // Params
  PARAMETER_VECTOR(W); // W = (U, V, beta); eta = design * W
  PARAMETER(logkappa); // Transformed Matern parameters
  PARAMETER(logtau);
  
  // Constants
  int n = y.size();
  double pi = 3.141592653589793115998;
  
  // Transformations
  vector<Type> eta = design * W; // dim(eta) = nrow(design) = 2 * n
  
  // Prob of nonstructural zero
  vector<Type> zeroprob(n);
  // Prob of infection
  vector<Type> infectionprob(n);
  for (int i = 0;i<n;i++) {
    zeroprob(i) = 1.0 / (1.0 + exp(-1.0 * eta(i)));
    infectionprob(i) = 1.0 / (1.0 + exp(-1.0 * eta(i + n)));
  }
  
  // rho and sigma
  Type kappa = exp(logkappa);
  Type tau = exp(logtau);
  Type rho = sqrt(8.0*nu) / kappa;
  Type sigma = tau / (pow(kappa,nu) * sqrt(exp(lgamma(nu + 1.0)) * (4.0 * pi) / exp(lgamma(nu))));

  // Log posterior
  Type lp = 0;
  // Prior for sigma, rho with dimension 2 fixed, formula simplifies
  // the logkappa + logtau at the end is the Jacobian
  
  Type lambda1 = -1.0 * (rho_u / sqrt(8.0 * nu)) * log(rho_alpha);
  Type lambda2 = (-1.0 * pow(kappa,-1.0 * nu) * sqrt(exp(lgamma(nu))  / (exp(lgamma(nu + 1.0)) * (4.0*pi)))) * log(sigma_alpha) / sigma_u;
    
  Type lpt = log(lambda1) + log(lambda2) - lambda1 * kappa - lambda2 * tau + logkappa + logtau;
  lp += lpt;
  
  // Prior for W
  // N(0, Matern())
  // From documentation: https://kaskr.github.io/adcomp/matern_8cpp-example.html
  // Incorporate the sqrt(8nu) difference...
  matrix<Type> C(D);
  for(int i=0; i<C.rows(); i++)
    for(int j=0; j<C.cols(); j++)
      C(i,j) = pow(sigma,2.0) * matern(D(i,j), rho / sqrt(8.0 * nu), nu);
  // Now split up the parameter vector
  // W[1:n] = zeroprob part
  // W[n+1] = beta for zeroprob part
  // W[(n+2):(2*n+1)] = infectprob part
  // W[2*n+2] = beta for infectprob part
  // UPDATE: this is true for the full model but only because ncol(design) = 2*n + 2.
  // Need to use ncol(design)
  int d = W.size() - 2;
  int Wd = d/2;
  vector<Type> W1(Wd);
  for (int i=0;i<Wd;i++) W1(i) = W(i);
  Type beta1 = W(Wd);
  vector<Type> W2(Wd);
  for (int i=0;i<Wd;i++) W2(i) = W(i + Wd + 1);
  Type beta2 = W(2*Wd+1);
  
  // Now add the priors
  Type nll1 = density::MVNORM_t<Type>(C)(W1);
  lp -= nll1; // Negative prior (?)
  Type nll2 = density::MVNORM_t<Type>(C)(W2);
  lp -= nll2;
  // Part for beta: soft fixed to MLEs
  lp += dnorm(beta1, Type(2.95), Type(0.001), true);
  lp += dnorm(beta2, Type(-1.99), Type(0.001), true);
  
  // Log likelihood
  Type L = 0;
  for (int i = 0;i < n;i++) {
    // Zero probability
    if (y(i) == 0) {
      L += 1.0 - zeroprob(i);
    }
    // Binomial variability
    L += zeroprob(i) * exp(lfactorial(N(i)) - lfactorial(y(i)) - lfactorial(N(i) - y(i))) *
      pow(infectionprob(i),y(i)) *
      pow(1.0 - infectionprob(i), N(i) - y(i));
    
    lp += log(L);
    L = 0;
  }
  
  REPORT(rho);
  REPORT(sigma);
  REPORT(kappa);
  REPORT(tau);
  REPORT(eta);
  REPORT(C);
  REPORT(nll1);
  REPORT(nll2);
  REPORT(betapart);
  REPORT(W1);
  REPORT(W2);
  REPORT(beta1);
  REPORT(beta2);
  REPORT(lpt);
  
  // Return NEGATED log posterior
  return -1.0 * lp;
}
