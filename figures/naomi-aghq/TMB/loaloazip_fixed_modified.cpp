// loalaozip_fixed_modified.cpp

#include <TMB.hpp>                                // Links in the TMB libraries
//#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(y); // Binomial response
  DATA_VECTOR(N); // Binomial sample size
  // DATA_SPARSE_MATRIX(design); // eta = design * W
  // UPDATE: design = (A,X,A,X)
  DATA_SPARSE_MATRIX(A);
  DATA_MATRIX(X);
  DATA_SCALAR(nu); // Matern shape
  DATA_SCALAR(rho_u);
  DATA_SCALAR(rho_alpha);
  DATA_SCALAR(sigma_u);
  DATA_SCALAR(sigma_alpha);
  DATA_SCALAR(betaprec);
  
  DATA_IVECTOR(W_starts);
  DATA_IVECTOR(W_lengths);
  DATA_INTEGER(i);
  
  DATA_MATRIX(D); // Distance matrix for calculating matern
  
  // Params
  PARAMETER_VECTOR(betarisk);
  PARAMETER_VECTOR(betazi);
  
  PARAMETER(W_i);
  PARAMETER_VECTOR(W_minus_i);
  
  vector<Type> W(380);
  int k = 0;
  for (int j = 0; j < 380; j++) {
    if (j + 1 == i) { // +1 because C++ does zero-indexing
      W(j) = W_i;
    } else {
      W(j) = W_minus_i(k);
      k++;
    }
  }
  
  vector<Type> Urisk = W.segment(W_starts(0), W_lengths(1));
  vector<Type> Uzi = W.segment(W_starts(1), W_lengths(2));

  PARAMETER(logkappa); // Transformed matern params
  PARAMETER(logtau);
  
  // Constants
  int n = y.size();
  double pi = 3.141592653589793115998;
  
  // Transformations
  // UPDATE: linear predictor with new matrices
  // vector<Type> eta = design * W; // dim(eta) = nrow(design) = 2*n;
  vector<Type> eta1 = A*Uzi + X*betazi;
  vector<Type> eta2 = A*Urisk + X*betarisk;
  vector<Type> eta(eta1.size()+eta2.size());
  eta << eta1,eta2;
  // Prob of nonstructural zero
  vector<Type> zeroprob(n);
  // Prob of infection
  vector<Type> infectionprob(n);
  for (int i = 0;i<n;i++) {
    zeroprob(i) = 1.0 / (1.0 + exp(-1.0 * eta(i)));
    infectionprob(i) = 1.0 / (1.0 + exp(-1.0 * eta(i+n)));
  }
  
  // rho and sigma
  Type kappa = exp(logkappa);
  Type tau = exp(logtau);
  Type rho = sqrt(8.0*nu) / kappa;
  Type sigma = tau / ( pow(kappa,nu) * sqrt( exp(lgamma(nu + 1.0)) * (4.0*pi) / exp(lgamma(nu)) ) );
  
  // Log posterior
  Type lp = 0;
  // Prior for sigma,rho. with dimension 2 fixed, formula simplifies
  // the logkappa + logtau at the end is the jacobian
  
  Type lambda1 = -1.0 * (rho_u / sqrt(8.0*nu)) * log(rho_alpha);
  Type lambda2 = ( -1.0 * pow(kappa,-1.0 * nu) * sqrt( exp(lgamma(nu))  / ( exp(lgamma(nu + 1.0)) * (4.0*pi) ) ) ) * log(sigma_alpha) / sigma_u;
  
  Type lpt = log(lambda1) + log(lambda2) - lambda1 * kappa - lambda2 * tau + logkappa + logtau;
  lp += lpt;
  
  // Prior for W
  // N(0,Matern())
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
  // UPDATE again: this is all obsolete now.
  // int d = W.size() - 2;
  // int Wd = d/2;
  // vector<Type> W1(Wd);
  // for (int i=0;i<Wd;i++) W1(i) = W(i);
  // Type beta1 = W(Wd);
  // vector<Type> W2(Wd);
  // for (int i=0;i<Wd;i++) W2(i) = W(i + Wd + 1);
  // Type beta2 = W(2*Wd+1);
  
  // // Now add the priors
  // Type nll1 = density::MVNORM_t<Type>(C)(W1);
  // lp -= nll1; // Negative prior (?)
  // Type nll2 = density::MVNORM_t<Type>(C)(W2);
  // lp -= nll2;
  // // Part for beta
  // Type betapart = -1.0 * log(2.0 * pi) + log(betaprec) - betaprec * 0.5 * (pow(beta1,2.0) + pow(beta2,2.0));
  // lp += betapart;
  // UPDATE: add the priors using the new param names
  Type nll1 = density::MVNORM_t<Type>(C)(Urisk);
  lp -= nll1; // Negative prior (?)
  Type nll2 = density::MVNORM_t<Type>(C)(Uzi);
  lp -= nll2;
  // Part for beta. Note both betas passed as vectors, but will have only one element
  Type betapart = -1.0 * log(2.0 * pi) + log(betaprec) - betaprec * 0.5 * (pow(betarisk(0),2.0) + pow(betazi(0),2.0));
  lp += betapart;
  
  // Log likelihood
  Type L = 0;
  for (int i = 0;i < n;i++) {
    // Zero prob
    if (y(i) == 0) {
      L += 1.0 - zeroprob(i);
    }
    // Binomial variability
    L += zeroprob(i) * exp(lfactorial(N(i)) - lfactorial(y(i)) - lfactorial(N(i) - y(i))) *
      pow(infectionprob(i),y(i)) *
      pow(1.0 - infectionprob(i),N(i) - y(i));
    
    lp += log(L);
    L = 0;
  }
  
  REPORT(rho);
  REPORT(sigma);
  REPORT(kappa);
  REPORT(tau);
  REPORT(eta);
  REPORT(eta1);
  REPORT(eta2);
  REPORT(C);
  REPORT(nll1);
  REPORT(nll2);
  REPORT(betapart);
  REPORT(Urisk);
  REPORT(Uzi);
  REPORT(betarisk);
  REPORT(betazi);
  REPORT(lpt);
  
  // Return NEGATED log posterior
  return -1.0 * lp;
}