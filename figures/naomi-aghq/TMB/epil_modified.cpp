// epil_modified.cpp

#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_INTEGER(N);
  DATA_INTEGER(J);
  DATA_INTEGER(K);
  DATA_MATRIX(X);
  DATA_VECTOR(y);
  DATA_MATRIX(E); // Epsilon matrix
  
  DATA_IVECTOR(x_starts);  // Start index of each subvector of x
  DATA_IVECTOR(x_lengths); // Length of each subvector of x
  DATA_INTEGER(i);         // Index i
  
  PARAMETER(x_i);
  PARAMETER_VECTOR(x_minus_i);
  
  vector<Type> x(301);
  int k = 0;
  for (int j = 0; j < 301; j++) {
    if (j + 1 == i) { // +1 because C++ does zero-indexing
      x(j) = x_i;
    } else {
      x(j) = x_minus_i(k);
      k++;
    }
  }
  
  vector<Type> beta = x.segment(x_starts(0), x_lengths(0));
  vector<Type> epsilon = x.segment(x_starts(1), x_lengths(1));
  vector<Type> nu = x.segment(x_starts(2), x_lengths(2));

  PARAMETER(l_tau_epsilon);
  PARAMETER(l_tau_nu);
  
  Type tau_epsilon = exp(l_tau_epsilon);
  Type tau_nu = exp(l_tau_nu);
  Type sigma_epsilon = sqrt(1 / tau_epsilon);
  Type sigma_nu = sqrt(1 / tau_nu);
  vector<Type> eta(X * beta + nu + E * epsilon);
  vector<Type> lambda(exp(eta));
  
  Type nll;
  nll = Type(0.0);
  
  // Note: dgamma() is parameterised as (shape, scale)
  // R-INLA is parameterised as (shape, rate)
  nll -= dlgamma(l_tau_epsilon, Type(0.001),
                 Type(1.0 / 0.001), true);
  nll -= dlgamma(l_tau_nu, Type(0.001), Type(1.0 / 0.001), true);
  nll -= dnorm(epsilon, Type(0), sigma_epsilon, true).sum();
  nll -= dnorm(nu, Type(0), sigma_nu, true).sum();
  nll -= dnorm(beta, Type(0), Type(100), true).sum();
  
  nll -= dpois(y, lambda, true).sum();
  
  ADREPORT(tau_epsilon);
  ADREPORT(tau_nu);
  
  return(nll);
}