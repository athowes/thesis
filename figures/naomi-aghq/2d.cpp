// 2d.cpp

#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  PARAMETER(theta1);
  PARAMETER(theta2);

  Type nll;
  nll = Type(0.0);

  Type phi1 = 0.5 * theta1;
  Type phi2 = 0.8 * theta1 - 0.5 * theta2;

  nll -= dsn(phi1, Type(2.0), true);
  nll -= dsn(phi2, Type(-2.0), true);

  return(nll);
}
