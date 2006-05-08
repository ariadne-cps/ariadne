#ifndef __BLAS_ASUM_HPP__
#define __BLAS_ASUM_HPP__

#include "blas.hpp"

template<typename real>
real
BLAS::asum (const int N, const real *X, const int incX)
{
{
  real r = 0.0;
  int i;
  int ix = 0;
  if (incX <= 0) {
    return 0;
  }
  for (i = 0; i < N; i++) {
    r += abs(X[ix]);
    ix += incX;
  }
  return r;
}
}

template<typename real>
real
BLAS::asum (const int N, const complex<real> *X, const int incX)
{
{
  real r = 0.0;
  int i;
  int ix = 0;
  if (incX <= 0) {
    return 0;
  }
  for (i = 0; i < N; i++) {
    r += abs(X[ix]);
    ix += incX;
  }
  return r;
}
}

#endif // __BLAS_ASUM_HPP__
