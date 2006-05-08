#ifndef __BLAS_AMAX_HPP__
#define __BLAS_AMAX_HPP__

#include "blas.hpp"

template<typename real>
real
BLAS::amax (const int N, const real *X, const int incX)
{
{
  real mx = 0.0;
  int ix = 0;
  int i;

  if (incX <= 0) {
    return 0;
  }

  for (i = 0; i < N; i++) {
    if (abs(X[ix]) > mx) {
      mx = abs(X[ix]);
    }
    ix += incX;
  }

  return mx;
}
}

template<typename real>
real
BLAS::amax (const int N, const complex<real> *X, const int incX)
{
{
  real mx = 0.0;
  int ix = 0;
  int i;

  if (incX <= 0) {
    return 0;
  }

  for (i = 0; i < N; i++) {
    if (abs(X[ix]) > mx) {
      mx = abs(X[ix]);
    }
    ix += incX;
  }

  return mx;
}
}


#endif // __BLAS_AMAX_HPP__
