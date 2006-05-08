#ifndef __BLAS_NRM2_HPP__
#define __BLAS_NRM2_HPP__

#include "blas.hpp"

template<typename real>
real
BLAS::nrm2 (const int N, const real *X, const int incX)
{
{
  real scale = 0.0;
  real ssq = 1.0;
  int i;
  int ix = 0;
  if (N <= 0 || incX <= 0) {
    return 0;
  } else if (N == 1) {
    return abs(X[0]);
  }
  for (i = 0; i < N; i++) {
    const real x = X[ix];
    if (x != 0.0) {
      const real ax = abs(x);
      if (scale < ax) {
        ssq = 1.0 + ssq * (scale / ax) * (scale / ax);
        scale = ax;
      } else {
        ssq += (ax / scale) * (ax / scale);
      }
    }
    ix += incX;
  }
  return scale * sqrt(ssq);
}
}

template<typename real>
real
BLAS::nrm2 (const int N, const complex<real> *X, const int incX)
{
{
  real scale = 0.0;
  real ssq = 1.0;
  int i;
  int ix = 0;
  if (N <= 0 || incX <= 0) {
    return 0;
  } else if (N == 1) {
    return abs(X[0]);
  }
  for (i = 0; i < N; i++) {
    const complex<real> x = X[ix];
    if (x != 0.0) {
      const real ax = abs(x);
      if (scale < ax) {
        ssq = 1.0 + ssq * (scale / ax) * (scale / ax);
        scale = ax;
      } else {
        ssq += (ax / scale) * (ax / scale);
      }
    }
    ix += incX;
  }
  return scale * sqrt(ssq);
}
}

#endif // __BLAS_NRM2_HPP__
