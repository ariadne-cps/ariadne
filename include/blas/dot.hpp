#ifndef __BLAS_DOT_HPP__
#define __BLAS_DOT_HPP__

#include "blas.hpp"

template<typename Real>
Real
BLAS::dot (const int N, const Real *X, const int incX, const Real *Y,
            const int incY)
{
{
  Real r = 0.0;
  int i;
  int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
  int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
  for (i = 0; i < N; i++) {
    r += X[ix] * Y[iy];
    ix += incX;
    iy += incY;
  }
  return r;
}
}

#endif // __BLAS_DOT_HPP__
