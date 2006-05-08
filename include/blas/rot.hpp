#ifndef __BLAS_ROT_HPP__
#define __BLAS_ROT_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::rot (const int N, Real *X, const int incX, Real *Y, const int incY,
            const Real c, const Real s)
{
{
  int i;
  int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
  int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
  for (i = 0; i < N; i++) {
    const Real x = X[ix];
    const Real y = Y[iy];
    X[ix] = c * x + s * y;
    Y[iy] = -s * x + c * y;
    ix += incX;
    iy += incY;
  }
}
}

#endif // __BLAS_ROT_HPP__
