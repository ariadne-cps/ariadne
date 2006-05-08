#ifndef __BLAS_SWAP_HPP__
#define __BLAS_SWAP_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::swap (const int N, Real *X, const int incX, Real *Y,
             const int incY)
{
{
  int i;
  int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
  int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
  for (i = 0; i < N; i++) {
    const Real tmp = X[ix];
    X[ix] = Y[iy];
    Y[iy] = tmp;
    ix += incX;
    iy += incY;
  }
}
}

#endif // __BLAS_SWAP_HPP__
