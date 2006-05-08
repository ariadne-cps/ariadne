#ifndef __BLAS_SETV_HPP__
#define __BLAS_SETV_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::setv (const int N, const Real alpha, Real *X, const int incX)
{
{
  int i;
  int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
  for (i = 0; i < N; i++) {
    X[ix] = alpha;
    ix += incX;
  }
}
}

#endif // __BLAS_SETV_HPP__
