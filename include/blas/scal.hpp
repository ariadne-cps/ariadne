#ifndef __BLAS_SCAL_HPP__
#define __BLAS_SCAL_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::scal (const int N, const Real alpha, Real *X, const int incX)
{
{
  int i;
  int ix;
  if (incX <= 0) {
    return;
  }
  ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
  for (i = 0; i < N; i++) {
    X[ix] *= alpha;
    ix += incX;
  }
}
}

#endif // __BLAS_SCAL_HPP__
