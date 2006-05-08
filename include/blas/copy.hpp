#ifndef __BLAS_COPY_HPP__
#define __BLAS_COPY_HPP__

#include "blas.hpp"

template<typename scalarX, typename scalarY>
void
BLAS::copy (const int N, const scalarX *X, const int incX, scalarY *Y,
            const int incY)
{
  int i;
  int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
  int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
  for (i = 0; i < N; i++) {
    Y[iy] = X[ix];
    ix += incX;
    iy += incY;
  }
}

#endif // __BLAS_COPY_HPP__
