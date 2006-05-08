#ifndef __BLAS_SET_HPP__
#define __BLAS_SET_HPP__

#include "blas.hpp"

template<typename scalarA, typename scalarX>
void
BLAS::set (const int N, const scalarA alpha, scalarX *X, const int incX)
{
  int i;
  int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
  for (i = 0; i < N; i++) {
    X[ix] = alpha;
    ix += incX;
  }
}

#endif // __BLAS_SET_HPP__
