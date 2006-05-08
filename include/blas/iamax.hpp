#ifndef __BLAS_IAMAX_HPP__
#define __BLAS_IAMAX_HPP__

#include "blas.hpp"

template<typename scalar>
int
BLAS::iamax (const int N, const scalar *X, const int incX)
{
{
  scalar mx = 0;
  int ix = 0;
  int i;
  int result = 0;

  if (incX <= 0) {
    return 0;
  }

  for (i = 0; i < N; i++) {
    if (abs(X[ix]) > mx) {
      mx = abs(X[ix]);
      result = i;
    }
    ix += incX;
  }

  return result;
}
}

#endif // __BLAS_IAMAX_HPP__
