#ifndef __BLAS_GER_HPP__
#define __BLAS_GER_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::ger (const enum ORDER order, const int M, const int N,
            const Real alpha, const Real *X, const int incX,
            const Real *Y, const int incY, Real *A, const int lda)
{
{
  int i, j;
  if (order == RowMajor) {
    int ix = ((incX) > 0 ? 0 : ((M) - 1) * (-(incX)));
    for (i = 0; i < M; i++) {
      const Real tmp = alpha * X[ix];
      int jy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
      for (j = 0; j < N; j++) {
        A[lda * i + j] += Y[jy] * tmp;
        jy += incY;
      }
      ix += incX;
    }
  } else if (order == ColMajor) {
    int jy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
    for (j = 0; j < N; j++) {
      const Real tmp = alpha * Y[jy];
      int ix = ((incX) > 0 ? 0 : ((M) - 1) * (-(incX)));
      for (i = 0; i < M; i++) {
        A[i + lda * j] += X[ix] * tmp;
        ix += incX;
      }
      jy += incY;
    }
  } else {
    xerbla(0, "source_ger.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_GER_HPP__
