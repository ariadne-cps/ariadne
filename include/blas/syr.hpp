#ifndef __BLAS_SYR_HPP__
#define __BLAS_SYR_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::syr (const enum ORDER order, const enum UPLO Uplo,
            const int N, const Real alpha, const Real *X, const int incX,
            Real *A, const int lda)
{
{
  int i, j;
  if (N == 0)
    return;
  if (alpha == 0.0)
    return;
  if ((order == RowMajor && Uplo == Upper)
      || (order == ColMajor && Uplo == Lower)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
    for (i = 0; i < N; i++) {
      const Real tmp = alpha * X[ix];
      int jx = ix;
      for (j = i; j < N; j++) {
        A[lda * i + j] += X[jx] * tmp;
        jx += incX;
      }
      ix += incX;
    }
  } else if ((order == RowMajor && Uplo == Lower)
             || (order == ColMajor && Uplo == Upper)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
    for (i = 0; i < N; i++) {
      const Real tmp = alpha * X[ix];
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
      for (j = 0; j <= i; j++) {
        A[lda * i + j] += X[jx] * tmp;
        jx += incX;
      }
      ix += incX;
    }
  } else {
    xerbla(0, "source_syr.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_SYR_HPP__
