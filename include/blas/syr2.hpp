#ifndef __BLAS_SYR2_HPP__
#define __BLAS_SYR2_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::syr2 (const enum ORDER order, const enum UPLO Uplo,
             const int N, const Real alpha, const Real *X, const int incX,
             const Real *Y, const int incY, Real *A, const int lda)
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
    int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
    for (i = 0; i < N; i++) {
      const Real tmp1 = alpha * X[ix];
      const Real tmp2 = alpha * Y[iy];
      int jx = ix;
      int jy = iy;
      for (j = i; j < N; j++) {
        A[lda * i + j] += tmp1 * Y[jy] + tmp2 * X[jx];
        jx += incX;
        jy += incY;
      }
      ix += incX;
      iy += incY;
    }
  } else if ((order == RowMajor && Uplo == Lower)
             || (order == ColMajor && Uplo == Upper)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
    int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
    for (i = 0; i < N; i++) {
      const Real tmp1 = alpha * X[ix];
      const Real tmp2 = alpha * Y[iy];
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
      int jy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
      for (j = 0; j <= i; j++) {
        A[lda * i + j] += tmp1 * Y[jy] + tmp2 * X[jx];
        jx += incX;
        jy += incY;
      }
      ix += incX;
      iy += incY;
    }
  } else {
    xerbla(0, "source_syr2.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_SYR2_HPP__
