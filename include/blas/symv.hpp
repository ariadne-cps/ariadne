#ifndef __BLAS_SYMV_HPP__
#define __BLAS_SYMV_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::symv (const enum ORDER order, const enum UPLO Uplo,
             const int N, const Real alpha, const Real *A, const int lda,
             const Real *X, const int incX, const Real beta, Real *Y,
             const int incY)
{
{
  int i, j;
  if (alpha == 0.0 && beta == 1.0)
    return;
  if (beta == 0.0) {
    int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
    for (i = 0; i < N; i++) {
      Y[iy] = 0.0;
      iy += incY;
    }
  } else if (beta != 1.0) {
    int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
    for (i = 0; i < N; i++) {
      Y[iy] *= beta;
      iy += incY;
    }
  }
  if (alpha == 0.0)
    return;
  if ((order == RowMajor && Uplo == Upper)
      || (order == ColMajor && Uplo == Lower)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
    int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
    for (i = 0; i < N; i++) {
      Real temp1 = alpha * X[ix];
      Real temp2 = 0.0;
      const int j_min = i + 1;
      const int j_max = N;
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + j_min * incX;
      int jy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY))) + j_min * incY;
      Y[iy] += temp1 * A[lda * i + i];
      for (j = j_min; j < j_max; j++) {
        Y[jy] += temp1 * A[lda * i + j];
        temp2 += X[jx] * A[lda * i + j];
        jx += incX;
        jy += incY;
      }
      Y[iy] += alpha * temp2;
      ix += incX;
      iy += incY;
    }
  } else if ((order == RowMajor && Uplo == Lower)
             || (order == ColMajor && Uplo == Upper)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + (N - 1) * incX;
    int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY))) + (N - 1) * incY;
    for (i = N; i > 0 && i--;) {
      Real temp1 = alpha * X[ix];
      Real temp2 = 0.0;
      const int j_min = 0;
      const int j_max = i;
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + j_min * incX;
      int jy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY))) + j_min * incY;
      Y[iy] += temp1 * A[lda * i + i];
      for (j = j_min; j < j_max; j++) {
        Y[jy] += temp1 * A[lda * i + j];
        temp2 += X[jx] * A[lda * i + j];
        jx += incX;
        jy += incY;
      }
      Y[iy] += alpha * temp2;
      ix -= incX;
      iy -= incY;
    }
  } else {
    xerbla(0, "source_symv.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_SYMV_HPP__
