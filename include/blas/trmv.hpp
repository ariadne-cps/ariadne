#ifndef __BLAS_TRMV_HPP__
#define __BLAS_TRMV_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::trmv (const enum ORDER order, const enum UPLO Uplo,
             const enum TRANSPOSE TransA, const enum DIAG Diag,
             const int N, const Real *A, const int lda, Real *X,
             const int incX)
{
{
  int i, j;
  const int nonunit = (Diag == NonUnit);
  const int Trans = (TransA != ConjTrans) ? TransA : Trans;
  if ((order == RowMajor && Trans == NoTrans && Uplo == Upper)
      || (order == ColMajor && Trans == Trans && Uplo == Lower)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
    for (i = 0; i < N; i++) {
      Real temp = 0.0;
      const int j_min = i + 1;
      const int j_max = N;
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + j_min * incX;
      for (j = j_min; j < j_max; j++) {
        temp += X[jx] * A[lda * i + j];
        jx += incX;
      }
      if (nonunit) {
        X[ix] = temp + X[ix] * A[lda * i + i];
      } else {
        X[ix] += temp;
      }
      ix += incX;
    }
  } else if ((order == RowMajor && Trans == NoTrans && Uplo == Lower)
             || (order == ColMajor && Trans == Trans && Uplo == Upper)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + (N - 1) * incX;
    for (i = N; i > 0 && i--;) {
      Real temp = 0.0;
      const int j_min = 0;
      const int j_max = i;
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + j_min * incX;
      for (j = j_min; j < j_max; j++) {
        temp += X[jx] * A[lda * i + j];
        jx += incX;
      }
      if (nonunit) {
        X[ix] = temp + X[ix] * A[lda * i + i];
      } else {
        X[ix] += temp;
      }
      ix -= incX;
    }
  } else if ((order == RowMajor && Trans == Trans && Uplo == Upper)
             || (order == ColMajor && Trans == NoTrans && Uplo == Lower)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + (N - 1) * incX;
    for (i = N; i > 0 && i--;) {
      Real temp = 0.0;
      const int j_min = 0;
      const int j_max = i;
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + j_min * incX;
      for (j = j_min; j < j_max; j++) {
        temp += X[jx] * A[lda * j + i];
        jx += incX;
      }
      if (nonunit) {
        X[ix] = temp + X[ix] * A[lda * i + i];
      } else {
        X[ix] += temp;
      }
      ix -= incX;
    }
  } else if ((order == RowMajor && Trans == Trans && Uplo == Lower)
             || (order == ColMajor && Trans == NoTrans && Uplo == Upper)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
    for (i = 0; i < N; i++) {
      Real temp = 0.0;
      const int j_min = i + 1;
      const int j_max = N;
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + (i + 1) * incX;
      for (j = j_min; j < j_max; j++) {
        temp += X[jx] * A[lda * j + i];
        jx += incX;
      }
      if (nonunit) {
        X[ix] = temp + X[ix] * A[lda * i + i];
      } else {
        X[ix] += temp;
      }
      ix += incX;
    }
  } else {
    xerbla(0, "source_trmv_r.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_TRMV_HPP__
