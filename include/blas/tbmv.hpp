#ifndef __BLAS_TBMV_HPP__
#define __BLAS_TBMV_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::tbmv (const enum ORDER order, const enum UPLO Uplo,
             const enum TRANSPOSE TransA, const enum DIAG Diag,
             const int N, const int K, const Real *A, const int lda,
             Real *X, const int incX)
{
{
  int i, j;
  const int nonunit = (Diag == NonUnit);
  const int Trans = (TransA != ConjTrans) ? TransA : Trans;
  if (N == 0)
    return;
  if ((order == RowMajor && Trans == NoTrans && Uplo == Upper)
      || (order == ColMajor && Trans == Trans && Uplo == Lower)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
    for (i = 0; i < N; i++) {
      Real temp = (nonunit ? A[lda * i + 0] : 1.0) * X[ix];
      const int j_min = i + 1;
      const int j_max = min(N, i + K + 1);
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + j_min * incX;
      for (j = j_min; j < j_max; j++) {
        temp += X[jx] * A[lda * i + (j - i)];
        jx += incX;
      }
      X[ix] = temp;
      ix += incX;
    }
  } else if ((order == RowMajor && Trans == NoTrans && Uplo == Lower)
             || (order == ColMajor && Trans == Trans && Uplo == Upper)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + (N - 1) * incX;
    for (i = N; i > 0 && i--;) {
      Real temp = (nonunit ? A[lda * i + K] : 1.0) * X[ix];
      const int j_min = (i > K ? i - K : 0);
      const int j_max = i;
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + j_min * incX;
      for (j = j_min; j < j_max; j++) {
        temp += X[jx] * A[lda * i + (K - i + j)];
        jx += incX;
      }
      X[ix] = temp;
      ix -= incX;
    }
  } else if ((order == RowMajor && Trans == Trans && Uplo == Upper)
             || (order == ColMajor && Trans == NoTrans && Uplo == Lower)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + (N - 1) * incX;
    for (i = N; i > 0 && i--;) {
      Real temp = 0.0;
      const int j_min = (K > i ? 0 : i - K);
      const int j_max = i;
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + j_min * incX;
      for (j = j_min; j < j_max; j++) {
        temp += X[jx] * A[lda * j + (i - j)];
        jx += incX;
      }
      if (nonunit) {
        X[ix] = temp + X[ix] * A[lda * i + 0];
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
      const int j_max = min(N, i + K + 1);
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + j_min * incX;
      for (j = j_min; j < j_max; j++) {
        temp += X[jx] * A[lda * j + (K - j + i)];
        jx += incX;
      }
      if (nonunit) {
        X[ix] = temp + X[ix] * A[lda * i + K];
      } else {
        X[ix] += temp;
      }
      ix += incX;
    }
  }
}
}

#endif // __BLAS_TBMV_HPP__
