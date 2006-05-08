#ifndef __BLAS_TBSV_HPP__
#define __BLAS_TBSV_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::tbsv (const enum ORDER order, const enum UPLO Uplo,
             const enum TRANSPOSE TransA, const enum DIAG Diag,
             const int N, const int K, const Real *A, const int lda,
             Real *X, const int incX)
{
{
  const int nonunit = (Diag == NonUnit);
  int i, j;
  const int Trans = (TransA != ConjTrans) ? TransA : Trans;
  if (N == 0)
    return;
  if ((order == RowMajor && Trans == NoTrans && Uplo == Upper)
      || (order == ColMajor && Trans == Trans && Uplo == Lower)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + incX * (N - 1);
    for (i = N; i > 0 && i--;) {
      Real tmp = X[ix];
      const int j_min = i + 1;
      const int j_max = min(N, i + K + 1);
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + j_min * incX;
      for (j = j_min; j < j_max; j++) {
        const Real Aij = A[lda * i + (j - i)];
        tmp -= Aij * X[jx];
        jx += incX;
      }
      if (nonunit) {
        X[ix] = tmp / A[lda * i + 0];
      } else {
        X[ix] = tmp;
      }
      ix -= incX;
    }
  } else if ((order == RowMajor && Trans == NoTrans && Uplo == Lower)
             || (order == ColMajor && Trans == Trans && Uplo == Upper)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
    for (i = 0; i < N; i++) {
      Real tmp = X[ix];
      const int j_min = (i > K ? i - K : 0);
      const int j_max = i;
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + j_min * incX;
      for (j = j_min; j < j_max; j++) {
        const Real Aij = A[lda * i + (K + j - i)];
        tmp -= Aij * X[jx];
        jx += incX;
      }
      if (nonunit) {
        X[ix] = tmp / A[lda * i + K];
      } else {
        X[ix] = tmp;
      }
      ix += incX;
    }
  } else if ((order == RowMajor && Trans == Trans && Uplo == Upper)
             || (order == ColMajor && Trans == NoTrans && Uplo == Lower)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
    for (i = 0; i < N; i++) {
      Real tmp = X[ix];
      const int j_min = (K > i ? 0 : i - K);
      const int j_max = i;
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + j_min * incX;
      for (j = j_min; j < j_max; j++) {
        const Real Aji = A[(i - j) + lda * j];
        tmp -= Aji * X[jx];
        jx += incX;
      }
      if (nonunit) {
        X[ix] = tmp / A[0 + lda * i];
      } else {
        X[ix] = tmp;
      }
      ix += incX;
    }
  } else if ((order == RowMajor && Trans == Trans && Uplo == Lower)
             || (order == ColMajor && Trans == NoTrans && Uplo == Upper)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + (N - 1) * incX;
    for (i = N; i > 0 && i--;) {
      Real tmp = X[ix];
      const int j_min = i + 1;
      const int j_max = min(N, i + K + 1);
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + j_min * incX;
      for (j = j_min; j < j_max; j++) {
        const Real Aji = A[(K + i - j) + lda * j];
        tmp -= Aji * X[jx];
        jx += incX;
      }
      if (nonunit) {
        X[ix] = tmp / A[K + lda * i];
      } else {
        X[ix] = tmp;
      }
      ix -= incX;
    }
  } else {
    xerbla(0, "source_tbsv_r.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_TBSV_HPP__
