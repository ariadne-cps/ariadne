#ifndef __BLAS_TRSV_HPP__
#define __BLAS_TRSV_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::trsv (const enum ORDER order, const enum UPLO Uplo,
             const enum TRANSPOSE TransA, const enum DIAG Diag,
             const int N, const Real *A, const int lda, Real *X,
             const int incX)
{
{
  const int nonunit = (Diag == NonUnit);
  int ix, jx;
  int i, j;
  const int Trans = (TransA != ConjTrans) ? TransA : Trans;
  if (N == 0)
    return;
  if ((order == RowMajor && Trans == NoTrans && Uplo == Upper)
      || (order == ColMajor && Trans == Trans && Uplo == Lower)) {
    ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + incX * (N - 1);
    if (nonunit) {
      X[ix] = X[ix] / A[lda * (N - 1) + (N - 1)];
    }
    ix -= incX;
    for (i = N - 1; i > 0 && i--;) {
      Real tmp = X[ix];
      jx = ix + incX;
      for (j = i + 1; j < N; j++) {
        const Real Aij = A[lda * i + j];
        tmp -= Aij * X[jx];
        jx += incX;
      }
      if (nonunit) {
        X[ix] = tmp / A[lda * i + i];
      } else {
        X[ix] = tmp;
      }
      ix -= incX;
    }
  } else if ((order == RowMajor && Trans == NoTrans && Uplo == Lower)
             || (order == ColMajor && Trans == Trans && Uplo == Upper)) {
    ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
    if (nonunit) {
      X[ix] = X[ix] / A[lda * 0 + 0];
    }
    ix += incX;
    for (i = 1; i < N; i++) {
      Real tmp = X[ix];
      jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
      for (j = 0; j < i; j++) {
        const Real Aij = A[lda * i + j];
        tmp -= Aij * X[jx];
        jx += incX;
      }
      if (nonunit) {
        X[ix] = tmp / A[lda * i + i];
      } else {
        X[ix] = tmp;
      }
      ix += incX;
    }
  } else if ((order == RowMajor && Trans == Trans && Uplo == Upper)
             || (order == ColMajor && Trans == NoTrans && Uplo == Lower)) {
    ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
    if (nonunit) {
      X[ix] = X[ix] / A[lda * 0 + 0];
    }
    ix += incX;
    for (i = 1; i < N; i++) {
      Real tmp = X[ix];
      jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
      for (j = 0; j < i; j++) {
        const Real Aji = A[lda * j + i];
        tmp -= Aji * X[jx];
        jx += incX;
      }
      if (nonunit) {
        X[ix] = tmp / A[lda * i + i];
      } else {
        X[ix] = tmp;
      }
      ix += incX;
    }
  } else if ((order == RowMajor && Trans == Trans && Uplo == Lower)
             || (order == ColMajor && Trans == NoTrans && Uplo == Upper)) {
    ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + (N - 1) * incX;
    if (nonunit) {
      X[ix] = X[ix] / A[lda * (N - 1) + (N - 1)];
    }
    ix -= incX;
    for (i = N - 1; i > 0 && i--;) {
      Real tmp = X[ix];
      jx = ix + incX;
      for (j = i + 1; j < N; j++) {
        const Real Aji = A[lda * j + i];
        tmp -= Aji * X[jx];
        jx += incX;
      }
      if (nonunit) {
        X[ix] = tmp / A[lda * i + i];
      } else {
        X[ix] = tmp;
      }
      ix -= incX;
    }
  } else {
    xerbla(0, "source_trsv_r.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_TRSV_HPP__
