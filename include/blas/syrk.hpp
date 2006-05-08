#ifndef __BLAS_SYRK_HPP__
#define __BLAS_SYRK_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::syrk (const enum ORDER Order, const enum UPLO Uplo,
             const enum TRANSPOSE Trans, const int N, const int K,
             const Real alpha, const Real *A, const int lda,
             const Real beta, Real *C, const int ldc)
{
{
  int i, j, k;
  int uplo, trans;
  if (alpha == 0.0 && beta == 1.0)
    return;
  if (Order == RowMajor) {
    uplo = Uplo;
    trans = (Trans == ConjTrans) ? Trans : Trans;
  } else {
    uplo = (Uplo == Upper) ? Lower : Upper;
    if (Trans == Trans || Trans == ConjTrans) {
      trans = NoTrans;
    } else {
      trans = Trans;
    }
  }
  if (beta == 0.0) {
    if (uplo == Upper) {
      for (i = 0; i < N; i++) {
        for (j = i; j < N; j++) {
          C[ldc * i + j] = 0.0;
        }
      }
    } else {
      for (i = 0; i < N; i++) {
        for (j = 0; j <= i; j++) {
          C[ldc * i + j] = 0.0;
        }
      }
    }
  } else if (beta != 1.0) {
    if (uplo == Upper) {
      for (i = 0; i < N; i++) {
        for (j = i; j < N; j++) {
          C[ldc * i + j] *= beta;
        }
      }
    } else {
      for (i = 0; i < N; i++) {
        for (j = 0; j <= i; j++) {
          C[ldc * i + j] *= beta;
        }
      }
    }
  }
  if (alpha == 0.0)
    return;
  if (uplo == Upper && trans == NoTrans) {
    for (i = 0; i < N; i++) {
      for (j = i; j < N; j++) {
        Real temp = 0.0;
        for (k = 0; k < K; k++) {
          temp += A[i * lda + k] * A[j * lda + k];
        }
        C[i * ldc + j] += alpha * temp;
      }
    }
  } else if (uplo == Upper && trans == Trans) {
    for (i = 0; i < N; i++) {
      for (j = i; j < N; j++) {
        Real temp = 0.0;
        for (k = 0; k < K; k++) {
          temp += A[k * lda + i] * A[k * lda + j];
        }
        C[i * ldc + j] += alpha * temp;
      }
    }
  } else if (uplo == Lower && trans == NoTrans) {
    for (i = 0; i < N; i++) {
      for (j = 0; j <= i; j++) {
        Real temp = 0.0;
        for (k = 0; k < K; k++) {
          temp += A[i * lda + k] * A[j * lda + k];
        }
        C[i * ldc + j] += alpha * temp;
      }
    }
  } else if (uplo == Lower && trans == Trans) {
    for (i = 0; i < N; i++) {
      for (j = 0; j <= i; j++) {
        Real temp = 0.0;
        for (k = 0; k < K; k++) {
          temp += A[k * lda + i] * A[k * lda + j];
        }
        C[i * ldc + j] += alpha * temp;
      }
    }
  } else {
    xerbla(0, "source_syrk_r.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_SYRK_HPP__
