#ifndef __BLAS_SYR2K_HPP__
#define __BLAS_SYR2K_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::syr2k (const enum ORDER Order, const enum UPLO Uplo,
              const enum TRANSPOSE Trans, const int N, const int K,
              const Real alpha, const Real *A, const int lda,
              const Real *B, const int ldb, const Real beta, Real *C,
              const int ldc)
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
          temp += (A[i * lda + k] * B[j * ldb + k]
                   + B[i * ldb + k] * A[j * lda + k]);
        }
        C[i * ldc + j] += alpha * temp;
      }
    }
  } else if (uplo == Upper && trans == Trans) {
    for (k = 0; k < K; k++) {
      for (i = 0; i < N; i++) {
        Real temp1 = alpha * A[k * lda + i];
        Real temp2 = alpha * B[k * ldb + i];
        for (j = i; j < N; j++) {
          C[i * lda + j] += temp1 * B[k * ldb + j] + temp2 * A[k * lda + j];
        }
      }
    }
  } else if (uplo == Lower && trans == NoTrans) {
    for (i = 0; i < N; i++) {
      for (j = 0; j <= i; j++) {
        Real temp = 0.0;
        for (k = 0; k < K; k++) {
          temp += (A[i * lda + k] * B[j * ldb + k]
                   + B[i * ldb + k] * A[j * lda + k]);
        }
        C[i * ldc + j] += alpha * temp;
      }
    }
  } else if (uplo == Lower && trans == Trans) {
    for (k = 0; k < K; k++) {
      for (i = 0; i < N; i++) {
        Real temp1 = alpha * A[k * lda + i];
        Real temp2 = alpha * B[k * ldb + i];
        for (j = 0; j <= i; j++) {
          C[i * lda + j] += temp1 * B[k * ldb + j] + temp2 * A[k * lda + j];
        }
      }
    }
  } else {
    xerbla(0, "source_syr2k_r.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_SYR2K_HPP__
