#ifndef __BLAS_TRMM_HPP__
#define __BLAS_TRMM_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::trmm (const enum ORDER Order, const enum SIDE Side,
             const enum UPLO Uplo, const enum TRANSPOSE TransA,
             const enum DIAG Diag, const int M, const int N,
             const Real alpha, const Real *A, const int lda, Real *B,
             const int ldb)
{
{
  int i, j, k;
  int n1, n2;
  const int nonunit = (Diag == NonUnit);
  int side, uplo, trans;
  if (Order == RowMajor) {
    n1 = M;
    n2 = N;
    side = Side;
    uplo = Uplo;
    trans = (TransA == ConjTrans) ? Trans : TransA;
  } else {
    n1 = N;
    n2 = M;
    side = (Side == Left) ? Right : Left;
    uplo = (Uplo == Upper) ? Lower : Upper;
    trans = (TransA == ConjTrans) ? Trans : TransA;
  }
  if (side == Left && uplo == Upper && trans == NoTrans) {
    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        Real temp = 0.0;
        if (nonunit) {
          temp = A[i * lda + i] * B[i * ldb + j];
        } else {
          temp = B[i * ldb + j];
        }
        for (k = i + 1; k < n1; k++) {
          temp += A[lda * i + k] * B[k * ldb + j];
        }
        B[ldb * i + j] = alpha * temp;
      }
    }
  } else if (side == Left && uplo == Upper && trans == Trans) {
    for (i = n1; i > 0 && i--;) {
      for (j = 0; j < n2; j++) {
        Real temp = 0.0;
        for (k = 0; k < i; k++) {
          temp += A[lda * k + i] * B[k * ldb + j];
        }
        if (nonunit) {
          temp += A[i * lda + i] * B[i * ldb + j];
        } else {
          temp += B[i * ldb + j];
        }
        B[ldb * i + j] = alpha * temp;
      }
    }
  } else if (side == Left && uplo == Lower && trans == NoTrans) {
    for (i = n1; i > 0 && i--;) {
      for (j = 0; j < n2; j++) {
        Real temp = 0.0;
        for (k = 0; k < i; k++) {
          temp += A[lda * i + k] * B[k * ldb + j];
        }
        if (nonunit) {
          temp += A[i * lda + i] * B[i * ldb + j];
        } else {
          temp += B[i * ldb + j];
        }
        B[ldb * i + j] = alpha * temp;
      }
    }
  } else if (side == Left && uplo == Lower && trans == Trans) {
    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        Real temp = 0.0;
        if (nonunit) {
          temp = A[i * lda + i] * B[i * ldb + j];
        } else {
          temp = B[i * ldb + j];
        }
        for (k = i + 1; k < n1; k++) {
          temp += A[lda * k + i] * B[k * ldb + j];
        }
        B[ldb * i + j] = alpha * temp;
      }
    }
  } else if (side == Right && uplo == Upper && trans == NoTrans) {
    for (i = 0; i < n1; i++) {
      for (j = n2; j > 0 && j--;) {
        Real temp = 0.0;
        for (k = 0; k < j; k++) {
          temp += A[lda * k + j] * B[i * ldb + k];
        }
        if (nonunit) {
          temp += A[j * lda + j] * B[i * ldb + j];
        } else {
          temp += B[i * ldb + j];
        }
        B[ldb * i + j] = alpha * temp;
      }
    }
  } else if (side == Right && uplo == Upper && trans == Trans) {
    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        Real temp = 0.0;
        if (nonunit) {
          temp = A[j * lda + j] * B[i * ldb + j];
        } else {
          temp = B[i * ldb + j];
        }
        for (k = j + 1; k < n2; k++) {
          temp += A[lda * j + k] * B[i * ldb + k];
        }
        B[ldb * i + j] = alpha * temp;
      }
    }
  } else if (side == Right && uplo == Lower && trans == NoTrans) {
    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        Real temp = 0.0;
        if (nonunit) {
          temp = A[j * lda + j] * B[i * ldb + j];
        } else {
          temp = B[i * ldb + j];
        }
        for (k = j + 1; k < n2; k++) {
          temp += A[lda * k + j] * B[i * ldb + k];
        }
        B[ldb * i + j] = alpha * temp;
      }
    }
  } else if (side == Right && uplo == Lower && trans == Trans) {
    for (i = 0; i < n1; i++) {
      for (j = n2; j > 0 && j--;) {
        Real temp = 0.0;
        for (k = 0; k < j; k++) {
          temp += A[lda * j + k] * B[i * ldb + k];
        }
        if (nonunit) {
          temp += A[j * lda + j] * B[i * ldb + j];
        } else {
          temp += B[i * ldb + j];
        }
        B[ldb * i + j] = alpha * temp;
      }
    }
  } else {
    xerbla(0, "source_trmm_r.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_TRMM_HPP__
