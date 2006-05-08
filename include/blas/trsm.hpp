#ifndef __BLAS_TRSM_HPP__
#define __BLAS_TRSM_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::trsm (const enum ORDER Order, const enum SIDE Side,
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
    if (alpha != 1.0) {
      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          B[ldb * i + j] *= alpha;
        }
      }
    }
    for (i = n1; i > 0 && i--;) {
      if (nonunit) {
        Real Aii = A[lda * i + i];
        for (j = 0; j < n2; j++) {
          B[ldb * i + j] /= Aii;
        }
      }
      for (k = 0; k < i; k++) {
        const Real Aki = A[k * lda + i];
        for (j = 0; j < n2; j++) {
          B[ldb * k + j] -= Aki * B[ldb * i + j];
        }
      }
    }
  } else if (side == Left && uplo == Upper && trans == Trans) {
    if (alpha != 1.0) {
      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          B[ldb * i + j] *= alpha;
        }
      }
    }
    for (i = 0; i < n1; i++) {
      if (nonunit) {
        Real Aii = A[lda * i + i];
        for (j = 0; j < n2; j++) {
          B[ldb * i + j] /= Aii;
        }
      }
      for (k = i + 1; k < n1; k++) {
        const Real Aik = A[i * lda + k];
        for (j = 0; j < n2; j++) {
          B[ldb * k + j] -= Aik * B[ldb * i + j];
        }
      }
    }
  } else if (side == Left && uplo == Lower && trans == NoTrans) {
    if (alpha != 1.0) {
      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          B[ldb * i + j] *= alpha;
        }
      }
    }
    for (i = 0; i < n1; i++) {
      if (nonunit) {
        Real Aii = A[lda * i + i];
        for (j = 0; j < n2; j++) {
          B[ldb * i + j] /= Aii;
        }
      }
      for (k = i + 1; k < n1; k++) {
        const Real Aki = A[k * lda + i];
        for (j = 0; j < n2; j++) {
          B[ldb * k + j] -= Aki * B[ldb * i + j];
        }
      }
    }
  } else if (side == Left && uplo == Lower && trans == Trans) {
    if (alpha != 1.0) {
      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          B[ldb * i + j] *= alpha;
        }
      }
    }
    for (i = n1; i > 0 && i--;) {
      if (nonunit) {
        Real Aii = A[lda * i + i];
        for (j = 0; j < n2; j++) {
          B[ldb * i + j] /= Aii;
        }
      }
      for (k = 0; k < i; k++) {
        const Real Aik = A[i * lda + k];
        for (j = 0; j < n2; j++) {
          B[ldb * k + j] -= Aik * B[ldb * i + j];
        }
      }
    }
  } else if (side == Right && uplo == Upper && trans == NoTrans) {
    if (alpha != 1.0) {
      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          B[ldb * i + j] *= alpha;
        }
      }
    }
    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        if (nonunit) {
          Real Ajj = A[lda * j + j];
          B[ldb * i + j] /= Ajj;
        }
        {
          Real Bij = B[ldb * i + j];
          for (k = j + 1; k < n2; k++) {
            B[ldb * i + k] -= A[j * lda + k] * Bij;
          }
        }
      }
    }
  } else if (side == Right && uplo == Upper && trans == Trans) {
    if (alpha != 1.0) {
      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          B[ldb * i + j] *= alpha;
        }
      }
    }
    for (i = 0; i < n1; i++) {
      for (j = n2; j > 0 && j--;) {
        if (nonunit) {
          Real Ajj = A[lda * j + j];
          B[ldb * i + j] /= Ajj;
        }
        {
          Real Bij = B[ldb * i + j];
          for (k = 0; k < j; k++) {
            B[ldb * i + k] -= A[k * lda + j] * Bij;
          }
        }
      }
    }
  } else if (side == Right && uplo == Lower && trans == NoTrans) {
    if (alpha != 1.0) {
      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          B[ldb * i + j] *= alpha;
        }
      }
    }
    for (i = 0; i < n1; i++) {
      for (j = n2; j > 0 && j--;) {
        if (nonunit) {
          Real Ajj = A[lda * j + j];
          B[ldb * i + j] /= Ajj;
        }
        {
          Real Bij = B[ldb * i + j];
          for (k = 0; k < j; k++) {
            B[ldb * i + k] -= A[j * lda + k] * Bij;
          }
        }
      }
    }
  } else if (side == Right && uplo == Lower && trans == Trans) {
    if (alpha != 1.0) {
      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          B[ldb * i + j] *= alpha;
        }
      }
    }
    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        if (nonunit) {
          Real Ajj = A[lda * j + j];
          B[ldb * i + j] /= Ajj;
        }
        {
          Real Bij = B[ldb * i + j];
          for (k = j + 1; k < n2; k++) {
            B[ldb * i + k] -= A[k * lda + j] * Bij;
          }
        }
      }
    }
  } else {
    xerbla(0, "source_trsm_r.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_TRSM_HPP__
