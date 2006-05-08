#ifndef __BLAS_GEMV_HPP__
#define __BLAS_GEMV_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::gemv (const enum ORDER order, const enum TRANSPOSE TransA,
             const int M, const int N, const Real alpha, const Real *A,
             const int lda, const Real *X, const int incX,
             const Real beta, Real *Y, const int incY)
{
{
  int i, j;
  int lenX, lenY;
  const int Trans = (TransA != ConjTrans) ? TransA : Trans;
  if (M == 0 || N == 0)
    return;
  if (alpha == 0.0 && beta == 1.0)
    return;
  if (Trans == NoTrans) {
    lenX = N;
    lenY = M;
  } else {
    lenX = M;
    lenY = N;
  }
  if (beta == 0.0) {
    int iy = offset(lenY,incY);
    for (i = 0; i < lenY; i++) {
      Y[iy] = 0.0;
      iy += incY;
    }
  } else if (beta != 1.0) {
    int iy = offset(lenY,incY);
    for (i = 0; i < lenY; i++) {
      Y[iy] *= beta;
      iy += incY;
    }
  }
  if (alpha == 0.0)
    return;
  if ((order == RowMajor && Trans == NoTrans)
      || (order == ColMajor && Trans == Trans)) {
    int iy = offset(lenY,incY);
    for (i = 0; i < lenY; i++) {
      Real temp = 0.0;
      int ix = offset(lenX,incX);
      for (j = 0; j < lenX; j++) {
        temp += X[ix] * A[lda * i + j];
        ix += incX;
      }
      Y[iy] += alpha * temp;
      iy += incY;
    }
  } else if ((order == RowMajor && Trans == Trans)
             || (order == ColMajor && Trans == NoTrans)) {
    int ix = offset(lenX,incX);
    for (j = 0; j < lenX; j++) {
      const Real temp = alpha * X[ix];
      if (temp != 0.0) {
        int iy = offset(lenY,incY);
        for (i = 0; i < lenY; i++) {
          Y[iy] += temp * A[lda * j + i];
          iy += incY;
        }
      }
      ix += incX;
    }
  } else {
    xerbla(0, "source_gemv_r.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_GEMV_HPP__
