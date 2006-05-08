#ifndef __BLAS_ZHER_HPP__
#define __BLAS_ZHER_HPP__

#include "blas.hpp"

template<typename real>
void
BLAS::her (const enum ORDER order, const enum UPLO Uplo,
            const int N, const real alpha, const complex<real> *X, const int incX,
            complex<real> *A, const int lda)
{
{
  const int conj = (order == ColMajor) ? -1 : 1;
  int i, j;
  if (alpha == 0.0)
    return;
  if ((order == RowMajor && Uplo == Upper)
      || (order == ColMajor && Uplo == Lower)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
    for (i = 0; i < N; i++) {
      const real tmp_real = alpha * (((const real *) X)[2*(ix)]);
      const real tmp_imag = alpha * conj * (((const real *) X)[2*(ix)+1]);
      int jx = ix;
      {
        const real X_real = (((const real *) X)[2*(jx)]);
        const real X_imag = -conj * (((const real *) X)[2*(jx)+1]);
        (((real *) A)[2*(lda * i + i)]) += X_real * tmp_real - X_imag * tmp_imag;
        (((real *) A)[2*(lda * i + i)+1]) = 0;
        jx += incX;
      }
      for (j = i + 1; j < N; j++) {
        const real X_real = (((const real *) X)[2*(jx)]);
        const real X_imag = -conj * (((const real *) X)[2*(jx)+1]);
        (((real *) A)[2*(lda * i + j)]) += X_real * tmp_real - X_imag * tmp_imag;
        (((real *) A)[2*(lda * i + j)+1]) += X_imag * tmp_real + X_real * tmp_imag;
        jx += incX;
      }
      ix += incX;
    }
  } else if ((order == RowMajor && Uplo == Lower)
             || (order == ColMajor && Uplo == Upper)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
    for (i = 0; i < N; i++) {
      const real tmp_real = alpha * (((const real *) X)[2*(ix)]);
      const real tmp_imag = alpha * conj * (((const real *) X)[2*(ix)+1]);
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
      for (j = 0; j < i; j++) {
        const real X_real = (((const real *) X)[2*(jx)]);
        const real X_imag = -conj * (((const real *) X)[2*(jx)+1]);
        (((real *) A)[2*(lda * i + j)]) += X_real * tmp_real - X_imag * tmp_imag;
        (((real *) A)[2*(lda * i + j)+1]) += X_imag * tmp_real + X_real * tmp_imag;
        jx += incX;
      }
      {
        const real X_real = (((const real *) X)[2*(jx)]);
        const real X_imag = -conj * (((const real *) X)[2*(jx)+1]);
        (((real *) A)[2*(lda * i + i)]) += X_real * tmp_real - X_imag * tmp_imag;
        (((real *) A)[2*(lda * i + i)+1]) = 0;
        jx += incX;
      }
      ix += incX;
    }
  } else {
    xerbla(0, "source_her.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_ZHER_HPP__
