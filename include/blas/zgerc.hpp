#ifndef __BLAS_ZGERC_HPP__
#define __BLAS_ZGERC_HPP__

#include "blas.hpp"

template<typename real>
void
BLAS::gerc (const enum ORDER order, const int M, const int N,
             const complex<real> alpha, const complex<real> *X, const int incX, const complex<real> *Y,
             const int incY, complex<real> *A, const int lda)
{
{
  int i, j;
  const real alpha_real = alpha.real(); 
  const real alpha_imag = alpha.imag(); 
  if (order == RowMajor) {
    int ix = ((incX) > 0 ? 0 : ((M) - 1) * (-(incX)));
    for (i = 0; i < M; i++) {
      const real X_real = (((const real *) X)[2*(ix)]);
      const real X_imag = (((const real *) X)[2*(ix)+1]);
      const real tmp_real = alpha_real * X_real - alpha_imag * X_imag;
      const real tmp_imag = alpha_imag * X_real + alpha_real * X_imag;
      int jy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
      for (j = 0; j < N; j++) {
        const real Y_real = (((const real *) Y)[2*(jy)]);
        const real Y_imag = -(((const real *) Y)[2*(jy)+1]);
        (((real *) A)[2*(lda * i + j)]) += Y_real * tmp_real - Y_imag * tmp_imag;
        (((real *) A)[2*(lda * i + j)+1]) += Y_imag * tmp_real + Y_real * tmp_imag;
        jy += incY;
      }
      ix += incX;
    }
  } else if (order == ColMajor) {
    int jy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
    for (j = 0; j < N; j++) {
      const real Y_real = (((const real *) Y)[2*(jy)]);
      const real Y_imag = -(((const real *) Y)[2*(jy)+1]);
      const real tmp_real = alpha_real * Y_real - alpha_imag * Y_imag;
      const real tmp_imag = alpha_imag * Y_real + alpha_real * Y_imag;
      int ix = ((incX) > 0 ? 0 : ((M) - 1) * (-(incX)));
      for (i = 0; i < M; i++) {
        const real X_real = (((const real *) X)[2*(ix)]);
        const real X_imag = (((const real *) X)[2*(ix)+1]);
        (((real *) A)[2*(i + lda * j)]) += X_real * tmp_real - X_imag * tmp_imag;
        (((real *) A)[2*(i + lda * j)+1]) += X_imag * tmp_real + X_real * tmp_imag;
        ix += incX;
      }
      jy += incY;
    }
  } else {
    xerbla(0, "source_gerc.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_ZGERC_HPP__
