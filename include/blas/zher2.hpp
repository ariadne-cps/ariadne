#ifndef __BLAS_ZHER2_HPP__
#define __BLAS_ZHER2_HPP__

#include "blas.hpp"

template<typename real>
void
BLAS::her2 (const enum ORDER order, const enum UPLO Uplo,
             const int N, const complex<real> alpha, const complex<real> *X, const int incX,
             const complex<real> *Y, const int incY, complex<real> *A, const int lda)
{
{
  const int conj = (order == ColMajor) ? -1 : 1;
  int i, j;
  const real alpha_real = alpha.real(); 
  const real alpha_imag = alpha.imag(); 
  if (alpha_real == 0.0 && alpha_imag == 0.0)
    return;
  if ((order == RowMajor && Uplo == Upper)
      || (order == ColMajor && Uplo == Lower)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
    int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
    for (i = 0; i < N; i++) {
      const real Xi_real = (((const real *) X)[2*(ix)]);
      const real Xi_imag = (((const real *) X)[2*(ix)+1]);
      const real tmp1_real = alpha_real * Xi_real - alpha_imag * Xi_imag;
      const real tmp1_imag = alpha_imag * Xi_real + alpha_real * Xi_imag;
      const real Yi_real = (((const real *) Y)[2*(iy)]);
      const real Yi_imag = (((const real *) Y)[2*(iy)+1]);
      const real tmp2_real = alpha_real * Yi_real + alpha_imag * Yi_imag;
      const real tmp2_imag = -alpha_imag * Yi_real + alpha_real * Yi_imag;
      int jx = ix + incX;
      int jy = iy + incY;
      (((real *) A)[2*(lda * i + i)]) += 2 * (tmp1_real * Yi_real + tmp1_imag * Yi_imag);
      (((real *) A)[2*(lda * i + i)+1]) = 0;
      for (j = i + 1; j < N; j++) {
        const real Xj_real = (((const real *) X)[2*(jx)]);
        const real Xj_imag = (((const real *) X)[2*(jx)+1]);
        const real Yj_real = (((const real *) Y)[2*(jy)]);
        const real Yj_imag = (((const real *) Y)[2*(jy)+1]);
        (((real *) A)[2*(lda * i + j)]) += ((tmp1_real * Yj_real + tmp1_imag * Yj_imag)
                                 + (tmp2_real * Xj_real + tmp2_imag * Xj_imag));
        (((real *) A)[2*(lda * i + j)+1]) +=
            conj * ((tmp1_imag * Yj_real - tmp1_real * Yj_imag) +
                    (tmp2_imag * Xj_real - tmp2_real * Xj_imag));
        jx += incX;
        jy += incY;
      }
      ix += incX;
      iy += incY;
    }
  } else if ((order == RowMajor && Uplo == Lower)
             || (order == ColMajor && Uplo == Upper)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
    int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
    for (i = 0; i < N; i++) {
      const real Xi_real = (((const real *) X)[2*(ix)]);
      const real Xi_imag = (((const real *) X)[2*(ix)+1]);
      const real tmp1_real = alpha_real * Xi_real - alpha_imag * Xi_imag;
      const real tmp1_imag = alpha_imag * Xi_real + alpha_real * Xi_imag;
      const real Yi_real = (((const real *) Y)[2*(iy)]);
      const real Yi_imag = (((const real *) Y)[2*(iy)+1]);
      const real tmp2_real = alpha_real * Yi_real + alpha_imag * Yi_imag;
      const real tmp2_imag = -alpha_imag * Yi_real + alpha_real * Yi_imag;
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
      int jy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
      for (j = 0; j < i; j++) {
        const real Xj_real = (((const real *) X)[2*(jx)]);
        const real Xj_imag = (((const real *) X)[2*(jx)+1]);
        const real Yj_real = (((const real *) Y)[2*(jy)]);
        const real Yj_imag = (((const real *) Y)[2*(jy)+1]);
        (((real *) A)[2*(lda * i + j)]) += ((tmp1_real * Yj_real + tmp1_imag * Yj_imag)
                                 + (tmp2_real * Xj_real + tmp2_imag * Xj_imag));
        (((real *) A)[2*(lda * i + j)+1]) +=
            conj * ((tmp1_imag * Yj_real - tmp1_real * Yj_imag) +
                    (tmp2_imag * Xj_real - tmp2_real * Xj_imag));
        jx += incX;
        jy += incY;
      }
      (((real *) A)[2*(lda * i + i)]) += 2 * (tmp1_real * Yi_real + tmp1_imag * Yi_imag);
      (((real *) A)[2*(lda * i + i)+1]) = 0;
      ix += incX;
      iy += incY;
    }
  } else {
    xerbla(0, "source_her2.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_ZHER2_HPP__
