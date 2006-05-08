#ifndef __BLAS_ZHPMV_HPP__
#define __BLAS_ZHPMV_HPP__

#include "blas.hpp"

template<typename real>
void
BLAS::hpmv (const enum ORDER order, const enum UPLO Uplo,
             const int N, const complex<real> alpha, const complex<real> *Ap, const complex<real> *X,
             const int incX, const complex<real> beta, complex<real> *Y, const int incY)
{
{
  const int conj = (order == ColMajor) ? -1 : 1;
  int i, j;
  const real alpha_real = alpha.real(); 
  const real alpha_imag = alpha.imag(); 
  const real beta_real = beta.real();
  const real beta_imag = beta.imag(); 
  if ((alpha_real == 0.0 && alpha_imag == 0.0)
      && (beta_real == 1.0 && beta_imag == 0.0))
    return;
  if (beta_real == 0.0 && beta_imag == 0.0) {
    int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
    for (i = 0; i < N; i++) {
      (((real *) Y)[2*(iy)]) = 0.0;
      (((real *) Y)[2*(iy)+1]) = 0.0;
      iy += incY;
    }
  } else if (!(beta_real == 1.0 && beta_imag == 0.0)) {
    int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
    for (i = 0; i < N; i++) {
      const real y_real = (((real *) Y)[2*(iy)]);
      const real y_imag = (((real *) Y)[2*(iy)+1]);
      const real tmpR = y_real * beta_real - y_imag * beta_imag;
      const real tmpI = y_real * beta_imag + y_imag * beta_real;
      (((real *) Y)[2*(iy)]) = tmpR;
      (((real *) Y)[2*(iy)+1]) = tmpI;
      iy += incY;
    }
  }
  if (alpha_real == 0.0 && alpha_imag == 0.0)
    return;
  if ((order == RowMajor && Uplo == Upper)
      || (order == ColMajor && Uplo == Lower)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
    int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
    for (i = 0; i < N; i++) {
      real x_real = (((const real *) X)[2*(ix)]);
      real x_imag = (((const real *) X)[2*(ix)+1]);
      real temp1_real = alpha_real * x_real - alpha_imag * x_imag;
      real temp1_imag = alpha_real * x_imag + alpha_imag * x_real;
      real temp2_real = 0.0;
      real temp2_imag = 0.0;
      const int j_min = i + 1;
      const int j_max = N;
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + j_min * incX;
      int jy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY))) + j_min * incY;
      real Aii_real = (((const real *) Ap)[2*(((((((i)-1)+1)*(2*(N)-((i)-1)))/2)+(i)-(i)))]);
      (((real *) Y)[2*(iy)]) += temp1_real * Aii_real;
      (((real *) Y)[2*(iy)+1]) += temp1_imag * Aii_real;
      for (j = j_min; j < j_max; j++) {
        real Aij_real = (((const real *) Ap)[2*(((((((i)-1)+1)*(2*(N)-((i)-1)))/2)+(j)-(i)))]);
        real Aij_imag = conj * (((const real *) Ap)[2*(((((((i)-1)+1)*(2*(N)-((i)-1)))/2)+(j)-(i)))+1]);
        (((real *) Y)[2*(jy)]) += temp1_real * Aij_real - temp1_imag * (-Aij_imag);
        (((real *) Y)[2*(jy)+1]) += temp1_real * (-Aij_imag) + temp1_imag * Aij_real;
        x_real = (((const real *) X)[2*(jx)]);
        x_imag = (((const real *) X)[2*(jx)+1]);
        temp2_real += x_real * Aij_real - x_imag * Aij_imag;
        temp2_imag += x_real * Aij_imag + x_imag * Aij_real;
        jx += incX;
        jy += incY;
      }
      (((real *) Y)[2*(iy)]) += alpha_real * temp2_real - alpha_imag * temp2_imag;
      (((real *) Y)[2*(iy)+1]) += alpha_real * temp2_imag + alpha_imag * temp2_real;
      ix += incX;
      iy += incY;
    }
  } else if ((order == RowMajor && Uplo == Lower)
             || (order == ColMajor && Uplo == Upper)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
    int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
    for (i = 0; i < N; i++) {
      real x_real = (((const real *) X)[2*(ix)]);
      real x_imag = (((const real *) X)[2*(ix)+1]);
      real temp1_real = alpha_real * x_real - alpha_imag * x_imag;
      real temp1_imag = alpha_real * x_imag + alpha_imag * x_real;
      real temp2_real = 0.0;
      real temp2_imag = 0.0;
      const int j_min = 0;
      const int j_max = i;
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + j_min * incX;
      int jy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY))) + j_min * incY;
      real Aii_real = (((const real *) Ap)[2*((((i)*((i)+1))/2 + (i)))]);
      (((real *) Y)[2*(iy)]) += temp1_real * Aii_real;
      (((real *) Y)[2*(iy)+1]) += temp1_imag * Aii_real;
      for (j = j_min; j < j_max; j++) {
        real Aij_real = (((const real *) Ap)[2*((((i)*((i)+1))/2 + (j)))]);
        real Aij_imag = conj * (((const real *) Ap)[2*((((i)*((i)+1))/2 + (j)))+1]);
        (((real *) Y)[2*(jy)]) += temp1_real * Aij_real - temp1_imag * (-Aij_imag);
        (((real *) Y)[2*(jy)+1]) += temp1_real * (-Aij_imag) + temp1_imag * Aij_real;
        x_real = (((const real *) X)[2*(jx)]);
        x_imag = (((const real *) X)[2*(jx)+1]);
        temp2_real += x_real * Aij_real - x_imag * Aij_imag;
        temp2_imag += x_real * Aij_imag + x_imag * Aij_real;
        jx += incX;
        jy += incY;
      }
      (((real *) Y)[2*(iy)]) += alpha_real * temp2_real - alpha_imag * temp2_imag;
      (((real *) Y)[2*(iy)+1]) += alpha_real * temp2_imag + alpha_imag * temp2_real;
      ix += incX;
      iy += incY;
    }
  } else {
    xerbla(0, "source_hpmv.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_ZHPMV_HPP__
