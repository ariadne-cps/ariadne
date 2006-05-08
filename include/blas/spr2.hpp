#ifndef __BLAS_SPR2_HPP__
#define __BLAS_SPR2_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::spr2 (const enum ORDER order, const enum UPLO Uplo,
             const int N, const Real alpha, const Real *X, const int incX,
             const Real *Y, const int incY, Real *Ap)
{
{
  int i, j;
  if (N == 0)
    return;
  if (alpha == 0.0)
    return;
  if ((order == RowMajor && Uplo == Upper)
      || (order == ColMajor && Uplo == Lower)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
    int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
    for (i = 0; i < N; i++) {
      const Real tmp1 = alpha * X[ix];
      const Real tmp2 = alpha * Y[iy];
      int jx = ix;
      int jy = iy;
      for (j = i; j < N; j++) {
        Ap[((((((i)-1)+1)*(2*(N)-((i)-1)))/2)+(j)-(i))] += tmp1 * Y[jy] + tmp2 * X[jx];
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
      const Real tmp1 = alpha * X[ix];
      const Real tmp2 = alpha * Y[iy];
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
      int jy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
      for (j = 0; j <= i; j++) {
        Ap[(((i)*((i)+1))/2 + (j))] += tmp1 * Y[jy] + tmp2 * X[jx];
        jx += incX;
        jy += incY;
      }
      ix += incX;
      iy += incY;
    }
  } else {
    xerbla(0, "source_spr2.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_SPR2_HPP__
