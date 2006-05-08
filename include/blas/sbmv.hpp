/*
 * Copyright (C) 2006 Pieter Collins <Pieter.Collins@cwi.nl>
 *
 * Based on the BLAS implementation in Gnu Scientific Library 1.8
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 */

/*  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef __BLAS_SBMV_HPP__
#define __BLAS_SBMV_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::sbmv (const enum ORDER order, const enum UPLO Uplo,
             const int N, const int K, const Real alpha, const Real *A,
             const int lda, const Real *X, const int incX,
             const Real beta, Real *Y, const int incY)
{
{
  int i, j;
  if (N == 0)
    return;
  if (alpha == 0.0 && beta == 1.0)
    return;
  if (beta == 0.0) {
    int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
    for (i = 0; i < N; i++) {
      Y[iy] = 0.0;
      iy += incY;
    }
  } else if (beta != 1.0) {
    int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
    for (i = 0; i < N; i++) {
      Y[iy] *= beta;
      iy += incY;
    }
  }
  if (alpha == 0.0)
    return;
  if ((order == RowMajor && Uplo == Upper)
      || (order == ColMajor && Uplo == Lower)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
    int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
    for (i = 0; i < N; i++) {
      Real tmp1 = alpha * X[ix];
      Real tmp2 = 0.0;
      const int j_min = i + 1;
      const int j_max = min(N, i + K + 1);
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + j_min * incX;
      int jy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY))) + j_min * incY;
      Y[iy] += tmp1 * A[0 + i * lda];
      for (j = j_min; j < j_max; j++) {
        Real Aij = A[(j - i) + i * lda];
        Y[jy] += tmp1 * Aij;
        tmp2 += Aij * X[jx];
        jx += incX;
        jy += incY;
      }
      Y[iy] += alpha * tmp2;
      ix += incX;
      iy += incY;
    }
  } else if ((order == RowMajor && Uplo == Lower)
             || (order == ColMajor && Uplo == Upper)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
    int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
    for (i = 0; i < N; i++) {
      Real tmp1 = alpha * X[ix];
      Real tmp2 = 0.0;
      const int j_min = (i > K) ? i - K : 0;
      const int j_max = i;
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + j_min * incX;
      int jy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY))) + j_min * incY;
      for (j = j_min; j < j_max; j++) {
        Real Aij = A[(K - i + j) + i * lda];
        Y[jy] += tmp1 * Aij;
        tmp2 += Aij * X[jx];
        jx += incX;
        jy += incY;
      }
      Y[iy] += tmp1 * A[K + i * lda] + alpha * tmp2;
      ix += incX;
      iy += incY;
    }
  } else {
    xerbla(0, "source_sbmv.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_SBMV_HPP__
