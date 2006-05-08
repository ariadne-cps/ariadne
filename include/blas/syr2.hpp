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

#ifndef __BLAS_SYR2_HPP__
#define __BLAS_SYR2_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::syr2 (const enum ORDER order, const enum UPLO Uplo,
             const int N, const Real alpha, const Real *X, const int incX,
             const Real *Y, const int incY, Real *A, const int lda)
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
        A[lda * i + j] += tmp1 * Y[jy] + tmp2 * X[jx];
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
        A[lda * i + j] += tmp1 * Y[jy] + tmp2 * X[jx];
        jx += incX;
        jy += incY;
      }
      ix += incX;
      iy += incY;
    }
  } else {
    xerbla(0, "source_syr2.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_SYR2_HPP__
