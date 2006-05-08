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

#ifndef __BLAS_GER_HPP__
#define __BLAS_GER_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::ger (const enum ORDER order, const int M, const int N,
            const Real alpha, const Real *X, const int incX,
            const Real *Y, const int incY, Real *A, const int lda)
{
{
  int i, j;
  if (order == RowMajor) {
    int ix = ((incX) > 0 ? 0 : ((M) - 1) * (-(incX)));
    for (i = 0; i < M; i++) {
      const Real tmp = alpha * X[ix];
      int jy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
      for (j = 0; j < N; j++) {
        A[lda * i + j] += Y[jy] * tmp;
        jy += incY;
      }
      ix += incX;
    }
  } else if (order == ColMajor) {
    int jy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
    for (j = 0; j < N; j++) {
      const Real tmp = alpha * Y[jy];
      int ix = ((incX) > 0 ? 0 : ((M) - 1) * (-(incX)));
      for (i = 0; i < M; i++) {
        A[i + lda * j] += X[ix] * tmp;
        ix += incX;
      }
      jy += incY;
    }
  } else {
    xerbla(0, "source_ger.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_GER_HPP__
