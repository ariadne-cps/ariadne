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

#ifndef __BLAS_GBMV_HPP__
#define __BLAS_GBMV_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::gbmv (const enum ORDER order, const enum TRANSPOSE TransA,
             const int M, const int N, const int KL, const int KU,
             const Real alpha, const Real *A, const int lda,
             const Real *X, const int incX, const Real beta, Real *Y,
             const int incY)
{
{
  int i, j;
  int lenX, lenY, L, U;
  const int Trans = (TransA != ConjTrans) ? TransA : Trans;
  if (M == 0 || N == 0)
    return;
  if (alpha == 0.0 && beta == 1.0)
    return;
  if (Trans == NoTrans) {
    lenX = N;
    lenY = M;
    L = KL;
    U = KU;
  } else {
    lenX = M;
    lenY = N;
    L = KU;
    U = KL;
  }
  if (beta == 0.0) {
    int iy = offset(lenY,incY);
    for (i = 0; i < lenY; i++) {
      Y[iy] = 0;
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
      const int j_min = (i > L ? i - L : 0);
      const int j_max = min(lenX, i + U + 1);
      int jx = offset(lenX,incX) + j_min * incX;
      for (j = j_min; j < j_max; j++) {
        temp += X[jx] * A[(L - i + j) + i * lda];
        jx += incX;
      }
      Y[iy] += alpha * temp;
      iy += incY;
    }
  } else if ((order == RowMajor && Trans == Trans)
             || (order == ColMajor && Trans == NoTrans)) {
    int jx = offset(lenX,incX);
    for (j = 0; j < lenX; j++) {
      const Real temp = alpha * X[jx];
      if (temp != 0.0) {
        const int i_min = (j > U ? j - U : 0);
        const int i_max = min(lenY, j + L + 1);
        int iy = offset(lenY,incY) + i_min * incY;
        for (i = i_min; i < i_max; i++) {
          Y[iy] += temp * A[lda * j + (U + i - j)];
          iy += incY;
        }
      }
      jx += incX;
    }
  } else {
    xerbla(0, "source_gbmv_r.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_GBMV_HPP__
