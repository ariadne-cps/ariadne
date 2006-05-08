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

#ifndef __BLAS_TPSV_HPP__
#define __BLAS_TPSV_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::tpsv (const enum ORDER order, const enum UPLO Uplo,
             const enum TRANSPOSE TransA, const enum DIAG Diag,
             const int N, const Real *Ap, Real *X, const int incX)
{
{
  int i, j;
  const int nonunit = (Diag == NonUnit);
  const int Trans = (TransA != ConjTrans) ? TransA : Trans;
  if (N == 0)
    return;
  if ((order == RowMajor && Trans == NoTrans && Uplo == Upper)
      || (order == ColMajor && Trans == Trans && Uplo == Lower)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + incX * (N - 1);
    if (nonunit) {
      X[ix] = X[ix] / Ap[(((((((N - 1))-1)+1)*(2*(N)-(((N - 1))-1)))/2)+((N - 1))-((N - 1)))];
    }
    ix -= incX;
    for (i = N - 1; i > 0 && i--;) {
      Real tmp = X[ix];
      int jx = ix + incX;
      for (j = i + 1; j < N; j++) {
        const Real Aij = Ap[((((((i)-1)+1)*(2*(N)-((i)-1)))/2)+(j)-(i))];
        tmp -= Aij * X[jx];
        jx += incX;
      }
      if (nonunit) {
        X[ix] = tmp / Ap[((((((i)-1)+1)*(2*(N)-((i)-1)))/2)+(i)-(i))];
      } else {
        X[ix] = tmp;
      }
      ix -= incX;
    }
  } else if ((order == RowMajor && Trans == NoTrans && Uplo == Lower)
             || (order == ColMajor && Trans == Trans && Uplo == Upper)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
    if (nonunit) {
      X[ix] = X[ix] / Ap[(((0)*((0)+1))/2 + (0))];
    }
    ix += incX;
    for (i = 1; i < N; i++) {
      Real tmp = X[ix];
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
      for (j = 0; j < i; j++) {
        const Real Aij = Ap[(((i)*((i)+1))/2 + (j))];
        tmp -= Aij * X[jx];
        jx += incX;
      }
      if (nonunit) {
        X[ix] = tmp / Ap[(((i)*((i)+1))/2 + (j))];
      } else {
        X[ix] = tmp;
      }
      ix += incX;
    }
  } else if ((order == RowMajor && Trans == Trans && Uplo == Upper)
             || (order == ColMajor && Trans == NoTrans && Uplo == Lower)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
    if (nonunit) {
      X[ix] = X[ix] / Ap[((((((0)-1)+1)*(2*(N)-((0)-1)))/2)+(0)-(0))];
    }
    ix += incX;
    for (i = 1; i < N; i++) {
      Real tmp = X[ix];
      int jx = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
      for (j = 0; j < i; j++) {
        const Real Aji = Ap[((((((j)-1)+1)*(2*(N)-((j)-1)))/2)+(i)-(j))];
        tmp -= Aji * X[jx];
        jx += incX;
      }
      if (nonunit) {
        X[ix] = tmp / Ap[((((((i)-1)+1)*(2*(N)-((i)-1)))/2)+(i)-(i))];
      } else {
        X[ix] = tmp;
      }
      ix += incX;
    }
  } else if ((order == RowMajor && Trans == Trans && Uplo == Lower)
             || (order == ColMajor && Trans == NoTrans && Uplo == Upper)) {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX))) + (N - 1) * incX;
    if (nonunit) {
      X[ix] = X[ix] / Ap[((((N - 1))*(((N - 1))+1))/2 + ((N - 1)))];
    }
    ix -= incX;
    for (i = N - 1; i > 0 && i--;) {
      Real tmp = X[ix];
      int jx = ix + incX;
      for (j = i + 1; j < N; j++) {
        const Real Aji = Ap[(((j)*((j)+1))/2 + (i))];
        tmp -= Aji * X[jx];
        jx += incX;
      }
      if (nonunit) {
        X[ix] = tmp / Ap[(((i)*((i)+1))/2 + (i))];
      } else {
        X[ix] = tmp;
      }
      ix -= incX;
    }
  } else {
    xerbla(0, "source_tpsv_r.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_TPSV_HPP__
