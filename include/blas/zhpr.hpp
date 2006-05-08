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

#ifndef __BLAS_ZHPR_HPP__
#define __BLAS_ZHPR_HPP__

#include "blas.hpp"

template<typename real>
void
BLAS::hpr (const enum ORDER order, const enum UPLO Uplo,
            const int N, const real alpha, const complex<real> *X, const int incX,
            complex<real> *Ap)
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
        (((real *) Ap)[2*(((((((i)-1)+1)*(2*(N)-((i)-1)))/2)+(i)-(i)))]) += X_real * tmp_real - X_imag * tmp_imag;
        (((real *) Ap)[2*(((((((i)-1)+1)*(2*(N)-((i)-1)))/2)+(i)-(i)))+1]) = 0;
        jx += incX;
      }
      for (j = i + 1; j < N; j++) {
        const real X_real = (((const real *) X)[2*(jx)]);
        const real X_imag = -conj * (((const real *) X)[2*(jx)+1]);
        (((real *) Ap)[2*(((((((i)-1)+1)*(2*(N)-((i)-1)))/2)+(j)-(i)))]) += X_real * tmp_real - X_imag * tmp_imag;
        (((real *) Ap)[2*(((((((i)-1)+1)*(2*(N)-((i)-1)))/2)+(j)-(i)))+1]) += X_imag * tmp_real + X_real * tmp_imag;
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
        (((real *) Ap)[2*((((i)*((i)+1))/2 + (j)))]) += X_real * tmp_real - X_imag * tmp_imag;
        (((real *) Ap)[2*((((i)*((i)+1))/2 + (j)))+1]) += X_imag * tmp_real + X_real * tmp_imag;
        jx += incX;
      }
      {
        const real X_real = (((const real *) X)[2*(jx)]);
        const real X_imag = -conj * (((const real *) X)[2*(jx)+1]);
        (((real *) Ap)[2*((((i)*((i)+1))/2 + (i)))]) += X_real * tmp_real - X_imag * tmp_imag;
        (((real *) Ap)[2*((((i)*((i)+1))/2 + (i)))+1]) = 0;
        jx += incX;
      }
      ix += incX;
    }
  } else {
    xerbla(0, "source_hpr.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_ZHPR_HPP__
