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

#ifndef __BLAS_ZGERU_HPP__
#define __BLAS_ZGERU_HPP__

#include "blas.hpp"

template<typename real>
void
BLAS::geru (const enum ORDER order, const int M, const int N,
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
        const real Y_imag = (((const real *) Y)[2*(jy)+1]);
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
      const real Y_imag = (((const real *) Y)[2*(jy)+1]);
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
    xerbla(0, "source_geru.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_ZGERU_HPP__
