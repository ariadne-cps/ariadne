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

#ifndef __BLAS_ZHERK_HPP__
#define __BLAS_ZHERK_HPP__

#include "blas.hpp"

template<typename real>
void
BLAS::herk (const enum ORDER Order, const enum UPLO Uplo,
             const enum TRANSPOSE Trans, const int N, const int K,
             const real alpha, const complex<real> *A, const int lda,
             const real beta, complex<real> *C, const int ldc)
{
{
  int i, j, k;
  int uplo, trans;
  if (beta == 1.0 && (alpha == 0.0 || K == 0))
    return;
  if (Order == RowMajor) {
    uplo = Uplo;
    trans = Trans;
  } else {
    uplo = (Uplo == Upper) ? Lower : Upper;
    trans = (Trans == NoTrans) ? ConjTrans : NoTrans;
  }
  if (beta == 0.0) {
    if (uplo == Upper) {
      for (i = 0; i < N; i++) {
        for (j = i; j < N; j++) {
          (((real *) C)[2*(ldc * i + j)]) = 0.0;
          (((real *) C)[2*(ldc * i + j)+1]) = 0.0;
        }
      }
    } else {
      for (i = 0; i < N; i++) {
        for (j = 0; j <= i; j++) {
          (((real *) C)[2*(ldc * i + j)]) = 0.0;
          (((real *) C)[2*(ldc * i + j)+1]) = 0.0;
        }
      }
    }
  } else if (beta != 1.0) {
    if (uplo == Upper) {
      for (i = 0; i < N; i++) {
        (((real *) C)[2*(ldc * i + i)]) *= beta;
        (((real *) C)[2*(ldc * i + i)+1]) = 0;
        for (j = i + 1; j < N; j++) {
          (((real *) C)[2*(ldc * i + j)]) *= beta;
          (((real *) C)[2*(ldc * i + j)+1]) *= beta;
        }
      }
    } else {
      for (i = 0; i < N; i++) {
        for (j = 0; j < i; j++) {
          (((real *) C)[2*(ldc * i + j)]) *= beta;
          (((real *) C)[2*(ldc * i + j)+1]) *= beta;
        }
        (((real *) C)[2*(ldc * i + i)]) *= beta;
        (((real *) C)[2*(ldc * i + i)+1]) = 0;
      }
    }
  } else {
    for (i = 0; i < N; i++) {
      (((real *) C)[2*(ldc * i + i)+1]) = 0.0;
    }
  }
  if (alpha == 0.0)
    return;
  if (uplo == Upper && trans == NoTrans) {
    for (i = 0; i < N; i++) {
      for (j = i; j < N; j++) {
        real temp_real = 0.0;
        real temp_imag = 0.0;
        for (k = 0; k < K; k++) {
          const real Aik_real = (((const real *) A)[2*(i * lda + k)]);
          const real Aik_imag = (((const real *) A)[2*(i * lda + k)+1]);
          const real Ajk_real = (((const real *) A)[2*(j * lda + k)]);
          const real Ajk_imag = -(((const real *) A)[2*(j * lda + k)+1]);
          temp_real += Aik_real * Ajk_real - Aik_imag * Ajk_imag;
          temp_imag += Aik_real * Ajk_imag + Aik_imag * Ajk_real;
        }
        (((real *) C)[2*(i * ldc + j)]) += alpha * temp_real;
        (((real *) C)[2*(i * ldc + j)+1]) += alpha * temp_imag;
      }
    }
  } else if (uplo == Upper && trans == ConjTrans) {
    for (i = 0; i < N; i++) {
      for (j = i; j < N; j++) {
        real temp_real = 0.0;
        real temp_imag = 0.0;
        for (k = 0; k < K; k++) {
          const real Aki_real = (((const real *) A)[2*(k * lda + i)]);
          const real Aki_imag = -(((const real *) A)[2*(k * lda + i)+1]);
          const real Akj_real = (((const real *) A)[2*(k * lda + j)]);
          const real Akj_imag = (((const real *) A)[2*(k * lda + j)+1]);
          temp_real += Aki_real * Akj_real - Aki_imag * Akj_imag;
          temp_imag += Aki_real * Akj_imag + Aki_imag * Akj_real;
        }
        (((real *) C)[2*(i * ldc + j)]) += alpha * temp_real;
        (((real *) C)[2*(i * ldc + j)+1]) += alpha * temp_imag;
      }
    }
  } else if (uplo == Lower && trans == NoTrans) {
    for (i = 0; i < N; i++) {
      for (j = 0; j <= i; j++) {
        real temp_real = 0.0;
        real temp_imag = 0.0;
        for (k = 0; k < K; k++) {
          const real Aik_real = (((const real *) A)[2*(i * lda + k)]);
          const real Aik_imag = (((const real *) A)[2*(i * lda + k)+1]);
          const real Ajk_real = (((const real *) A)[2*(j * lda + k)]);
          const real Ajk_imag = -(((const real *) A)[2*(j * lda + k)+1]);
          temp_real += Aik_real * Ajk_real - Aik_imag * Ajk_imag;
          temp_imag += Aik_real * Ajk_imag + Aik_imag * Ajk_real;
        }
        (((real *) C)[2*(i * ldc + j)]) += alpha * temp_real;
        (((real *) C)[2*(i * ldc + j)+1]) += alpha * temp_imag;
      }
    }
  } else if (uplo == Lower && trans == ConjTrans) {
    for (i = 0; i < N; i++) {
      for (j = 0; j <= i; j++) {
        real temp_real = 0.0;
        real temp_imag = 0.0;
        for (k = 0; k < K; k++) {
          const real Aki_real = (((const real *) A)[2*(k * lda + i)]);
          const real Aki_imag = -(((const real *) A)[2*(k * lda + i)+1]);
          const real Akj_real = (((const real *) A)[2*(k * lda + j)]);
          const real Akj_imag = (((const real *) A)[2*(k * lda + j)+1]);
          temp_real += Aki_real * Akj_real - Aki_imag * Akj_imag;
          temp_imag += Aki_real * Akj_imag + Aki_imag * Akj_real;
        }
        (((real *) C)[2*(i * ldc + j)]) += alpha * temp_real;
        (((real *) C)[2*(i * ldc + j)+1]) += alpha * temp_imag;
      }
    }
  } else {
    xerbla(0, "source_herk.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_ZHERK_HPP__
