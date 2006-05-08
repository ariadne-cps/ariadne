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

#ifndef __BLAS_ZHEMM_HPP__
#define __BLAS_ZHEMM_HPP__

#include "blas.hpp"

template<typename real>
void
BLAS::hemm (const enum ORDER Order, const enum SIDE Side,
             const enum UPLO Uplo, const int M, const int N,
             const complex<real> alpha, const complex<real> *A, const int lda, const complex<real> *B,
             const int ldb, const complex<real> beta, complex<real> *C, const int ldc)
{
{
  int i, j, k;
  int n1, n2;
  int uplo, side;
  const real alpha_real = alpha.real();
  const real alpha_imag = alpha.imag();
  const real beta_real = beta.real();
  const real beta_imag = beta.imag();
  if ((alpha_real == 0.0 && alpha_imag == 0.0)
      && (beta_real == 1.0 && beta_imag == 0.0))
    return;
  if (Order == RowMajor) {
    n1 = M;
    n2 = N;
    uplo = Uplo;
    side = Side;
  } else {
    n1 = N;
    n2 = M;
    uplo = (Uplo == Upper) ? Lower : Upper;
    side = (Side == Left) ? Right : Left;
  }
  if (beta_real == 0.0 && beta_imag == 0.0) {
    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        (((real *) C)[2*(ldc * i + j)]) = 0.0;
        (((real *) C)[2*(ldc * i + j)+1]) = 0.0;
      }
    }
  } else if (!(beta_real == 1.0 && beta_imag == 0.0)) {
    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        const real Cij_real = (((real *) C)[2*(ldc * i + j)]);
        const real Cij_imag = (((real *) C)[2*(ldc * i + j)+1]);
        (((real *) C)[2*(ldc * i + j)]) = beta_real * Cij_real - beta_imag * Cij_imag;
        (((real *) C)[2*(ldc * i + j)+1]) = beta_real * Cij_imag + beta_imag * Cij_real;
      }
    }
  }
  if (alpha_real == 0.0 && alpha_imag == 0.0)
    return;
  if (side == Left && uplo == Upper) {
    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        const real Bij_real = (((const real *) B)[2*(ldb * i + j)]);
        const real Bij_imag = (((const real *) B)[2*(ldb * i + j)+1]);
        const real temp1_real = alpha_real * Bij_real - alpha_imag * Bij_imag;
        const real temp1_imag = alpha_real * Bij_imag + alpha_imag * Bij_real;
        real temp2_real = 0.0;
        real temp2_imag = 0.0;
        {
          const real Aii_real = (((const real *) A)[2*(i * lda + i)]);
          (((real *) C)[2*(i * ldc + j)]) += temp1_real * Aii_real;
          (((real *) C)[2*(i * ldc + j)+1]) += temp1_imag * Aii_real;
        }
        for (k = i + 1; k < n1; k++) {
          const real Aik_real = (((const real *) A)[2*(i * lda + k)]);
          const real Aik_imag = (((const real *) A)[2*(i * lda + k)+1]);
          const real Bkj_real = (((const real *) B)[2*(ldb * k + j)]);
          const real Bkj_imag = (((const real *) B)[2*(ldb * k + j)+1]);
          (((real *) C)[2*(k * ldc + j)]) += Aik_real * temp1_real - (-Aik_imag) * temp1_imag;
          (((real *) C)[2*(k * ldc + j)+1]) += Aik_real * temp1_imag + (-Aik_imag) * temp1_real;
          temp2_real += Aik_real * Bkj_real - Aik_imag * Bkj_imag;
          temp2_imag += Aik_real * Bkj_imag + Aik_imag * Bkj_real;
        }
        (((real *) C)[2*(i * ldc + j)]) += alpha_real * temp2_real - alpha_imag * temp2_imag;
        (((real *) C)[2*(i * ldc + j)+1]) += alpha_real * temp2_imag + alpha_imag * temp2_real;
      }
    }
  } else if (side == Left && uplo == Lower) {
    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        const real Bij_real = (((const real *) B)[2*(ldb * i + j)]);
        const real Bij_imag = (((const real *) B)[2*(ldb * i + j)+1]);
        const real temp1_real = alpha_real * Bij_real - alpha_imag * Bij_imag;
        const real temp1_imag = alpha_real * Bij_imag + alpha_imag * Bij_real;
        real temp2_real = 0.0;
        real temp2_imag = 0.0;
        for (k = 0; k < i; k++) {
          const real Aik_real = (((const real *) A)[2*(i * lda + k)]);
          const real Aik_imag = (((const real *) A)[2*(i * lda + k)+1]);
          const real Bkj_real = (((const real *) B)[2*(ldb * k + j)]);
          const real Bkj_imag = (((const real *) B)[2*(ldb * k + j)+1]);
          (((real *) C)[2*(k * ldc + j)]) += Aik_real * temp1_real - (-Aik_imag) * temp1_imag;
          (((real *) C)[2*(k * ldc + j)+1]) += Aik_real * temp1_imag + (-Aik_imag) * temp1_real;
          temp2_real += Aik_real * Bkj_real - Aik_imag * Bkj_imag;
          temp2_imag += Aik_real * Bkj_imag + Aik_imag * Bkj_real;
        }
        {
          const real Aii_real = (((const real *) A)[2*(i * lda + i)]);
          (((real *) C)[2*(i * ldc + j)]) += temp1_real * Aii_real;
          (((real *) C)[2*(i * ldc + j)+1]) += temp1_imag * Aii_real;
        }
        (((real *) C)[2*(i * ldc + j)]) += alpha_real * temp2_real - alpha_imag * temp2_imag;
        (((real *) C)[2*(i * ldc + j)+1]) += alpha_real * temp2_imag + alpha_imag * temp2_real;
      }
    }
  } else if (side == Right && uplo == Upper) {
    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        const real Bij_real = (((const real *) B)[2*(ldb * i + j)]);
        const real Bij_imag = (((const real *) B)[2*(ldb * i + j)+1]);
        const real temp1_real = alpha_real * Bij_real - alpha_imag * Bij_imag;
        const real temp1_imag = alpha_real * Bij_imag + alpha_imag * Bij_real;
        real temp2_real = 0.0;
        real temp2_imag = 0.0;
        {
          const real Ajj_real = (((const real *) A)[2*(j * lda + j)]);
          (((real *) C)[2*(i * ldc + j)]) += temp1_real * Ajj_real;
          (((real *) C)[2*(i * ldc + j)+1]) += temp1_imag * Ajj_real;
        }
        for (k = j + 1; k < n2; k++) {
          const real Ajk_real = (((const real *) A)[2*(j * lda + k)]);
          const real Ajk_imag = (((const real *) A)[2*(j * lda + k)+1]);
          const real Bik_real = (((const real *) B)[2*(ldb * i + k)]);
          const real Bik_imag = (((const real *) B)[2*(ldb * i + k)+1]);
          (((real *) C)[2*(i * ldc + k)]) += temp1_real * Ajk_real - temp1_imag * Ajk_imag;
          (((real *) C)[2*(i * ldc + k)+1]) += temp1_real * Ajk_imag + temp1_imag * Ajk_real;
          temp2_real += Bik_real * Ajk_real - Bik_imag * (-Ajk_imag);
          temp2_imag += Bik_real * (-Ajk_imag) + Bik_imag * Ajk_real;
        }
        (((real *) C)[2*(i * ldc + j)]) += alpha_real * temp2_real - alpha_imag * temp2_imag;
        (((real *) C)[2*(i * ldc + j)+1]) += alpha_real * temp2_imag + alpha_imag * temp2_real;
      }
    }
  } else if (side == Right && uplo == Lower) {
    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        const real Bij_real = (((const real *) B)[2*(ldb * i + j)]);
        const real Bij_imag = (((const real *) B)[2*(ldb * i + j)+1]);
        const real temp1_real = alpha_real * Bij_real - alpha_imag * Bij_imag;
        const real temp1_imag = alpha_real * Bij_imag + alpha_imag * Bij_real;
        real temp2_real = 0.0;
        real temp2_imag = 0.0;
        for (k = 0; k < j; k++) {
          const real Ajk_real = (((const real *) A)[2*(j * lda + k)]);
          const real Ajk_imag = (((const real *) A)[2*(j * lda + k)+1]);
          const real Bik_real = (((const real *) B)[2*(ldb * i + k)]);
          const real Bik_imag = (((const real *) B)[2*(ldb * i + k)+1]);
          (((real *) C)[2*(i * ldc + k)]) += temp1_real * Ajk_real - temp1_imag * Ajk_imag;
          (((real *) C)[2*(i * ldc + k)+1]) += temp1_real * Ajk_imag + temp1_imag * Ajk_real;
          temp2_real += Bik_real * Ajk_real - Bik_imag * (-Ajk_imag);
          temp2_imag += Bik_real * (-Ajk_imag) + Bik_imag * Ajk_real;
        }
        {
          const real Ajj_real = (((const real *) A)[2*(j * lda + j)]);
          (((real *) C)[2*(i * ldc + j)]) += temp1_real * Ajj_real;
          (((real *) C)[2*(i * ldc + j)+1]) += temp1_imag * Ajj_real;
        }
        (((real *) C)[2*(i * ldc + j)]) += alpha_real * temp2_real - alpha_imag * temp2_imag;
        (((real *) C)[2*(i * ldc + j)+1]) += alpha_real * temp2_imag + alpha_imag * temp2_real;
      }
    }
  } else {
    xerbla(0, "source_hemm.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_ZHEMM_HPP__
