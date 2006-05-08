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

#ifndef __BLAS_SYRK_HPP__
#define __BLAS_SYRK_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::syrk (const enum ORDER Order, const enum UPLO Uplo,
             const enum TRANSPOSE Trans, const int N, const int K,
             const Real alpha, const Real *A, const int lda,
             const Real beta, Real *C, const int ldc)
{
{
  int i, j, k;
  int uplo, trans;
  if (alpha == 0.0 && beta == 1.0)
    return;
  if (Order == RowMajor) {
    uplo = Uplo;
    trans = (Trans == ConjTrans) ? Trans : Trans;
  } else {
    uplo = (Uplo == Upper) ? Lower : Upper;
    if (Trans == Trans || Trans == ConjTrans) {
      trans = NoTrans;
    } else {
      trans = Trans;
    }
  }
  if (beta == 0.0) {
    if (uplo == Upper) {
      for (i = 0; i < N; i++) {
        for (j = i; j < N; j++) {
          C[ldc * i + j] = 0.0;
        }
      }
    } else {
      for (i = 0; i < N; i++) {
        for (j = 0; j <= i; j++) {
          C[ldc * i + j] = 0.0;
        }
      }
    }
  } else if (beta != 1.0) {
    if (uplo == Upper) {
      for (i = 0; i < N; i++) {
        for (j = i; j < N; j++) {
          C[ldc * i + j] *= beta;
        }
      }
    } else {
      for (i = 0; i < N; i++) {
        for (j = 0; j <= i; j++) {
          C[ldc * i + j] *= beta;
        }
      }
    }
  }
  if (alpha == 0.0)
    return;
  if (uplo == Upper && trans == NoTrans) {
    for (i = 0; i < N; i++) {
      for (j = i; j < N; j++) {
        Real temp = 0.0;
        for (k = 0; k < K; k++) {
          temp += A[i * lda + k] * A[j * lda + k];
        }
        C[i * ldc + j] += alpha * temp;
      }
    }
  } else if (uplo == Upper && trans == Trans) {
    for (i = 0; i < N; i++) {
      for (j = i; j < N; j++) {
        Real temp = 0.0;
        for (k = 0; k < K; k++) {
          temp += A[k * lda + i] * A[k * lda + j];
        }
        C[i * ldc + j] += alpha * temp;
      }
    }
  } else if (uplo == Lower && trans == NoTrans) {
    for (i = 0; i < N; i++) {
      for (j = 0; j <= i; j++) {
        Real temp = 0.0;
        for (k = 0; k < K; k++) {
          temp += A[i * lda + k] * A[j * lda + k];
        }
        C[i * ldc + j] += alpha * temp;
      }
    }
  } else if (uplo == Lower && trans == Trans) {
    for (i = 0; i < N; i++) {
      for (j = 0; j <= i; j++) {
        Real temp = 0.0;
        for (k = 0; k < K; k++) {
          temp += A[k * lda + i] * A[k * lda + j];
        }
        C[i * ldc + j] += alpha * temp;
      }
    }
  } else {
    xerbla(0, "source_syrk_r.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_SYRK_HPP__
