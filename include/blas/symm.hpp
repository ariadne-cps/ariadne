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

#ifndef __BLAS_SYMM_HPP__
#define __BLAS_SYMM_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::symm (const enum ORDER Order, const enum SIDE Side,
             const enum UPLO Uplo, const int M, const int N,
             const Real alpha, const Real *A, const int lda,
             const Real *B, const int ldb, const Real beta, Real *C,
             const int ldc)
{
{
  int i, j, k;
  int n1, n2;
  int uplo, side;
  if (alpha == 0.0 && beta == 1.0)
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
  if (beta == 0.0) {
    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        C[ldc * i + j] = 0.0;
      }
    }
  } else if (beta != 1.0) {
    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        C[ldc * i + j] *= beta;
      }
    }
  }
  if (alpha == 0.0)
    return;
  if (side == Left && uplo == Upper) {
    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        const Real temp1 = alpha * B[ldb * i + j];
        Real temp2 = 0.0;
        C[i * ldc + j] += temp1 * A[i * lda + i];
        for (k = i + 1; k < n1; k++) {
          const Real Aik = A[i * lda + k];
          C[k * ldc + j] += Aik * temp1;
          temp2 += Aik * B[ldb * k + j];
        }
        C[i * ldc + j] += alpha * temp2;
      }
    }
  } else if (side == Left && uplo == Lower) {
    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        const Real temp1 = alpha * B[ldb * i + j];
        Real temp2 = 0.0;
        for (k = 0; k < i; k++) {
          const Real Aik = A[i * lda + k];
          C[k * ldc + j] += Aik * temp1;
          temp2 += Aik * B[ldb * k + j];
        }
        C[i * ldc + j] += temp1 * A[i * lda + i] + alpha * temp2;
      }
    }
  } else if (side == Right && uplo == Upper) {
    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        const Real temp1 = alpha * B[ldb * i + j];
        Real temp2 = 0.0;
        C[i * ldc + j] += temp1 * A[j * lda + j];
        for (k = j + 1; k < n2; k++) {
          const Real Ajk = A[j * lda + k];
          C[i * ldc + k] += temp1 * Ajk;
          temp2 += B[ldb * i + k] * Ajk;
        }
        C[i * ldc + j] += alpha * temp2;
      }
    }
  } else if (side == Right && uplo == Lower) {
    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        const Real temp1 = alpha * B[ldb * i + j];
        Real temp2 = 0.0;
        for (k = 0; k < j; k++) {
          const Real Ajk = A[j * lda + k];
          C[i * ldc + k] += temp1 * Ajk;
          temp2 += B[ldb * i + k] * Ajk;
        }
        C[i * ldc + j] += temp1 * A[j * lda + j] + alpha * temp2;
      }
    }
  } else {
    xerbla(0, "source_symm_r.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_SYMM_HPP__
