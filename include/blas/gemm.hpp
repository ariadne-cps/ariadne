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

#ifndef __BLAS_GEMM_HPP__
#define __BLAS_GEMM_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::gemm (const enum ORDER Order, const enum TRANSPOSE TransA,
             const enum TRANSPOSE TransB, const int M, const int N,
             const int K, const Real alpha, const Real *A, const int lda,
             const Real *B, const int ldb, const Real beta, Real *C,
             const int ldc)
{
{
  int i, j, k;
  int n1, n2;
  int ldf, ldg;
  int TransF, TransG;
  const Real *F, *G;
  if (alpha == 0.0 && beta == 1.0)
    return;
  if (Order == RowMajor) {
    n1 = M;
    n2 = N;
    F = A;
    ldf = lda;
    TransF = (TransA == ConjTrans) ? Trans : TransA;
    G = B;
    ldg = ldb;
    TransG = (TransB == ConjTrans) ? Trans : TransB;
  } else {
    n1 = N;
    n2 = M;
    F = B;
    ldf = ldb;
    TransF = (TransB == ConjTrans) ? Trans : TransB;
    G = A;
    ldg = lda;
    TransG = (TransA == ConjTrans) ? Trans : TransA;
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
  if (TransF == NoTrans && TransG == NoTrans) {
    for (k = 0; k < K; k++) {
      for (i = 0; i < n1; i++) {
        const Real temp = alpha * F[ldf * i + k];
        if (temp != 0.0) {
          for (j = 0; j < n2; j++) {
            C[ldc * i + j] += temp * G[ldg * k + j];
          }
        }
      }
    }
  } else if (TransF == NoTrans && TransG == Trans) {
    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        Real temp = 0.0;
        for (k = 0; k < K; k++) {
          temp += F[ldf * i + k] * G[ldg * j + k];
        }
        C[ldc * i + j] += alpha * temp;
      }
    }
  } else if (TransF == Trans && TransG == NoTrans) {
    for (k = 0; k < K; k++) {
      for (i = 0; i < n1; i++) {
        const Real temp = alpha * F[ldf * k + i];
        if (temp != 0.0) {
          for (j = 0; j < n2; j++) {
            C[ldc * i + j] += temp * G[ldg * k + j];
          }
        }
      }
    }
  } else if (TransF == Trans && TransG == Trans) {
    for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        Real temp = 0.0;
        for (k = 0; k < K; k++) {
          temp += F[ldf * k + i] * G[ldg * j + k];
        }
        C[ldc * i + j] += alpha * temp;
      }
    }
  } else {
    xerbla(0, "source_gemm_r.h", "unrecognized operation");;
  }
}
}

#endif // __BLAS_GEMM_HPP__
