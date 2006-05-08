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

#ifndef __BLAS_FBLAS_HPP__
#define __BLAS_FBLAS_HPP__

/*! \file fortran_blas.h
 *  \brief Interface to Fortran BLAS operations.
 */

#include "blas.hpp"

extern "C" {
#include <cblas.h>
}

namespace BLAS {

#ifndef NO_CBLAS

/* Specialisations for T=double use cblas */

/* Level 1 BLAS functions */

template<>
inline 
double dot(const int N, const double *X, const int incX,
           const double *Y, const int incY)
{
  // std::cerr << "cblas_ddot\n";
  return cblas_ddot(N,X,incX,Y,incY);
}

template<>
inline
double nrm2(const int N, const double *X, const int incX)
{
  // std::cerr << "cblas_dnrm2\n";
  return cblas_dnrm2(N,X,incX);
}

template<>
inline
double asum(const int N, const double *X, const int incX)
{
  // std::cerr << "cblas_dasum\n";
  return cblas_dasum(N,X,incX);
}

template<>
inline
int iamax(const int N, const double *X, const int incX)
{
  // std::cerr << "cblas_idamax\n";
  return cblas_idamax(N,X,incX);
}


/* Level 1 BLAS routines */


template<>
inline
void
swap<double>(const int N, double *X, const int incX,
             double* Y, const int incY)
{
  // std::cerr << "cblas_dswap\n";
  cblas_dswap(N,X,incX,Y,incY);
}

template<>
inline
void copy<double,double>(const int N, const double *X,
                         const int incX, double *Y, const int incY)
{
  // std::cerr << "cblas_dcopy!\n";
  cblas_dcopy(N,X,incX,Y,incY);
}

template<>
inline
void axpy<double>(const int N, const double alpha, const double *X,
                  const int incX, double *Y, const int incY)
{
  // std::cerr << "cblas_daxpy!\n";
  cblas_daxpy(N,alpha,X,incX,Y,incY);
}

template<>
inline
void
scal<double>(const int N, const double alpha,
             double* X, const int incX)
{
  // std::cerr << "cblas_dscal\n";
  cblas_dscal(N,alpha,X,incX);
}


/* Level 2 */

template<>
inline
void
gemv<double>(const ORDER order,
             const TRANSPOSE transA, const int M, const int N,
             const double alpha, const double *A, const int ldA,
             const double *X, const int incX, const double beta,
             double *Y, const int incY)
{
  // std::cerr << "cblas_dgemv\n";
  cblas_dgemv(CBLAS_ORDER(order),CBLAS_TRANSPOSE(transA),
              M,N,alpha,A,ldA,X,incX,beta,Y,incY);
}

template<>
inline 
void trmv(const enum ORDER order, const enum UPLO uplo,
          const enum TRANSPOSE transA, const enum DIAG diag,
          const int N, const double *A, const int ldA, 
          double *X, const int incX)
{
  // std::cerr << "cblas_dtrmv\n";
  cblas_dtrmv(CBLAS_ORDER(order),CBLAS_UPLO(uplo),
              CBLAS_TRANSPOSE(transA), CBLAS_DIAG(diag),
              N,A,ldA,X,incX);
}

template<>
inline 
void trsv(const enum ORDER order, const enum UPLO uplo,
          const enum TRANSPOSE transA, const enum DIAG diag,
          const int N, const double *A, const int ldA, double *X,
          const int incX)
{
  // std::cerr << "cblas_dtrsv\n";
  cblas_dtrsv(CBLAS_ORDER(order),CBLAS_UPLO(uplo),
              CBLAS_TRANSPOSE(transA), CBLAS_DIAG(diag),
              N,A,ldA,X,incX);
}

template<>
inline 
void ger(const enum ORDER order, const int M, const int N,
         const double alpha, const double *X, const int incX,
         const double *Y, const int incY, double *A, const int ldA)
{
  // std::cerr << "cblas_dger\n";
  cblas_dger(CBLAS_ORDER(order),M,N,alpha,X,incX,Y,incY,A,ldA);
}



/* Level 3 */

template<>
inline
void
gemm<double>(const  ORDER order, const  TRANSPOSE transA,
                     const  TRANSPOSE transB, const int M, const int N,
                     const int K, const double alpha, const double *A,
                     const int ldA, const double *B, const int ldB,
                     const double beta, double *C, const int ldC)
{
  // std::cerr << "cblas_dgemm\n";
  cblas_dgemm(CBLAS_ORDER(order),
              CBLAS_TRANSPOSE(transA),CBLAS_TRANSPOSE(transB),
              M,N,K,alpha,A,ldA,B,ldB,beta,C,ldC);
}

template<>
inline 
void trmm(const enum ORDER order, const enum SIDE side,
          const enum UPLO uplo, const enum TRANSPOSE transA,
          const enum DIAG diag, const int M, const int N,
          const double alpha, const double *A, const int ldA,
          double *B, const int ldB)
{
  // std::cerr << "cblas_dtrmm\n";
  cblas_dtrmm(CBLAS_ORDER(order),CBLAS_SIDE(side),CBLAS_UPLO(uplo),
              CBLAS_TRANSPOSE(transA),CBLAS_DIAG(diag),
              M,N,alpha,A,ldA,B,ldB);
}


template<>
inline 
void trsm(const enum ORDER order, const enum SIDE side,
          const enum UPLO uplo, const enum TRANSPOSE transA,
          const enum DIAG diag, const int M, const int N,
          const double alpha, const double *A, const int ldA,
          double *B, const int ldB)
{
  // std::cerr << "cblas_dtrsm\n";
  cblas_dtrsm(CBLAS_ORDER(order),CBLAS_SIDE(side),CBLAS_UPLO(uplo),
              CBLAS_TRANSPOSE(transA),CBLAS_DIAG(diag),
              M,N,alpha,A,ldA,B,ldB);
}


/* Level 1 extensions */

template<>
inline
double amax(const int N, const double *X, const int incX)
{
  // std::cerr << "cblas_idamax\n";
  return X[iamax(N,X,incX)*incX];
}

template<>
inline
void set(const int N, const double alpha, double *X, const int incX)
{
  int i;
  int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
  for (i = 0; i < N; i++) {
    X[ix] = alpha;
    ix += incX;
  }
}

#endif // NO_CBLAS

}


#endif // __BLAS_FBLAS_HPP__
