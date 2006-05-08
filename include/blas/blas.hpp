/*  blas.hpp
 *
 *  Copyright (C) 2004 Pieter Collins <Pieter.Collins@cwi.nl>
 *
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

#ifndef __BLAS_HPP__
#define __BLAS_HPP__

/*! \file blas.hpp
 *  \brief Templated BLAS operations in C++.
 *
 *  Matrix types 
 *    -- GE General
 *    -- SY Symmetric
 *    -- HE Hermitian
 *    -- TR Triangular
 *    -- GB General Band
 *    -- SB Symmetric Band
 *    -- HB Hermitian Band
 *    -- TB Triangular Band
 *    -- SP Symmetric Packed
 *    -- HP Hermitian Packed
 *    -- TB Triangular Packed
 */

#include <cmath>
#include <complex>
#include <iostream>

namespace BLAS {

using std::complex;

enum ORDER {RowMajor=101, ColMajor=102};
enum TRANSPOSE {NoTrans=111, Trans=112, ConjTrans=113};
enum UPLO {Upper=121, Lower=122};
enum DIAG {NonUnit=131, Unit=132};
enum SIDE {Left=141, Right=142};

void xerbla(int p, const char *rout, const char *form, ...);

inline int offset(const int& N, const int& incX) { return incX > 0 ?  0 : (N - 1) * (-incX); }

template<typename real> inline real max(const real& a, const real& b) { return a > b ? a : b; }
template<typename real> inline real min(const real& a, const real& b) { return a < b ? a : b; }
template<typename real> inline real sign(const real& x) { return x >= 0 ? 1 : -1; }
template<typename real> inline real abs(const real& x) { return x >= 0 ? x : -x; }

template<typename real> 
real amax(const int N, const real *X, const int incX);

template<typename real> 
real amax(const int N, const complex<real> *X, const int incX);

template<typename scalar, typename scalarX> 
void set(const int N, const scalar alpha, 
         scalarX *X, const int incX);

/*
 * ===========================================================================
 * Prototypes for level 1 BLAS functions
 * ===========================================================================
 */

/*
 * Functions for real and complex scalars
 */
template<typename real> 
real nrm2(const int N, const real *X, const int incX);

template<typename real> 
real nrm2(const int N, const complex<real> *X, const int incX);

template<typename real> 
real asum(const int N, const real *X, const int incX);

template<typename real> 
real asum(const int N, const complex<real> *X, const int incX);

template<typename scalar> 
int iamax(const int N, const scalar *X, const int incX);

/*
 * Functions for real scalars only
 */
template<typename scalar> 
scalar dot(const int N, const scalar *X, const int incX,
           const scalar *Y, const int incY);

/* 
 * Functions for complex scalars only
 */
template<typename real> 
complex<real> dotu(const int N, const complex<real> *X, const int incX,
                   const complex<real> *Y, const int incY);

template<typename real> 
complex<real> dotc(const int N, const complex<real> *X, const int incX,
                   const complex<real> *Y, const int incY);


/*
 * ===========================================================================
 * Prototypes for level 1 BLAS routines
 * ===========================================================================
 */

/* 
 * Routines for real and complex scalars
 */
template<typename scalar> 
void swap(const int N, scalar *X, const int incX, 
          scalar *Y, const int incY);

template<typename scalarX, typename scalarY> 
void copy(const int N, const scalarX *X, const int incX, 
          scalarY *Y, const int incY);

template<typename scalar> 
void axpy(const int N, const scalar alpha, const scalar *X,
          const int incX, scalar *Y, const int incY);

template<typename scalar> 
void scal(const int N, const scalar alpha, scalar *X, const int incX);
  

/* 
 * Routines for real scalars only
 */
template<typename real> 
void rotg(real *a, real *b, real *c, real *s);

template<typename real> 
void rotmg(real *d1, real *d2, real *b1, const real b2, real *P);
  
template<typename real> 
void rot(const int N, real *X, const int incX,
         real *Y, const int incY, const real c, const real  s);
template<typename real> 
void rotm(const int N, real *X, const int incX,
          real *Y, const int incY, const real *P);



/*
 * ===========================================================================
 * Prototypes for level 2 BLAS
 * ===========================================================================
 */

/* 
 * Non-standard extensions
 */
template<typename scalarA, typename scalarB>
void gecpy(const enum ORDER order, int m, int n, 
           const scalarA *A, int ldA, scalarB *B, int ldB);
  
template<typename scalarA, typename scalarB>
void trcpy(const enum ORDER order, const enum UPLO uplo, const enum DIAG diag, 
           int m, int n, const scalarA *A, int ldA, scalarB *B, int ldB);

template<typename scalar, typename scalarA>
void geset(const enum ORDER order, int m, int n, 
           scalar alpha, scalar beta, scalarA *A, int ldA);
  
template<typename scalar, typename scalarA>
void trset(const enum ORDER order, const enum UPLO uplo, const enum DIAG diag, 
           int m, int n, scalar alpha, scalarA *A, int ldA);


/* 
 * Routines for real and complex scalars
 */
template<typename scalar> 
void gemv(const enum ORDER order,
          const enum TRANSPOSE TransA, const int M, const int N,
          const scalar alpha, const scalar *A, const int lda,
          const scalar *X, const int incX, const scalar beta,
          scalar *Y, const int incY);

template<typename scalar> 
void gbmv(const enum ORDER order,
          const enum TRANSPOSE TransA, const int M, const int N,
          const int KL, const int KU, const scalar alpha,
          const scalar *A, const int lda, const scalar *X,
          const int incX, const scalar beta, scalar *Y, const int incY);

template<typename scalar> 
void trmv(const enum ORDER order, const enum UPLO Uplo,
          const enum TRANSPOSE TransA, const enum DIAG Diag,
          const int N, const scalar *A, const int lda, 
          scalar *X, const int incX);

template<typename scalar> 
void tbmv(const enum ORDER order, const enum UPLO Uplo,
          const enum TRANSPOSE TransA, const enum DIAG Diag,
          const int N, const int K, const scalar *A, const int lda, 
          scalar *X, const int incX);

template<typename scalar> 
void tpmv(const enum ORDER order, const enum UPLO Uplo,
          const enum TRANSPOSE TransA, const enum DIAG Diag,
          const int N, const scalar *Ap, scalar *X, const int incX);

template<typename scalar> 
void trsv(const enum ORDER order, const enum UPLO Uplo,
          const enum TRANSPOSE TransA, const enum DIAG Diag,
          const int N, const scalar *A, const int lda, scalar *X,
          const int incX);

template<typename scalar> 
void tbsv(const enum ORDER order, const enum UPLO Uplo,
          const enum TRANSPOSE TransA, const enum DIAG Diag,
          const int N, const int K, const scalar *A, const int lda,
          scalar *X, const int incX);

template<typename scalar> 
void tpsv(const enum ORDER order, const enum UPLO Uplo,
          const enum TRANSPOSE TransA, const enum DIAG Diag,
          const int N, const scalar *Ap, scalar *X, const int incX);


/* 
 * Routines for real scalars only
 */
template<typename real> 
void symv(const enum ORDER order, const enum UPLO Uplo,
          const int N, const real alpha, const real *A,
          const int lda, const real *X, const int incX,
          const real beta, real *Y, const int incY);

template<typename real> 
void sbmv(const enum ORDER order, const enum UPLO Uplo,
          const int N, const int K, const real alpha, const real *A,
          const int lda, const real *X, const int incX,
          const real beta, real *Y, const int incY);

template<typename real> 
void spmv(const enum ORDER order, const enum UPLO Uplo,
          const int N, const real alpha, const real *Ap,
          const real *X, const int incX,
          const real beta, real *Y, const int incY);

template<typename real> 
void ger(const enum ORDER order, const int M, const int N,
         const real alpha, const real *X, const int incX,
         const real *Y, const int incY, real *A, const int lda);

template<typename real> 
void syr(const enum ORDER order, const enum UPLO Uplo,
         const int N, const real alpha, const real *X,
         const int incX, real *A, const int lda);
  
template<typename real> 
void spr(const enum ORDER order, const enum UPLO Uplo,
         const int N, const real alpha, const real *X,
         const int incX, real *Ap);

template<typename real> 
void syr2(const enum ORDER order, const enum UPLO Uplo,
          const int N, const real alpha, const real *X,
          const int incX, const real *Y, const int incY, real *A,
          const int lda);

template<typename real> 
void spr2(const enum ORDER order, const enum UPLO Uplo,
          const int N, const real alpha, const real *X,
          const int incX, const real *Y, const int incY, real *A);


/* 
 * Routines for complex scalars only
 */
template<typename real> 
void hemv(const enum ORDER order, const enum UPLO Uplo,
          const int N, const complex<real> alpha, const complex<real> *A,
          const int lda, const complex<real> *X, const int incX,
          const complex<real> beta, complex<real> *Y, const int incY);

template<typename real> 
void hbmv(const enum ORDER order, const enum UPLO Uplo,
          const int N, const int K, const complex<real> alpha, 
          const complex<real> *A, const int lda, 
          const complex<real> *X, const int incX,
          const complex<real> beta, complex<real> *Y, const int incY);

template<typename real> 
void hpmv(const enum ORDER order, const enum UPLO Uplo,
          const int N, const complex<real> alpha, const complex<real> *Ap,
          const complex<real> *X, const int incX,
          const complex<real> beta, complex<real> *Y, const int incY);

template<typename real> 
void geru(const enum ORDER order, const int M, const int N,
          const complex<real> alpha, const complex<real> *X, const int incX,
          const complex<real> *Y, const int incY, 
          complex<real> *A, const int lda);

template<typename real> 
void gerc(const enum ORDER order, const int M, const int N,
          const complex<real> alpha, const complex<real> *X, const int incX,
          const complex<real> *Y, const int incY, 
          complex<real> *A, const int lda);

template<typename real> 
void her(const enum ORDER order, const enum UPLO Uplo,
         const int N, const real alpha, const complex<real> *X, const int incX,
         complex<real> *A, const int lda);

template<typename real> 
void hpr(const enum ORDER order, const enum UPLO Uplo,
         const int N, const real alpha, const complex<real> *X,
         const int incX, complex<real> *A);
  
template<typename real> 
void her2(const enum ORDER order, const enum UPLO Uplo, const int N,
          const complex<real> alpha, const complex<real> *X, const int incX,
          const complex<real> *Y, const int incY, 
          complex<real> *A, const int lda);

template<typename real> 
void hpr2(const enum ORDER order, const enum UPLO Uplo, const int N,
          const complex<real> alpha, const complex<real> *X, const int incX,
          const complex<real> *Y, const int incY, complex<real> *Ap);

/*
 * ===========================================================================
 * Prototypes for level 3 BLAS
 * ===========================================================================
 */

/* 
 * Routines for real and complex scalars
 */
template<typename scalar> 
void gemm(const enum ORDER Order, const enum TRANSPOSE TransA,
          const enum TRANSPOSE TransB, const int M, const int N,
          const int K, const scalar alpha, const scalar *A,
          const int lda, const scalar *B, const int ldb,
          const scalar beta, scalar *C, const int ldc);

template<typename scalar> 
void symm(const enum ORDER Order, const enum SIDE Side,
          const enum UPLO Uplo, const int M, const int N,
          const scalar alpha, const scalar *A, const int lda,
          const scalar *B, const int ldb, const scalar beta,
          scalar *C, const int ldc);

template<typename scalar> 
void syrk(const enum ORDER Order, const enum UPLO Uplo,
          const enum TRANSPOSE Trans, const int N, const int K,
          const scalar alpha, const scalar *A, const int lda,
          const scalar beta, scalar *C, const int ldc);

template<typename scalar> 
void syr2k(const enum ORDER Order, const enum UPLO Uplo,
           const enum TRANSPOSE Trans, const int N, const int K,
           const scalar alpha, const scalar *A, const int lda,
           const scalar *B, const int ldb, const scalar beta,
           scalar *C, const int ldc);

template<typename scalar> 
void trmm(const enum ORDER Order, const enum SIDE Side,
          const enum UPLO Uplo, const enum TRANSPOSE TransA,
          const enum DIAG Diag, const int M, const int N,
          const scalar alpha, const scalar *A, const int lda,
          scalar *B, const int ldb);

template<typename scalar> 
void trsm(const enum ORDER Order, const enum SIDE Side,
          const enum UPLO Uplo, const enum TRANSPOSE TransA,
          const enum DIAG Diag, const int M, const int N,
          const scalar alpha, const scalar *A, const int lda,
          scalar *B, const int ldb);



/* 
 * Routines for complex scalars only
 */
template<typename real> 
void hemm(const enum ORDER Order, const enum SIDE Side,
          const enum UPLO Uplo, const int M, const int N,
          const complex<real> alpha, const complex<real> *A, const int lda,
          const complex<real> *B, const int ldb, const complex<real> beta,
          complex<real> *C, const int ldc);

template<typename real> 
void herk(const enum ORDER Order, const enum UPLO Uplo,
          const enum TRANSPOSE Trans, const int N, const int K,
          const real alpha, const complex<real> *A, const int lda,
          const real beta, complex<real> *C, const int ldc);

template<typename real> 
void her2k(const enum ORDER Order, const enum UPLO Uplo,
           const enum TRANSPOSE Trans, const int N, const int K,
           const complex<real> alpha, const complex<real> *A, const int lda,
           const complex<real> *B, const int ldb, const real beta,
           complex<real> *C, const int ldc);

}

#include "fblas.hpp"
#include "xerbla.hpp"

#endif // __BLAS_HPP__
