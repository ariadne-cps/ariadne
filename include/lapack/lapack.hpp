/*  lapack.hpp
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

#ifndef __LAPACK_HPP__
#define __LAPACK_HPP__

/*! \file lapack.hpp
 *  \brief Templated LAPACK operations in C++.
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

#include <blas/blas.hpp>

/*
 * Recall cblas enum's.
 *
 * enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
 * enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
 * enum CBLAS_UPLO {CblasUpper=121, CblasLower=122};
 * enum CBLAS_DIAG {CblasNonUnit=131, CblasUnit=132};
 * enum CBLAS_SIDE {CblasLeft=141, CblasRight=142};
 */

#include <iostream>
#include <algorithm>

/*! \brief Templated LAPACK routines. */
namespace LAPACK {

enum DIRECTION { Forward, Backward };
enum PIVOT { Variable, Top, Bottom };

using std::max;
using std::min;
using std::sqrt;

template<typename S>
inline
S 
abs(const S& s) 
{ 
  return (s<0) ? -s :s;
}

template<typename S>
inline
S
sign(const S& a, const S& b) 
{ 
  S x = (a >= 0 ? a : - a);
  return( b >= 0 ? x : -x);
}

template<typename S>
inline
S
pow(const S& x, int n) 
{
  S y=1; while(n>0) { y*=x; --n; } while(n<0) { y/=x; ++n; } return y;
}


inline
bool
lsame(char* s1, const char* s2) 
{
  return s1[0]==s2[0];
}

/*
 * ===========================================================================
 * Prototypes for linear equations
 * ===========================================================================
 */
 
template<typename Scalar>
void 
gesv(BLAS::ORDER order, int n, int nrhs, Scalar *A, int ldA, int *P, Scalar *B, int ldB);

/*! \brief Compute LU factorization of \c m x \c n matrix \c A with pivots \c piv. 
 *
 *  \c A is overwritten so that \f$ A_{ij}\f$ with \f$i>j\f$ contains the 
 *  nonzero, non-diagonal elements of \f$L\f$, and \f$A_{ij}\f$ with \f$i<=j\f$ contains
 *  the nonzero elements of \f$ U \f$.
 */
template<typename Scalar>
void 
getrf(BLAS::ORDER order, int m, int n, Scalar *A, int ldA, int *piv);


template<typename Scalar>
void 
getrs(BLAS::ORDER order, BLAS::TRANSPOSE trans, 
      int n, int nrhs, Scalar *a, int ldA, int *piv, Scalar *B, int ldB);

/*
 * ===========================================================================
 * Prototypes for QR factorization
 * ===========================================================================
 */
 
/*! \brief Compute QR factorization of general matrix A. */
template<typename Real>
void 
geqrf(BLAS::ORDER order, int m, int n, Real *A, int ldA, 
      Real *tau, Real *work);

/*! \brief Generate orthogonal matrix Q. */
template<typename Real>
void 
orgqr(BLAS::ORDER order, int m, int n, int k, Real *A, int ldA, 
      const Real *tau, Real *work);

/*
 * ===========================================================================
 * Prototypes for singular value decomposition
 * ===========================================================================
 */
 
/*! \brief Simple (non-standard) driver routine. */
template<typename real>
void
gesvd(BLAS::ORDER order, int m, int n, real *A, int ldA,
      real *S, real *U, int ldU, real *V, int ldV);

/*
template<typename real>
void
gesvd(BLAS::ORDER order, int m, int n, 
      const real *A, int ldA, real *S, real *U, int ldU, real *Vt, int ldVt, 
      real *work, int *lwork);
 
template<typename real>
void
gebrd(BLAS::ORDER order, int m, int n, real *A, int ldA, real *D, real *E, 
      real *tauq, real *taup, real *work);
*/

/*
 * ===========================================================================
 * Prototypes for BLAS extensions
 * ===========================================================================
 */
 
/*! \brief Sets diagonal elements to \a beta and off-diagonal elements to \a alpha. */
template<typename scalar>
void
laset(BLAS::ORDER order, int m, int n, scalar alpha, scalar beta, scalar *A, int ldA);

/*! \brief Swap rows \c k1 and \c k2 of matrix \c A with \c n columns and leading dimension \c ldA. */
template<typename Scalar>
void 
laswp(BLAS::ORDER order, int n, Scalar *A, int ldA, 
      int k1, int k2, const int *piv, int incPiv);
  
 
/*! \brief Apply elementary Householder transfomation. */
template<typename Real>
void 
larf(BLAS::ORDER order, BLAS::SIDE, 
     int m, int n, 
     Real *V, int incV, 
     Real tau, 
     Real *C, int ldC, 
     Real *work);

/*! \brief Generate elementary Householder transfomation. */
template<typename Real>
void 
larfg(BLAS::ORDER order, int n, Real &alpha, Real *X, 
      int incX, Real &tau);

/*! \brief Transforms A by a product of Givens rotations. */
template<typename real>
void
lasr(BLAS::ORDER order, BLAS::SIDE side, LAPACK::PIVOT pivot, LAPACK::DIRECTION direct,
     int m, int n, real *C, real *S, real *A, int ldA);

/*
 * ===========================================================================
 * Prototypes for environment detection routines
 * ===========================================================================
 */


}

#endif // __LAPACK_HPP__
