/*
 * Copyright (C) 2006 Pieter Collins <Pieter.Collins@cwi.nl>
 *
 * Based on the routine in LAPACK (version 3.0)
 *   Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
 *   Courant Institute, Argonne National Lab, and Rice University
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

#ifndef __LAPACK_GETRF_HPP__
#define __LAPACK_GETRF_HPP__

#include "lapack.hpp"

#include <blas/iamax.hpp>
#include <blas/scal.hpp>
#include <blas/swap.hpp>
#include <blas/ger.hpp>
  
template<typename Scalar>
void 
LAPACK::getrf(BLAS::ORDER order, 
              int m, int n, Scalar *A, int ldA, 
              int *piv)
{
    std::cerr << "LAPACK::getrf\n";
  
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1992   


    Purpose   
    =======   

    DGETRF computes an LU factorization of a general m-by-n matrix A   
    using partial pivoting with row interchanges.   

    The factorization has the form   
       A = P * L * U   
    where P is a permutation matrix, L is lower triangular with unit   
    diagonal elements (lower trapezoidal if m > n), and U is upper   
    triangular (upper trapezoidal if m < n).   

    This is the right-looking Level 2 BLAS version of the algorithm.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the m by n matrix to be factored.   
            On exit, the factors L and U from the factorization   
            A = P*L*U; the unit diagonal elements of L are not stored.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    IPIV    (output) INTEGER array, dimension (min(M,N))   
            The pivot indices; for 1 <= i <= min(M,N), row i of the   
            matrix was interchanged with row IPIV(i).   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -k, the k-th argument had an illegal value   
            > 0: if INFO = k, U(k,k) is exactly zero. The factorization   
                 has been completed, but the factor U is exactly   
                 singular, and division by zero will occur if it is used   
                 to solve a system of equations.   

    =====================================================================   

*/
    assert(order==BLAS::RowMajor);
   
    /* System generated locals */
    int  i__1, i__2, i__3;
    Scalar d__1;
    
    /* Local variables */
    static int j;
    static int jp;

    #define a_ref(a_1,a_2) A[(a_1)*ldA + a_2]


    /* Function Body */
    int info = 0;
    if (m < 0) {
        info = -1;
    } else if (n < 0) {
        info = -2;
   } else if (ldA < max(1,m)) {
        info = -4;
    }
    if (info != 0) {
        std::cerr << "DGETF2" << -info << std::endl; 
        return;
    }

/*     Quick return if possible */

    if (m == 0 || n == 0) {
        return;
    }

    i__1 = min(m,n);
    for (j = 0; j != i__1; ++j) {

      // Find pivot and test for singularity.
        i__2 = m - j;
        jp = j + BLAS::iamax(i__2, &a_ref(j, j), ldA);
        piv[j] = jp;
        
        if (a_ref(jp, j) != 0) {
            // Apply the interchange of rows j, jp to columns 0:N.
            if (jp != j) {
                BLAS::swap(n, &a_ref(j, 0), 1, &a_ref(jp, 0), 1);
            }
            // Compute elements J+1:M of J-th column. */
            if (j < m) {
                i__2 = m - j - 1;
                d__1 = 1. / a_ref(j, j);
                BLAS::scal(i__2, d__1, &a_ref(j+1,j), ldA);
            }
        } else if (info == 0) {
            // Singular matrix
            info = j;
        }

        if (j < min(m,n)) {
            // Update trailing submatrix.
            i__2 = m - j - 1;
            i__3 = n - j - 1;
            BLAS::ger(BLAS::RowMajor, i__2, i__3, Scalar(-1), 
                         &a_ref(j+1, j), ldA, &a_ref(j, j+1), 1, &a_ref(j+1, j+1), ldA);
        }
    }
    return;
}

#undef a_ref

#endif // __LAPACK_GETRF_HPP__
