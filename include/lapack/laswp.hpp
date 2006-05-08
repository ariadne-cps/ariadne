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

#ifndef __LAPACK_LASWP_HPP__
#define __LAPACK_LASWP_HPP__

#include "lapack.hpp"

template<typename Scalar>
void
LAPACK::laswp(BLAS::ORDER order, int n, Scalar *A, int ldA, 
              int k1, int k2, const int *iPiv, int incP)
{
  std::cerr << "LAPACK::laswp\n";
  
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    DLASWP performs a series of row interchanges on the matrix A.   
    One row interchange is initiated for each of rows K1 through K2 of A.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of columns of the matrix A.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the matrix of column dimension N to which the row   
            interchanges will be applied.   
            On exit, the permuted matrix.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.   

    K1      (input) INTEGER   
            The first element of IPIV for which a row interchange will   
            be done.   

    K2      (input) INTEGER   
            The last element of IPIV for which a row interchange will   
            be done.   

    IPIV    (input) INTEGER array, dimension (M*abs(INCX))   
            The vector of pivot indices.  Only the elements in positions   
            K1 through K2 of IPIV are accessed.   
            IPIV(K) = L implies rows K and L are to be interchanged.   

    INCP    (input) INTEGER   
            The increment between successive values of IPIV.  If IPIV   
            is negative, the pivots are applied in reverse order.   

    Further Details   
    ===============   

    Modified by   
     R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA   

   =====================================================================   


    Interchange row I with row IPIV(I) for each of rows K1 through K2. */

    assert(order==BLAS::RowMajor);
    assert(incP==1 || incP==-1);
    
    /* Local variables */
    static Scalar temp;
    static int ip, i1, i2, inci;

#define a_ref(a_1,a_2) A[(a_1)*ldA + a_2]

    /* Function Body */
    if (incP > 0) {
        i1 = k1;
        i2 = k2;
        inci = 1;
    } else if (incP < 0) {
        i1 = k2-1;
        i2 = k1-1;
        inci = -1;
      } else {
        return;
    }

    for (int i = i1; i != i2; i+=inci) {
        ip = iPiv[i];
        if (ip != i) {
            for (int k = 0; k < n; ++k) {
                temp = a_ref(i, k);
                a_ref(i, k) = a_ref(ip, k);
                a_ref(ip, k) = temp;
            }
        }
    }

    return;
}

#undef a_ref

#endif // __LAPACK_LASWP_HPP__
