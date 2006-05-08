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

#ifndef __LAPACK_LARF_HPP__
#define __LAPACK_LARF_HPP__

#include "lapack.hpp"

#include <blas/ger.hpp>
#include <blas/gemm.hpp>

/* Subroutine */ 
template<typename real>
void 
LAPACK::larf(BLAS::ORDER order, BLAS::SIDE side, 
             int m, int n, 
             real *V, int incV, 
             real tau, 
             real *C, int ldC, 
             real *work)
{
  std::cerr << "LAPACK::larf\n";
  
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLARF applies a real elementary reflector H to a real m by n matrix   
    C, from either the left or the right. H is represented in the form   

          H = I - tau * v * v'   

    where tau is a real scalar and v is a real vector.   

    If tau = 0, then H is taken to be the unit matrix.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': form  H * C   
            = 'R': form  C * H   

    M       (input) INTEGER   
            The number of rows of the matrix C.   

    N       (input) INTEGER   
            The number of columns of the matrix C.   

    V       (input) DOUBLE PRECISION array, dimension   
                       (1 + (M-1)*abs(INCV)) if SIDE = 'L'   
                    or (1 + (N-1)*abs(INCV)) if SIDE = 'R'   
            The vector v in the representation of H. V is not used if   
            TAU = 0.   

    INCV    (input) INTEGER   
            The increment between elements of v. INCV <> 0.   

    TAU     (input) DOUBLE PRECISION   
            The value tau in the representation of H.   

    C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the m by n matrix C.   
            On exit, C is overwritten by the matrix H * C if SIDE = 'L',   
            or C * H if SIDE = 'R'.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >= max(1,M).   

    WORK    (workspace) DOUBLE PRECISION array, dimension   
                           (N) if SIDE = 'L'   
                        or (M) if SIDE = 'R'   

    =====================================================================   


       Parameter adjustments */
    /* Table of constant values */
    const real one = 1;
    const real zero = 0;
    
    /* Function Body */
    if (side==BLAS::Left) {
        // Form  H * C 
        if (tau != 0) {
            // w := C' * v
            BLAS::gemv(BLAS::RowMajor, BLAS::Trans, m, n, one, C, ldC, V, incV,
                          zero, work, 1);
            // C := C - v * w'
            BLAS::ger(BLAS::RowMajor, m, n, -tau, V, incV, work, 1, C, ldC);
        }
    } else {
        // Form  C * H
        if (tau != 0.) {
            // w := C * v
            BLAS::gemv(BLAS::RowMajor, BLAS::NoTrans, m, n, one, C, ldC, V, 
                          incV, zero, work, 1);
            // C := C - w * v'
            BLAS::ger(BLAS::RowMajor, m, n, -tau, work, 1, V, incV, C, ldC);
        }
    }
    return;

}

#endif // __LAPACK_LARF_HPP__
