#ifndef __LAPACK_LARFG_HPP__
#define __LAPACK_LARFG_HPP__

#include "lapack.hpp"

#include <blas/nrm2.hpp>
#include <blas/scal.hpp>

/* Subroutine */ 
template<typename real>
void 
LAPACK::larfg(BLAS::ORDER order, int n, real &alpha, real *X, 
                 int incX, real &tau)
{
  std::cerr  << "LAPACK::larfg\n";
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLARFG generates a real elementary reflector H of order n, such   
    that   

          H * ( alpha ) = ( beta ),   H' * H = I.   
              (   x   )   (   0  )   

    where alpha and beta are scalars, and x is an (n-1)-element real   
    vector. H is represented in the form   

          H = I - tau * ( 1 ) * ( 1 v' ) ,   
                        ( v )   

    where tau is a real scalar and v is a real (n-1)-element   
    vector.   

    If the elements of x are all zero, then tau = 0 and H is taken to be   
    the unit matrix.   

    Otherwise  1 <= tau <= 2.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the elementary reflector.   

    ALPHA   (input/output) DOUBLE PRECISION   
            On entry, the value alpha.   
            On exit, it is overwritten with the value beta.   

    X       (input/output) DOUBLE PRECISION array, dimension   
                           (1+(N-2)*abs(INCX))   
            On entry, the vector x.   
            On exit, it is overwritten with the vector v.   

    INCX    (input) INTEGER   
            The increment between elements of X. INCX > 0.   

    TAU     (output) DOUBLE PRECISION   
            The value tau.   

    =====================================================================   

*/
    assert(order==BLAS::RowMajor);
    
    /* Function Body */
    if (n <= 1) {
        tau = 0;
        return;
    }
    
    // std::cerr << "alpha=" << alpha << "  x=" << BLAS::vector(n-1,X,incX) << "\n";
    
    real xnorm = BLAS::nrm2(n-1, X, incX);
    
    if (xnorm == 0) {
        // H  =  I
        tau = 0;
    } else {
        /* general case */
        real norm = sqrt(alpha*alpha+xnorm*xnorm);
        real beta = (alpha >= 0) ? -norm : norm;

        tau = (beta - alpha) / beta;
        real sf = 1. / (alpha - beta);
        BLAS::scal(n-1, sf, X, incX);
        alpha = beta;
    }
    
    // std::cerr << "xnorm=" << xnorm << "  norm=" << norm << "  tau=" << tau << "  beta=" << beta << "  sf=" << sf << std::endl;
 
    return;

} /* dlarfg_ */

#endif // __LAPACK_LARFG_HPP__
