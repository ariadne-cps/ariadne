#ifndef __LAPACK_ORGQR_HPP__
#define __LAPACK_ORGQR_HPP__

#include "lapack.hpp"

#include "larf.hpp"
#include <blas/scal.hpp>

template<typename real>
void 
LAPACK::orgqr(BLAS::ORDER order, int m, int n, int k, real *A, int ldA, 
                 const real *tau, real *work)
{
  std::cerr << "LAPACK::orgqr\n";
  
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DORGQR generates an m by n real matrix Q with orthonormal columns,   
    which is defined as the first n columns of a product of k elementary   
    reflectors of order m   

          Q  =  H(1) H(2) . . . H(k)   

    as returned by DGEQRF.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix Q. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix Q. M >= N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines the   
            matrix Q. N >= K >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the i-th column must contain the vector which   
            defines the elementary reflector H(i), for i = 1,2,...,k, as   
            returned by DGEQRF in the first k columns of its array   
            argument A.   
            On exit, the m-by-n matrix Q.   

    LDA     (input) INTEGER   
            The first dimension of the array A. LDA >= max(1,M).   

    TAU     (input) DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGEQRF.   

    WORK    (workspace) DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument has an illegal value   

    =====================================================================   

*/
   
    #define a_ref(a_1,a_2) A[(a_1)*ldA + a_2]

    /* Function Body */
    int info = 0;
    if (m < 0) {
        info = -1;
    } else if (n < 0 || n > m) {
        info = -2;
    } else if (k < 0 || k > n) {
        info = -3;
    } else if (ldA < max(1,m)) {
        info = -5;
    }
    if (info != 0) {
        std::cerr << "DORGQR" << -info << std::endl;
        return;
    }

    // Quick return if possible */

    if (n <= 0) {
        return;
    }

    // Initialise columns k+1:n to columns of the unit matrix */

   for (int j = k; j != n; ++j) {
        for (int i = 0; i != m; ++i) {
            a_ref(i, j) = 0;
        }
        a_ref(j, j) = 1;
    }

    // std::cerr << matrix(3,3,A) << std::endl;
    for (int i = k-1; i != -1; --i) {
        // std::cerr << "i=" << i << "\n";

        // Apply H(i) to A(i:m,i+1:n) from the left */

        if (i < n) {
            a_ref(i, i) = 1;
            int i1 = m - i;
            int i2 = n - i -1;
            LAPACK::larf(order, BLAS::Left, i1, i2, &a_ref(i, i), ldA, tau[i], &a_ref(i, i+1), ldA, work);
        }
        if (i < m) {
            int i1 = m - i - 1;
            real d1 = -tau[i];
            BLAS::scal(i1, d1, &a_ref(i+1, i), ldA);
        }
        a_ref(i, i) = 1 - tau[i];

        // Set A(0:i,i) to zero

        for(int l = 0; l < i; ++l) {
            a_ref(l, i) = 0.;
        }
        // std::cerr << matrix(3,3,A) << std::endl;
    }
    return;

} 

#undef a_ref

#endif // __LAPACK_ORGQR_HPP__
