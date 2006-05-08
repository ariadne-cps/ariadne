#ifndef __LAPACK_GEQRF_HPP__
#define __LAPACK_GEQRF_HPP__

#include "lapack.hpp"

#include "larf.hpp"
#include "larfg.hpp"


template<typename real>
void 
LAPACK::geqrf(BLAS::ORDER order, int m, int n, real *A, int ldA, 
              real *tau, real *work)
{
  std::cerr << "LAPACK::geqrf\n";
  
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DGEQRF computes a QR factorization of a real m by n matrix A:   
    A = Q * R.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the m by n matrix A.   
            On exit, the elements on and above the diagonal of the array   
            contain the min(m,n) by n upper trapezoidal matrix R (R is   
            upper triangular if m >= n); the elements below the diagonal,   
            with the array TAU, represent the orthogonal matrix Q as a   
            product of elementary reflectors (see Further Details).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))   
            The scalar factors of the elementary reflectors (see Further   
            Details).   

    WORK    (workspace) DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    The matrix Q is represented as a product of elementary reflectors   

       Q = H(1) H(2) . . . H(k), where k = min(m,n).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),   
    and tau in TAU(i).   

    =====================================================================   
*/

   assert(order==BLAS::RowMajor);
    
    /* Local variables */
#define a_ref(a_1,a_2) A[a_1*ldA+a_2]


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
        std::cerr << "DGEQR2" << -info << "\n";
        return;
    }

    int k = min(m,n);

    for (int i = 0; i < k; ++i) {
        // Generate elementary reflector H(i) to annihilate A(i+1:m,i)

        // Computing MIN
        int i2 = i + 1;
        int i3 = m - i;
        larfg(BLAS::RowMajor, m-i, a_ref(i, i), &a_ref(min(i2,m), i), ldA, tau[i]);
        if (i < n) {
            // Apply H(i) to A(i+1:m,i+1:n) from the left
            real aii = a_ref(i, i);
            a_ref(i, i) = 1;
            i2 = m - i;
            i3 = n - i - 1;
            larf(BLAS::RowMajor, BLAS::Left, i2, i3, &a_ref(i, i), ldA, tau[i],
                 &a_ref(i, i + 1), ldA, &work[1]);
            a_ref(i, i) = aii;
        }
    }
    return;

}

#undef a_ref

#endif // __LAPACK_GEQRF_HPP__
