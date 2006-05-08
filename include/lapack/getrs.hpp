#ifndef __LAPACK_GETRS_HPP__
#define __LAPACK_GETRS_HPP__

#include <blas/blas.hpp>

template<typename Scalar>
void 
LAPACK::getrs(BLAS::ORDER order, BLAS::TRANSPOSE trans, 
              int n, int nrhs, Scalar *A, int ldA, int *piv, Scalar *B, int ldB)
{
  std::cerr << "LAPACK::getrs\n";
  assert(order==BLAS::RowMajor);
  
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DGETRS solves a system of linear equations   
       A * X = B  or  A' * X = B   
    with a general N-by-N matrix A using the LU factorization computed   
    by DGETRF.   

    Arguments   
    =========   

    TRANS   (input) CHARACTER*1   
            Specifies the form of the system of equations:   
            = 'N':  A * X = B  (No transpose)   
            = 'T':  A'* X = B  (Transpose)   
            = 'C':  A'* X = B  (Conjugate transpose = Transpose)   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,N)   
            The factors L and U from the factorization A = P*L*U   
            as computed by DGETRF.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    IPIV    (input) INTEGER array, dimension (N)   
            The pivot indices from DGETRF; for 1<=i<=N, row i of the   
            matrix was interchanged with row IPIV(i).   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the right hand side matrix B.   
            On exit, the solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    /* Table of constant values */
    const Scalar one = 1.;
    static int c__1 = 1;
    static int c_n1 = -1;
    
    /* System generated locals */
    int a_dim1, a_offset, b_dim1, b_offset, i__1;
    /* Local variables */
    static bool notran;


    /* Function Body */
    int info = 0;
    if (trans!=BLAS::Trans && trans!=BLAS::NoTrans 
                               && trans!=BLAS::ConjTrans) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (nrhs < 0) {
        info = -3;
    } else if (ldA < max(1,n)) {
        info = -5;
    } else if (ldB < nrhs) {
        info = -8;
    }
    if (info != 0) {
        std::cerr << "DGETRS " << -info;
        return;
    }

    // Quick return if possible
    if (n == 0 || nrhs == 0) {
        return;
    }

    if(trans==BLAS::NoTrans) {
        // Solve A * X = B.   

        // Apply row interchanges to the right hand sides. 
        LAPACK::laswp(BLAS::RowMajor, nrhs, B, ldB, 0, n, piv, 1);

       // Solve L*X = B, overwriting B with X.
        BLAS::trsm(BLAS::RowMajor, BLAS::Left, BLAS::Lower, BLAS::NoTrans,BLAS::Unit, n, nrhs, one, A, ldA, B, ldB);
                
        // Solve U*X = B, overwriting B with X.
        BLAS::trsm(BLAS::RowMajor, BLAS::Left, BLAS::Upper, BLAS::NoTrans, BLAS::NonUnit, n, nrhs, one, A, ldA, B, ldB);
    } else {
        // Solve A' * X = B. 

        // Solve U'*X = B, overwriting B with X.
        BLAS::trsm(BLAS::RowMajor, BLAS::Left, BLAS::Upper, BLAS::Trans, BLAS::NonUnit, n, nrhs, one, A, ldA, B, ldB);

        // Solve L'*X = B, overwriting B with X.
        BLAS::trsm(BLAS::RowMajor, BLAS::Left, BLAS::Lower, BLAS::Trans, BLAS::Unit, n, nrhs, one, A, ldA, B, ldB);

        // Apply row interchanges to the solution vectors.
        LAPACK::laswp(BLAS::RowMajor, nrhs, B, ldB, 1, n, piv, -1);
    }

    return;

}

#endif // __LAPACK_GETRS_HPP__
