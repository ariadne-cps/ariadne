#ifndef __LAPACK_GESV_HPP__
#define __LAPACK_GESV_HPP__

#include "lapack.hpp"

template<typename Scalar>
void 
LAPACK::gesv(BLAS::ORDER order, int n, int nrhs, Scalar *A, int ldA, int *P, Scalar *B, int ldB)
{
  std::cerr << "LAPACK::gesv\n";
/*  -- LAPACK driver routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DGESV computes the solution to a real system of linear equations   
       A * X = B,   
    where A is an N-by-N matrix and X and B are N-by-NRHS matrices.   

    The LU decomposition with partial pivoting and row interchanges is   
    used to factor A as   
       A = P * L * U,   
    where P is a permutation matrix, L is unit lower triangular, and U is   
    upper triangular.  The factored form of A is then used to solve the   
    system of equations A * X = B.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of linear equations, i.e., the order of the   
            matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the N-by-N coefficient matrix A.   
            On exit, the factors L and U from the factorization   
            A = P*L*U; the unit diagonal elements of L are not stored.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    IPIV    (output) INTEGER array, dimension (N)   
            The pivot indices that define the permutation matrix P;   
            row i of the matrix was interchanged with row IPIV(i).   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the N-by-NRHS matrix of right hand side matrix B.   
            On exit, if INFO = 0, the N-by-NRHS solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization   
                  has been completed, but the factor U is exactly   
                  singular, so the solution could not be computed.   

    =====================================================================   

    /* Function Body */
    assert(order==BLAS::RowMajor);
    
    // Test the input parameters.
    int info = 0;
    if (n < 0) {
        info = -1;
    } else if (nrhs < 0) {
        info = -2;
    } else if (ldA < max(1,n)) {
        info = -4;
    } else if (ldB < nrhs) {
        info = -7;
    }
    if (info != 0) {
        std::cerr << "DGESV " << -info << std::endl;
        return;
    }

    // Compute the LU factorization of A.
    getrf(order, n, n, A, ldA, P);
    // Solve the system A*X = B, overwriting B with X.
    getrs(order, BLAS::NoTrans, n, nrhs, A, ldA, P, B, ldB);

    return;

}

#endif // __LAPACK_GESV_HPP__
