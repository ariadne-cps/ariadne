#ifndef __LAPACK_LASET_HPP__
#define __LAPACK_LASET_HPP__

#include "lapack.hpppp"

template<typename scaler>
void
laset(BLAS::ORDER order, int m, int n, scalar alpha, scalar beta, scalar *A, int ldA)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLASET initializes an m-by-n matrix A to BETA on the diagonal and   
    ALPHA on the offdiagonals.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    ALPHA   (input) DOUBLE PRECISION   
            The constant to which the offdiagonal elements are to be set.   

    BETA    (input) DOUBLE PRECISION   
            The constant to which the diagonal elements are to be set.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On exit, the leading m-by-n submatrix of A is set as follows:   

            if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,   
            if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,   
            otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,   

            and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

   =====================================================================   
*/
  assert(order==RowMajor);
  
  if(order==RowMajor) {
    for(int i=0; i!=m; ++i) {
      for(int j=0; j!=n; ++j) {
        A[i*ldA+j]=alpha;
      }
    }
  } else {
    for(int i=0; i!=m; ++i) {
      for(int j=0; j!=n; ++j) {
        A[i+j*ldA]=alpha;
      }
    }
  }

  for(int i=0; i!=min(m,n); ++i) {
    A[i*ldA+i]=beta;
  }

}

#endif // __LAPACK_LASET_HPP__
