#ifndef __BLAS_TRSET_HPP__
#define __BLAS_TRSET_HPP__

#include "blas.hpp"
 
template<typename scalar, typename scalarA>
void
BLAS::trset(ORDER order, UPLO uplo, DIAG diag, int m, int n, 
            scalar alpha, scalarA *A, int ldA)
{
  if ( (order==RowMajor && uplo==Lower) || (order==ColMajor && uplo==Upper)) {
    for(int i=0; i!=m; ++i) {
      for(int j=0; j!=i; ++j) {
        A[i*ldA+j]=alpha;
      }
    }
  } else {  
    for(int i=0; i!=m; ++i) {
      for(int j=i+1; j!=n; ++j) {
        A[i*ldA+j]=alpha;
      }
    }
  }
  
  if(diag==NonUnit) {
    for(int i=0; i!=min(m,n); ++i) {
      A[i*ldA+i]=alpha;
    }
  }
  return;
}

#endif // __BLAS_TRSET_HPP__
