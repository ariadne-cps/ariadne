#ifndef __BLAS_TRCPY_HPP__
#define __BLAS_TRCPY_HPP__

#include "blas.hpp"
 
template<typename scalarA, typename scalarB>
void
BLAS::trcpy(ORDER order, UPLO uplo, DIAG diag, int m, int n, 
            const scalarA *A, int ldA, scalarB *B, int ldB)
{
  if ( (order==RowMajor && uplo==Lower) || (order==ColMajor && uplo==Upper)) {
    for(int i=0; i!=m; ++i) {
      for(int j=0; j!=i; ++j) {
        B[i*ldB+j]=A[i*ldA+j];
      }
    }
  } else {  
    for(int i=0; i!=m; ++i) {
      for(int j=i+1; j!=n; ++j) {
        B[i*ldB+j]=A[i*ldA+j];
      }
    }
  }
  
  if(diag==Unit) {
    for(int i=0; i!=min(m,n); ++i) {
      B[i*ldB+i]=1;
    }
  } else {
    for(int i=0; i!=min(m,n); ++i) {
      B[i*ldB+i]=A[i*ldA+i];
    }
  }
  
  return;
}

#endif // __BLAS_TRCPY_HPP__
