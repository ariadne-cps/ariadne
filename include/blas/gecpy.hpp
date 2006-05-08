#ifndef __BLAS_GECPY_HPP__
#define __BLAS_GECPY_HPP__

#include "blas.hpp"
 
template<typename scalarA, typename scalarB>
void
BLAS::gecpy(ORDER order, int m, int n, 
            const scalarA *A, int ldA, scalarB *B, int ldB)
{
  if (order==RowMajor) {
    for(int i=0; i!=m; ++i) {
      for(int j=0; j!=n; ++j) {
        B[i*ldB+j]=A[i*ldA+j];
      }
    }
  } else {  
    for(int i=0; i!=m; ++i) {
      for(int j=0; j!=n; ++j) {
        B[j*ldB+i]=A[j*ldA+i];
      }
    }
  }
  return;
}

#endif // __BLAS_GECPY_HPP__
