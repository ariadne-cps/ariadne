#ifndef __BLAS_GESET_HPP__
#define __BLAS_GESET_HPP__

#include "blas.hpp"
 
template<typename scalar, typename scalarA>
void
BLAS::geset(ORDER order, int m, int n, 
            scalar alpha, scalar beta, scalarA *A, int ldA)
{
  if (order==RowMajor) {
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
  return;
}

#endif // __BLAS_GESET_HPP__
