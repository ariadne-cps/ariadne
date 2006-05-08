#ifndef __BLAS_SETM_HPP__
#define __BLAS_SETM_HPP__

#include "blas.hpp"

template<typename scalar>
void
BLAS::setm (const enum ORDER Order, const int M, const int N, 
            const scalar alpha, scalar *A, const int ldA)
{

  int i,j;
  if(Order==RowMajor) {
    for (i = 0; i < M; i++) {
      for (j = 0; j < N; j++) {
        A[ldA*i+j] = alpha;
      }
    }
  } else {
    for (i = 0; i < M; i++) {
      for (j = 0; j < N; j++) {
        A[i+ldA*j] = alpha;
      }
    }
  }

}

  

#endif // __BLAS_SETM_HPP__
