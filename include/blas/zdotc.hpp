#ifndef __BLAS_ZDOTC_HPP__
#define __BLAS_ZDOTC_HPP__

#include "blas.hpp"

template<typename real>
BLAS::complex<real>
BLAS::dotc (const int N, const complex<real> *X, const int incX, 
            const complex<real> *Y, const int incY)
{
{
  complex<real> result=0;
  int i;
  int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
  int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
  for (i = 0; i < N; i++) {
    result += conj(*X) * *Y;
    ix += incX;
    iy += incY;
  }
  return result;
}
}

#endif // __BLAS_ZDOTC_HPP__
