#ifndef __BLAS_ROTM_HPP__
#define __BLAS_ROTM_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::rotm (const int N, Real *X, const int incX, Real *Y,
             const int incY, const Real *P)
{
{
  int n;
  int i = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
  int j = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
  Real h11, h21, h12, h22;
  if (P[0] == -1.0) {
    h11 = P[1];
    h21 = P[2];
    h12 = P[3];
    h22 = P[4];
  } else if (P[0] == 0.0) {
    h11 = 1.0;
    h21 = P[2];
    h12 = P[3];
    h22 = 1.0;
  } else if (P[0] == 1.0) {
    h11 = P[1];
    h21 = -1.0;
    h12 = 1.0;
    h22 = P[4];
  } else if (P[0] == -2.0) {
    return;
  } else {
    xerbla(0, "source_rotm.h", "unrecognized value of P[0]");;
    return;
  }
  for (n = 0; n < N; n++) {
    const Real w = X[i];
    const Real z = Y[j];
    X[i] = h11 * w + h12 * z;
    Y[j] = h21 * w + h22 * z;
    i += incX;
    j += incY;
  }
}
}

#endif // __BLAS_ROTM_HPP__
