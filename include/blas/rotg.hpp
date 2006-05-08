#ifndef __BLAS_ROTG_HPP__
#define __BLAS_ROTG_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::rotg (Real *a, Real *b, Real *c, Real *s)
{
{
  const Real roe = (abs(*a) > abs(*b) ? *a : *b);
  const Real scale = abs(*a) + abs(*b);
  Real r, z;
  if (scale != 0.0) {
    const Real aos = *a / scale;
    const Real bos = *b / scale;
    r = scale * sqrt(aos * aos + bos * bos);
    r = sign(roe) * r;
    *c = *a / r;
    *s = *b / r;
    z = 1.0;
    if (abs(*a) > abs(*b))
      z = *s;
    if (abs(*b) >= abs(*a) && *c != 0.0)
      z = 1.0 / (*c);
  } else {
    *c = 1.0;
    *s = 0.0;
    r = 0.0;
    z = 0.0;
  }
  *a = r;
  *b = z;
}
}

#endif // __BLAS_ROTG_HPP__
