/*
 * Copyright (C) 2006 Pieter Collins <Pieter.Collins@cwi.nl>
 *
 * Based on the BLAS implementation in Gnu Scientific Library 1.8
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 */

/*  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef __BLAS_NRM2_HPP__
#define __BLAS_NRM2_HPP__

#include "blas.hpp"

template<typename real>
real
BLAS::nrm2 (const int N, const real *X, const int incX)
{
{
  real scale = 0.0;
  real ssq = 1.0;
  int i;
  int ix = 0;
  if (N <= 0 || incX <= 0) {
    return 0;
  } else if (N == 1) {
    return abs(X[0]);
  }
  for (i = 0; i < N; i++) {
    const real x = X[ix];
    if (x != 0.0) {
      const real ax = abs(x);
      if (scale < ax) {
        ssq = 1.0 + ssq * (scale / ax) * (scale / ax);
        scale = ax;
      } else {
        ssq += (ax / scale) * (ax / scale);
      }
    }
    ix += incX;
  }
  return scale * sqrt(ssq);
}
}

template<typename real>
real
BLAS::nrm2 (const int N, const complex<real> *X, const int incX)
{
{
  real scale = 0.0;
  real ssq = 1.0;
  int i;
  int ix = 0;
  if (N <= 0 || incX <= 0) {
    return 0;
  } else if (N == 1) {
    return abs(X[0]);
  }
  for (i = 0; i < N; i++) {
    const complex<real> x = X[ix];
    if (x != 0.0) {
      const real ax = abs(x);
      if (scale < ax) {
        ssq = 1.0 + ssq * (scale / ax) * (scale / ax);
        scale = ax;
      } else {
        ssq += (ax / scale) * (ax / scale);
      }
    }
    ix += incX;
  }
  return scale * sqrt(ssq);
}
}

#endif // __BLAS_NRM2_HPP__
