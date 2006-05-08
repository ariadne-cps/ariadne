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

#ifndef __BLAS_AMAX_HPP__
#define __BLAS_AMAX_HPP__

#include "blas.hpp"

template<typename real>
real
BLAS::amax (const int N, const real *X, const int incX)
{
{
  real mx = 0.0;
  int ix = 0;
  int i;

  if (incX <= 0) {
    return 0;
  }

  for (i = 0; i < N; i++) {
    if (abs(X[ix]) > mx) {
      mx = abs(X[ix]);
    }
    ix += incX;
  }

  return mx;
}
}

template<typename real>
real
BLAS::amax (const int N, const complex<real> *X, const int incX)
{
{
  real mx = 0.0;
  int ix = 0;
  int i;

  if (incX <= 0) {
    return 0;
  }

  for (i = 0; i < N; i++) {
    if (abs(X[ix]) > mx) {
      mx = abs(X[ix]);
    }
    ix += incX;
  }

  return mx;
}
}


#endif // __BLAS_AMAX_HPP__
