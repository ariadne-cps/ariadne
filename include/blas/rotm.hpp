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
