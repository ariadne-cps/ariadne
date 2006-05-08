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

#ifndef __BLAS_AXPY_HPP__
#define __BLAS_AXPY_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::axpy (const int N, const Real alpha, const Real *X, const int incX,
             Real *Y, const int incY)
{
{
  int i;
  if (alpha == 0.0) {
    return;
  }
  if (incX == 1 && incY == 1) {
    const int m = N % 4;
    for (i = 0; i < m; i++) {
      Y[i] += alpha * X[i];
    }
    for (i = m; i + 3 < N; i += 4) {
      Y[i] += alpha * X[i];
      Y[i + 1] += alpha * X[i + 1];
      Y[i + 2] += alpha * X[i + 2];
      Y[i + 3] += alpha * X[i + 3];
    }
  } else {
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));
    int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));
    for (i = 0; i < N; i++) {
      Y[iy] += alpha * X[ix];
      ix += incX;
      iy += incY;
    }
  }
}
}

#endif // __BLAS_AXPY_HPP__
