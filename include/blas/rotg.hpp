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
