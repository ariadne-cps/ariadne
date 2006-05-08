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

#ifndef __BLAS_ROTMG_HPP__
#define __BLAS_ROTMG_HPP__

#include "blas.hpp"

template<typename Real>
void
BLAS::rotmg (Real *d1, Real *d2, Real *b1, const Real b2, Real *P)
{
{
  const Real G = 4096.0, G2 = G * G;
  Real D1 = *d1, D2 = *d2, x = *b1, y = b2;
  Real h11, h12, h21, h22, u;
  Real c, s;
  if (D1 < 0.0) {
    P[0] = -1;
    P[1] = 0;
    P[2] = 0;
    P[3] = 0;
    P[4] = 0;
    *d1 = 0;
    *d2 = 0;
    *b1 = 0;
    return;
  }
  if (D2 * y == 0.0) {
    P[0] = -2;
    return;
  }
  c = abs(D1 * x * x);
  s = abs(D2 * y * y);
  if (c > s) {
    P[0] = 0.0;
    h11 = 1;
    h12 = (D2 * y) / (D1 * x);
    h21 = -y / x;
    h22 = 1;
    u = 1 - h21 * h12;
    if (u <= 0.0) {
      P[0] = -1;
      P[1] = 0;
      P[2] = 0;
      P[3] = 0;
      P[4] = 0;
      *d1 = 0;
      *d2 = 0;
      *b1 = 0;
      return;
    }
    D1 /= u;
    D2 /= u;
    x *= u;
  } else {
    if (D2 * y * y < 0.0) {
      P[0] = -1;
      P[1] = 0;
      P[2] = 0;
      P[3] = 0;
      P[4] = 0;
      *d1 = 0;
      *d2 = 0;
      *b1 = 0;
      return;
    }
    P[0] = 1;
    h11 = (D1 * x) / (D2 * y);
    h12 = 1;
    h21 = -1;
    h22 = x / y;
    u = 1 + h11 * h22;
    D1 /= u;
    D2 /= u;
    {
      Real tmp = D2;
      D2 = D1;
      D1 = tmp;
    }
    x = y * u;
  }
  while (D1 <= 1.0 / G2 && D1 != 0.0) {
    P[0] = -1;
    D1 *= G2;
    x /= G;
    h11 /= G;
    h12 /= G;
  }
  while (D1 >= G2) {
    P[0] = -1;
    D1 /= G2;
    x *= G;
    h11 *= G;
    h12 *= G;
  }
  while (abs(D2) <= 1.0 / G2 && D2 != 0.0) {
    P[0] = -1;
    D2 *= G2;
    h21 /= G;
    h22 /= G;
  }
  while (abs(D2) >= G2) {
    P[0] = -1;
    D2 /= G2;
    h21 *= G;
    h22 *= G;
  }
  *d1 = D1;
  *d2 = D2;
  *b1 = x;
  if (P[0] == -1.0) {
    P[1] = h11;
    P[2] = h21;
    P[3] = h12;
    P[4] = h22;
  } else if (P[0] == 0.0) {
    P[2] = h21;
    P[3] = h12;
  } else if (P[0] == 1.0) {
    P[1] = h11;
    P[4] = h22;
  }
}
}

#endif // __BLAS_ROTMG_HPP__
