/***************************************************************************
 *            arithmetic.h
 *
 *  Wed 18 May 2005
 *  Copyright 2005  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
#ifndef _ARITHMETIC_H
#define _ARITHMETIC_H

/// Arithmetic for double, dyadic and rational types and intervals.

#include <cmath>
#include "numerical_type.h"

namespace Ariadne {
  uint factorial(const uint& n) {
    uint result=1;
    for(uint i=1; i<=n; ++i) {
      result*=i;
    }
    return result;
  }

  Rational pow(const Rational& x, const uint& n) {
    Rational result=1;
    for(uint i=0; i!=n; ++i) {
      result*=x;
    }
    return result;
  }

  int pow(int x, uint n) {
    int result=1;
    for(uint i=0; i!=n; ++i) {
      result*=x;
    }
    return result;
  }
}



#endif /* _ARITHMETIC_H */
