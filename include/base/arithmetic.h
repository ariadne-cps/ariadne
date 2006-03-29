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
#include "../base/numerical_type.h"

namespace Ariadne {
  inline
  uint 
  factorial(const uint& n) {
    uint result=1;
    for(uint i=1; i<=n; ++i) {
      result*=i;
    }
    return result;
  }

  inline
  Rational 
  pow(const Rational& x, const uint& n) {
    Rational result=1;
    for(uint i=0; i!=n; ++i) {
      result*=x;
    }
    return result;
  }

  inline
  int 
  pow(int x, uint n) {
    int result=1;
    for(uint i=0; i!=n; ++i) {
      result*=x;
    }
    return result;
  }
  
  inline
  uint 
  pow(uint x, uint n) {
    uint result=1;
    for(uint i=0; i!=n; ++i) {
      result*=x;
    }
    return result;
  }
  
  inline
  uint 
  log_floor(uint n, uint x) {
    assert(n>1 && x>0);
    uint result=0;
    while(x>=n) {
      x/=n;
      result+=1;
    }
    return result;
  }
  
  inline
  uint 
  log_ceil(uint n, uint x) {
    assert(n>1 && x>0);
    if(x==1) {
      return 0;
    }
    return log_floor(n,x-1)+1;
  }
  
  template<typename R>
  inline
  R 
  abs(const R& x) {
    return (x>=R(0)) ? x : R(-x);
  }


}



#endif /* _ARITHMETIC_H */
