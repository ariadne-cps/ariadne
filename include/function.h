/***************************************************************************
 *            function.h
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
  
#ifndef _ARIADNE_FUNCTION_H
#define _ARIADNE_FUNCTION_H

/// Elementary functions for double, dyadic and rational types and intervals.

#include <cmath>
#include "numerical_type.h"
#include "interval.h"

namespace Ariadne {
  
  template<typename R> 
  R approx_exp(R x, R e) {
    R result=0;
    R term=1;
    R error=term;
    R n=0;    
    while(error>e) {
      result+=term;
      n+=1;
      term*=x/n;
      error=term;
    }
    return result;
  }
  
  template<typename R> 
  R approx_cos(R x, R e) {
    R result=0;
    R term=1;
    R error=term;
    R n=0;
    
    while(error>e) {
      result+=term;
      n+=2;
      term*=-x*x/(n*(n-1));
      error=term;
    }
    
    return result;
  }
  
  template<typename R> 
  R approx_sin(R x, R e) {
    R result=0;
    R term=x;
    R error=term;
    R n=1;
    
    while(error>e) {
      result+=term;
      n+=2;
      term*=-x*x/(n*(n-1));
      error=term;
    }
    
    return result;
  }

}

#endif /* _ARIADNE_FUNCTION_H */
