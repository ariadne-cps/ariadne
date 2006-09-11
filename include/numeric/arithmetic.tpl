/***************************************************************************
 *            arithmetic.tpl
 *
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
 

#include "arithmetic.h"
#include "float64.h"
#include "mpfloat.h"
#include "dyadic.h"
#include "rational.h"

namespace Ariadne {
  namespace Numeric {
    
    template<> int quot(const Float64& x, const Float64& y)
    {
      assert(y>=0.0);
      int q = int(x/y+0.5);
      if(x < q*y) { q=q-1; }
      assert(q*y <= x && x < (q+1)*y);
      return q;
    }
  
    template<> int quot(const Rational& x, const Rational& y) 
    {
      assert(y>=0);

      Rational d = x/y;
      Integer qz= d.get_num() / d.get_den();
      if(d<0 && qz*y!=x) {
        qz-=1;
      }
      int q = qz.get_si();
      assert(q*y <= x && x < (q+1)*y);
      return q;
    }
  
    template<> int quot(const MPFloat& x, const MPFloat& y)
    {
      return quot<int>(Rational(x),Rational(y));
    }
  
    template<> int quot(const Dyadic& x, const Dyadic& y) 
    {
      return quot<int>(Rational(x),Rational(y));
    }
  
  }
}
