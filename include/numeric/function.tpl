/***************************************************************************
 *            function.tpl
 *
 *  Copyright 2005-6  Alberto Casagrande, Pieter Collins
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
  

#include "function.h"

#include "../numeric/integer.h"
#include "../numeric/numerical_types.h"
#include "../numeric/arithmetic.h"
#include "../numeric/approximation.h"

namespace Ariadne {
  namespace Numeric {
    
    template<>
    MPFloat 
    div_prec(const MPFloat& x1, const MPFloat& x2, const uint& n) 
    {
      MPFloat r(0,n);
      mpf_div(r.get_mpf_t(),x1.get_mpf_t(),x2.get_mpf_t());
      return r;
    }
 

    template<>
    Float64 
    div_approx(const Float64& x1, const Float64& x2, const Float64& e) 
    {
      return x1/x2;
    }
    
    template<>
    MPFloat div_approx(const MPFloat& x1, const MPFloat& x2, const MPFloat& e) 
    {
      return div_prec(x1,x2,convert_to<int>(precision(e)));
    }
    
    template<>
    Dyadic 
    div_approx(const Dyadic& x1, const Dyadic& x2, const Dyadic& e) 
    {
      Rational q=Numeric::div(Rational(x1),Rational(x2));
      return approximate<Dyadic>(q,e);
    }
    
    
    template<>
    Rational 
    div_approx(const Rational& x1, const Rational& x2, const Rational& e) 
    {
      return x1/x2;
    }
    
    template<typename R> 
    R 
    sqrt_approx(const R& x, const R& e) 
    {
      assert(false);
    }

    template<typename R> 
    R 
    exp_approx(const R& x, const R& e) 
    {
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
    R 
    exp_down(const R& x, const R& e) 
    {
      R result=0;
      R term=1;
      R error=term;
      R n=0;    
      while(error>e && x>n) {
        result+=term;
        n+=1;
        term*=x/n;
        error=term;
      }
      if(term<0) {
        result+=term;
      }
      return result;
    }
  
    template<typename R> 
    R 
    exp_up(const R& x, const R& e) 
    {
      R result=0;
      R term=1;
      R error=term;
      R n=0;    
      while(error>e && x>n) {
        result+=term;
        n+=1;
        term*=x/n;
        error=term;
      }
      if(term>0) {
        result+=term;
      }
      return result;
    }

    template<typename R> 
    R 
    cos_approx(const R& x, const R& e) 
    {
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
    R sin_approx(const R& x, const R& e) 
    {
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
}
