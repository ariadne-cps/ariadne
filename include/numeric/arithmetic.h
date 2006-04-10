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
 
/*! \file arithmetic.h
 *  Simple arithmetic functions for double, dyadic and rational types and intervals.
 */
 
#ifndef _ARIADNE_ARITHMETIC_H
#define _ARIADNE_ARITHMETIC_H

#include <algorithm>

#include "../numeric/numerical_types.h"
#include "../numeric/approximation.h"

namespace Ariadne {
  namespace Numeric {
    
    /*! \brief Minimum. */
    template<typename R>
    inline
    R 
    min(const R& x1, const R& x2) {
      return (x1<x2) ? x1 : x2;
    }
  
    /*! \brief Maximum. */
    template<typename R>
    inline
    R 
    max(const R& x1, const R& x2) {
      return (x1>x2) ? x1 : x2;
    }
  
    /*! \brief Absolute value. */
    template<typename R>
    inline
    R 
    abs(const R& x) {
      return (x>=R(0)) ? x : R(-x);
    }
  
    /*! \brief Unary negation. */
    template<typename R> 
    inline
    R neg(const R& x) {
      return -x;
    }
    
    /*! \brief Addition. */
    template<typename R> 
    inline
    R add(const R& x1,const R& x2) {
      return x1+x2;
    }
    
    /*! \brief Subtraction. */
    template<typename R> 
    inline
    R sub(const R& x1,const R& x2) {
      return x1-x2;
    }
    
    /*! \brief Multiplication. */
    template<typename R> 
    inline
    R mul(const R& x1,const R& x2) {
      return x1*x2;
    }
    
    /*! \brief Division. */
    template<typename R> 
    inline
    typename numerical_traits<R>::field_extension_type
    div(const R& x1, const R& x2)
    {
      return _div(x1,x2, typename numerical_traits<R>::algebraic_category());
    }
    
    template<typename R> 
    inline
    R
    _div(const R& x1, const R& x2, field_tag)
    {
      return x1/x2;
    }
    
    template<typename R> 
    inline
    typename numerical_traits<R>::field_extension_type
    _div(const R& x1, const R& x2, ring_tag)
    {
      typedef typename numerical_traits<R>::field_extension_type F;
      return F(x1)/F(x2);
    }
  
    /*! \brief An integer n such that \f$n\leq x/y < n+1\f$. */
    template<typename R> int quotient(const R& x, const R& y);
    
    /*! \brief Factorial. */
    template<typename N>
    inline
    N 
    factorial(const N& n) {
      uint result=1;
      for(N i=1; i!=n; ++i) {
        result*=i;
      }
      return result;
    }
  
  
    /*! \brief Integer power. */
    template<typename R, typename N>
    inline
    R 
    pow(const R& x, const N& n) {
      R result=1;
      for(N i=0; i!=n; ++i) {
        result*=x;
      }
      return result;
    }
  
  
    /*! \brief The floor of the logarithm of \a x in base \a n. */
    template<typename R, typename N>
    inline
    N 
    log_floor(const N& n, const R& x) {
      assert(n>1 && x>=1);
      N result=0;
      R y=x;
      while(y>=n) {
        y/=n;
        result+=1;
      }
      return result;
    }
      
    
    /*! \brief The ceiling of the logarithm of \a x in base \a n. */
    template<typename R, typename N>
    inline
    N 
    log_ceil(const N& n, const R& x) {
      assert(n>1 && x>=1);
      N result=0;
      R y=x;
      while(y>1) {
        y/=n;
        result+=1;
      }
      return result;
    }
    
    
    /*! \brief Greatest common divisor. */
    template <typename N>
    inline 
    N 
    gcd(const N &a, const N &b)
    {
      N c=a%b;
      if (c==0) { 
        return b;
      }
      return (gcd(b, c));
    }
  
    
    /*! \brief Least common multiple. */
    template <typename N>
    inline 
    N 
    lcm(const N &a, const N &b) {
      return ((a*b)/gcd(a,b));
    }
  
  }
}



#endif /* _ARITHMETIC_H */
