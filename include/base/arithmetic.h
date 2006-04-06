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
 
#ifndef _ARITHMETIC_H
#define _ARITHMETIC_H

#include <cmath>
#include "../base/numerical_type.h"
#include "../base/approximation.h"

namespace Ariadne {
  
  /*! \brief Unary negation. */
  template<typename R> 
  inline
  R neg(const R& x) {
    return -x;
  }
  
  /*! \brief Absolute value. */
  template<typename R>
  inline
  R 
  abs(const R& x) {
    return (x>=R(0)) ? x : R(-x);
  }

  using std::min;
  using std::max;
  
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
    typedef typename numerical_traits<R>::field_extension_type Fld;
    return Fld(x1)/Fld(x2);
  }

  inline MPFloat div_approx(const MPFloat& x1, const MPFloat& x2, const MPFloat& e) {
    return div_approx(x1,x2,precision(e));
  }
    
  inline MPFloat div_approx(const MPFloat& x1, const MPFloat& x2, const uint& n) {
    MPFloat r(0,n);
    mpf_div(r.get_mpf_t(),x1.get_mpf_t(),x2.get_mpf_t());
    return r;
  }
 
  inline Dyadic div_approx(const Dyadic& x1, const Dyadic& x2, const Dyadic& e) {
    Rational q=div(Rational(x1),Rational(x2));
    return approximate<Dyadic>(q,e);
  }
    
  /*! \brief An integer n such that \f$n\leq x/y < n+1\f$. */
  template<typename R> int quotient(R x, R y);

  template<> inline int quotient(Float64 x, Float64 y)
  {
    assert(y!=0.0);
    int q = int(x/y+0.5);
    if(x < q*y) { q=q-1; }
    assert(q*y <= x && x < (q+1)*y);
    return q;
  }

  template<> inline int quotient(Rational x, Rational y) 
  {
    assert(y!=Rational(0));
    Rational d = x/y;
    Integer qz = d.get_num() / d.get_den();
    int q = int(qz.get_si());
    assert(q*y <= x && x < (q+1)*y);
    return q;
  }

  template<> inline int quotient(MPFloat x, MPFloat y)
  {
    return quotient(Rational(x),Rational(y));
  }

  template<> inline int quotient(Dyadic x, Dyadic y) 
  {
    return quotient(Rational(x),Rational(y));
  }
  
  
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
  log_floor(N n, R x) {
    assert(n>1 && x>=1);
    N result=0;
    while(x>=n) {
      x/=n;
      result+=1;
    }
    return result;
  }
    
  
  /*! \brief The ceiling of the logarithm of \a x in base \a n. */
  template<typename R, typename N>
  inline
  N 
  log_ceil(N n, R x) {
    assert(n>1 && x>=1);
    N result=0;
    while(x>1) {
      x/=n;
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



#endif /* _ARITHMETIC_H */
