/***************************************************************************
 *            interval_rational.h
 *
 *  Copyright 2005-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file interval_rational.h
 *  \brief The rational representing intervals of rational number types. This header can be included to obtain declarations without needing to include definitions.
 */
 
#ifndef ARIADNE_INTERVAL_RATIONAL_H
#define ARIADNE_INTERVAL_RATIONAL_H

#include <iostream>
#include <stdexcept>

#include "base/tribool.h"
#include "numeric/rational.h"
#include "numeric/traits.h"

namespace Ariadne {
  namespace Numeric {

    template<class R> class Interval;
  
    template<>
    class Interval<Rational>
      : public Value< Interval<Rational> >
    {
      typedef Rational Q;
     public:
      Q _lower; Q _upper;
     public:
      //@{
      //! \name Constructors and assignment operators
      /*! \brief Default constructer constructs zero interval. */
      Interval() :  _lower(0), _upper(0) { }
      /*! \brief Construct a one-point interval. */
      template<class X> Interval(const X& x) : _lower(x), _upper(x) { }
      Interval(const Rational& x) : _lower(x), _upper(x) { }
      /*! \brief Construct from lower and upper bounds. */
      template<class XL, class XU> Interval(const XL& l, const XU& u) : _lower(l), _upper(u) { }
      /*! \brief Copy constructor. */
      template<class X> Interval(const Interval<X>& ivl) : _lower(ivl.lower()), _upper(ivl.upper()) { }

      /*! \brief Assign from a number. */
      template<class X> Interval<Q>& operator=(const X& x) { _lower=x; _upper=x; return *this; }
      /*! \brief Copy assignment operator. */
      template<class X> Interval<Q>& operator=(const Interval<X>& ivl) { _lower=ivl.lower(); _upper=ivl.upper(); return *this; }

      /*! \brief Construct from an expression. */
      template<class E> Interval(const Expression<E>& e) { e.assign_to(*this); }
      /*! \brief Assign from an expression. */
      template<class E> Interval<Q>& operator=(const Expression<E>& e) { e.assign_to(*this); return *this;  }

      //@}
      
      //@{
      //! \name Data access
      /*! \brief The lower bound. */
      const Q& lower() const { return this->_lower; }
      /*! \brief The upper bound. */
      const Q& upper() const { return this->_upper; }
      /*! \brief The midpoint of the interval, given by \f$(a+b)/2\f$. */
      Q midpoint() const { return (this->_lower+this->_upper)/2; }
      /*! \brief The radius of the interval, given by \f$(b-a)/2\f$. */
      Q radius() const { return (this->_upper-this->_lower)/2; }
      /*! \brief The width of the interval, given by \f$b-a\f$. */
      Q width() const { return this->_upper-this->_lower; };
      //@}

      //@{
      //! \name Geometric operations
      /*! \brief Tests if the interval is empty. */
      bool empty() const { return this->_lower > this->_upper; }
      /*! \brief Tests if the interval consists of a single point. */
      bool singleton() const { return this->_lower == this->_upper; };
      /*! \brief Tests if the interval contains \a x. */
      template<class RX> bool encloses(const RX& x) const { 
        return this->_lower <= x && x <= this->_upper; }
      /*! \brief Tests if the interval contains \a r. */
      template<class RX> bool refines(const Interval<RX>& ivl) const { 
        return ivl._lower <= this->_lower && this->_upper <= ivl._upper; }
      //}
    };      

    template<> inline std::string name<Interval<Rational> >() {
      return "QInterval"; 
    }

    inline void pos_(Interval<Rational>& r, const Interval<Rational>& x) {
      r=x; }

    inline void neg_noalias_(Interval<Rational>& r, const Interval<Rational>& x) {
      r._lower=-x._upper; r._upper=-x._lower; } 

    inline void neg_(Interval<Rational>& r, const Interval<Rational>& x) {
      if(&r==&x) { Interval<Rational> t=x; neg_noalias_(r,t); } else { neg_noalias_(r,x); } }

    inline void add_(Interval<Rational>& r, const Interval<Rational>& x, const Interval<Rational>& y) {
      r._lower=x._lower+y._lower; r._upper=x._upper+y._upper; }

    inline void sub_noalias_(Interval<Rational>& r, const Interval<Rational>& x, const Interval<Rational>& y) {
      r._lower=x._lower-y._upper; r._upper=x._upper-y._lower; }
 
    inline void sub_(Interval<Rational>& r, const Interval<Rational>& x, const Interval<Rational>& y) {
      if(&r==&y) { Interval<Rational> t=y; sub_noalias_(r,x,t); } else { sub_noalias_(r,x,y); } }

    inline void mul_(Interval<Rational>& r, const Interval<Rational>& x, const Interval<Rational>& y) {
      Rational b[4];
      b[0]=x.lower()*y.lower();
      b[1]=x.lower()*y.upper();
      b[2]=x.upper()*y.lower();
      b[3]=x.upper()*y.upper();
      r._lower=min(min(b[0],b[1]),min(b[2],b[3]));
      r._upper=max(max(b[0],b[1]),max(b[2],b[3]));
    }
  
    inline void div_(Interval<Rational>& r, const Interval<Rational>& x, const Interval<Rational>& y) {
      if(y._lower<=0 && 0 <= y._upper) {
        inf_(r._upper); neg_(r._lower,r._upper); return; }
      Interval<Rational> t;
      div_(t._lower,1,y._upper);
      div_(t._upper,1,y._lower);
      mul_(r,x,t);
    }
  
    

  /*
    inline Interval<Rational> operator+(const Interval<Rational>& ivl) {
      return Interval<Rational>(ivl.upper(),ivl.lower());
    }
       
    inline Interval<Rational> operator-(const Interval<Rational>& ivl) {
      return Interval<Rational>(ivl.upper(),ivl.lower());
    }
       
    inline Interval<Rational> operator+(const Interval<Rational>& ivl1, const Interval<Rational>& ivl2) {
      return Interval<Rational>(ivl1.lower()+ivl2.lower(),ivl1.upper()+ivl2.upper()); 
    }
       
    inline Interval<Rational> operator-(const Interval<Rational>& ivl1, const Interval<Rational>& ivl2) {
      return Interval<Rational>(ivl1.lower()-ivl2.upper(),ivl1.upper()-ivl2.lower()); 
    }
       
    inline Interval<Rational> operator*(const Interval<Rational>& ivl1, const Interval<Rational>& ivl2) {
      Rational b[4];
      b[0]=ivl1.lower()*ivl2.lower();
      b[1]=ivl1.lower()*ivl2.upper();
      b[2]=ivl1.upper()*ivl2.lower();
      b[3]=ivl1.upper()*ivl2.upper();
      const Rational& l=min(min(b[0],b[1]),min(b[2],b[3]));
      const Rational& u=max(max(b[0],b[1]),max(b[2],b[3]));
      return Interval<Rational>(l,u);
    }

    inline Interval<Rational> operator/(const Interval<Rational>& ivl1, const Interval<Rational>& ivl2) {
      if(ivl2.lower()<=0 && 0<=ivl2.upper()) {
        Rational inf=Numeric::inf<Rational>(); 
        return Interval<Rational>(-inf,inf); 
      }
      Rational b[4];
      b[0]=ivl1.lower()/ivl2.lower();
      b[1]=ivl1.lower()/ivl2.upper();
      b[2]=ivl1.upper()/ivl2.lower();
      b[3]=ivl1.upper()/ivl2.upper();
      const Rational& l=min(min(b[0],b[1]),min(b[2],b[3]));
      const Rational& u=max(max(b[0],b[1]),max(b[2],b[3]));
      return Interval<Rational>(l,u);
    }

    inline Interval<Rational> operator/(const Interval<Rational>& ivl, const Rational& x) {
      return ivl/Interval<Rational>(x); }
     
    inline Interval<Rational> operator/(const Rational& x,const Interval<Rational>& ivl) {
      return Interval<Rational>(x)/ivl; }
     
  */

  } 
}



#endif /* ARIADNE_INTERVAL_RATIONAL_H */
