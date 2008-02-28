/***************************************************************************
 *            interval_integer.h
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
 
/*! \file interval_integer.h
 *  \brief The integer representing intervals of integer number types. This header can be included to obtain declarations without needing to include definitions.
 */
 
#ifndef ARIADNE_INTERVAL_INTEGER_H
#define ARIADNE_INTERVAL_INTEGER_H

#include <iostream>
#include <stdexcept>

#include "base/tribool.h"
#include "numeric/expression.h"
#include "numeric/integer.h"

namespace Ariadne {
  namespace Numeric {

    template<class R> class Interval;
  
    template<>
    class Interval<Integer>
      : public Value< Interval<Integer> >
    {
      typedef Integer Z;
     public:
      Z _lower; Z _upper;
     public:
      //@{
      //! \name Constructors and assignment operators
      /*! \brief Default constructer constructs zero interval. */
      Interval() :  _lower(0), _upper(0) { }
      /*! \brief Construct a one-point interval. */
      template<class X> Interval(const X& x) : _lower(x), _upper(x) { }
      Interval(const Integer& x) : _lower(x), _upper(x) { }
      /*! \brief Construct from lower and upper bounds. */
      template<class XL, class XU> Interval(const XL& l, const XU& u) : _lower(l), _upper(u) { }
      /*! \brief Copy constructor. */
      template<class X> Interval(const Interval<X>& ivl) : _lower(ivl.lower()), _upper(ivl.upper()) { }

      /*! \brief Assign from a number. */
      template<class X> Interval<Z>& operator=(const X& x) { _lower=x; _upper=x; return *this; }
      /*! \brief Copy assignment operator. */
      template<class X> Interval<Z>& operator=(const Interval<X>& ivl) { _lower=ivl.lower(); _upper=ivl.upper(); return *this; }

      //@}
      
      //@{
      //! \name Data access
      /*! \brief The lower bound. */
      const Z& lower() const { return this->_lower; }
      /*! \brief The upper bound. */
      const Z& upper() const { return this->_upper; }
      /*! \brief The width of the interval, given by \f$b-a\f$. */
      Z width() const { return this->_upper-this->_lower; };
      //@}
    };      

    template<> inline std::string name<Interval<Integer> >() {
      return "ZInterval"; 
    }

    inline void pos_(Interval<Integer>& r, const Interval<Integer>& x, const Interval<Integer>& y) {
      r._lower=x._lower; r._upper=x._upper; 
    }
       
    inline void neg_noalias_(Interval<Integer>& r, const Interval<Integer>& x) {
      r._lower=-x._upper; r._upper=-x._lower; 
    }
       
    inline void neg_(Interval<Integer>& r, const Interval<Integer>& x) {
      if(&r==&x) { Interval<Integer> t(x); neg_noalias_(r,x); } else { neg_noalias_(r,x); }
    }
       
    inline void add_(Interval<Integer>& r, const Interval<Integer>& x, const Interval<Integer>& y) {
      r._lower=x._lower+y._lower; r._upper=x._upper+y._upper; 
    }
       
    inline void sub_noalias_(Interval<Integer>& r, const Interval<Integer>& x, const Interval<Integer>& y) {
      r._lower=x._lower-y._upper; r._upper=x._upper-y._lower; 
    }
       
    inline void sub_(Interval<Integer>& r, const Interval<Integer>& x, const Interval<Integer>& y) {
      if(&r==&y) { Interval<Integer> t(y); sub_noalias_(r,x,y); } else { sub_noalias_(r,x,y); }
    }
       
    inline void mul_(Interval<Integer>& r, const Interval<Integer>& x, const Interval<Integer>& y) {
      Integer b[4];
      b[0]=x.lower()*y.lower();
      b[1]=x.lower()*y.upper();
      b[2]=x.upper()*y.lower();
      b[3]=x.upper()*y.upper();
      r._lower=min(min(b[0],b[1]),min(b[2],b[3]));
      r._upper=max(max(b[0],b[1]),max(b[2],b[3]));
    }
       
  /*
    Interval<Integer> inline operator+(const Interval<Integer>& ivl1, const Interval<Integer>& ivl2) {
      return Interval<Integer>(ivl1.lower()+ivl2.lower(),ivl1.upper()+ivl2.upper()); 
    }
       
    Interval<Integer> inline operator-(const Interval<Integer>& ivl1, const Interval<Integer>& ivl2) {
      return Interval<Integer>(ivl1.lower()-ivl2.upper(),ivl1.upper()-ivl2.lower()); 
    }
       
    Interval<Integer> inline operator*(const Interval<Integer>& ivl1, const Interval<Integer>& ivl2) {
      Integer b[4];
      b[0]=ivl1.lower()*ivl2.lower();
      b[1]=ivl1.lower()*ivl2.upper();
      b[2]=ivl1.upper()*ivl2.lower();
      b[3]=ivl1.upper()*ivl2.upper();
      const Integer& l=min(min(b[0],b[1]),min(b[2],b[3]));
      const Integer& u=max(max(b[0],b[1]),max(b[2],b[3]));
      return Interval<Integer>(l,u);
    }

  */

  } 
}


#endif /* ARIADNE_INTERVAL_INTEGER_H */
