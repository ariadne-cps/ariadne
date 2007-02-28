/***************************************************************************
 *            python/python_utilities.h
 *
 *  16 November 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

/*! \file python_utilities.h
 *  Commonly used inline methods for the Python interface.
 */
 
#ifndef ARIADNE_PYTHON_UTILITIES_H
#define ARIADNE_PYTHON_UTILITIES_H

#include <cstring>
#include <functional>


#include "numeric/float64.h"
#include "numeric/mpfloat.h"
#include "numeric/rational.h"

template<class R> inline std::string python_name(const std::string& bn);

template<> inline std::string python_name<bool>(const std::string& bn) { return "Boolean"+bn; }
template<> inline std::string python_name<Ariadne::index_type>(const std::string& bn) { return "Index"+bn; }
template<> inline std::string python_name<Ariadne::size_type>(const std::string& bn) { return "Size"+bn; }
template<> inline std::string python_name<Ariadne::Numeric::Integer>(const std::string& bn) { return "Z"+bn; }
template<> inline std::string python_name<Ariadne::Numeric::Rational>(const std::string& bn) { return "Q"+bn; }

template<> inline std::string python_name<Ariadne::Numeric::Float64>(const std::string& bn) { return "F64"+bn; }
//template<> inline std::string python_name<Ariadne::Numeric::MPFloat>(const std::string& bn) { return "MPF"+bn; }
template<> inline std::string python_name<Ariadne::Numeric::MPFloat>(const std::string& bn) { return ""+bn; }



template<class C> 
inline
typename C::value_type 
get_item(const C& c, int n) {
  if(n<0) {
    n+=c.size();
  }
  if(n<0) { throw std::out_of_range("Index out-of-range"); }
  size_t m=size_t(n);
  if(c.size()<=m) { throw std::out_of_range("Index out-of-range"); }
  return c[m];
}


template<class C, class T> 
inline
void
set_item_from(C& c, int n, const T& x) {
  if(n<0) {
    n+=c.size();
  }
  if(n<0) { throw std::out_of_range("Index out-of-range"); }
  size_t m=size_t(n);
  if(c.size()<=m) { throw std::out_of_range("Index out-of-range"); }
  c[n]=x;
}


template<class C> 
inline
void
set_item(C& c, int n, const typename C::value_type& x) {
  if(n<0) {
    n+=c.size();
  }
  if(n<0) { throw std::out_of_range("Index out-of-range"); }
  size_t m=size_t(n);
  if(c.size()<=m) { throw std::out_of_range("Index out-of-range"); }
  c[n]=x;
}



template<class Res, class Arg1, class Arg2, class Op>
inline
Res
evaluate(const Arg1& a1, const Arg2& a2)
{
  return Res(Op()(a1,a2));
}


template<class Res, class Arg> inline
Res neg(const Arg& a) {
  return Res(-a);
}

template<class Res, class Arg1, class Arg2, class Tmp1, class Tmp2> inline
Res add(const Arg1& a1, const Arg2& a2) {
  return Res(Tmp1(a1)+Tmp2(a2));
}

template<class Res, class Arg1, class Arg2> inline
Res add(const Arg1& a1, const Arg2& a2) {
  return Res(a1+a2);
}


template<class Res, class Arg1, class Arg2, class Tmp1, class Tmp2> inline
Res sub(const Arg1& a1, const Arg2& a2) {
  return Res(Tmp1(a1)-Tmp2(a2));
}

template<class Res, class Arg1, class Arg2> inline
Res sub(const Arg1& a1, const Arg2& a2) {
  return Res(a1-a2);
}

template<class Res, class Arg1, class Arg2, class Tmp1, class Tmp2> inline
Res rsub(const Arg1& a1, const Arg2& a2) {
  return Res(Tmp2(a2)-Tmp1(a1));
}

template<class Res, class Arg1, class Arg2> inline
Res rsub(const Arg1& a1, const Arg2& a2) {
  return Res(a2-a1);
}


template<class Res, class Arg1, class Arg2, class Tmp1, class Tmp2> inline
Res mul(const Arg1& a1, const Arg2& a2) {
 return Res(Tmp1(a1)*Tmp2(a2));
}

template<class Res, class Arg1, class Arg2> inline
Res mul(const Arg1& a1, const Arg2& a2) {
  return Res(a1*a2);
}


template<class Res, class Arg1, class Arg2, class Tmp1, class Tmp2> inline
Res div(const Arg1& a1, const Arg2& a2) {
  return Res(Tmp1(a1)/Tmp2(a2));
}

template<class Res, class Arg1, class Arg2> inline
Res div(const Arg1& a1, const Arg2& a2) {
  return Res(a1/a2);
}

template<class Res, class Arg1, class Arg2, class Tmp1, class Tmp2> inline
Res rdiv(const Arg1& a1, const Arg2& a2) {
  return Res(Tmp2(a2)/Tmp1(a1));
}

template<class Res, class Arg1, class Arg2> inline
Res rdiv(const Arg1& a1, const Arg2& a2) {
  return Res(a2/a1);
}



template<class Res,class Arg1,class Arg2> inline
Res eq(const Arg1& a1, const Arg2& a2) {
  return a1==a2;
}

template<class Res,class Arg1,class Arg2> inline
Res ne(const Arg1& a1, const Arg2& a2) {
  return a1!=a2;
}

template<class Res,class Arg1,class Arg2> inline
Res lt(const Arg1& a1, const Arg2& a2) {
  return a1< a2;
}

template<class Res,class Arg1,class Arg2> inline
Res le(const Arg1& a1, const Arg2& a2) {
  return a1<=a2;
}

template<class Res,class Arg1,class Arg2> inline
Res gt(const Arg1& a1, const Arg2& a2) {
  return a1> a2;
}

template<class Res,class Arg1,class Arg2> inline
Res ge(const Arg1& a1, const Arg2& a2) {
  return a1>=a2;
}


#endif /* ARIADNE_PYTHON_UTILITIES_H */
