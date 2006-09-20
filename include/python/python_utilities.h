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
 
#ifndef _ARIADNE_PYTHON_UTILITIES_H
#define _ARIADNE_PYTHON_UTILITIES_H

#include <cstring>
#include <functional>


#include "numeric/float64.h"
#include "numeric/mpfloat.h"
#include "numeric/dyadic.h"
#include "numeric/rational.h"

template<typename R> inline std::string python_name(const std::string& bn);
template<> inline std::string python_name<Ariadne::Numeric::Float64>(const std::string& bn) { return "F"+bn; }
//template<> inline const char* python_name<Ariadne::Numeric::MPFloat>(const std::string& bn) { return "R"+bn; }
template<> inline std::string python_name<Ariadne::Numeric::MPFloat>(const std::string& bn) { return "R"+bn; }
template<> inline std::string python_name<Ariadne::Numeric::Dyadic>(const std::string& bn) { return "D"+bn; }
template<> inline std::string python_name<Ariadne::Numeric::Rational>(const std::string& bn) { return "Q"+bn; }


template<typename R> inline const char* prefix(); 
template<> inline const char* prefix<Ariadne::Numeric::Float64>() { return "F"; }
//template<> inline std::string prefix<Ariadne::Numeric::MPFloat>() { return "R"; }
template<> inline const char* prefix<Ariadne::Numeric::MPFloat>() { return ""; }
template<> inline const char*  prefix<Ariadne::Numeric::Dyadic>() { return "D"; }
template<> inline const char*  prefix<Ariadne::Numeric::Rational>() { return "Q"; }


template<class C> 
inline
typename C::value_type 
get_item(const C& c, int n) {
  if(n<0) {
    n+=c.size();
  }
  assert(0<=n);
  size_t m=size_t(n);
  assert(m<c.size());
  return c[m];
}

template<class C, class T> 
inline
void
set_item_from(C& c, int n, const T& x) {
  if(n<0) {
    n+=c.size();
  }
  assert(0<=n);
  size_t m=size_t(n);
  assert(m<c.size());
  c[n]=x;
}

template<class C> 
inline
void
set_item(C& c, int n, const typename C::value_type& x) {
  if(n<0) {
    n+=c.size();
  }
  assert(0<=n);
  size_t m=size_t(n);
  assert(m<c.size());
  c[n]=x;
}

template<typename Res, typename Arg1, typename Arg2, typename Op>
inline
Res
evaluate(const Arg1& a1, const Arg2& a2)
{
  Op op;
  return Res(op(a1,a2));
}

template<typename Res, typename Arg1, typename Arg2>
inline
Res
add(const Arg1& a1, const Arg2& a2)
{
  return Res(a1+a2);
}

template<typename Res, typename Arg1, typename Arg2>
inline
Res
sub(const Arg1& a1, const Arg2& a2)
{
  return Res(a1-a2); 
}

template<typename Res, typename Arg1, typename Arg2>
inline
Res
mul(const Arg1& a1, const Arg2& a2)
{
  return Res(a1*a2);
}


template<typename Res, typename Arg1, typename Arg2>
inline
Res
div(const Arg1& a1, const Arg2& a2)
{
  return Res(a1/a2);
}


#endif /* _ARIADNE_PYTHON_UTILITIES_H */
