/***************************************************************************
 *            geometry/interval_set.inline.h
 *
 *  Copyright 2007  Pieter Collins
 *  
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

#include "numeric/arithmetic.h"
 
namespace Ariadne {
  
template<class X> inline 
Geometry::IntervalSet<X>::IntervalSet()
  : _lower(1), _upper(0) 
{
}

template<class X> inline 
Geometry::IntervalSet<X>::IntervalSet(const X& x)
  : _lower(x), _upper(x) 
{
}

template<class X> inline 
Geometry::IntervalSet<X>::IntervalSet(const X& l, const X& u)
  : _lower(l), _upper(u) 
{
}

template<class X> inline 
Geometry::IntervalSet<X>::IntervalSet(const Geometry::IntervalSet<X>& ivl)
  : _lower(ivl._lower), _upper(ivl._upper) 
{
}

template<class X> inline 
Geometry::IntervalSet<X>& 
Geometry::IntervalSet<X>::operator=(const X& x) 
{
  this->_lower=x; this->_upper=x; return *this;
}

template<class X> inline 
Geometry::IntervalSet<X>& 
Geometry::IntervalSet<X>::operator=(const Geometry::IntervalSet<X>& ivl) 
{
  this->_lower=ivl._lower; this->_upper=ivl._upper; return *this;
}

template<class X> template<class XL,class XU> inline 
Geometry::IntervalSet<X>::IntervalSet(const XL& l, const XU& u)
  : _lower(l), _upper(u)
{
}

template<class X> template<class XX> inline 
Geometry::IntervalSet<X>::IntervalSet(const IntervalSet<XX>& ivl)
  : _lower(ivl.lower()), _upper(ivl.upper()) 
{
}

template<class X> template<class XX> inline 
Geometry::IntervalSet<X>::IntervalSet(const XX& x)
  : _lower(x), _upper(x) 
{
}

template<class X> template<class XX> inline 
Geometry::IntervalSet<X>& 
Geometry::IntervalSet<X>::operator=(const XX& x) 
{
  this->_lower=x; this->_upper=x; return *this;
}



template<class X> inline 
bool
Geometry::IntervalSet<X>::operator==(const IntervalSet<X>& ivl) {
  return this->_lower==ivl._lower && this->_upper==ivl._upper;
}


template<class X> inline 
bool
Geometry::IntervalSet<X>::operator!=(const IntervalSet<X>& ivl) {
  return !(*this==ivl);
}


template<class X> inline 
const X& 
Geometry::IntervalSet<X>::lower() const 
{ 
  return this->_lower; 
}

template<class X> inline 
const X& 
Geometry::IntervalSet<X>::lower_bound() const 
{ 
  return this->_lower; 
}

template<class X> inline 
const X& 
Geometry::IntervalSet<X>::upper_bound() const 
{ 
  return this->_upper; 
}

template<class X> inline 
const X& 
Geometry::IntervalSet<X>::upper() const 
{ 
  return this->_upper; 
}





template<class X> 
std::ostream& 
Geometry::IntervalSet<X>::write(std::ostream& os) const {
  if(this->_lower > this->_upper) {
    return os << "[1,0]";
  }
  else {
    return os << "[" << this->lower() << "," << this->upper() << "]";
  }
}


template<class X> 
std::istream& 
Geometry::IntervalSet<X>::read(std::istream& is) {
  char c;
  X l;
  X u;
  is >> c;
  if(c=='[') {
    is >> l >> c;
    if(c!=',' && c!=':') {
      is.setstate(std::ios_base::failbit);
    }
    is >> u >> c;
    if(c!=']') {
      is.setstate(std::ios_base::failbit);
    }
    (*this)=Geometry::IntervalSet<X>(l,u);
  } else {
    is.putback(c);
    is >> l;
    (*this)=Geometry::IntervalSet<X>(l);
  }
  return is;
}





template<class X>
X
Geometry::centre(const IntervalSet<X>& ivl) 
{
  return Numeric::med_approx(ivl.lower(),ivl.upper());
}

template<class X>
X
Geometry::radius(const IntervalSet<X>& ivl) 
{
  return Numeric::rad_up(ivl.lower(),ivl.upper());
}

template<class X>
tribool
Geometry::empty(const IntervalSet<X>& ivl) 
{
  if(ivl.lower() < ivl.upper()) {
    return false;
  } else if(ivl.lower() > ivl.upper()) {
    return true;
  } else {
    return indeterminate;
  }
}

template<class X>
tribool
Geometry::bounded(const IntervalSet<X>& ivl) 
{
  X inf=Numeric::inf<X>();
  X minf=-inf;
  if(ivl.lower() > minf && ivl.upper() < inf) {
    return true;
  } else if(ivl.lower() == minf || ivl.upper() == inf) {
    return false;
  } else { // May occur if one of the bounds is \a nan or an unbounded interval
    return indeterminate;
  }
}





template<class X1, class X2> inline
tribool 
Geometry::contains(const IntervalSet<X1>& ivl, const X2& x)  
{ 
  if(ivl.lower() < x && x < ivl.upper()) {
    return true;
  } else if(ivl.lower() > x || x > ivl.upper()) {
    return false;
  } else {
    return indeterminate;
  }
}

template<class X1, class X2> inline
tribool 
Geometry::disjoint(const IntervalSet<X1>& ivl1, const IntervalSet<X2>& ivl2)  
{ 
  if(ivl1.lower() > ivl2.upper() || ivl1.upper() < ivl2.lower()) {
    return true;
  } else if(ivl1.lower() < ivl2.upper() && ivl1.upper() > ivl2.lower()) {
    return false;
  } else {
    return indeterminate;
  }
}

template<class X1, class X2> inline
tribool 
Geometry::intersect(const IntervalSet<X1>& ivl1, const IntervalSet<X2>& ivl2)  
{ 
  return !Geometry::disjoint(ivl1,ivl2);
}

template<class X1, class X2> inline
tribool 
Geometry::subset(const IntervalSet<X1>& ivl1, const IntervalSet<X2>& ivl2)  
{ 
  if(ivl1.lower() > ivl2.lower() && ivl1.upper() < ivl2.upper()) {
    return true;
  } else if(ivl1.lower() < ivl2.lower() || ivl1.upper() > ivl2.upper()) {
    return false;
  } else {
    return indeterminate;
  }
}

template<class X1, class X2> inline
tribool 
Geometry::superset(const IntervalSet<X1>& ivl1, const IntervalSet<X2>& ivl2)  
{ 
  return Geometry::subset(ivl2,ivl1);
}

template<class X> inline
Geometry::IntervalSet<X>
Geometry::closed_intersection(const IntervalSet<X>& ivl1, const IntervalSet<X>& ivl2)  
{ 
  if(disjoint(ivl1,ivl2)) {
    return IntervalSet<X>(1,0);
  }
  return IntervalSet<X>(Numeric::max(ivl1.lower(),ivl2.lower()),Numeric::min(ivl1.upper(),ivl2.upper()));
}

template<class X> inline
Geometry::IntervalSet<X>
Geometry::open_intersection(const IntervalSet<X>& ivl1, const IntervalSet<X>& ivl2)  
{ 
  if(possibly(disjoint(ivl1,ivl2))) {
    return IntervalSet<X>(1,0);
  }
  return IntervalSet<X>(Numeric::max(ivl1.lower(),ivl2.lower()),Numeric::min(ivl1.upper(),ivl2.upper()));
}

template<class X> inline
Geometry::IntervalSet<X>
Geometry::interval_hull(const IntervalSet<X>& ivl1, const IntervalSet<X>& ivl2)  
{ 
  return IntervalSet<X>(Numeric::min(ivl1.lower(),ivl2.lower()),Numeric::max(ivl1.upper(),ivl2.upper()));
}



template<class X> inline
std::ostream& 
Geometry::operator<<(std::ostream& os, const IntervalSet<X>& x) {
  return x.write(os);
}


template<class X> inline
std::istream& 
Geometry::operator>>(std::istream& is, IntervalSet<X>& x) {
  return x.read(is);
}

} // namespace Ariadne
