/***************************************************************************
 *            interval.inline.h
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
 
namespace Ariadne {
  
template<class X> inline 
Geometry::Interval<X>::Interval()
  : _lower(1), _upper(0) 
{
}

template<class X> inline 
Geometry::Interval<X>::Interval(const X& x)
  : _lower(x), _upper(x) 
{
}

template<class X> inline 
Geometry::Interval<X>::Interval(const X& l, const X& u)
  : _lower(l), _upper(u) 
{
}

template<class X> inline 
Geometry::Interval<X>::Interval(const Geometry::Interval<X>& ivl)
  : _lower(ivl._lower), _upper(ivl._upper) 
{
}

template<class X> inline 
Geometry::Interval<X>& 
Geometry::Interval<X>::operator=(const X& x) 
{
  this->_lower=x; this->_upper=x; return *this;
}

template<class X> inline 
Geometry::Interval<X>& 
Geometry::Interval<X>::operator=(const Geometry::Interval<X>& ivl) 
{
  this->_lower=ivl._lower; this->_upper=ivl._upper; return *this;
}

template<class X> template<class XL,class XU> inline 
Geometry::Interval<X>::Interval(const XL& l, const XU& u)
  : _lower(Numeric::conv_down<X>(l)), _upper(Numeric::conv_up<X>(u)) 
{
}

template<class X> template<class XX> inline 
Geometry::Interval<X>::Interval(const Interval<XX>& ivl)
  : _lower(Numeric::conv_down<X>(ivl.lower())), _upper(Numeric::conv_up<X>(ivl.upper())) 
{
}

template<class X> template<class XX> inline 
Geometry::Interval<X>::Interval(const XX& x)
  : _lower(Numeric::conv_down<X>(x)), _upper(Numeric::conv_up<X>(x)) 
{
}

template<class X> template<class XX> inline 
Geometry::Interval<X>& 
Geometry::Interval<X>::operator=(const XX& x) 
{
  this->_lower=Numeric::conv_down<X>(x); this->_upper=Numeric::conv_up<X>(x); return *this;
}



template<class X> inline 
bool
Geometry::Interval<X>::operator==(const Interval<X>& ivl) {
  return this->_lower==ivl._lower && this->_upper==ivl._upper;
}


template<class X> inline 
bool
Geometry::Interval<X>::operator!=(const Interval<X>& ivl) {
  return !(*this==ivl);
}


template<class X> inline 
const X& 
Geometry::Interval<X>::lower() const 
{ 
  return this->_lower; 
}

template<class X> inline 
const X& 
Geometry::Interval<X>::lower_bound() const 
{ 
  return this->_lower; 
}

template<class X> inline 
const X& 
Geometry::Interval<X>::upper_bound() const 
{ 
  return this->_upper; 
}

template<class X> inline 
const X& 
Geometry::Interval<X>::upper() const 
{ 
  return this->_upper; 
}





template<class X> 
std::ostream& 
Geometry::Interval<X>::write(std::ostream& os) const {
  if(this->_lower > this->_upper) {
    return os << "[1,0]";
  }
  else {
    return os << "[" << this->lower() << "," << this->upper() << "]";
  }
}


template<class X> 
std::istream& 
Geometry::Interval<X>::read(std::istream& is) {
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
    (*this)=Geometry::Interval<X>(l,u);
  } else {
    is.putback(c);
    is >> l;
    (*this)=Geometry::Interval<X>(l);
  }
  return is;
}





template<class X>
X
Geometry::centre(const Interval<X>& ivl) 
{
  return Numeric::div_approx(Numeric::add_approx(ivl.lower(),ivl.upper()),X(2));
}

template<class X>
X
Geometry::radius(const Interval<X>& ivl) 
{
  return Numeric::div_up(Numeric::sub_up(ivl.upper(),ivl.lower()),X(2)); 
}

template<class X>
tribool
Geometry::empty(const Interval<X>& ivl) 
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
Geometry::bounded(const Interval<X>& ivl) 
{
  X inf=Numeric::inf<X>();
  if(ivl.lower() > -inf && ivl.upper() < +inf) {
    return true;
  } else if(ivl.lower() == -inf || ivl.upper() == +inf) {
    return false;
  } else { // May occur if one of the bounds is \a nan or an unbounded interval
    return indeterminate;
  }
}





template<class X1, class X2> inline
tribool 
Geometry::contains(const Interval<X1>& ivl, const X2& x)  
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
Geometry::disjoint(const Interval<X1>& ivl1, const Interval<X2>& ivl2)  
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
Geometry::intersect(const Interval<X1>& ivl1, const Interval<X2>& ivl2)  
{ 
  return !Geometry::disjoint(ivl1,ivl2);
}

template<class X1, class X2> inline
tribool 
Geometry::subset(const Interval<X1>& ivl1, const Interval<X2>& ivl2)  
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
Geometry::superset(const Interval<X1>& ivl1, const Interval<X2>& ivl2)  
{ 
  return Geometry::subset(ivl2,ivl1);
}

template<class X> inline
Geometry::Interval<X>
Geometry::closed_intersection(const Interval<X>& ivl1, const Interval<X>& ivl2)  
{ 
  if(disjoint(ivl1,ivl2)) {
    return Interval<X>(1,0);
  }
  return Interval<X>(Numeric::max(ivl1.lower(),ivl2.lower()),Numeric::min(ivl1.upper(),ivl2.upper()));
}

template<class X> inline
Geometry::Interval<X>
Geometry::open_intersection(const Interval<X>& ivl1, const Interval<X>& ivl2)  
{ 
  if(possibly(disjoint(ivl1,ivl2))) {
    return Interval<X>(1,0);
  }
  return Interval<X>(Numeric::max(ivl1.lower(),ivl2.lower()),Numeric::min(ivl1.upper(),ivl2.upper()));
}

template<class X> inline
Geometry::Interval<X>
Geometry::interval_hull(const Interval<X>& ivl1, const Interval<X>& ivl2)  
{ 
  return Interval<X>(Numeric::min(ivl1.lower(),ivl2.lower()),Numeric::max(ivl1.upper(),ivl2.upper()));
}



template<class X> inline
std::ostream& 
Geometry::operator<<(std::ostream& os, const Interval<X>& x) {
  return x.write(os);
}


template<class X> inline
std::istream& 
Geometry::operator>>(std::istream& is, Interval<X>& x) {
  return x.read(is);
}

} // namespace Ariadne
