/***************************************************************************
 *            interval.inline.h
 *
 *  Copyright 2005-8  Alberto Casagrande, Pieter Collins
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

#include "numeric/expression.h"
#include "numeric/rounding.h"

namespace Ariadne {

  
template<class R> inline
bool empty(const Interval<R>& x) {
  return x.empty();
}

template<class R> inline
R lower(const Interval<R>& x) {
  return x.lower();
}

template<class R> inline
R upper(const Interval<R>& x) {
  return x.upper();
}

template<class R> inline
R midpoint(const Interval<R>& x) {
  return x.midpoint();
}

template<class R> inline
R radius(const Interval<R>& x) { 
  return x.radius();
}

template<class R> inline
R width(const Interval<R>& x) { 
  return x.width();
}


template<class R> inline
R mignitude(const Interval<R>& x) { 
  if(x.lower()>0) { return x.lower(); }
  else if(x.upper()<0) { return -x.upper(); }
  else { return static_cast<R>(0); }
}


template<class X, class Y> inline
bool encloses(const Interval<X>& x, const Interval<Y>& y) { 
  return x.lower()<=y.lower() && x.upper()>=y.upper();
}

template<class X, class Y> inline
bool encloses(const Interval<X>& ivl, const Y& y) { 
  return ivl.encloses(y);
}

template<class X, class Y> inline
bool refines(const Interval<X>& x, const Interval<Y>& y) { 
  return x.refines(y);
}


template<class X, class Y> inline
tribool operator==(const Interval<X>& x, const Interval<Y>& y) { 
  if(x.lower()==y.upper() && x.upper()==y.lower()) { return true; }
  if(x.lower()>y.upper() || x.upper()<y.lower()) { return false; }
  return indeterminate;
}


template<class X, class Y> inline
tribool operator!=(const Interval<X>& x, const Interval<Y>& y) { 
  if(x.lower()>y.upper() || x.upper()<y.lower()) { return true; }
  if(x.lower()==y.upper() && x.upper()==y.lower()) { return false; }
  return indeterminate;
}


template<class X, class Y> inline
tribool operator<(const Interval<X>& x, const Interval<Y>& y) { 
  if(x.upper()<y.lower()) { return true; } 
  if(x.lower()>=y.upper()) { return false; }
  return indeterminate;
}


template<class X, class Y> inline
tribool operator>(const Interval<X>& x, const Interval<Y>& y) { 
  if(x.lower()>y.upper()) { return true; }
  if(x.upper()<=y.lower()) { return false; }
  return indeterminate;
}


template<class X, class Y> inline
tribool operator<=(const Interval<X>& x, const Interval<Y>& y) { 
  if(x.upper()<=y.lower()) { return true; } 
  if(x.lower()>y.upper()) { return false; }
  return indeterminate;
}


template<class X, class Y> inline
tribool operator>=(const Interval<X>& x, const Interval<Y>& y) { 
  if(x.lower()>=y.upper()) { return true; }
  if(x.upper()<y.lower()) { return false; }
  return indeterminate;
}




template<class X, class Y> inline
tribool operator==(const Interval<X>& ivl, const Y& x) { 
  if(ivl.lower()==x && ivl.upper()==x) { return true; }
  if(ivl.lower()>x || ivl.upper()<x) { return false; }
  return indeterminate;
}


template<class X, class Y> inline
tribool operator!=(const Interval<X>& ivl, const Y& x) { 
  if(ivl.lower()>x || ivl.upper()<x) { return true; }
  if(ivl.lower()==x && ivl.upper()==x) { return false; }
  return indeterminate;
}


template<class X, class Y> inline
tribool operator<(const Interval<X>& ivl, const Y& x) { 
  if(ivl.upper()<x) { return true; }
  if(ivl.lower()>=x) { return false; }
  return indeterminate;
}


template<class X, class Y> inline
tribool operator>(const Interval<X>& ivl, const Y& x) { 
  if(ivl.lower()>x) { return true; }
  if(ivl.upper()<=x) { return false; }
  return indeterminate;
}


template<class X, class Y> inline
tribool operator<=(const Interval<X>& ivl, const Y& x) { 
  if(ivl.upper()<=x) { return true; }
  if(ivl.lower()>x) { return false; }
  return indeterminate;
}


template<class X, class Y> inline
tribool operator>=(const Interval<X>& ivl, const Y& x) { 
  if(ivl.lower()>=x) { return true; }
  if(ivl.upper()<x) { return false; }
  return indeterminate;
}


template<class X, class Y> inline
tribool operator==(const X& x, const Interval<Y>& ivl) { 
  return ivl==x; 
}



template<class X, class Y> inline
tribool operator!=(const X& x, const Interval<Y>& ivl) {
  return ivl!=x;
}


template<class X, class Y> inline
tribool operator<(const X& x, const Interval<Y>& ivl) {
  return ivl>x;
}


template<class X, class Y> inline
tribool operator>(const X& x, const Interval<Y>& ivl) {
  return ivl<x;
}


template<class X, class Y> inline
tribool operator<=(const X& x, const Interval<Y>& ivl) {
  return ivl>=x;
}


template<class X, class Y> inline
tribool operator>=(const X& x, const Interval<Y>& ivl) {
  return ivl<=x;
}





template<class R> inline
void min_(Interval<R>& r, const Interval<R>& x1, const Interval<R>& x2) {
  r._lower=min(x1.lower(),x2.lower()); r._upper=min(x1.upper(),x2.upper());
}

template<class R> inline
void max_(Interval<R>& r, const Interval<R>& x1, const Interval<R>& x2) {
  r._lower=max(x1.lower(),x2.lower()); r._upper=max(x1.upper(),x2.upper());
}


template<class R, class X> inline
void abs_(Interval<R>& r, const Interval<X>& x) {
  if(x.lower()>=0) { r=x; return; } 
  if(x.upper() < 0) { r=-x; return; } 
  R nxl; neg_(nxl,x.lower()); 
  const R& xu=x.upper();
  set_(r._lower,0); max_(r._upper,nxl,xu);
}

template<class R, class X> inline
void abs_(Interval<R>& r, const X& x) {
  if(x>=0) { r=x; } else { r=-x; }
}




template<class R> inline
Interval<R> abs(const Interval<R>& x) {
  Interval<R> r; abs_(r,x); return r; }

template<class R> inline
R mag(const Interval<R>& x) {
  return max(abs(x._lower),abs(x._upper));
}

template<class R> inline
R mig(const Interval<R>& x) {
  if(x._lower>0) {
    return x._lower;
  } else if(x._upper<0) {
    return -x._upper;
  } else { 
    return 0;
  }
}


template<class R> inline
bool equal(const Interval<R>& x1, const Interval<R>& x2) {
  return (x1.lower()==x2.lower() && x1.upper()==x2.upper());
}

template<class R> inline
bool disjoint(const Interval<R>& x1, const Interval<R>& x2) {
  return (x1.upper()<x2.lower() || x1.lower()>x2.upper());
}

template<class R> inline
bool overlap(const Interval<R>& x1, const Interval<R>& x2) {
  return (x1.upper()>x2.lower() && x1.lower()<x2.upper());
}

template<class R> inline
bool subset(const Interval<R>& x1, const Interval<R>& x2) {
  return (x1.lower()>=x2.lower() && x1.upper()<=x2.upper());
}

template<class R> inline
bool inside(const Interval<R>& x1, const Interval<R>& x2) {
  return (x1.lower()>x2.lower() && x1.upper()<x2.upper());
}

template<class R> inline
Interval<R> intersection(const Interval<R>& x1, const Interval<R>& x2) {
  return Interval<R>(max(x1.lower(),x2.lower()),
                     min(x1.upper(),x2.upper()));
}

template<class R> inline
Interval<R> hull(const Interval<R>& x, const Interval<R>& y) {
  if(x.empty()) { return y; }
  if(y.empty()) { return x; }
  return Interval<R>(min(x.lower(),y.lower()),
                     max(x.upper(),y.upper()));
}




} // namespace Ariadne






namespace TBLAS {

template<class R>
int
iamax_ (const int N, const Ariadne::Interval<R> *X, const int incX)
{
#ifdef DEBUG
  std::cerr << "TBLAS::iamax_(const int N, const interval<real> *X, const int incX)\n";
#endif
  
  R mx = 0;
  int ix = 0;
  int i;
  int result = 0;
  
  if (incX <= 0) {
    return 0;
  }
  
  for (i = 0; i < N; i++) {
    R av=Ariadne::min_(
                                   Ariadne::abs_(X[ix].lower()),
                                   Ariadne::abs_(X[ix].upper()));
    if (av > mx) {
      mx = av;
      result = i;
    }
    ix += incX;
  }
  
  return result;
}

} // namespace TBLAS
