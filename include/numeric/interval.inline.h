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

#include "numeric/expression.h"
#include "numeric/rounding.h"
 
namespace Ariadne {
namespace Numeric {
  

using Base::tribool;
using Base::indeterminate;

using std::min;
using std::max;


template<class R> inline 
Interval<R>::Interval()
  : _lower(0), _upper(0) { }

template<class R> inline 
Interval<R>::Interval(const R& x)
  : _lower(x), _upper(x) { }

template<class R> inline 
Interval<R>::Interval(const R& l, const R& u)
  : _lower(l), _upper(u) { }

template<class R> inline 
Interval<R>::Interval(const Interval<R>& ivl)
  : _lower(ivl._lower), _upper(ivl._upper) { }

template<class R> template<class RL,class RU> inline 
Interval<R>::Interval(const RL& l, const RU& u)
  : _lower(l,round_down), _upper(u,round_up) { }

template<class R> template<class RX> inline 
Interval<R>::Interval(const Interval<RX>& ivl)
  : _lower(ivl.lower(),round_down), _upper(ivl.upper(),round_up) { }

template<class R> template<class RX> inline 
Interval<R>::Interval(const RX& x)
  : _lower(x,round_down), _upper(x,round_up) { }


template<class R> inline 
Interval<R>& Interval<R>::operator=(const R& x) {
  this->_lower=x; 
  this->_upper=x; 
  return *this;
}

template<class R> inline 
Interval<R>& Interval<R>::operator=(const Interval<R>& ivl) {
  this->_lower=ivl._lower; 
  this->_upper=ivl._upper; 
  return *this;
}

template<class R> template<class RX> inline 
Interval<R>& Interval<R>::operator=(const RX& x) {
  set_(this->_lower,x,round_down); 
  set_(this->_upper,x,round_up); 
  return *this;
}

template<class R> template<class RX> inline 
Interval<R>& Interval<R>::operator=(const Interval<RX>& ivl) {
  set_(this->_lower,ivl.lower(),round_down); 
  set_(this->_upper,ivl.upper(),round_up); 
  return *this;
}



template<class R> inline 
const R& Interval<R>::lower() const { 
  return this->_lower; 
}

template<class R> inline 
const R& Interval<R>::upper() const { 
  return this->_upper; 
}


template<class R> inline 
R Interval<R>::midpoint() const { 
  R r; 
  med_(r,this->_lower,this->_upper,round_approx); 
  return r;
}

template<class R> inline 
R Interval<R>::radius() const { 
  R r; 
  rad_(r,this->_lower,this->_upper,round_up); 
  return r;
}

template<class R> inline
R Interval<R>::width() const { 
  R r;
  sub_(r,this->_upper,this->_lower,round_up); 
  return r;
}


template<class R> inline 
bool Interval<R>::empty() const { 
  return this->lower()>this->upper(); 
}

template<class R> inline 
bool Interval<R>::singleton() const { 
  return this->lower()==this->upper();
}

template<class R> template<class RX> inline 
bool Interval<R>::encloses(const RX& x) const { 
  return this->lower()<=x && x<=this->upper();
}

template<class R> template<class RX> inline 
bool Interval<R>::refines(const Interval<RX>& ivl) const { 
  return ivl.lower()<=this->lower() && this->upper()<=ivl.upper(); 
}







template<class R> inline 
IntervalReference<R>::IntervalReference(R& l, R& u)
  : _lower(&l), _upper(&u) { }

template<class R> inline 
IntervalReference<R>::IntervalReference(Interval<R>& ivl)
  : _lower(&const_cast<R&>(ivl.lower())), _upper(&const_cast<R&>(ivl.upper())) { }

template<class R> inline
void IntervalReference<R>::operator=(const Interval<R>& ivl) {
  *_lower=ivl.lower(); *_upper=ivl.upper(); }

template<class R> inline
IntervalReference<R>::operator Interval<R> () const {
  return Interval<R>(*_lower,*_upper); }

template<class R> inline 
R IntervalReference<R>::lower() const { 
  return *_lower; }

template<class R> inline 
R IntervalReference<R>::upper() const { 
  return *_upper; }

template<class R> inline 
R IntervalReference<R>::midpoint() const { 
  R r; med_(r,*_lower,*_upper,round_approx); return r; 
}

template<class R> inline 
R IntervalReference<R>::radius() const { 
  R r; rad_(r,*_lower,*_upper,round_approx); return r; 
}






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



template<class R, class X> inline
void pos_(Interval<R>& r, const Interval<X>& x) {
  pos_(r._lower,x._lower); 
  pos_(r._upper,x._upper);
}

template<class R, class X> inline
void pos_(Interval<R>& r, const X& x) {
  pos_(r._lower,x); 
  pos_(r._upper,x);
}


template<class R, class X> inline
void neg_(Interval<R>& r, const Interval<X>& x) {
  neg_(r._lower,x._upper);
  neg_(r._upper,x._lower);
}

template<class R, class X> inline
void neg_(Interval<R>& r, const X& x) {
  neg_(r._lower,x);
  neg_(r._upper,x);
}



template<class R, class X, class Y> inline
void add_(Interval<R>& r, const Interval<X>& x, const Interval<Y>& y) {
  add_(r._lower,x._lower,y._lower,round_down); 
  add_(r._upper,x._upper,y._upper,round_up);
}

template<class R, class X, class Y> inline
void add_(Interval<R>& r, const Interval<X>& x, const Y& y) {
  add_(r._lower,x._lower,y,round_down); 
  add_(r._upper,x._upper,y,round_up);
}

template<class R, class X, class Y> inline
void add_(Interval<R>& r, const X& x, const Interval<Y>& y) {
  add_(r._lower,x,y._lower,round_down); 
  add_(r._upper,x,y._upper,round_up);
}

template<class R, class X, class Y> inline
void add_(Interval<R>& r, const X& x, const Y& y) {
  add_(r._lower,x,y,round_down); 
  add_(r._upper,x,y,round_up);
}



template<class R, class X, class Y> inline
void sub_noalias_(Interval<R>& r, const Interval<X>& x, const Interval<Y>& y) {
  sub_(r._lower,x._lower,y._upper,round_down); 
  sub_(r._upper,x._upper,y._lower,round_up);
}

template<class R, class X, class Y> inline
void sub_noalias_(Interval<R>& r, const X& x, const Interval<Y>& y) {
  sub_(r._lower,x,y._upper,round_down); 
  sub_(r._upper,x,y._lower,round_up);
}

template<class R, class X, class Y> inline
void sub_(Interval<R>& r, const Interval<X>& x, const Interval<Y>& y) {
  if(&r==&y) { Interval<R> t; sub_noalias_(t,x,y); r=t; }
  else { sub_noalias_(r,x,y); }
}

template<class R, class X, class Y> inline
void sub_(Interval<R>& r, const X& x, const Interval<Y>&  y) {
  if(&r==&y) { Interval<R> t; sub_noalias_(t,x,y); r=t; }
  else { sub_noalias_(r,x,y); }
}

template<class R, class X, class Y> inline
void sub_(Interval<R>& r, const Interval<X>& x, const Y& y) {
  sub_(r._lower,x._lower,y,round_down); 
  sub_(r._upper,x._upper,y,round_up);
}

template<class R, class X, class Y> inline
void sub_(Interval<R>& r, const X& x, const Y& y) {
  sub_(r._lower,x,y,round_down); 
  sub_(r._upper,x,y,round_up);
}








template<class R, class X, class Y> 
void mul_noalias_(Interval<R>& r, const Interval<X>& x, const Interval<Y>& y) {
  R& rl = r._lower;
  R& ru = r._upper;
  const X& xl = x.lower();
  const X& xu = x.upper();
  const Y& yl = y.lower();
  const Y& yu = y.upper();
  RoundDown d;
  RoundUp u;
  
  //FIXME: Case where &r==&x or &r==&y;
  if (xl>=0) {
    if (yl>=0) {
      mul_(rl,xl,yl,d); mul_(ru,xu,yu,u);
    } else if(yu<=0) {
      mul_(rl,xu,yl,d); mul_(ru,xl,yu,u);
    } else {
      mul_(rl,xu,yl,d); mul_(ru,xu,yu,u);
    }
  } else if (xu<=0) {
    if (yl>=0) {
      mul_(rl,xl,yu,d); mul_(ru,xu,yl,u);
    } else if(yu<=0) {
      mul_(rl,xu,yu,d); mul_(ru,xl,yl,u); 
    } else {
      mul_(rl,xl,yu,d); mul_(ru,xl,yl,u);
    }
  } else {
    if (yl>=0) {
      mul_(rl,xl,yu,d); mul_(ru,xu,yu,u);
    } else if(yu<=0) {
      mul_(rl,xu,yl,d); mul_(ru,xl,yl,u);
    } else {
      R t1; R t2; 
      mul_(t1,xl,yu,d); mul_(t2,xu,yl,d); min_(rl,t1,t2);
      mul_(t1,xl,yl,u); mul_(t2,xu,yu,u); max_(ru,t1,t2);
    }
  }
}    



template<class R, class X, class Y> 
void mul_noalias_(Interval<R>& r, const Interval<X>& x, const Y& y) {
  typedef Interval<R> I;
  R& rl = r._lower;
  R& ru = r._upper;
  const X& xl = x.lower();
  const X& xu = x.upper();
  RoundDown d;
  RoundUp u;
  
  if (y>=0) {
    mul_(rl,xl,y,d); mul_(ru,xu,y,u);
  } else {
    mul_(rl,xu,y,d); mul_(ru,xl,y,u);
  }
}

template<class R, class X, class Y> 
void mul_(Interval<R>& r, const Interval<X>& x, const Interval<Y>& y) {
  if(&r==&x || &r==&y) {
    Interval<R> t; mul_noalias_(t,x,y); r=t;
  } else {
    mul_noalias_(r,x,y);
  }
}

template<class R, class X, class Y> 
void mul_(Interval<R>& r, const Interval<X>& x, const Y& y) {
  if(&r==&x) {
    Interval<R> t; mul_noalias_(t,x,y); r=t;
  } else {
    mul_noalias_(r,x,y);
  }
}

template<class R, class X, class Y> inline
void mul_(Interval<R>& r, const X& x, const Interval<Y>& y) {
  mul_(r,y,x);
}

template<class R, class X, class Y> inline
void mul_(Interval<R>& r, const X& x, const Y& y) {
  mul_(r._lower,x,y,round_down); mul_(r._upper,x,y,round_up);
}


template<class R, class X, class Y> 
void div_noalias_(Interval<R>& r, const Interval<X>& x, const Interval<Y>& y) {
  typedef Interval<R> I;
  R& rl = r._lower;
  R& ru = r._upper;
  const X& xl = x.lower();
  const X& xu = x.upper();
  const Y& yl = y.lower();
  const Y& yu = y.upper();
  RoundDown d;
  RoundUp u;
  
  if (yl>0) {
    if (xl>=0) {
      div_(rl,xl,yu,d); div_(ru,xu,yl,u);
    } else if(xu<=0) {
      div_(rl,xl,yl,d); div_(ru,xu,yu,u);
    } else {
      div_(rl,xl,yl,d); div_(ru,xu,yl,u);
    }
  } else if (yu<0) {
    if (xl>=0) {
      div_(rl,xu,yu,d); div_(ru,xl,yl,u);
    } else if(xu<=0) {
      div_(rl,xu,yl,d); div_(ru,xl,yu,u);
    } else {
      div_(rl,xu,yu,d); div_(ru,xl,yu,u);
    }
  } else {
    inf_(ru); neg_(rl,ru);
  }
}

template<class R, class X, class Y> 
void div_(Interval<R>& r, const Interval<X>& x, const Interval<Y>& y) {
  if(&r==&x || &r==&y) { Interval<R> t; div_noalias_(t,x,y); r=t; }
  else { div_noalias_(r,x,y); }
}

template<class R, class X, class Y> inline
void div_(Interval<R>& r, const Interval<X>& x, const Y& y) {
  div_(r,x,Interval<R>(y));
}

template<class R, class X, class Y> inline
void div_(Interval<R>& r, const X& x, const Interval<Y>& y) {
  div_(r,Interval<R>(x),y);
}

template<class R, class X, class Y> inline
void div_(Interval<R>& r, const X& x, const Y& y) {
  div_(r._lower,x,y,round_down);
  div_(r._upper,x,y,round_up);
}




template<class R> inline
void min_(Interval<R>& r, const Interval<R>& x1, const Interval<R>& x2) {
  r._lower=min(x1.lower(),x2.lower()); r._upper=min(x1.upper(),x2.upper());
}

template<class R> inline
void max_(Interval<R>& r, const Interval<R>& x1, const Interval<R>& x2) {
  r._lower=max(x1.lower(),x2.lower()); r._upper=max(x1.upper(),x2.upper());
}


template<class R> inline
void abs_(Interval<R>& r, const Interval<R>& x) {
  if(x.lower()>=0) { r=x; return; } 
  if(x.upper() < 0) { r=-x; return; } 
  R nxl; neg_(nxl,x.lower()); 
  const R& xu=x.upper();
  set_(r._lower,0); max_(r._upper,nxl,xu);
}

template<class R, class X> inline
void abs_(Interval<R>& r, const Interval<X>& x) {
  if(x.lower()>=0) { r=x; return; } 
  if(x.upper() < 0) { r=-x; return; } 
  R nxl; neg_(nxl,x.lower(),round_up); 
  R xu; pos_(xu,x.upper(),round_up);
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
void pow_(Interval<R>& r, const Interval<R>& x, const uint& n) {
  Interval<R> y=x;
  r=1;
  for(uint i=0; i!=n; ++i) {
    r=r*y;
  }
}

template<class R> inline
void pow_(Interval<R>& r, const Interval<R>& x, const int& n) {
  if(n>=0) { pow_(r,x,uint(n)); }
  else { pow_(r,Interval<R>(1/x),uint(-n)); }
}

template<class R> inline
void pow_(Interval<R>& r, const Interval<R>& x, const Integer& n) {
  Integer m=abs(n);
  Interval<R> a = (n>=0) ? x : 1/x;
  r=1;
  for(Integer i=0; i!=m; ++i) {
    r*=a;
  }
}








// Algebraic and transcendental functions 
template<class R> 
inline void sqrt_(Interval<R>& r, const R& x) {
  sqrt_(r._lower,x,round_down); 
  sqrt_(r._upper,x,round_up); } 

template<class R> 
inline void sqrt_(Interval<R>& r, const Interval<R>& x) {
  sqrt_(r._lower,x._lower,round_down); 
  sqrt_(r._upper,x._upper,round_up); } 


template<class R> inline 
void hypot_(Interval<R>& r, const R& x, const R& y) {
  hypot_(r._lower,x,y,round_down); 
  hypot_(r._upper,x,y,round_up); } 

template<class R, class X, class Y> inline 
void hypot_(Interval<R> r, const Interval<X>& x, const Interval<Y>& y) {
  hypot_(r._lower,x._lower,y._lower,round_down); 
  hypot_(r._upper,x._lower,y._lower,round_up); } 


template<class R> inline
void exp_(Interval<R>& r, const R& x) {
  exp_(r._lower,x,round_down); 
  exp_(r._upper,x,round_up); } 

template<class R, class X> inline 
void exp_(Interval<R>& r, const Interval<X>& x) {
  exp_(r._lower,x._lower,round_down); 
  exp_(r._upper,x._lower,round_up); } 


template<class R> inline
void log_(Interval<R>& r, const R& x) {
  log_(r._lower,x,round_down); 
  log_(r._upper,x,round_up); } 

template<class R, class X> inline 
void log_(Interval<R>& r, const Interval<X>& x) {
  log_(r._lower,x._lower,round_down); 
  log_(r._upper,x._lower,round_up); } 


// Trigonometric functions for non-interval classes.
template<class R, class X> inline
void sin_(Interval<R>& r, const X& x) {
  sin_(r,Interval<R>(x)); }

template<class R, class X> inline
void cos_(Interval<R>& r, const X& x) {
  cos_(r,Interval<R>(x)); }

template<class R, class X> inline
void tan_(Interval<R>& r, const X& x) {
  tan_(r,Interval<R>(x)); }

template<class R, class X> inline
void asin_(Interval<R>& r, const X& x) {
  asin_(r,Interval<R>(x)); }

template<class R, class X> inline
void acos_(Interval<R>& r, const X& x) {
  acos_(r,Interval<R>(x)); }

template<class R, class X> inline
void atan_(Interval<R>& r, const X& x) {
  atan_(r,Interval<R>(x)); }



// Directly evaluated transcendental functions
template<class R> 
inline Interval<R> sqrt(const Interval<R>& x) { 
  Interval<R> r; sqrt_(r,x); return r; }
template<class R> 
inline Interval<R> exp(const Interval<R>& x) { 
  Interval<R> r; exp_(r,x); return r; }
template<class R> 
inline Interval<R> log(const Interval<R>& x) { 
  Interval<R> r; log_(r,x); return r; }
template<class R> 
inline Interval<R> sin(const Interval<R>& x) { 
  Interval<R> r; sin_(r,x); return r; }
template<class R> 
inline Interval<R> cos(const Interval<R>& x) { 
  Interval<R> r; cos_(r,x); return r; }
template<class R> 
inline Interval<R> tan(const Interval<R>& x) { 
  Interval<R> r; tan_(r,x); return r; }
template<class R> 
inline Interval<R> asin(const Interval<R>& x) { 
  Interval<R> r; asin_(r,x); return r; }
template<class R> 
inline Interval<R> acos(const Interval<R>& x) { 
  Interval<R> r; acos_(r,x); return r; }
template<class R> 
inline Interval<R> atan(const Interval<R>& x) { 
  Interval<R> r; atan_(r,x); return r; }




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








template<class R> template<class E> inline
Interval<R>::Interval(const Expression<E>& e) { e.assign_to(*this); }

template<class R> template<class E> inline
Interval<R>& Interval<R>::operator=(const Expression<E>& e) { 
  e.assign_to(*this); return *this; }







} // namespace Numeric
} // namespace Ariadne






namespace TBLAS {

template<class real>
int
iamax_ (const int N, const Ariadne::Numeric::Interval<real> *X, const int incX)
{
#ifdef DEBUG
  std::cerr << "TBLAS::iamax_(const int N, const interval<real> *X, const int incX)\n";
#endif
  
  real mx = 0;
  int ix = 0;
  int i;
  int result = 0;
  
  if (incX <= 0) {
    return 0;
  }
  
  for (i = 0; i < N; i++) {
    real av=Ariadne::Numeric::min_(
                                   Ariadne::Numeric::abs_(X[ix].lower()),
                                   Ariadne::Numeric::abs_(X[ix].upper()));
    if (av > mx) {
      mx = av;
      result = i;
    }
    ix += incX;
  }
  
  return result;
}

} // namespace TBLAS
