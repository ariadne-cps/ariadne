/***************************************************************************
 *            interval-float.inline.h
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

  

using std::min;
using std::max;


template<class T> inline 
Interval< Float<T> >::Interval()
  : _lower(0), _upper(0) { }

template<class T> inline 
Interval< Float<T> >::Interval(const R& x)
  : _lower(x), _upper(x) { }

template<class T> inline 
Interval< Float<T> >::Interval(const R& l, const R& u)
  : _lower(l), _upper(u) { }

template<class T> inline 
Interval< Float<T> >::Interval(const Interval< Float<T> >& ivl)
  : _lower(ivl._lower), _upper(ivl._upper) { }

template<class T> template<class RL,class RU> inline 
Interval< Float<T> >::Interval(const RL& l, const RU& u)
  : _lower(l,round_down), _upper(u,round_up) { }

template<class T> template<class RX> inline 
Interval< Float<T> >::Interval(const Interval<RX>& ivl)
  : _lower(ivl.lower(),round_down), _upper(ivl.upper(),round_up) { }

template<class T> template<class RX> inline 
Interval< Float<T> >::Interval(const RX& x)
  : _lower(x,round_down), _upper(x,round_up) { }


template<class T> inline 
Interval< Float<T> >& Interval< Float<T> >::operator=(const R& x) {
  this->_lower=x; 
  this->_upper=x; 
  return *this;
}

template<class T> inline 
Interval< Float<T> >& Interval< Float<T> >::operator=(const Interval< Float<T> >& ivl) {
  this->_lower=ivl._lower; 
  this->_upper=ivl._upper; 
  return *this;
}

template<class T> template<class RX> inline 
Interval< Float<T> >& Interval< Float<T> >::operator=(const RX& x) {
  set_(this->_lower,x,round_down); 
  set_(this->_upper,x,round_up); 
  return *this;
}

template<class T> template<class RX> inline 
Interval< Float<T> >& Interval< Float<T> >::operator=(const Interval<RX>& ivl) {
  set_(this->_lower,ivl.lower(),round_down); 
  set_(this->_upper,ivl.upper(),round_up); 
  return *this;
}



template<class T> inline 
const Float<T>& Interval< Float<T> >::lower() const { 
  return this->_lower; 
}

template<class T> inline 
const Float<T>& Interval< Float<T> >::upper() const { 
  return this->_upper; 
}


template<class T> inline 
Float<T> Interval< Float<T> >::midpoint() const { 
  Float<T> r; 
  med_(r,this->_lower,this->_upper,round_approx); 
  return r;
}

template<class T> inline 
Float<T> Interval< Float<T> >::radius() const { 
  Float<T> r; 
  rad_(r,this->_lower,this->_upper,round_up); 
  return r;
}

template<class T> inline
Float<T> Interval< Float<T> >::width() const { 
  Float<T> r;
  sub_(r,this->_upper,this->_lower,round_up); 
  return r;
}

template<class T> inline
Interval< Float<T> >::Interval(const std::string& str) {
  std::stringstream ss(str); ss>>*this; }


template<class T> inline 
bool Interval< Float<T> >::empty() const { 
  return this->lower()>this->upper(); 
}

template<class T> inline 
bool Interval< Float<T> >::singleton() const { 
  return this->lower()==this->upper();
}

template<class T> template<class RX> inline 
bool Interval< Float<T> >::encloses(const RX& x) const { 
  return this->lower()<=x && x<=this->upper();
}

template<class T> template<class RX> inline 
bool Interval< Float<T> >::refines(const Interval<RX>& ivl) const { 
  return ivl.lower()<=this->lower() && this->upper()<=ivl.upper(); 
}







template<class T, class X> inline
void pos_(Interval< Float<T> >& r, const Interval<X>& x) {
  pos_(r._lower,x._lower,round_down); 
  pos_(r._upper,x._upper,round_up);
}

template<class T, class X> inline
void pos_(Interval< Float<T> >& r, const X& x) {
  pos_(r._lower,x,round_down); 
  pos_(r._upper,x,round_up);
}


template<class T, class X> inline
void neg_(Interval< Float<T> >& r, const Interval<X>& x) {
  neg_(r._lower,x._upper,round_down);
  neg_(r._upper,x._lower,round_up);
}

template<class T, class X> inline
void neg_(Interval< Float<T> >& r, const X& x) {
  neg_(r._lower,x,round_down);
  neg_(r._upper,x,round_up);
}



template<class T, class X, class Y> inline
void add_(Interval< Float<T> >& r, const Interval<X>& x, const Interval<Y>& y) {
  add_(r._lower,x._lower,y._lower,round_down); 
  add_(r._upper,x._upper,y._upper,round_up);
}

template<class T, class X, class Y> inline
void add_(Interval< Float<T> >& r, const Interval<X>& x, const Y& y) {
  add_(r._lower,x._lower,y,round_down); 
  add_(r._upper,x._upper,y,round_up);
}

template<class T, class X, class Y> inline
void add_(Interval< Float<T> >& r, const X& x, const Interval<Y>& y) {
  add_(r._lower,x,y._lower,round_down); 
  add_(r._upper,x,y._upper,round_up);
}

template<class T, class X, class Y> inline
void add_(Interval< Float<T> >& r, const X& x, const Y& y) {
  add_(r._lower,x,y,round_down); 
  add_(r._upper,x,y,round_up);
}



template<class T, class X, class Y> inline
void sub_noalias_(Interval< Float<T> >& r, const Interval<X>& x, const Interval<Y>& y) {
  sub_(r._lower,x._lower,y._upper,round_down); 
  sub_(r._upper,x._upper,y._lower,round_up);
}

template<class T, class X, class Y> inline
void sub_noalias_(Interval< Float<T> >& r, const X& x, const Interval<Y>& y) {
  sub_(r._lower,x,y._upper,round_down); 
  sub_(r._upper,x,y._lower,round_up);
}

template<class T, class X, class Y> inline
void sub_(Interval< Float<T> >& r, const Interval<X>& x, const Interval<Y>& y) {
  if(&r==&y) { Interval< Float<T> > t; sub_noalias_(t,x,y); r=t; }
  else { sub_noalias_(r,x,y); }
}

template<class T, class X, class Y> inline
void sub_(Interval< Float<T> >& r, const X& x, const Interval<Y>&  y) {
  if(&r==&y) { Interval< Float<T> > t; sub_noalias_(t,x,y); r=t; }
  else { sub_noalias_(r,x,y); }
}

template<class T, class X, class Y> inline
void sub_(Interval< Float<T> >& r, const Interval<X>& x, const Y& y) {
  sub_(r._lower,x._lower,y,round_down); 
  sub_(r._upper,x._upper,y,round_up);
}

template<class T, class X, class Y> inline
void sub_(Interval< Float<T> >& r, const X& x, const Y& y) {
  sub_(r._lower,x,y,round_down); 
  sub_(r._upper,x,y,round_up);
}








template<class T, class X, class Y> 
void mul_noalias_(Interval< Float<T> >& r, const Interval<X>& x, const Interval<Y>& y) {
  Float<T>& rl = r._lower;
  Float<T>& ru = r._upper;
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
      Float<T> t1; Float<T> t2; 
      mul_(t1,xl,yu,d); mul_(t2,xu,yl,d); min_(rl,t1,t2);
      mul_(t1,xl,yl,u); mul_(t2,xu,yu,u); max_(ru,t1,t2);
    }
  }
}    



template<class T, class X, class Y> 
void mul_noalias_(Interval< Float<T> >& r, const Interval<X>& x, const Y& y) {
  typedef Interval< Float<T> > I;
  Float<T>& rl = r._lower;
  Float<T>& ru = r._upper;
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

template<class T, class X, class Y> 
void mul_(Interval< Float<T> >& r, const Interval<X>& x, const Interval<Y>& y) {
  if(&r==&x || &r==&y) {
    Interval< Float<T> > t; mul_noalias_(t,x,y); r=t;
  } else {
    mul_noalias_(r,x,y);
  }
}

template<class T, class X, class Y> 
void mul_(Interval< Float<T> >& r, const Interval<X>& x, const Y& y) {
  if(&r==&x) {
    Interval< Float<T> > t; mul_noalias_(t,x,y); r=t;
  } else {
    mul_noalias_(r,x,y);
  }
}

template<class T, class X, class Y> inline
void mul_(Interval< Float<T> >& r, const X& x, const Interval<Y>& y) {
  mul_(r,y,x);
}

template<class T, class X, class Y> inline
void mul_(Interval< Float<T> >& r, const X& x, const Y& y) {
  mul_(r._lower,x,y,round_down); mul_(r._upper,x,y,round_up);
}


template<class T, class X, class Y> 
void div_noalias_(Interval< Float<T> >& r, const Interval<X>& x, const Interval<Y>& y) {
  typedef Interval< Float<T> > I;
  Float<T>& rl = r._lower;
  Float<T>& ru = r._upper;
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

template<class T, class X, class Y> 
void div_(Interval< Float<T> >& r, const Interval<X>& x, const Interval<Y>& y) {
  if(&r==&x || &r==&y) { Interval< Float<T> > t; div_noalias_(t,x,y); r=t; }
  else { div_noalias_(r,x,y); }
}

template<class T, class X, class Y> inline
void div_(Interval< Float<T> >& r, const Interval<X>& x, const Y& y) {
  div_(r,x,Interval< Float<T> >(y));
}

template<class T, class X, class Y> inline
void div_(Interval< Float<T> >& r, const X& x, const Interval<Y>& y) {
  div_(r,Interval< Float<T> >(x),y);
}

template<class T, class X, class Y> inline
void div_(Interval< Float<T> >& r, const X& x, const Y& y) {
  div_(r._lower,x,y,round_down);
  div_(r._upper,x,y,round_up);
}


template<class T, class X> inline
void abs_(Interval< Float<T> >& r, const Interval<X>& x) {
  if(x.lower()>=0) { r=x; return; } 
  if(x.upper() < 0) { r=-x; return; } 
  Float<T> nxl; neg_(nxl,x.lower(),round_up); 
  Float<T> xu; pos_(xu,x.upper(),round_up);
  set_(r._lower,0); max_(r._upper,nxl,xu);
}




template<class T> inline
void pow_(Interval< Float<T> >& r, const Interval< Float<T> >& x, const uint& n) {
  if(n%2 || x._lower>=0) {
    pow_(r._lower,x._lower,n,round_down);
    pow_(r._upper,x._upper,n,round_up);
  } else if(x._upper<=0) {
    Interval< Float<T> > t(x);
    pow_(r._lower,t._upper,n,round_down);
    pow_(r._upper,t._lower,n,round_up);
  } else {
    Float<T> rl, ru;
    pow_(rl,x._lower,n,round_up);
    pow_(ru,x._upper,n,round_up);
    r._lower=0;
    r._upper=max(rl,ru);
  }
}

template<class T> inline
void pow_(Interval< Float<T> >& r, const Interval< Float<T> >& x, const int& n) {
  if(n>=0) { pow_(r,x,uint(n)); }
  else { pow_(r,Interval< Float<T> >(1/x),uint(-n)); }
}

template<class T> inline
void pow_(Interval< Float<T> >& r, const Interval< Float<T> >& x, const Integer& n) {
  //TODO: Check this code for bugs; improve efficiency
  Integer m=abs(n);
  Interval< Float<T> > a = (n>=0) ? x : 1/x;
  r=1;
  for(Integer i=0; i!=m; ++i) {
    r*=a;
  }
  if(m%2==0) {
    if(r._lower<0) {
      r._lower=0;
    }
  }
}








// Algebraic and transcendental functions 
template<class T> 
inline void sqrt_(Interval< Float<T> >& r, const Float<T>& x) {
  sqrt_(r._lower,x,round_down); 
  sqrt_(r._upper,x,round_up); } 

template<class T> 
inline void sqrt_(Interval< Float<T> >& r, const Interval< Float<T> >& x) {
  sqrt_(r._lower,x._lower,round_down); 
  sqrt_(r._upper,x._upper,round_up); } 


template<class T> inline 
void hypot_(Interval< Float<T> >& r, const Float<T>& x, const Float<T>& y) {
  hypot_(r._lower,x,y,round_down); 
  hypot_(r._upper,x,y,round_up); } 

template<class T> inline 
void hypot_(Interval< Float<T> > r, const Interval< Float<T> >& x, const Interval< Float<T> >& y) {
  hypot_(r._lower,x._lower,y._lower,round_down); 
  hypot_(r._upper,x._upper,y._upper,round_up); } 


template<class T> inline
void exp_(Interval< Float<T> >& r, const Float<T>& x) {
  exp_(r._lower,x,round_down); 
  exp_(r._upper,x,round_up); } 

template<class T> inline 
void exp_(Interval< Float<T> >& r, const Interval< Float<T> >& x) {
  exp_(r._lower,x._lower,round_down); 
  exp_(r._upper,x._upper,round_up); } 


template<class T> inline
void log_(Interval< Float<T> >& r, const Float<T>& x) {
  log_(r._lower,x,round_down); 
  log_(r._upper,x,round_up); } 

template<class T> inline 
void log_(Interval< Float<T> >& r, const Interval< Float<T> >& x) {
  log_(r._lower,x._lower,round_down); 
  log_(r._upper,x._upper,round_up); } 


// Trigonometric functions for non-interval classes.
template<class T> inline
void sin_(Interval< Float<T> >& r, const Float<T>& x) {
  sin_(r,Interval< Float<T> >(x)); }

template<class T> inline
void cos_(Interval< Float<T> >& r, const Float<T>& x) {
  cos_(r,Interval< Float<T> >(x)); }

template<class T> inline
void tan_(Interval< Float<T> >& r, const Float<T>& x) {
  tan_(r,Interval< Float<T> >(x)); }

template<class T> inline
void asin_(Interval< Float<T> >& r, const Float<T>& x) {
  asin_(r,Interval< Float<T> >(x)); }

template<class T> inline
void acos_(Interval< Float<T> >& r, const Float<T>& x) {
  acos_(r,Interval< Float<T> >(x)); }

template<class T> inline
void atan_(Interval< Float<T> >& r, const Float<T>& x) {
  atan_(r,Interval< Float<T> >(x)); }



// Directly evaluated transcendental functions
template<class T> 
inline Interval< Float<T> > sqrt(const Interval< Float<T> >& x) { 
  Interval< Float<T> > r; sqrt_(r,x); return r; }
template<class T> 
inline Interval< Float<T> > exp(const Interval< Float<T> >& x) { 
  Interval< Float<T> > r; exp_(r,x); return r; }
template<class T> 
inline Interval< Float<T> > log(const Interval< Float<T> >& x) { 
  Interval< Float<T> > r; log_(r,x); return r; }
template<class T> 
inline Interval< Float<T> > sin(const Interval< Float<T> >& x) { 
  Interval< Float<T> > r; sin_(r,x); return r; }
template<class T> 
inline Interval< Float<T> > cos(const Interval< Float<T> >& x) { 
  Interval< Float<T> > r; cos_(r,x); return r; }
template<class T> 
inline Interval< Float<T> > tan(const Interval< Float<T> >& x) { 
  Interval< Float<T> > r; tan_(r,x); return r; }
template<class T> 
inline Interval< Float<T> > asin(const Interval< Float<T> >& x) { 
  Interval< Float<T> > r; asin_(r,x); return r; }
template<class T> 
inline Interval< Float<T> > acos(const Interval< Float<T> >& x) { 
  Interval< Float<T> > r; acos_(r,x); return r; }
template<class T> 
inline Interval< Float<T> > atan(const Interval< Float<T> >& x) { 
  Interval< Float<T> > r; atan_(r,x); return r; }






template<class T> template<class E> inline
Interval< Float<T> >::Interval(const Expression<E>& e) { e.assign_to(*this); }

template<class T> template<class E> inline
Interval< Float<T> >& Interval< Float<T> >::operator=(const Expression<E>& e) { 
  e.assign_to(*this); return *this; }



} // namespace Ariadne

