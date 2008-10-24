/***************************************************************************
 *            numeric.cc
 *
 *  Copyright 2008  Pieter Collins
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
 
#include <iostream>

#include "config.h"
#include "macros.h"
#include "numeric.h"

namespace Ariadne {


uint 
fac(uint n) 
{ 
  ARIADNE_ASSERT(n<13); 
  uint r=1; 
  for(uint i=1; i<=n; ++i) { 
    r*=i; 
  } 
  return r; 
}


uint 
bin(uint n, uint k) 
{ 
  ARIADNE_ASSERT(n<13); 
  if(k>n+1) { std::cerr << "bin("<<n<<","<<k<<")\n"; }
  if(k==n+1) { return 0; }
  ARIADNE_ASSERT(k>=0 && k<=n);
  uint r=1; 
  for(uint i=1; i<=k; ++i) { 
    r*=(n+1-i); 
    r/=i; 
  } 
  return r; 
}


Interval neg(Interval i) 
{
  return Interval(-i.u,-i.l);
}

Interval rec(Interval i) 
{
  if(i.l>0 || i.u<0) {
    return Interval(down(1/i.u),up(1/i.l)); 
  } else {
    ARIADNE_THROW(DivideByZeroException,"rec(Interval i)","i="<<i<<")");
  }
}

Interval add(Interval i1, Interval i2) 
{
  return Interval(down(i1.l+i2.l),up(i1.u+i2.u));
}

Interval sub(Interval i1, Interval i2) 
{
  return Interval(down(i1.l-i2.u),up(i1.u-i2.l));
}

Interval mul(Interval i1, Interval i2) 
{
  if(i1.l>=0) {
    if(i2.l>=0) {
      return Interval(down(i1.l*i2.l),up(i1.u*i2.u));
    } else if(i2.u<=0) {
      return Interval(down(i1.u*i2.l),up(i1.l*i2.u));
    } else {
      return Interval(down(i1.u*i2.l),up(i1.u*i2.u));
    }
  }
  else if(i1.u<=0) {
    if(i2.l>=0) {
      return Interval(down(i1.l*i2.u),up(i1.u*i2.l));
    } else if(i2.u<=0) {
      return Interval(down(i1.u*i2.u),up(i1.l*i2.l));
    } else {
      return Interval(down(i1.l*i2.u),up(i1.l*i2.l));
    }
  } else {
    if(i2.l>=0) {
      return Interval(down(i1.l*i2.u),up(i1.u*i2.u));
    } else if(i2.u<=0) {
      return Interval(down(i1.u*i2.l),up(i1.l*i2.l));;
    } else {
      return Interval(down(min(i1.u*i2.l,i1.l*i2.u)),up(max(i1.l*i2.l,i1.u*i2.u)));
    }
  }
}


Interval div(Interval i1, Interval i2) 
{
  if(i2.l>=0) {
    if(i1.l>=0) {
      return Interval(down(i1.l/i2.u),up(i1.u/i2.l));
    } else if(i1.u<=0) {
      return Interval(down(i1.l/i2.l),up(i1.u/i2.u));
    } else {
      return Interval(down(i1.l/i2.l),up(i1.u/i2.l));
    }
  }
  else if(i2.u<=0) {
    if(i1.l>=0) {
      return Interval(down(i1.l/i2.l),up(i1.u/i2.u));
    } else if(i1.u<=0) {
      return Interval(down(i1.u/i2.l),up(i1.l/i2.u));
    } else {
      return Interval(down(i1.u/i2.u),up(i1.l/i2.u));;
    } 
  }
  else {
    return Interval(-inf(),+inf());
  }
}



Interval mul(Interval i, Float x) 
{
  if(x>=0) {
    return Interval(down(i.l*x),up(i.u*x)); 
  } else {
    return Interval(down(i.u*x),up(i.l*x));
  }
}


Interval div(Interval i, Float x) 
{
  if(x>0) {
    return Interval(down(i.l/x),up(i.u/x)); 
  } else if(x<0) {
    return Interval(down(i.u/x),up(i.l/x));
  } else {
    return Interval(-inf(),+inf());
  }
}


Interval trunc(Interval x, uint n) 
{
  Interval e=Interval(pow(2.0,52-n));
  Interval y=x+e;
  return y-e;
}

Interval abs(Interval i) 
{
  if(i.lower()>=0) {
    return i;
  } else if(i.upper()<=0) {
    return -i;
  } else {
    return Interval(0,max(-i.lower(),i.upper()));
  }
}


Interval sqr(Interval i) 
{
  if(i.l >=0) {
    return Interval(down(i.l*i.l),up(i.u*i.u));
  } else if(i.u<=0) {
    return Interval(down(i.u*i.u),up(i.l*i.l));
  } else {
    return Interval(0.0,up(max(i.l*i.l,i.u*i.u)));
  }
}

Interval pow(Interval i, int n) 
{
  if(n<0) { return pow(rec(i),-n); }
  Interval r=1; Interval p=i;
  while(n>0) { if(n%2==1) { r*=p; } p=sqr(p); n/=2; }
  return r;
}

Interval pow(Interval i, uint m) 
{
  Interval r=1; Interval p=i;
  while(m>0) { if(m%2==1) { r*=p; } p=sqr(p); m/=2; }
  return r;
}

#warning Interval transcendental functions not correct

Interval sqrt(Interval i)
{
  return Interval(down(sqrt(i.l)),up(sqrt(i.u)));
}

Interval exp(Interval i)
{
  return Interval(down(exp(i.l)),up(exp(i.u)));
}

Interval log(Interval i)
{
  return Interval(down(log(i.l)),up(log(i.u)));
}

template<> Interval pi<Interval>()
{
  return Interval(3.1415926535897927,3.1415926535897936);
}


Interval sin(Interval i)
{
  ARIADNE_NOT_IMPLEMENTED;
}

Interval cos(Interval i)
{
  ARIADNE_NOT_IMPLEMENTED;
}

Interval tan(Interval i)
{
  ARIADNE_NOT_IMPLEMENTED;
}

Interval asin(Interval i)
{
  ARIADNE_NOT_IMPLEMENTED;
}

Interval acos(Interval i)
{
  ARIADNE_NOT_IMPLEMENTED;
}

Interval atan(Interval i)
{
  ARIADNE_NOT_IMPLEMENTED;
}



#ifdef HAVE_GMPXX_H 

typedef mpq_class Rational;

Rational sqr(const Rational& q) { return q*q; }

Rational pow(const Rational& q, uint n) {
  if(n==0) { return 1; }
  Rational r=1; Rational p=q; uint m=n;
  while(m>=1) { if(m%2) { r*=p; } m/=2; p=p*p; }
  return r;
}

Rational pow(const Rational& q, int n) {
  if(n>=0) { return pow(q,uint(n)); }
  else { return pow(1/q,uint(-n)); }
}


Interval::Interval(Rational q)
  : l(down(q.get_d())), u(up(q.get_d())) { }

Interval::Interval(Rational lower, Rational upper)
  : l(down(lower.get_d())), u(up(upper.get_d())) { }

#endif // HAVE_GMPXX_H 


std::ostream& 
operator<<(std::ostream& os, const Interval& ivl)
{
  return os << '[' << ivl.l << ':' << ivl.u << ']';
}

std::istream& 
operator>>(std::istream& is, Interval& ivl)
{
  char cl,cm,cr;
  is >> cl >> ivl.l >> cm >> ivl.u >> cr;
  ARIADNE_ASSERT(is);
  ARIADNE_ASSERT(cl=='[' || cl=='(');
  ARIADNE_ASSERT(cm==':' || cm==',' || cm==';');
  ARIADNE_ASSERT(cr==']' || cr==')');
  return is;
}

} // namespace Ariadne

