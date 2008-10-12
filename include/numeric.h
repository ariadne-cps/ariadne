/***************************************************************************
 *            numeric.h
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
 
/*! \file numeric.h
 *  \brief Numerical classes.
 */
#ifndef ARIADNE_NUMERIC_H
#define ARIADNE_NUMERIC_H

#ifdef HAVE_GMPXX_H
#include <gmpxx.h>
#endif // HAVE_GMPXX_H

#include <cmath>
#include <limits>

#include "tribool.h"

namespace Ariadne {

using std::min;
using std::max;

uint fac(uint n);
uint bin(uint n, uint k);

#ifdef HAVE_GMPXX_H
typedef mpz_class Integer;
typedef mpq_class Rational;
#endif // HAVE_GMPXX_H 


#ifdef DOXYGEN
//! \brief Floating point numbers.
class Float { };
#endif

typedef double Float;
//class Float : public double { template<class T> Float(const T& t) : double(t) { } };

inline Float inf() { return std::numeric_limits<double>::max(); }
inline Float eps() { return std::numeric_limits<double>::epsilon(); }

inline Float down(Float x) { return x>0 ? x*(1-2e-16) : x*(1+2e-16); }
inline Float up(Float x) { return x>0 ? x*(1+2e-16) : x*(1-2e-16); }

inline Float abs(Float x) { return std::fabs(x); }
inline Float sqr(Float x) { return x*x; }
inline Float pow(Float x, int n) { return std::pow(x,Float(n)); }
inline Float pow(Float x, uint n) { return std::pow(x,Float(n)); }

inline Float sqrt(Float x) { return std::sqrt(x); }
inline Float exp(Float x) { return std::exp(x); }
inline Float log(Float x) { return std::log(x); }

inline Float add_up(Float x, Float y) { return up(x+y); }
inline Float sub_up(Float x, Float y) { return up(x-y); }
inline Float mul_up(Float x, Float y) { return up(x*y); }
inline Float div_up(Float x, Float y) { return up(x/y); }
inline Float rad_up(Float x, Float y) { return up((y-x)/2); }

inline Float add_approx(Float x, Float y) { return x+y; }
inline Float sub_approx(Float x, Float y) { return x-y; }
inline Float mul_approx(Float x, Float y) { return x*y; }
inline Float div_approx(Float x, Float y) { return x/y; }
inline Float med_approx(Float x, Float y) { return (x+y)/2; }



//! \brief Intervals supporting interval arithmetic.
class Interval {
 public:
  Interval() : l(0.0), u(0.0) { }
  Interval(uint m) : l(m), u(m) { }
  Interval(int n) : l(n), u(n) { }
  Interval(Float x) : l(x), u(x) { }
  Interval(const Interval& i) : l(i.l), u(i.u) { }
  
  Interval(Float lower, Float upper) : l(lower), u(upper) { }
#ifdef HAVE_GMPXX_H
  Interval(Rational lower, Rational upper) : l(down(lower.get_d())), u(up(upper.get_d())) { }
#endif // HAVE_GMPXX_H 

  const Float& lower() const { return l; }
  const Float& upper() const { return u; }
  const Float midpoint() const { return (l+u)/2; }
  const Float radius() const { return up((u-l)/2); }
  const Float width() const { return up(u-l); }
 public:
  double l, u;
};

inline Float midpoint(Interval i) { 
  return (i.l+i.u)/2; 
}

inline Interval intersection(Interval i1, Interval i2) { 
  if(i1.l>i2.u || i1.u<i2.l) { return Interval(1,-1); }
  return Interval(max(i1.l,i2.l),min(i1.u,i2.u));
}

inline Interval hull(Interval i1, Interval i2) { 
  return Interval(min(i1.l,i2.l),max(i1.u,i2.u));
}

Interval trunc(Interval, uint eps);

inline Float med(Interval i) { return (i.l+i.u)/2; }
inline Float rad(Interval i) { return up((i.u-i.l)/2); }
inline Float diam(Interval i) { return up(i.u-i.l); }

Interval neg(Interval);
Interval add(Interval, Interval);
Interval sub(Interval, Interval);
Interval mul(Interval, Interval);
Interval mul(Float, Interval);
Interval abs(Interval);
Interval sqr(Interval);
Interval pow(Interval, int);
Interval pow(Interval, uint);
Interval rec(Interval);

Interval sqrt(Interval);
Interval exp(Interval);
Interval log(Interval);

inline bool operator==(const Interval& i1, const Interval& i2) { return i1.l==i2.l && i1.u==i2.u; }
inline bool operator!=(const Interval& i1, const Interval& i2) { return i1.l!=i2.l || i1.u!=i2.u; }

inline Interval operator+(Interval i) { return Interval(i.l,i.u); }
inline Interval operator-(Interval i) { return Interval(-i.u,-i.l); }
inline Interval operator+(Interval i1, Interval i2) { return Interval(down(i1.l+i2.l),up(i1.u+i2.u)); }
inline Interval operator-(Interval i1, Interval i2) { return Interval(down(i1.l-i2.u),up(i1.u-i2.l)); };
inline Interval operator*(Interval i1, Interval i2) { return mul(i1,i2); }
inline Interval operator/(Interval i1, Interval i2) { return mul(i1,rec(i2)); };

inline Interval& operator+=(Interval& i1, Interval i2) { i1=add(i1,i2); return i1; }
inline Interval& operator-=(Interval& i1, Interval i2) { i1=sub(i1,i2); return i1; }
inline Interval& operator*=(Interval& i1, Interval i2) { i1=mul(i1,i2); return i1; }
inline Interval& operator/=(Interval& i1, Interval i2) { i1=mul(i1,rec(i2)); return i1; }

inline Interval operator*(Interval i1, Float x2) { return mul(x2,i1); }

inline tribool operator>(Interval i1, Interval i2) { 
  if(i1.lower()>i2.upper()) { return true; }
  else if(i1.upper()<=i2.lower()) { return false; }
  else { return indeterminate; }
}

inline tribool operator<(Interval i1, Interval i2) { 
  if(i1.upper()<i2.lower()) { return true; }
  else if(i1.lower()>=i2.upper()) { return false; }
  else { return indeterminate; }
}


inline Float mag(Interval i) { return max(abs(i.l),abs(i.u)); }
inline Float mig(Interval i) { return min(abs(i.l),abs(i.u)); }

inline bool subset(Float x, Interval i) { return i.l<=x && x<=i.u; }
inline bool subset(Interval i1, Interval i2) { return i2.l<=i1.l && i1.u<=i2.u; }

std::ostream& operator<<(std::ostream&, const Interval&);

} // namespace Ariadne 

#endif
