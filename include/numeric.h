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

// Simplifying typedefs for unsigned types
// These may be inclused in other headers,
// but repeating a typedef is not an error
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;

namespace Ariadne {

#ifdef DOXYGEN
//! \brief Integers of arbitrary size. 
//! (Only available if the Gnu Multiple Precision library (GMP) is installed.)
class Integer { };
//! \brief Rationals numbers.
//! (Only available if the Gnu Multiple Precision library (GMP) is installed.)
class Rational { };
//! \brief Floating point numbers (double precision).
class Float { };
#endif // DOXYGEN


#ifdef HAVE_GMPXX_H
typedef mpq_class Rational;
Rational sqr(const Rational& q);
Rational pow(const Rational& q, int n);
Rational pow(const Rational& q, uint n);
#else
#endif // HAVE_GMPXX_H 


typedef double Float;


using std::min;
using std::max;
using std::abs;

uint16_t fac(uint16_t n);
uint32_t fac(uint32_t n);
uint64_t fac(uint64_t n);
uint16_t bin(uint16_t n, uint16_t k);
uint32_t bin(uint32_t n, uint32_t k);
uint64_t bin(uint64_t n, uint64_t k);


template<class X> X pi();

inline Float inf() { return std::numeric_limits<double>::max(); }
inline Float eps() { return std::numeric_limits<double>::epsilon(); }

inline Float down(Float x) { return x>0 ? x*(1-2e-16) : x*(1+2e-16); }
inline Float up(Float x) { return x>0 ? x*(1+2e-16) : x*(1-2e-16); }

inline Float neg(Float x) { return -x; }
inline Float rec(Float x) { return 1.0/x; }

inline Float add(Float x, Float y) { return x+y; }
inline Float sub(Float x, Float y) { return x-y; }
inline Float mul(Float x, Float y) { return x*y; }
inline Float div(Float x, Float y) { return x/y; }

inline Float sqr(Float x) { return x*x; }
inline Float pow(Float x, int n) { return std::pow(x,Float(n)); }
inline Float pow(Float x, uint n) { return std::pow(x,Float(n)); }

inline Float sqrt(Float x) { return std::sqrt(x); }
inline Float exp(Float x) { return std::exp(x); }
inline Float log(Float x) { return std::log(x); }

template<> inline Float pi<Float>() { return 3.1415926535897931; }
inline Float sin(Float x) { return std::sin(x); }
inline Float cos(Float x) { return std::cos(x); }
inline Float tan(Float x) { return std::tan(x); }
inline Float asin(Float x) { return std::asin(x); }
inline Float acos(Float x) { return std::acos(x); }
inline Float atan(Float x) { return std::atan(x); }

inline Float add_up(Float x, Float y) { return up(x+y); }
inline Float sub_up(Float x, Float y) { return up(x-y); }
inline Float mul_up(Float x, Float y) { return up(x*y); }
inline Float div_up(Float x, Float y) { return up(x/y); }
inline Float pow_up(Float x, int n) { return up(pow(x,n)); }

inline Float add_down(Float x, Float y) { return down(x+y); }
inline Float sub_down(Float x, Float y) { return down(x-y); }
inline Float mul_down(Float x, Float y) { return down(x*y); }
inline Float div_down(Float x, Float y) { return down(x/y); }
inline Float pow_down(Float x, int n) { return down(pow(x,n)); }

inline Float add_approx(Float x, Float y) { return x+y; }
inline Float sub_approx(Float x, Float y) { return x-y; }
inline Float mul_approx(Float x, Float y) { return x*y; }
inline Float div_approx(Float x, Float y) { return x/y; }
inline Float pow_approx(Float x, int n) { return pow(x,n); }

inline Float rad_up(Float x, Float y) { return up((y-x)/2); }
inline Float med_approx(Float x, Float y) { return (x+y)/2; }

inline Float add_rnd(Float x, Float y) { return x+y; }
inline Float sub_rnd(Float x, Float y) { return x-y; }
inline Float mul_rnd(Float x, Float y) { return x*y; }
inline Float div_rnd(Float x, Float y) { return x/y; }

inline Float add_opp(Float x, Float y) { volatile Float t=(-x)-y; return -t; }
inline Float sub_opp(Float x, Float y) { volatile Float t=(-x)+y; return -t; }
inline Float mul_opp(Float x, Float y) { volatile Float t=(-x)*y; return -t; }
inline Float div_opp(Float x, Float y) { volatile Float t=(-x)/y; return -t; }


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
    Interval(Rational q);
    Interval(Rational lower, Rational upper);
#endif // HAVE_GMPXX_H 

    const Float& lower() const { return l; }
    const Float& upper() const { return u; }
    const Float midpoint() const { return (l+u)/2; }
    const Float radius() const { return up((u-l)/2); }
    const Float width() const { return up(u-l); }

    bool empty() const { return l>u; }
    bool singleton() const { return l==u; }

    void set_lower(const Float& lower) { l=lower; }
    void set_upper(const Float& upper) { u=upper; }
  public:
    double l, u;
};

inline Float midpoint(Interval i) { 
    return (i.l+i.u)/2; 
}

inline Float radius(Interval i) { 
    return up((i.u-i.l)/2); 
}

inline Float width(Interval i) { 
    return up(i.u-i.l); 
}

inline bool equal(Interval i1, Interval i2) { 
    return i1.l==i2.l && i1.u==i2.u;
}

inline bool empty(Interval i) { 
    return i.l>i.u;
}

inline bool bounded(Interval i) { 
    return i.l!=-inf() && i.u!=+inf();
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



Interval abs(Interval);
Interval neg(Interval);
Interval rec(Interval);
Interval add(Interval, Interval);
Interval sub(Interval, Interval);
Interval mul(Interval, Interval);
Interval div(Interval, Interval);
Interval mul(Interval, Float);
Interval div(Interval, Float);
Interval sqr(Interval);
Interval pow(Interval, int);
Interval pow(Interval, uint);

Interval sqrt(Interval);
Interval exp(Interval);
Interval log(Interval);

template<> Interval pi<Interval>();
Interval sin(Interval);
Interval cos(Interval);
Interval tan(Interval);
Interval asin(Interval);
Interval acos(Interval);
Interval atan(Interval);

// Standard equality operators
inline bool operator==(const Interval& i1, const Interval& i2) { return i1.l==i2.l && i1.u==i2.u; }
inline bool operator!=(const Interval& i1, const Interval& i2) { return i1.l!=i2.l || i1.u!=i2.u; }

// Boost-style tribool (in)equality operators
//inline tribool operator==(const Interval& i1, const Interval& i2) { 
//  if(i1.l>i2.u || i1.u<i2.l) { return false; } else if(i1.l==i2.u && i1.u==i2.l) { return true; } else { return indeterminate; } }
//inline tribool operator!=(const Interval& i1, const Interval& i2) { return !(i1==i2); }


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

inline Interval operator+(Interval i1, Float x2) { return Interval(down(i1.l+x2),up(i1.u+x2)); };
inline Interval operator+(Float x1, Interval i2) { return Interval(down(x1+i2.l),up(x1+i2.u)); };
inline Interval operator-(Interval i1, Float x2) { return Interval(down(i1.l-x2),up(i1.u-x2)); };
inline Interval operator-(Float x1, Interval i2) { return Interval(down(x1-i2.u),up(x1-i2.l)); };
inline Interval operator*(Interval i1, Float x2) { return mul(i1,x2); }
inline Interval operator*(Float x1, Interval i2) { return mul(i2,x1); }
inline Interval operator/(Interval i1, Float x2) { return div(i1,x2); }
inline Interval operator/(Float x1, Interval i2) { return mul(rec(i2),x1); }

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

inline bool contains(Interval i, Float x) { return i.l<=x && x<=i.u; }
inline bool subset(Float x, Interval i) { return i.l<=x && x<=i.u; }
inline bool subset(Interval i1, Interval i2) { return i2.l<=i1.l && i1.u<=i2.u; }

template<class A> void serialize(A& a, Interval& ivl, const uint version) {
    a & ivl.l & ivl.u; }

std::ostream& operator<<(std::ostream&, const Interval&);
std::istream& operator>>(std::istream&, Interval&);

} // namespace Ariadne 

#endif
