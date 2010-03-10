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
#include <stdint.h>
#include <cstdlib>

#include "tribool.h"
#include "rounding.h"
#include "macros.h"

// Simplifying typedefs for unsigned types
// These may be inclused in other headers,
// but repeating a typedef is not an error
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;

namespace Ariadne {

#ifdef DOXYGEN
//! \brief Integers of arbitrary size with exact arithmetic.
//! (Only available if the Gnu Multiple Precision library (GMP) is installed.)
//! \details
//! Unlike C++ and the Python 2, integer division is performed exactly and returns a rational.
//! The operations \c quot(Integer,Integer) and \c rem(Integer,Integer) can be used to perform integer division.
class Integer { };
//! \brief %Rational numbers with exact arithmetic.
//! (Only available if the Gnu Multiple Precision library (GMP) is installed.)
class Rational { };
//! \brief Floating point numbers (double precision) using approxiamate arithmetic.
//! \details
//! The \c Float class represents floating-point numbers. Since most arithmetic operations on floating-point numbers can only be performed approximately, %Ariadne uses <em>interval arithmetic</em> to represent the results of floating-point computations. The result of any floating-point computation is represented as an interval \f$[l,u]\f$ enclosing the exact value of the result. In this way, round-off errors can be propagated automatically.
//!
//! Ariadne floating-point numbers can be constructed by conversion from built-in C++ types or from string literals. Note that a string literal representing a \c Float must be exacly representable on the machine. Hence <c>%Float(3.3)</c> and <c>%Float("3.25")</c> are both valid (the former has a value of \f$3.2999999999999998224\ldots\f$) but <c>%Float("3.3")</c> is an error.
//! \note Constructing a %Float from a string literal is currently not supported!
//! \sa Interval
class Float { };
#endif // DOXYGEN


#ifdef HAVE_GMPXX_H
class Integer : public mpz_class {
  public:
    Integer() : mpz_class() { }
    Integer(const int& n) : mpz_class(n) { }
    Integer(const std::string& s) : mpz_class(s) { }
};

#define HAVE_RATIONAL 1
typedef mpq_class Rational;
Rational sqr(const Rational& q);
Rational pow(const Rational& q, int n);
Rational pow(const Rational& q, uint n);
#else
class Integer {
  public:
    Integer() : _value() { }
    Integer(const int& n) : _value(n) { }
    Integer(const std::string& s) : _value(std::atoi(s.c_str())) { }
    operator int () const  { return _value; }
  private:
    int _value;
};
inline Integer operator+(const Integer& z) {
    return Integer(+int(z)); }
inline Integer operator-(const Integer& z) {
    return Integer(-int(z)); }
inline Integer operator+(const Integer& z1, const Integer& z2) {
    return Integer(int(z1)+int(z2)); }
inline Integer operator-(const Integer& z1, const Integer& z2) {
    return Integer(int(z1)-int(z2)); }
inline Integer operator*(const Integer& z1, const Integer& z2) {
    return Integer(int(z1)*int(z2)); }
inline bool operator==(const Integer& z1, const Integer& z2) {
    return int(z1)==int(z2); }
inline bool operator!=(const Integer& z1, const Integer& z2) {
    return int(z1)!=int(z2); }
inline bool operator<=(const Integer& z1, const Integer& z2) {
    return int(z1)<=int(z2); }
inline bool operator>=(const Integer& z1, const Integer& z2) {
    return int(z1)>=int(z2); }
inline bool operator< (const Integer& z1, const Integer& z2) {
    return int(z1)< int(z2); }
inline bool operator> (const Integer& z1, const Integer& z2) {
    return int(z1)> int(z2); }

#endif // HAVE_GMPXX_H


typedef double Float;


using std::min;
using std::max;
using std::abs;

uint32_t fac(uint8_t n);
uint16_t fac(uint16_t n);
uint32_t fac(uint32_t n);
uint64_t fac(uint64_t n);
uint32_t bin(uint8_t n, uint8_t k);
uint16_t bin(uint16_t n, uint16_t k);
uint32_t bin(uint32_t n, uint32_t k);
uint64_t bin(uint64_t n, uint64_t k);


template<class X> X pi();
template<class X> inline X inf();

template<> inline Float inf<Float>() { return std::numeric_limits<double>::infinity(); }
#ifdef HAVE_RATIONAL
template<> inline Rational inf<Rational>() { return Rational(1,-0); }
#endif


inline Float mx() { return std::numeric_limits<double>::max(); }
inline Float eps() { return std::numeric_limits<double>::epsilon(); }

inline Float down(Float x) { return x>0 ? x*(1-2e-16) : x*(1+2e-16); }
inline Float up(Float x) { return x>0 ? x*(1+2e-16) : x*(1-2e-16); }

inline Float min(Float x, Float y) { return std::min(x,y); }
inline Float max(Float x, Float y) { return std::max(x,y); }

// The magnitude of a quantity is the funciton norm; always a real
//   Provided for consistency with other types for which mag and abs are different
inline Float mag(Float x) { return abs(x); }

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

inline Float add_rnd(Float x, Float y) { return (volatile double&)x+(volatile double&)y; }
inline Float sub_rnd(Float x, Float y) { return (volatile double&)x-(volatile double&)y; }
inline Float mul_rnd(Float x, Float y) { return (volatile double&)x*(volatile double&)y; }
inline Float div_rnd(Float x, Float y) { return (volatile double&)x/(volatile double&)y; }

inline Float add_opp(Float x, Float y) { volatile double t=(-x)-y; return -t; }
inline Float sub_opp(Float x, Float y) { volatile double t=(-x)+y; return -t; }
inline Float mul_opp(Float x, Float y) { volatile double t=(-x)*y; return -t; }
inline Float div_opp(Float x, Float y) { volatile double t=(-x)/y; return -t; }

Float pow_rnd(Float x, int n);
Float pow_opp(Float x, int n);

Float sqrt_rnd(Float x);
Float exp_rnd(Float x);
Float log_rnd(Float x);
Float sin_rnd(Float x);
Float cos_rnd(Float x);
Float tan_rnd(Float x);

inline Float add_approx(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(to_nearest);
    volatile Float r=add_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float sub_approx(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(to_nearest);
    volatile Float r=sub_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float mul_approx(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(to_nearest);
    volatile Float r=mul_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float div_approx(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(to_nearest);
    volatile Float r=div_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float pow_approx(Float x, int n) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(to_nearest);
    volatile Float r=pow_rnd(x,n); set_rounding_mode(rounding_mode); return r; }

inline Float add_up(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(upward);
    volatile Float r=add_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float sub_up(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(upward);
    volatile Float r=sub_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float mul_up(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(upward);
    volatile Float r=mul_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float div_up(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(upward);
    volatile Float r=div_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float pow_up(Float x, int n) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(upward);
    volatile Float r=pow_rnd(x,n); set_rounding_mode(rounding_mode); return r; }

inline Float add_down(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(downward);
    volatile Float r=add_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float sub_down(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(downward);
    volatile Float r=sub_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float mul_down(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(downward);
    volatile Float r=mul_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float div_down(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(downward);
    volatile Float r=div_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float pow_down(Float x, int n) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(downward);
    volatile Float r=pow_rnd(x,n); set_rounding_mode(rounding_mode); return r; }

inline Float med_approx(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(to_nearest);
    volatile Float r=(x+y)/2; set_rounding_mode(rounding_mode); return r; }
inline Float rad_up(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(upward);
    volatile Float r=(y-x)/2; set_rounding_mode(rounding_mode); return r; }




//! \brief Intervals with floating-point endpoints supporting outwardly-rounded arithmetic.
//! \details
//! Note that <c>%Interval(3.3)</c> yields the singleton interval \f$[3.2999999999999998224,3.2999999999999998224]\f$ (the constant is first interpreted by the C++ compiler to give a C++ \c double, whereas <c>%Interval("3.3")</c> yields the interval \f$[3.2999999999999998224,3.3000000000000002665]\f$ enclosing \f$3.3\f$.
//!
//! Comparison tests on \c Interval use the idea that an interval represents a single number with an unknown value. Hence the result is of type \c tribool, which can take values { \c True, \c False, \c Indeterminate }.  Hence a test \f$[l_1,u_1]\leq [l_2,u_2]\f$ returns \c True if \f$u_1\leq u_2\f$, since in this case \f$x_1\leq x_2\f$ whenever \f$x_1\in[l_1,u_2]\f$ and \f$x_2\in[l_2,u_2]\f$, \c False if \f$l_1>u_2\f$, since in this case we know \f$x_1>x_2\f$, and \c Indeterminate otherwise, since in this case we can find \f$x_1,x_2\f$ making the result either true or false. In the case of equality, the comparison \f$[l_1,u_1]\f$==\f$[l_2,u_2]\f$ only returns \c True if both intervals are singletons, since otherwise we can find values making the result either true of false.
//!
//! To obtain the lower and upper bounds of an interval, use \c ivl.lower() and \c ivl.upper().
//! To obtain the midpoint and radius, use \c ivl.midpoint() and \c ivl.radius().
//! Alternatives \c midpoint(ivl) and \c radius(ivl) are also provided.
//! Note that \c midpoint and \c radius return approximations to the true midpoint and radius of the interval. If \f$m\f$ and \f$r\f$ are the returned midpoint and radius of the interval \f$[l,u]\f$, the using exact arithmetic, we guarentee \f$m-r\leq l\f$ and \f$m+r\geq u\f$
//!
//! To test if an interval contains a point or another interval, use \c encloses(Interval,Float) or \c encloses(Interval,Interval).
//! The test \c refines(Interval,Interval) can also be used.
//! \sa Float
//!
//! \par Python interface
//!
//! In the Python interface, %Ariadne intervals can be constructed from Python literals of the form \c {a:b} or (deprecated) \c [a,b] .
//! The former is preferred, as it cannot be confused with literals for other classes such as Vector and Array types.
//! Automatic conversion is used to convert Interval literals of the form \c {a,b} to an Interval in functions.
//! 
//! Care must be taken when defining intervals using floating-point coefficients, since values are first converted to the nearest
//! representable value by the Python interpreter. <br><br>
//! \code
//!   Interval({1.1:2.3}) # Create the interval [1.1000000000000001, 2.2999999999999998]
//!   Interval({2.5:4.25}) # Create the interval [2.5, 4.25], which can be represented exactly
//!   Interval([2.5,4.25]) # Alternative syntax for creating the interval [2.5, 4.25]
//! \endcode
class Interval {
  public:
    Interval() : l(0.0), u(0.0) { }
    Interval(uint m) : l(m), u(m) { }
    Interval(int n) : l(n), u(n) { }
    Interval(Float x) : l(x), u(x) { }
    Interval(const Interval& i) : l(i.l), u(i.u) { }

    Interval(Float lower, Float upper) : l(lower), u(upper) { ARIADNE_ASSERT_MSG(lower<=upper, "lower = "<<lower<<", upper ="<<upper); }
#ifdef HAVE_RATIONAL
    Interval(Rational q);
    Interval(Rational lower, Rational upper);
#endif // HAVE_RATIONAL

    const Float& lower() const { return l; }
    const Float& upper() const { return u; }
    const Float midpoint() const { return add_approx(l,u)/2; }
    const Float radius() const { return sub_up(u,l)/2; }
    const Float width() const { return sub_up(u,l); }

    bool empty() const { return l>u; }
    bool singleton() const { return l==u; }

    void set_lower(const Float& lower) { ARIADNE_ASSERT(lower<=this->u); l=lower; }
    void set_upper(const Float& upper) { ARIADNE_ASSERT(this->l<=upper); u=upper; }
    void set(const Float& lower, const Float& upper) { ARIADNE_ASSERT(lower<=upper); l=lower; u=upper; }
  public:
  private:
    double l, u;
};

std::ostream& operator<<(std::ostream& os, const Interval& ivl);

inline Float midpoint(Interval i) {
    return add_approx(i.lower(),i.upper())/2;
}

inline Float radius(Interval i) {
    return sub_up(i.upper(),i.lower())/2;
}

inline Float width(Interval i) {
    return sub_up(i.upper(),i.lower());
}

inline bool equal(Interval i1, Interval i2) {
    //std::cerr<<"equal(i1,i2) with i1="<<i1<<"; i2="<<i2<<std::endl;
    return i1.lower()==i2.lower() && i1.upper()==i2.upper();
}

inline bool empty(Interval i) {
    return i.lower()>i.upper();
}

inline bool bounded(Interval i) {
    return i.lower()!=-inf<Float>() && i.upper()!=+inf<Float>();
}

inline Interval intersection(Interval i1, Interval i2) {
    if(i1.lower()>i2.upper() || i1.upper()<i2.lower()) { return Interval(1,-1); }
    return Interval(max(i1.lower(),i2.lower()),min(i1.upper(),i2.upper()));
}

inline Interval hull(Interval i1, Interval i2) {
    return Interval(min(i1.lower(),i2.lower()),max(i1.upper(),i2.upper()));
}

Interval trunc(Interval, uint eps);

inline Float med(Interval i) { return (i.lower()+i.upper())/2; }
inline Float rad(Interval i) { return up((i.upper()-i.lower())/2); }
inline Float diam(Interval i) { return up(i.upper()-i.lower()); }

inline Interval max(Interval,Interval);
inline Interval min(Interval,Interval);
inline Interval abs(Interval);
inline Interval neg(Interval);
inline Interval add(Interval, Interval);
inline Interval add(Interval, Float);
inline Interval sub(Interval, Interval);
inline Interval sub(Interval, Float);
inline Interval sub(Float, Interval);

Interval rec(Interval);
Interval mul(Interval, Interval);
Interval div(Interval, Interval);
Interval mul(Interval, Float);
Interval div(Interval, Float);
Interval div(Float, Interval);

inline Interval add_ivl(Float, Float);
inline Interval sub_ivl(Float, Float);
inline Interval mul_ivl(Float, Float);
inline Interval div_ivl(Float, Float);

Interval sqr(Interval);
Interval pow(Interval, uint);
Interval pow(Interval, int);

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


inline Float mag(Interval i) { return max(abs(i.lower()),abs(i.upper())); }
inline Float mig(Interval i) { return min(abs(i.lower()),abs(i.upper())); }

inline bool contains(Interval i, Float x) { return i.lower()<=x && x<=i.upper(); }

inline bool subset(Interval i1, Interval i2) { return i1.lower()>=i2.lower() && i1.upper()<=i2.upper(); }
inline bool intersect(Interval i1, Interval i2) { return i1.lower()<=i2.upper() && i1.upper()>=i2.lower(); }
inline bool disjoint(Interval i1, Interval i2) { return i1.lower()>i2.upper() || i1.upper()<i2.lower(); }
inline bool overlap(Interval i1, Interval i2) { return i1.lower()<i2.upper() && i1.upper()>i2.lower(); }
inline bool inside(Interval i1, Interval i2) { return i1.lower()>i2.lower() && i1.upper()<i2.upper(); }
inline bool covers(Interval i1, Interval i2) { return i1.lower()<i2.lower() && i1.upper()>i2.upper(); }

inline Interval max(Interval i1, Interval i2)
{
    return Interval(max(i1.lower(),i2.lower()),max(i1.upper(),i2.upper()));
}

inline Interval min(Interval i1, Interval i2)
{
    return Interval(min(i1.lower(),i2.lower()),min(i1.upper(),i2.upper()));
}


inline Interval abs(Interval i)
{
    if(i.lower()>=0) {
        return Interval(i.lower(),i.upper());
    } else if(i.upper()<=0) {
        return Interval(-i.upper(),-i.lower());
    } else {
        return Interval(0.0,max(-i.lower(),i.upper()));
    }
}

inline Interval neg(Interval i)
{
    return Interval(-i.upper(),-i.lower());
}

inline Interval add(Interval i1, Interval i2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& i1l=const_cast<volatile double&>(i1.lower());
    volatile double& i1u=const_cast<volatile double&>(i1.upper());
    volatile double& i2l=const_cast<volatile double&>(i2.lower());
    volatile double& i2u=const_cast<volatile double&>(i2.upper());
    set_rounding_mode(downward);
    volatile double rl=i1l+i2l;
    set_rounding_mode(upward);
    volatile double ru=i1u+i2u;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

inline Interval add(Interval i1, Float x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& i1l=const_cast<volatile double&>(i1.lower());
    volatile double& i1u=const_cast<volatile double&>(i1.upper());
    volatile double& x2v=x2;
    set_rounding_mode(downward);
    volatile double rl=i1l+x2v;
    set_rounding_mode(upward);
    volatile double ru=i1u+x2v;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

inline Interval add_ivl(Float x1, Float x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& x1v=x1;
    volatile double& x2v=x2;
    set_rounding_mode(downward);
    volatile double rl=x1v+x2v;
    set_rounding_mode(upward);
    volatile double ru=x1v+x2v;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

inline Interval sub(Interval i1, Interval i2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& i1l=const_cast<volatile double&>(i1.lower());
    volatile double& i1u=const_cast<volatile double&>(i1.upper());
    volatile double& i2l=const_cast<volatile double&>(i2.lower());
    volatile double& i2u=const_cast<volatile double&>(i2.upper());
    set_rounding_mode(downward);
    volatile double rl=i1l-i2u;
    set_rounding_mode(upward);
    volatile double ru=i1u-i2l;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

inline Interval sub(Interval i1, Float x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& i1l=const_cast<volatile double&>(i1.lower());
    volatile double& i1u=const_cast<volatile double&>(i1.upper());
    volatile double& x2v=x2;
    set_rounding_mode(downward);
    volatile double rl=i1l-x2v;
    set_rounding_mode(upward);
    volatile double ru=i1u-x2v;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

inline Interval sub(Float x1, Interval i2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& x1v=x1;
    volatile double& i2l=const_cast<volatile double&>(i2.lower());
    volatile double& i2u=const_cast<volatile double&>(i2.upper());
    set_rounding_mode(downward);
    volatile double rl=x1v-i2u;
    set_rounding_mode(upward);
    volatile double ru=x1v-i2l;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

inline Interval sub_ivl(Float x1, Float x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& x1v=x1;
    volatile double& x2v=x2;
    set_rounding_mode(downward);
    volatile double rl=x1v-x2v;
    set_rounding_mode(upward);
    volatile double ru=x1v-x2v;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

inline Interval mul_ivl(Float x1, Float x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& x1v=x1;
    volatile double& x2v=x2;
    set_rounding_mode(downward);
    volatile double rl=x1v*x2v;
    set_rounding_mode(upward);
    volatile double ru=x1v*x2v;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

inline Interval div_ivl(Float x1, Float x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& x1v=x1;
    volatile double& x2v=x2;
    set_rounding_mode(downward);
    volatile double rl=x1v/x2v;
    set_rounding_mode(upward);
    volatile double ru=x1v/x2v;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

inline Interval med_ivl(Float x1, Float x2)
{
    return add_ivl(x1/2,x2/2);
}

inline Interval rad_ivl(Float x1, Float x2)
{
    return sub_ivl(x2/2,x1/2);
}

inline Interval med_ivl(Interval i) {
    return add_ivl(i.lower()/2,i.upper()/2);
}

inline Interval rad_ivl(Interval i) {
    return sub_ivl(i.upper()/2,i.lower()/2);
}

// Standard equality operators
inline bool operator==(const Interval& i1, const Interval& i2) { return i1.lower()==i2.lower() && i1.upper()==i2.upper(); }
inline bool operator!=(const Interval& i1, const Interval& i2) { return i1.lower()!=i2.lower() || i1.upper()!=i2.upper(); }

// Boost-style tribool (in)equality operators
//inline tribool operator==(const Interval& i1, const Interval& i2) {
//  if(i1.lower()>i2.upper() || i1.upper()<i2.lower()) { return false; } else if(i1.lower()==i2.upper() && i1.upper()==i2.lower()) { return true; } else { return indeterminate; } }
//inline tribool operator!=(const Interval& i1, const Interval& i2) { return !(i1==i2); }


inline Interval operator+(Interval i) { return Interval(i.lower(),i.upper()); }
inline Interval operator-(Interval i) { return Interval(-i.upper(),-i.lower()); }
inline Interval operator+(Interval i1, Interval i2) { return add(i1,i2); }
inline Interval operator-(Interval i1, Interval i2) { return sub(i1,i2); }
inline Interval operator*(Interval i1, Interval i2) { return mul(i1,i2); }
inline Interval operator/(Interval i1, Interval i2) { return div(i1,i2); };

inline Interval& operator+=(Interval& i1, Interval i2) { i1=add(i1,i2); return i1; }
inline Interval& operator-=(Interval& i1, Interval i2) { i1=sub(i1,i2); return i1; }
inline Interval& operator*=(Interval& i1, Interval i2) { i1=mul(i1,i2); return i1; }
inline Interval& operator/=(Interval& i1, Interval i2) { i1=div(i1,i2); return i1; }

inline Interval operator+(Interval i1, Float x2) { return add(i1,x2); }
inline Interval operator+(Float x1, Interval i2) { return add(i2,x1); }
inline Interval operator-(Interval i1, Float x2) { return sub(i1,x2); }
inline Interval operator-(Float x1, Interval i2) { return sub(x1,i2); }
inline Interval operator*(Interval i1, Float x2) { return mul(i1,x2); }
inline Interval operator*(Float x1, Interval i2) { return mul(i2,x1); }
inline Interval operator/(Interval i1, Float x2) { return div(i1,x2); }
inline Interval operator/(Float x1, Interval i2) { return div(x1,i2); }

inline Interval& operator+=(Interval& i1, Float x2) { i1=add(i1,x2); return i1; }
inline Interval& operator-=(Interval& i1, Float x2) { i1=sub(i1,x2); return i1; }
inline Interval& operator*=(Interval& i1, Float x2) { i1=mul(i1,x2); return i1; }
inline Interval& operator/=(Interval& i1, Float x2) { i1=div(i1,x2); return i1; }

inline tribool operator==(Interval i1, Float x2) {
    if(i1.upper()<x2 || i1.lower()>x2) { return false; }
    else if(i1.lower()==x2 && i1.upper()==x2) { return true; }
    else { return indeterminate; }
}

inline tribool operator!=(Interval i1, Float x2) {
    if(i1.upper()<x2 || i1.lower()>x2) { return true; }
    else if(i1.lower()==x2 && i1.upper()==x2) { return false; }
    else { return indeterminate; }
}

inline tribool operator> (Interval i1, Float x2) {
    if(i1.lower()> x2) { return true; }
    else if(i1.upper()<=x2) { return false; }
    else { return indeterminate; }
}

inline tribool operator< (Interval i1, Float x2) {
    if(i1.upper()< x2) { return true; }
    else if(i1.lower()>=x2) { return false; }
    else { return indeterminate; }
}

inline tribool operator>=(Interval i1, Float x2) {
    if(i1.lower()>=x2) { return true; }
    else if(i1.upper()< x2) { return false; }
    else { return indeterminate; }
}

inline tribool operator<=(Interval i1, Float x2) {
    if(i1.upper()<=x2) { return true; }
    else if(i1.lower()> x2) { return false; }
    else { return indeterminate; }
}



inline tribool operator> (Interval i1, Interval i2) {
    if(i1.lower()> i2.upper()) { return true; }
    else if(i1.upper()<=i2.lower()) { return false; }
    else { return indeterminate; }
}

inline tribool operator< (Interval i1, Interval i2) {
    if(i1.upper()< i2.lower()) { return true; }
    else if(i1.lower()>=i2.upper()) { return false; }
    else { return indeterminate; }
}

inline tribool operator>=(Interval i1, Interval i2) {
    if(i1.lower()>=i2.upper()) { return true; }
    else if(i1.upper()< i2.lower()) { return false; }
    else { return indeterminate; }
}

inline tribool operator<=(Interval i1, Interval i2) {
    if(i1.upper()<=i2.lower()) { return true; }
    else if(i1.lower()> i2.upper()) { return false; }
    else { return indeterminate; }
}

template<class A> void serialize(A& a, Interval& ivl, const uint version) {
    a & ivl.lower() & ivl.upper(); }

std::ostream& operator<<(std::ostream&, const Interval&);
std::istream& operator>>(std::istream&, Interval&);

} // namespace Ariadne

#endif
