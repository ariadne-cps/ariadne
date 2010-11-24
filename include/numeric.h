/***************************************************************************
 *      numeric.h
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
#include <cassert>
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

// Forward declarations
class Float;
class Interval;
class Real;

template<class X> X pi();
template<class X> X inf();

#ifdef DOXYGEN
//! \ingroup NumericModule
//! \brief Integers of arbitrary size with exact arithmetic.
//! (Only available if the Gnu Multiple Precision library (GMP) is installed.)
//! \details
//! Unlike C++ and the Python 2, integer division is performed exactly and returns a rational.
//! The operations \c quot(Integer,Integer) and \c rem(Integer,Integer) can be used to perform integer division.
class Integer { };
//! \ingroup NumericModule
//! \brief %Rational numbers with exact arithmetic.
//! (Only available if the Gnu Multiple Precision library (GMP) is installed.)
class Rational { };
#endif // DOXYGEN


#ifdef HAVE_GMPXX_H
class Integer : public mpz_class {
  public:
    Integer() : mpz_class() { }
    Integer(const int& n) : mpz_class(n) { }
    Integer(const std::string& s) : mpz_class(s) { }
};

typedef mpq_class Rational;
Rational sqr(const Rational& q);
Rational pow(const Rational& q, int n);
Rational pow(const Rational& q, uint n);
template<> inline Rational inf<Rational>() { return Rational(1,-0); }
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

uint32_t fac(uint8_t n);
uint16_t fac(uint16_t n);
uint32_t fac(uint32_t n);
uint64_t fac(uint64_t n);
uint32_t bin(uint8_t n, uint8_t k);
uint16_t bin(uint16_t n, uint16_t k);
uint32_t bin(uint32_t n, uint32_t k);
uint64_t bin(uint64_t n, uint64_t k);

using std::min;
using std::max;

//! \ingroup NumericModule
//! \brief Floating point numbers (double precision) using approxiamate arithmetic.
//! \details
//! The \c %Float class represents floating-point numbers.
//! Unless otherwise mentioned, operations on floating-point numbers are performed approximately, with no guarantees
//! on the output.
//!
//! To implement <em>interval arithmetic</em>, arithmetical operations of \c %Float can be performed with guaranteed rounding by
//! specifying \c _up and \c _down suffixes to arithmetical functions \c add, \c sub, \c mul and \c div.
//! Additionally, operations can be performed in the current <em>rounding mode</em> by using the \c _rnd suffix,
//! or with rounding reversed using the \c _opp suffix.
//! Operations can be specified to return an \c %Interval answer by using the \c _ivl suffix.
//! The \c _approx suffix is provided to specifically indicate that the operation is computed approximately.
//!
//! %Ariadne floating-point numbers can be constructed by conversion from built-in C++ types.
//! Note that the value of a built-in floating-point value may differ from the mathematical value of the literal.
//! For example, while <c>%Float(3.25)</c> is represented exactly, <c>%Float(3.3)</c> has a value of \f$3.2999999999999998224\ldots\f$.
//! \note In the future, the construction of a \c %Float from a string literal may be supported.
//! \sa Interval, Real
class Float {
  public:
    volatile double v;
  public:
    typedef Float NumericType;
  public:
    //! \brief Default constructor creates an uninitialised number.
    Float() : v() { }
    Float(unsigned int m) : v(m) { }
    Float(int n) : v(n) { }
    //! \brief Convert from a built-in double-precision floating-point number.
    Float(double x) : v(x) { }
    //! \brief Copy constructor.
    Float(const Float& x) : v(x.v) { }
    //! \brief Convert from a general real number by generating a representable approximation,
    //! not necessarily the nearest.
    Float(const Real& x);
    Float& operator=(double x) { v=x; return *this; }
    Float& operator=(const Float& x) { v=x.v; return *this; }
    Float& operator=(const Real& x);
    //! \brief An approximation by a built-in double-precision floating-point number.
    double get_d() const { return this->v; }
};

template<class R, class A> inline R internal_cast(const A& a) { return static_cast<R>(a); }
template<> inline const double& internal_cast(const Float& x) { return const_cast<const double&>(x.v); }
template<> inline double& internal_cast(const Float& x) { return const_cast<double&>(x.v); }
template<> inline volatile double& internal_cast(const Float& x) { return const_cast<volatile double&>(x.v); }
template<class R, class A> inline R internal_cast(A& a) { return static_cast<R>(a); }
template<> inline const double& internal_cast(Float& x) { return const_cast<const double&>(x.v); }
template<> inline double& internal_cast(Float& x) { return const_cast<double&>(x.v); }
template<> inline volatile double& internal_cast(Float& x) { return const_cast<volatile double&>(x.v); }

template<class R, class A> inline R integer_cast(const A& a);
template<> inline int integer_cast(const Float& a) { return static_cast<int>(a.v); }
template<> inline uint integer_cast(const Float& a) { return static_cast<uint>(a.v); }

template<class R, class A> inline R approx_cast(const A& a);
template<> inline double approx_cast(const Float& a) { return a.v; }


inline std::ostream& operator<<(std::ostream& os, const Float& x) { return os << x.v; }
inline std::istream& operator>>(std::istream& is, Float& x) { double v; is >> v; x=Float(v); return is; }

// Constants related to numerical limits
inline Float mx() { return std::numeric_limits<double>::max(); }


inline Float eps() { return std::numeric_limits<double>::epsilon(); }

// General constants
template<> inline Float pi<Float>() { return 3.1415926535897931; }
template<> inline Float inf<Float>() { return std::numeric_limits<double>::infinity(); }

// Operations for finding nearest representable values
inline Float down(Float x) { return x.v>0 ? x.v*(1-2e-16) : x.v*(1+2e-16); } // Deprecated
inline Float up(Float x) { return x.v>0 ? x.v*(1+2e-16) : x.v*(1-2e-16); } // Deprecated
//! \related Float \brief The next representable value above the given value.
inline Float above(Float x) { return x.v>0 ? x.v*(1-2e-16) : x.v*(1+2e-16); }
//! \related Float \brief The next representable value below the given value.
inline Float below(Float x) { return x.v>0 ? x.v*(1+2e-16) : x.v*(1-2e-16); }

// Discontinuous integer-valued functions
//! \related Float \brief The next lowest integer, represented as a floating-point type.
inline Float floor(Float x) { return std::floor(x.v); }
//! \related Float \brief The next highest integer, represented as a floating-point type.
inline Float ceil(Float x) { return std::ceil(x.v); }

// Non-smooth operations
//! \related Float \brief The magnitude of a floating-point number. Equal to the absolute value.
inline Float mag(Float x) { return std::abs(x.v); }
//! \related Float \brief The absolute value of a floating-point number. Can be computed exactly.
inline Float abs(Float x) { return std::abs(x.v); }
//! \related Float \brief The maximum of two floating-point numbers. Can be computed exactly.
inline Float max(Float x, Float y) { return std::max(x.v,y.v); }
//! \related Float \brief The minimum of two floating-point numbers. Can be computed exactly.
inline Float min(Float x, Float y) { return std::min(x.v,y.v); }

// Standard arithmetic functions
//! \related Float \brief The unary plus function \c +x. Guaranteed to be computed exactly in any rounding mode.
inline Float pos(Float x) { return +x.v; }
//! \related Float \brief The unary negation function \c -x. Guaranteed to be computed exactly in any rounding mode.
inline Float neg(Float x) { return -x.v; }
//! \related Float \brief The square function \c x*x. Also available with \c _rnd, \c _opp, \c _approx, \c _up, \c _down and \c ivl suffixes.
inline Float sqr(Float x) { return x.v*x.v; }
//! \related Float \brief The reciprocal function \c 1/x. Also available with \c _rnd, \c _opp, \c _approx, \c _up, \c _down and \c ivl suffixes.
inline Float rec(Float x) { return 1.0/x.v; }
//! \related Float \brief The binary addition function \c x+y. Also available with \c _rnd, \c _opp, \c _approx, \c _up, \c _down and \c _ivl suffixes.
inline Float add(Float x, Float y) { return x.v+y.v; }
//! \related Float \brief The subtraction function \c x+y. Also available with \c _rnd, \c _opp, \c _approx, \c _up, \c _down and \c _ivl suffixes.
inline Float sub(Float x, Float y) { return x.v-y.v; }
//! \related Float \brief The binary multiplication function \c x+y. Also available with \c _rnd, \c _opp, \c _approx, \c _up, \c _down and \c _ivl suffixes.
inline Float mul(Float x, Float y) { return x.v*y.v; }
//! \related Float \brief The division function \c x+y. Also available with \c _rnd, \c _opp, \c _approx, \c _up, \c _down and \c _ivl suffixes.
inline Float div(Float x, Float y) { return x.v/y.v; }

//! \related Float \brief The positive integer power operator \c x^m.
//! Note that there is no power operator in C++, so the named version must be used. In Python, the power operator is \c x**m.
inline Float pow(Float x, uint n) { return std::pow(x.v,double(n)); }
//! \related Float \brief The integer power operator \c x^n.
//! Also available with \c _rnd, \c _opp, \c _approx, \c _up, \c _down and \c _ivl suffixes.
//! Note that there is no power operator in C++, so the named version must be used. In Python, the power operator is \c x**n.
inline Float pow(Float x, int n) { return std::pow(x.v,double(n)); }

// Standard algebraic and transcendental functions
//! \related Float \brief The square-root function. Not guaranteed to be correctly rounded. Also available with \c _rnd suffix.
inline Float sqrt(Float x) { return std::sqrt(x.v); }
//! \related Float \brief The exponential function. Not guaranteed to be correctly rounded. Also available with \c _rnd suffix.
inline Float exp(Float x) { return std::exp(x.v); }
//! \related Float \brief The logarithm function. Not guaranteed to be correctly rounded. Also available with \c _rnd suffix.
inline Float log(Float x) { return std::log(x.v); }

//! \related Float \brief The sine function. Not guaranteed to be correctly rounded. Also available with \c _rnd suffix.
inline Float sin(Float x) { return std::sin(x.v); }
//! \related Float \brief The cosine function. Not guaranteed to be correctly rounded. Also available with \c _rnd suffix.
inline Float cos(Float x) { return std::cos(x.v); }
//! \related Float \brief The tangent function. Not guaranteed to be correctly rounded. Also available with \c _rnd suffix.
inline Float tan(Float x) { return std::tan(x.v); }
//! \related Float \brief The arcsine function. Not guaranteed to be correctly rounded.
inline Float asin(Float x) { return std::asin(x.v); }
//! \related Float \brief The arccosine function. Not guaranteed to be correctly rounded.
inline Float acos(Float x) { return std::acos(x.v); }
//! \related Float \brief The arctangent function. Not guaranteed to be correctly rounded.
inline Float atan(Float x) { return std::atan(x.v); }


// Correctly rounded arithmetic
inline Float pos_rnd(const Float& x) { return +x.v; }
inline Float neg_rnd(const Float& x) { return -x.v; }
inline Float sqr_rnd(const Float& x) { return x.v*x.v; }
inline Float rec_rnd(const Float& x) { return 1.0/x.v; }
inline Float add_rnd(const Float& x, const Float& y) { return x.v+y.v; }
inline Float sub_rnd(const Float& x, const Float& y) { return x.v-y.v; }
inline Float mul_rnd(const Float& x, const Float& y) { return x.v*y.v; }
inline Float div_rnd(const Float& x, const Float& y) { return x.v/y.v; }
Float pow_rnd(Float x, int n);

// Opposite rounded arithmetic
inline Float pos_opp(const Float& x) { volatile double t=-x.v; return -t; }
inline Float neg_opp(const Float& x) { volatile double t=x.v; return -t; }
inline Float sqr_opp(const Float& x) { volatile double t=(-x.v)*x.v; return -t; }
inline Float rec_opp(const Float& x) { volatile double t=-1.0/x.v; return -t; }
inline Float add_opp(Float x, Float y) { volatile double t=(-x.v)-y.v; return -t; }
inline Float sub_opp(Float x, Float y) { volatile double t=(-x.v)+y.v; return -t; }
inline Float mul_opp(Float x, Float y) { volatile double t=(-x.v)*y.v; return -t; }
inline Float div_opp(Float x, Float y) { volatile double t=(-x.v)/y.v; return -t; }
Float pow_opp(Float x, int n);

// Correctly rounded algebraic and transcendental functions
Float sqrt_rnd(Float x);
Float exp_rnd(Float x);
Float log_rnd(Float x);
Float sin_rnd(Float x);
Float cos_rnd(Float x);
Float tan_rnd(Float x);


// Arithmetic operators
//! \related Float \brief Unary plus (identity) operator. Guaranteed to be exact.
inline Float operator+(const Float& x) { return pos(x); }
//! \related Float \brief Unary negation operator. Guaranteed to be exact.
inline Float operator-(const Float& x) { return neg(x); }
//! \related Float \brief The addition operator. Guaranteed to respect the current rounding mode.
inline Float operator+(const Float& x1, const Float& x2) { return add_rnd(x1,x2); }
//! \related Float \brief The subtraction operator. Guaranteed to respect the current rounding mode.
inline Float operator-(const Float& x1, const Float& x2) { return sub_rnd(x1,x2); }
//! \related Float \brief The multiplication operator. Guaranteed to respect the current rounding mode.
inline Float operator*(const Float& x1, const Float& x2) { return mul_rnd(x1,x2); }
//! \related Float \brief The division operator. Guaranteed to respect the current rounding mode.
inline Float operator/(const Float& x1, const Float& x2) { return div_rnd(x1,x2); }

//! \related Float \brief The in-place addition operator. Guaranteed to respect the current rounding mode.
inline Float& operator+=(Float& x, const Float& y) { x.v+=y.v; return x; }
//! \related Float \brief The in-place subtraction operator. Guaranteed to respect the current rounding mode.
inline Float& operator-=(Float& x, const Float& y) { x.v-=y.v; return x; }
//! \related Float \brief The in-place multiplication operator. Guaranteed to respect the current rounding mode.
inline Float& operator*=(Float& x, const Float& y) { x.v*=y.v; return x; }
//! \related Float \brief The in-place division operator. Guaranteed to respect the current rounding mode.
inline Float& operator/=(Float& x, const Float& y) { x.v/=y.v; return x; }

inline Float operator+(const Float& x1, double x2) { return x1.v+x2; }
inline Float operator+(double x1, const Float& x2) { return x1+x2.v; }
inline Float operator-(const Float& x1, double x2) { return x1.v-x2; }
inline Float operator-(double x1, const Float& x2) { return x1-x2.v; }
inline Float operator*(const Float& x1, double x2) { return x1.v*x2; }
inline Float operator*(double x1, const Float& x2) { return x1*x2.v; }
inline Float operator/(const Float& x1, double x2) { return x1.v/x2; }
inline Float operator/(double x1, const Float& x2) { return x1/x2.v; }

inline Float& operator+=(Float& x1, double x2) { x1.v+=x2; return x1; }
inline Float& operator-=(Float& x1, double x2) { x1.v-=x2; return x1; }
inline Float& operator*=(Float& x1, double x2) { x1.v*=x2; return x1; }
inline Float& operator/=(Float& x1, double x2) { x1.v/=x2; return x1; }

// Comparison operators
//! \related Float \brief The equality operator.
inline bool operator==(const Float& x1, const Float& x2) { return x1.v==x2.v; }
//! \related Float \brief The inequality operator.
inline bool operator!=(const Float& x1, const Float& x2) { return x1.v!=x2.v; }
//! \related Float \brief The less-than-or-equal-to comparison operator.
inline bool operator<=(const Float& x1, const Float& x2) { return x1.v<=x2.v; }
//! \related Float \brief The greater-than-or-equal-to comparison operator.
inline bool operator>=(const Float& x1, const Float& x2) { return x1.v>=x2.v; }
//! \related Float \brief The less-than comparison operator.
inline bool operator< (const Float& x1, const Float& x2) { return x1.v< x2.v; }
//! \related Float \brief The greater-than comparison operator.
inline bool operator> (const Float& x1, const Float& x2) { return x1.v> x2.v; }

inline bool operator==(const Float& x1, double x2) { return x1.v==x2; }
inline bool operator!=(const Float& x1, double x2) { return x1.v!=x2; }
inline bool operator<=(const Float& x1, double x2) { return x1.v<=x2; }
inline bool operator>=(const Float& x1, double x2) { return x1.v>=x2; }
inline bool operator< (const Float& x1, double x2) { return x1.v< x2; }
inline bool operator> (const Float& x1, double x2) { return x1.v> x2; }

inline bool operator==(double x1, const Float& x2) { return x1==x2.v; }
inline bool operator!=(double x1, const Float& x2) { return x1!=x2.v; }
inline bool operator<=(double x1, const Float& x2) { return x1<=x2.v; }
inline bool operator>=(double x1, const Float& x2) { return x1>=x2.v; }
inline bool operator< (double x1, const Float& x2) { return x1< x2.v; }
inline bool operator> (double x1, const Float& x2) { return x1> x2.v; }

#ifdef HAVE_GMPXX_H
inline bool operator==(const Float& x, const Rational& q) { return x.v==q; }
inline bool operator!=(const Float& x, const Rational& q) { return x.v!=q; }
inline bool operator<=(const Float& x, const Rational& q) { return x.v<=q; }
inline bool operator>=(const Float& x, const Rational& q) { return x.v>=q; }
inline bool operator< (const Float& x, const Rational& q) { return x.v< q; }
inline bool operator> (const Float& x, const Rational& q) { return x.v> q; }

inline bool operator==(const Rational& q, const Float& x) { return q==x.v; }
inline bool operator!=(const Rational& q, const Float& x) { return q!=x.v; }
inline bool operator<=(const Rational& q, const Float& x) { return q<=x.v; }
inline bool operator>=(const Rational& q, const Float& x) { return q>=x.v; }
inline bool operator< (const Rational& q, const Float& x) { return q< x.v; }
inline bool operator> (const Rational& q, const Float& x) { return q> x.v; }
#endif



inline Float add_approx(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(to_nearest);
    Float r=add_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float sub_approx(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(to_nearest);
    Float r=sub_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float mul_approx(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(to_nearest);
    Float r=mul_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float div_approx(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(to_nearest);
    Float r=div_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float pow_approx(Float x, int n) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(to_nearest);
    Float r=pow_rnd(x,n); set_rounding_mode(rounding_mode); return r; }

inline Float add_up(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(upward);
    Float r=add_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float sub_up(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(upward);
    Float r=sub_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float mul_up(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(upward);
    Float r=mul_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float div_up(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(upward);
    Float r=div_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float pow_up(Float x, int n) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(upward);
    Float r=pow_rnd(x,n); set_rounding_mode(rounding_mode); return r; }

inline Float add_down(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(downward);
    Float r=add_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float sub_down(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(downward);
    Float r=sub_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float mul_down(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(downward);
    Float r=mul_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float div_down(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(downward);
    Float r=div_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float pow_down(Float x, int n) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(downward);
    Float r=pow_rnd(x,n); set_rounding_mode(rounding_mode); return r; }

//! \related Float \brief The average of two values, computed with nearest rounding. Also available with \c _ivl suffix.
inline Float med_approx(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(to_nearest);
    Float r=add_rnd(x,y)/2; set_rounding_mode(rounding_mode); return r; }
//! \related Float \brief Half of the difference of two values, computed with upward rounding. Also available with \c _ivl suffix.
inline Float rad_up(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(upward);
    Float r=sub_rnd(y,x)/2; set_rounding_mode(rounding_mode); return r; }

inline Interval sqr_ivl(Float x);
inline Interval rec_ivl(Float x);
inline Interval add_ivl(Float x, Float y);
inline Interval sub_ivl(Float x, Float y);
inline Interval mul_ivl(Float x, Float y);
inline Interval div_ivl(Float x, Float y);
inline Interval pow_ivl(Float x, int n);

inline Interval rad_ivl(Float x, Float y);
inline Interval med_ivl(Float x, Float y);



//! \ingroup NumericModule
//! \brief Intervals with floating-point endpoints supporting outwardly-rounded arithmetic.
//! \details
//! Note that <c>%Interval(3.3)</c> yields the singleton interval \f$[3.2999999999999998224,3.2999999999999998224]\f$ (the constant is first interpreted by the C++ compiler to give a C++ \c double, whereas <c>%Interval("3.3")</c> yields the interval \f$[3.2999999999999998224,3.3000000000000002665]\f$ enclosing \f$3.3\f$.
//!
//! Comparison tests on \c Interval use the idea that an interval represents a single number with an unknown value.
//! Hence the result is of type \c tribool, which can take values { \c True, \c False, \c Indeterminate }.
//! Hence a test \f$[l_1,u_1]\leq [l_2,u_2]\f$ returns \c True if \f$u_1\leq u_2\f$, since in this case \f$x_1\leq x_2\f$ whenever \f$x_1\in[l_1,u_2]\f$ and \f$x_2\in[l_2,u_2]\f$, \c False if \f$l_1>u_2\f$, since in this case we know \f$x_1>x_2\f$, and \c Indeterminate otherwise, since in this case we can find \f$x_1,x_2\f$ making the result either true or false.
//! In the case of equality, the comparison \f$[l_1,u_1]\f$==\f$[l_2,u_2]\f$ only returns \c True if both intervals are singletons, since otherwise we can find values making the result either true of false.
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
    typedef Interval NumericType;
  public:
    //! \brief Default constructor yields the singleton zero interval \a [0,0].
    Interval() : l(0.0), u(0.0) { }
    Interval(uint m) : l(m), u(m) { }
    Interval(int n) : l(n), u(n) { }
    //! \brief Convert from a builtin double-precision floating-point value. Yields the singleton interval \a [x,x].
    Interval(double x) : l(x), u(x) { }
    //! \brief Create from a floating-point value. Yields the singleton interval \a [x,x].
    //! Cannot be used in conversions since the \c %Interval class provides stronger accuracy guarantees than the \c %Float class.
    explicit Interval(const Float& x) : l(x), u(x) { }
    //! \brief Copy constructor.
    Interval(const Interval& i) : l(i.l), u(i.u) { }
    //! \brief Convert from a general real number. Yields an interval containing the exact value.
    Interval(const Real& x);

    //! \brief Create from explicitly given lower and upper bounds. Yields the interval \a [lower,upper].
    Interval(double lower, double upper) : l(lower), u(upper) { }
    //! \brief Create from explicitly given lower and upper bounds. Yields the interval \a [lower,upper].
    Interval(const Float& lower, const Float& upper) : l(lower), u(upper) { }
     // ARIADNE_ASSERT_MSG(lower<=upper, "lower = "<<lower<<", upper ="<<upper);
#ifdef HAVE_GMPXX_H
    Interval(const Rational& q);
    Interval& operator=(const Rational& q);
    Interval(const Rational& lower, const Rational& upper);
#endif // HAVE_GMPXX_H

    Interval& operator=(uint m) { l=m; u=m; return *this; }
    Interval& operator=(int n) { l=n; u=n; return *this; }
    Interval& operator=(double c) { l=c; u=c; return *this; }
    Interval& operator=(const Float& x) { l=x; u=x; return *this; }
    Interval& operator=(const Real& x);

    //! \brief The lower bound of the interval.
    const Float& lower() const { return l; }
    //! \brief The upper bound of the interval.
    const Float& upper() const { return u; }
    //! \brief An approximation to the midpoint of the interval.
    const Float midpoint() const { return add_approx(l,u)/2; }
    //! \brief An over-approximation to the radius of the interval.
    const Float radius() const { return sub_up(u,l)/2; }
    //! \brief An over-approximation to the width of the interval.
    const Float width() const { return sub_up(u,l); }

    //! \brief Tests if the interval is empty.
    bool empty() const { return l>u; }
    //! \brief Tests if the interval is a singleton.
    bool singleton() const { return l==u; }

    //! \brief Sets the interval to a "canonical" empty interval \a [1,0].
    void set_empty() { l=1.0; u=0.0; }
    void set_lower(const Float& lower) { l=lower; }
     // ARIADNE_ASSERT(lower<=this->u);
    void set_upper(const Float& upper) { u=upper; }
     // ARIADNE_ASSERT(this->l<=upper);
    void set(const Float& lower, const Float& upper) { l=lower; u=upper; }
     // ARIADNE_ASSERT(lower<=upper);
  public:
    //! \brief Extract a double-precision point approximation to the value represented by the interval.
    double get_d() const { return (this->l.get_d()+this->u.get_d())/2; }
  private:
    Float l, u;
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

//! \related Interval \brief Test if the intervals are equal (as sets).
inline bool equal(Interval i1, Interval i2) {
    //std::cerr<<"equal(i1,i2) with i1="<<i1<<"; i2="<<i2<<std::endl;
    return i1.lower()==i2.lower() && i1.upper()==i2.upper();
}

//! \related Interval \brief Test if the interval is empty.
inline bool empty(Interval i) {
    return i.lower()>i.upper();
}

//! \related Interval \brief Test if the interval is bounded.
inline bool bounded(Interval i) {
    return i.lower()!=-inf<Float>() && i.upper()!=+inf<Float>();
}

//! \related Interval \brief The intersection of two intervals.
inline Interval intersection(Interval i1, Interval i2) {
    return Interval(max(i1.lower(),i2.lower()),min(i1.upper(),i2.upper()));
}

//! \related Interval \brief The hull of two intervals, equal to the smallest interval containing both as subsets.
inline Interval hull(Interval i1, Interval i2) {
    assert(i1.lower()<=i1.upper() && i2.lower()<=i2.upper());
    return Interval(min(i1.lower(),i2.lower()),max(i1.upper(),i2.upper()));
}

//! \related Interval \brief The hull of an interval and a point, equal to the smallest interval containing both.
inline Interval hull(Interval i1, Float x2) {
    return Interval(min(i1.lower(),x2),max(i1.upper(),x2));
}

// An interval one ulp wider
//! \related Interval \brief An interval containing the given interval in its interior.
Interval widen(Interval i);

// Over-approximate by an interval with float coefficients
//! \related Interval \brief Over-approximate the interval by one using builtin single-precision floating-point values as endpoints.
Interval trunc(Interval);
Interval trunc(Interval, uint eps);

//! \related Interval \brief The midpoint of the interval.
inline Float med(Interval i) { return (i.lower()+i.upper())/2; }
//! \related Interval \brief An over-approximation to the radius of the interval.
inline Float rad(Interval i) { return up((i.upper()-i.lower())/2); }
//! \related Interval \brief An over-approximation to the width of the interval.
inline Float diam(Interval i) { return up(i.upper()-i.lower()); }

//! \related Interval \brief The interval of possible maximum values. Yields the interval between \c i1.upper() and \c i2.upper().
inline Interval max(Interval i1,Interval i2);
//! \related Interval \brief The interval of possible minimum values. Yields the interval between \c i1.lower() and \c i2.lower().
inline Interval min(Interval,Interval);
//! \related Interval \brief The interval of possible absolute values. Yields \f$\{ |x| \mid x\in I\}\f$.
inline Interval abs(Interval);

//! \related Interval \brief Unary plus function. Yields the identity \f$I=\{+x | x\in I\}\f$.
inline Interval pos(Interval i);
//! \related Interval \brief Unary negation function. Yields the exact interval \f$\{-x | x\in I\}\f$.
inline Interval neg(Interval i);
//! \related Interval \brief Unary square function. Yields an over-approximation to \f$\{ x^2 \mid x\in I\}\f$.
//! Note that if \a I contains positive and negative values, \c sqr(I) is tighter than \c I*I .
Interval sqr(Interval i);
//! \related Interval \brief Unary reciprocal function. Yields an over-approximation to \f$\{ 1/x \mid x\in I\}\f$.
//! Yields \f$[-\infty,+\infty]\f$ if \a I contains \a 0 in its interior, and an interval containing \f$[1/u,+\infty]\f$ if \a I=[0,u] .
Interval rec(Interval i);

//! \related Interval \brief Binary addition function. Yields an over-approximation to \f$\{ x_1+x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
inline Interval add(Interval, Interval);
//! \related Interval \brief Subtraction function. Yields an over-approximation to \f$\{ x_1-x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
inline Interval sub(Interval, Interval);
//! \related Interval \brief Binary multiplication function. Yields an over-approximation to \f$\{ x_1\times x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
Interval mul(Interval, Interval);
//! \related Interval \brief Division function. Yields an over-approximation to \f$\{ x_1 \div x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
Interval div(Interval, Interval);

inline Interval add(Interval, Float);
inline Interval add(Float, Interval);
inline Interval sub(Interval, Float);
inline Interval sub(Float, Interval);
Interval mul(Interval, Float);
Interval mul(Float,Interval);
Interval div(Interval, Float);
Interval div(Float, Interval);

inline Interval neg_ivl(Float);
inline Interval rec_ivl(Float);
inline Interval add_ivl(Float, Float);
inline Interval sub_ivl(Float, Float);
inline Interval mul_ivl(Float, Float);
inline Interval div_ivl(Float, Float);

//! \related Interval \brief Positive integer power function. Yields an over-approximation to \f$\{ x^m \mid x\in I\}\f$.
Interval pow(Interval i, uint m);
//! \related Interval \brief %Integer power function. Yields an over-approximation to \f$\{ x^n \mid x\in I\}\f$.
Interval pow(Interval i, int n);

//! \related Interval \brief Square-root function. Yields an over-approximation to \f$\{ \sqrt{x} \mid x\in I\}\f$.
//! Requires \c I.lower()>=0 .
Interval sqrt(Interval);
//! \related Interval \brief Exponential function. Yields an over-approximation to \f$\{ \exp{x} \mid x\in I\}\f$.
Interval exp(Interval);
//! \related Interval \brief Natural logarithm function. Yields an over-approximation to \f$\{ \log{x} \mid x\in I\}\f$.
//! Requires \c I.lower()>0 .
Interval log(Interval);

template<> Interval pi<Interval>();
//! \related Interval \brief Sine function. Yields an over-approximation to \f$\{ \sin{x} \mid x\in I\}\f$.
Interval sin(Interval);
//! \related Interval \brief Cosine function. Yields an over-approximation to \f$\{ \sin{x} \mid x\in I\}\f$.
Interval cos(Interval);
//! \related Interval \brief Tangent function. Yields an over-approximation to \f$\{ \sin{x} \mid x\in I\}\f$.
Interval tan(Interval);
Interval asin(Interval);
Interval acos(Interval);
Interval atan(Interval);


//! \related Interval \brief The magnitude of the interval \a I. Yields \f$ \max\{ |x|\,\mid\,x\in I \}\f$.
inline Float mag(Interval i) { return max(abs(i.lower()),abs(i.upper())); }
//! \related Interval \brief The mignitude of the interval \a I. Yields \f$ \min\{ |x|\,\mid\,x\in I \}\f$.
inline Float mig(Interval i) { return min(abs(i.lower()),abs(i.upper())); }

//! \related Interval \brief Test if the interval \a I contains the number \a x.
inline bool contains(Interval i, Float x) { return i.lower()<=x && x<=i.upper(); }

//! \related Interval \brief Test if the interval \a I1 is a subset of \a I2.
inline bool subset(Interval i1, Interval i2) { return i1.lower()>=i2.lower() && i1.upper()<=i2.upper(); }
//! \related Interval \brief Test if the interval \a I1 intersects \a I2. Returns \c true even if the two intervals only have an endpoint in common.
inline bool intersect(Interval i1, Interval i2) { return i1.lower()<=i2.upper() && i1.upper()>=i2.lower(); }
//! \related Interval \brief Test if the interval \a I1 is disjoint from \a I2. Returns \c false even if the two intervals only have an endpoint in common.
inline bool disjoint(Interval i1, Interval i2) { return i1.lower()>i2.upper() || i1.upper()<i2.lower(); }
//! \related Interval \brief Test if the interval \a I1 overlaps \a I2.
//! Returns \c false if the two intervals only have an endpoint in common.
//! Returns \c true if one of the intervals is a singleton in the interior of the other.
inline bool overlap(Interval i1, Interval i2) { return i1.lower()<i2.upper() && i1.upper()>i2.lower(); }
//! \related Interval \brief Test if the (closed) interval \a I1 is a subset of the interior of \a I2.
inline bool inside(Interval i1, Interval i2) { return i1.lower()>i2.lower() && i1.upper()<i2.upper(); }
//! \related Interval \brief Test if the interior of the interval \a I1 is a superset of the (closed) interval \a I2.
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
     return Interval(static_cast<Float>(0.0),max(-i.lower(),i.upper()));
    }
}

inline Interval pos(Interval i)
{
    return Interval(+i.lower(),+i.upper());
}

inline Interval pos_ivl(Float x)
{
    return Interval(+x,+x);
}

inline Interval neg(Interval i)
{
    return Interval(-i.upper(),-i.lower());
}

inline Interval neg_ivl(Float x)
{
    return Interval(-x,-x);
}

inline Interval sqr_ivl(Float x)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& xv=internal_cast<volatile double&>(x);
    set_rounding_mode(downward);
    volatile double rl=xv*xv;
    set_rounding_mode(upward);
    volatile double ru=xv*xv;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

inline Interval rec_ivl(Float x)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& xv=internal_cast<volatile double&>(x);
    set_rounding_mode(downward);
    volatile double rl=1.0/xv;
    set_rounding_mode(upward);
    volatile double ru=1.0/xv;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}



inline Interval add(Interval i1, Interval i2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& i1l=internal_cast<volatile double&>(i1.lower());
    volatile double& i1u=internal_cast<volatile double&>(i1.upper());
    volatile double& i2l=internal_cast<volatile double&>(i2.lower());
    volatile double& i2u=internal_cast<volatile double&>(i2.upper());
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
    volatile double& i1l=internal_cast<volatile double&>(i1.lower());
    volatile double& i1u=internal_cast<volatile double&>(i1.upper());
    volatile double& x2v=internal_cast<volatile double&>(x2);
    set_rounding_mode(downward);
    volatile double rl=i1l+x2v;
    set_rounding_mode(upward);
    volatile double ru=i1u+x2v;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

inline Interval add(Float x1, Interval i2)
{
    return add(i2,x1);
}

inline Interval add_ivl(Float x1, Float x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& x1v=internal_cast<volatile double&>(x1);
    volatile double& x2v=internal_cast<volatile double&>(x2);
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
    volatile double& i1l=internal_cast<volatile double&>(i1.lower());
    volatile double& i1u=internal_cast<volatile double&>(i1.upper());
    volatile double& i2l=internal_cast<volatile double&>(i2.lower());
    volatile double& i2u=internal_cast<volatile double&>(i2.upper());
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
    volatile double& i1l=internal_cast<volatile double&>(i1.lower());
    volatile double& i1u=internal_cast<volatile double&>(i1.upper());
    volatile double& x2v=internal_cast<volatile double&>(x2);
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
    volatile double& x1v=internal_cast<volatile double&>(x1);
    volatile double& i2l=internal_cast<volatile double&>(i2.lower());
    volatile double& i2u=internal_cast<volatile double&>(i2.upper());
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
    volatile double& x1v=internal_cast<volatile double&>(x1);
    volatile double& x2v=internal_cast<volatile double&>(x2);
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
    volatile double& x1v=internal_cast<volatile double&>(x1);
    volatile double& x2v=internal_cast<volatile double&>(x2);
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
    volatile double& x1v=internal_cast<volatile double&>(x1);
    volatile double& x2v=internal_cast<volatile double&>(x2);
    set_rounding_mode(downward);
    volatile double rl=x1v/x2v;
    set_rounding_mode(upward);
    volatile double ru=x1v/x2v;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

inline Interval pow_ivl(Float x1, int n2)
{
    return pow(Interval(x1),n2);
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


//! \related Interval \brief Unary plus operator. Should be implemented exactly and yield \f$\{ +x \mid x\in I\}\f$.
inline Interval operator+(Interval i) { return Interval(i.lower(),i.upper()); }
//! \related Interval \brief Unary negation operator. Should be implemented exactly and yield \f$\{ -x \mid x\in I\}\f$.
inline Interval operator-(Interval i) { return Interval(-i.upper(),-i.lower()); }
//! \related Interval \brief Binary addition operator. Guaranteed to yield an over approximation to \f$\{ x_1+x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
inline Interval operator+(Interval i1, Interval i2) { return add(i1,i2); }
//! \related Interval \brief Binary addition operator. Guaranteed to yield an over approximation to \f$\{ x_1-x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
inline Interval operator-(Interval i1, Interval i2) { return sub(i1,i2); }
//! \related Interval \brief Binary addition operator. Guaranteed to yield an over approximation to \f$\{ x_1*x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
inline Interval operator*(Interval i1, Interval i2) { return mul(i1,i2); }
//! \related Interval \brief Binary addition operator. Guaranteed to yield an over approximation to \f$\{ x_1/x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$. Yields \f$[-\infty,+\infty]\f$ if \f$0\in I_2\f$.
inline Interval operator/(Interval i1, Interval i2) { return div(i1,i2); };

//! \related Interval \brief Inplace addition operator.
inline Interval& operator+=(Interval& i1, Interval i2) { i1=add(i1,i2); return i1; }
//! \related Interval \brief Inplace subtraction operator.
inline Interval& operator-=(Interval& i1, Interval i2) { i1=sub(i1,i2); return i1; }
//! \related Interval \brief Inplace multiplication operator.
inline Interval& operator*=(Interval& i1, Interval i2) { i1=mul(i1,i2); return i1; }
//! \related Interval \brief Inplace division operator.
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

inline Interval operator+(Interval i1, double x2) { return add(i1,static_cast<Float>(x2)); }
inline Interval operator+(double x1, Interval i2) { return add(i2,static_cast<Float>(x1)); }
inline Interval operator-(Interval i1, double x2) { return sub(i1,static_cast<Float>(x2)); }
inline Interval operator-(double x1, Interval i2) { return sub(static_cast<Float>(x1),i2); }
inline Interval operator*(Interval i1, double x2) { return mul(i1,static_cast<Float>(x2)); }
inline Interval operator*(double x1, Interval i2) { return mul(i2,static_cast<Float>(x1)); }
inline Interval operator/(Interval i1, double x2) { return div(i1,static_cast<Float>(x2)); }
inline Interval operator/(double x1, Interval i2) { return div(static_cast<Float>(x1),i2); }

inline Interval& operator+=(Interval& i1, double x2) { i1=add(i1,static_cast<Float>(x2)); return i1; }
inline Interval& operator-=(Interval& i1, double x2) { i1=sub(i1,static_cast<Float>(x2)); return i1; }
inline Interval& operator*=(Interval& i1, double x2) { i1=mul(i1,static_cast<Float>(x2)); return i1; }
inline Interval& operator/=(Interval& i1, double x2) { i1=div(i1,static_cast<Float>(x2)); return i1; }

//inline Interval operator/(Interval i1, int n2) { return div(i1,Float(n2)); }
//inline Interval operator/(Interval i1, double x2) { return div(i1,Float(x2)); }

// Standard equality operators
//! \related Interval \brief Equality operator. Tests equality of intervals as geometric objects, so \c [0,1]==[0,1] returns \c true.
inline bool operator==(const Interval& i1, const Interval& i2) { return i1.lower()==i2.lower() && i1.upper()==i2.upper(); }
//! \related Interval \brief Inequality operator. Tests equality of intervals as geometric objects, so \c [0,2]!=[1,3] returns \c true,
//! even though the intervals possibly represent the same exact real value.
inline bool operator!=(const Interval& i1, const Interval& i2) { return i1.lower()!=i2.lower() || i1.upper()!=i2.upper(); }

// Boost-style tribool (in)equality operators
//inline tribool operator==(const Interval& i1, const Interval& i2) {
//  if(i1.lower()>i2.upper() || i1.upper()<i2.lower()) { return false; } else if(i1.lower()==i2.upper() && i1.upper()==i2.lower()) { return true; } else { return indeterminate; } }
//inline tribool operator!=(const Interval& i1, const Interval& i2) { return !(i1==i2); }

//! \related Interval \brief Equality operator. Tests equality of represented real-point value.
//! Hence \c [0.0,2.0]==1.0 yields \c indeterminate since the interval may represent a real number other than \c 1.0 .
inline tribool operator==(Interval i1, Float x2) {
    if(i1.upper()<x2 || i1.lower()>x2) { return false; }
    else if(i1.lower()==x2 && i1.upper()==x2) { return true; }
    else { return indeterminate; }
}

//! \related Interval \brief Equality operator. Tests equality of represented real-point value.
//! Hence \c [0.0,2.0]!=1.0 yields \c indeterminate since the interval may represent a real number equal to \c 1.0 .
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

inline tribool operator==(Interval i1, double x2) { return i1==static_cast<Float>(x2); }
inline tribool operator!=(Interval i1, double x2) { return i1!=static_cast<Float>(x2); }
inline tribool operator<=(Interval i1, double x2) { return i1<=static_cast<Float>(x2); }
inline tribool operator>=(Interval i1, double x2) { return i1>=static_cast<Float>(x2); }
inline tribool operator< (Interval i1, double x2) { return i1< static_cast<Float>(x2); }
inline tribool operator> (Interval i1, double x2) { return i1> static_cast<Float>(x2); }



//! \related Interval \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
//! Hence \c [1.0,3.0]>[0.0,2.0] yields \c indeterminate since the first interval could represent the number 1.25 and the second 1.75.
inline tribool operator> (Interval i1, Interval i2) {
    if(i1.lower()> i2.upper()) { return true; }
    else if(i1.upper()<=i2.lower()) { return false; }
    else { return indeterminate; }
}

//! \related Interval \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
inline tribool operator< (Interval i1, Interval i2) {
    if(i1.upper()< i2.lower()) { return true; }
    else if(i1.lower()>=i2.upper()) { return false; }
    else { return indeterminate; }
}

//! \related Interval \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
inline tribool operator>=(Interval i1, Interval i2) {
    if(i1.lower()>=i2.upper()) { return true; }
    else if(i1.upper()< i2.lower()) { return false; }
    else { return indeterminate; }
}

//! \related Interval \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
inline tribool operator<=(Interval i1, Interval i2) {
    if(i1.upper()<=i2.lower()) { return true; }
    else if(i1.lower()> i2.upper()) { return false; }
    else { return indeterminate; }
}

template<class A> void serialize(A& a, Interval& ivl, const uint version) {
    a & ivl.lower() & ivl.upper(); }

std::ostream& operator<<(std::ostream&, const Interval&);
std::istream& operator>>(std::istream&, Interval&);


//! \ingroup NumericModule
//! \brief Computable real numbers.
//!
//! Support over-approximation by an Interval and approximation by a Float.
class Real {
    Interval _ivl;
  public:
    typedef Real NumericType;
  public:
    //! Default constructor yields the exact value \a 0.
    Real() : _ivl() { }
    //! \brief Construct from a string literal.
    //! This can be used to create real numbers representing decimal values e.g. \c %Real("4.2") .
    explicit Real(const std::string& s);
    Real(unsigned int m) : _ivl(m) { }
    Real(int n) : _ivl(n) { }
    //! \brief Convert from a builtin double-precision floating-point value.
    //! A numeric literal is first processed by the language support, and the resulting %Real may not have
    //! the same value as the mathematical literal. e.g. \c %Real(4.2) has the value 4.2000000000000002 to 16 decimal places.
    Real(double x) : _ivl(x) { }
#ifdef HAVE_GMPXX_H
    //! \brief Construct from a rational number.
    Real(const Rational& q);
#endif
    //! \brief Construct from a floating-point value.
    explicit Real(const Float& x) : _ivl(x) { }
    //! \brief Construct from a interval. The resulting %Real object does not describe a number arbitrarily accurately.
    //! \deprecated This constructor should be avoided in user code.
    explicit Real(const Interval& ivl) : _ivl(ivl) { }
    explicit Real(double l, double u) : _ivl(l,u) { }
    explicit Real(double l, double x, double u) : _ivl(l,u) { }
    Real& operator=(const double& x) { this->_ivl=x; return *this; }
    Real& operator=(const Float& x) { this->_ivl=x; return *this; }
    Real& operator=(const Interval& x) { this->_ivl=x; return *this; }
    // Can't use conversion operators below in g++ since compiler complains
    // about ambiguous conversion to Interval through Interval(Real::operator Float())
    //operator Float() const { return this->_ivl.midpoint(); }
    //operator Interval() const { return this->_ivl; }
  public:
    //! \brief Get an approximation as a builtin double-precision floating-point number.
    double get_d() const { return this->_ivl.get_d(); }
  private:
    friend Float::Float(const Real&);
    friend Interval::Interval(const Real&);
};

inline Float::Float(const Real& x) : v(x._ivl.midpoint().v) { }
inline Interval::Interval(const Real& x) : l(x._ivl.l), u(x._ivl.u) { }
inline Float& Float::operator=(const Real& x) { *this=Float(x); return *this; }
inline Interval& Interval::operator=(const Real& x) { *this=Interval(x); return *this; }

//@{
//! \related Real \name Arithmetic Operators

//!  \brief Unary plus operator.
inline Real operator+(const Real& x) { return Real(+static_cast<Interval>(x)); }
//!  \brief Unary negation operator.
inline Real operator-(const Real& x) { return Real(-static_cast<Interval>(x)); }
//!  \brief Binary addition operator.
inline Real operator+(const Real& x, const Real& y) { return Real(static_cast<Interval>(x)+static_cast<Interval>(y)); }
//!  \brief Subtraction operator.
inline Real operator-(const Real& x, const Real& y) { return Real(static_cast<Interval>(x)-static_cast<Interval>(y)); }
//! \brief Binary multiplication operator.
inline Real operator*(const Real& x, const Real& y) { return Real(static_cast<Interval>(x)*static_cast<Interval>(y)); }
//! \brief Division operator.
inline Real operator/(const Real& x, const Real& y) { return Real(static_cast<Interval>(x)/static_cast<Interval>(y)); }
//!  \brief Inplace addition operator.
inline Real& operator+=(Real& x, const Real& y) { return x=x+y; }
//!  \brief Inplace subtraction operator.
inline Real& operator-=(Real& x, const Real& y) { return x=x-y; }
//!  \brief Inplace multiplication operator.
inline Real& operator*=(Real& x, const Real& y) { return x=x*y; }
//!  \brief Inplace division operator.
inline Real& operator/=(Real& x, const Real& y) { return x=x/y; }

//@}

inline Real operator+(const Real& x, double y) { return x+static_cast<Real>(y); }
inline Real operator-(const Real& x, double y) { return x-static_cast<Real>(y); }
inline Real operator*(const Real& x, double y) { return x*static_cast<Real>(y); }
inline Real operator/(const Real& x, double y) { return x/static_cast<Real>(y); }
inline Real operator+(double x, const Real& y) { return static_cast<Real>(x)+y; }
inline Real operator-(double x, const Real& y) { return static_cast<Real>(x)-y; }
inline Real operator*(double x, const Real& y) { return static_cast<Real>(x)*y; }
inline Real operator/(double x, const Real& y) { return static_cast<Real>(x)/y; }

inline Float operator+(const Real& x, const Float& y) { return static_cast<Float>(x)+y; }
inline Float operator-(const Real& x, const Float& y) { return static_cast<Float>(x)-y; }
inline Float operator*(const Real& x, const Float& y) { return static_cast<Float>(x)*y; }
inline Float operator/(const Real& x, const Float& y) { return static_cast<Float>(x)/y; }
inline Float operator+(const Float& x, const Real& y) { return x+static_cast<Float>(y); }
inline Float operator-(const Float& x, const Real& y) { return x-static_cast<Float>(y); }
inline Float operator*(const Float& x, const Real& y) { return x*static_cast<Float>(y); }
inline Float operator/(const Float& x, const Real& y) { return x/static_cast<Float>(y); }
inline Float& operator+=(Float& x, const Real& y) { return x+=static_cast<Float>(y); }
inline Float& operator-=(Float& x, const Real& y) { return x-=static_cast<Float>(y); }
inline Float& operator*=(Float& x, const Real& y) { return x*=static_cast<Float>(y); }
inline Float& operator/=(Float& x, const Real& y) { return x/=static_cast<Float>(y); }

inline Interval operator+(const Real& x, const Interval& y) { return static_cast<Interval>(x)+y; }
inline Interval operator-(const Real& x, const Interval& y) { return static_cast<Interval>(x)-y; }
inline Interval operator*(const Real& x, const Interval& y) { return static_cast<Interval>(x)*y; }
inline Interval operator/(const Real& x, const Interval& y) { return static_cast<Interval>(x)/y; }
inline Interval operator+(const Interval& x, const Real& y) { return x+static_cast<Interval>(y); }
inline Interval operator-(const Interval& x, const Real& y) { return x-static_cast<Interval>(y); }
inline Interval operator*(const Interval& x, const Real& y) { return x*static_cast<Interval>(y); }
inline Interval operator/(const Interval& x, const Real& y) { return x/static_cast<Interval>(y); }
inline Interval& operator+=(Interval& x, const Real& y) { return x+=static_cast<Interval>(y); }
inline Interval& operator-=(Interval& x, const Real& y) { return x-=static_cast<Interval>(y); }
inline Interval& operator*=(Interval& x, const Real& y) { return x*=static_cast<Interval>(y); }
inline Interval& operator/=(Interval& x, const Real& y) { return x/=static_cast<Interval>(y); }

//@{
//! \related Real \name Comparison Operators

//! \brief Equality operator.
//! Returns \c true if the numbers have the same representation or can be evaluated exactly to the same value.
//! Returns \c false if the numbers can be evalated to different values.
//! Returns \c indeterminate if the numbers cannot be shown to be the same or different. Implementation dependent.
inline tribool operator==(const Real& x, const Real& y) { return static_cast<Interval>(x)==static_cast<Interval>(y); }
//! \brief Inequality operator.
//! Returns \c true if the numbers have been proved to be different, such as by evaluation to different values.
//! Returns \c false if the numbers have the same representation or can be evaluated exactly to the same value.
//! Returns \c indeterminate if the numbers cannot be shown to be the same or different. Implementation dependent.
inline tribool operator!=(const Real& x, const Real& y) { return static_cast<Interval>(x)!=static_cast<Interval>(y); }
//! \brief Greater-than-or-equal-to comparison operator.
inline tribool operator>=(const Real& x, const Real& y) { return static_cast<Interval>(x)>=static_cast<Interval>(y); }
//! \brief Less-than-or-equal-to comparison operator.
inline tribool operator<=(const Real& x, const Real& y) { return static_cast<Interval>(x)<=static_cast<Interval>(y); }
//! \brief Strictly-greater-than comparison operator.
inline tribool operator> (const Real& x, const Real& y) { return static_cast<Interval>(x)> static_cast<Interval>(y); }
//! \brief Strictly-less-than comparison operator.
inline tribool operator< (const Real& x, const Real& y) { return static_cast<Interval>(x)< static_cast<Interval>(y); }
//@}

inline tribool operator==(const Real& x, double y) { return static_cast<Interval>(x)==static_cast<Interval>(y); }
inline tribool operator!=(const Real& x, double y) { return static_cast<Interval>(x)!=static_cast<Interval>(y); }
inline tribool operator>=(const Real& x, double y) { return static_cast<Interval>(x)>=static_cast<Interval>(y); }
inline tribool operator<=(const Real& x, double y) { return static_cast<Interval>(x)<=static_cast<Interval>(y); }
inline tribool operator> (const Real& x, double y) { return static_cast<Interval>(x)> static_cast<Interval>(y); }
inline tribool operator< (const Real& x, double y) { return static_cast<Interval>(x)< static_cast<Interval>(y); }

//@{
//! \related Real \name Arithmetical, algebraic and transcendental functions

//!  \brief The absolute value function \c |x|.
inline Real abs(const Real& x) { return Real(abs(static_cast<Interval>(x))); }
//!  \brief The unary plus function \c +x.
inline Real pos(const Real& x) { return Real(pos(static_cast<Interval>(x))); }
//!  \brief The unary negation function \c -x.
inline Real neg(const Real& x) { return Real(neg(static_cast<Interval>(x))); }
//!  \brief The reciprocal function \c 1/x.
inline Real sqr(const Real& x) { return Real(sqr(static_cast<Interval>(x))); }
//!  \brief The reciprocal function \c 1/x.
inline Real rec(const Real& x) { return Real(rec(static_cast<Interval>(x))); }
//!  \brief The binary addition function \c x+y.
inline Real add(const Real& x, const Real& y) { return Real(add(static_cast<Interval>(x),static_cast<Interval>(y))); }
//!  \brief The subtraction function \c x-y.
inline Real sub(const Real& x, const Real& y) { return Real(sub(static_cast<Interval>(x),static_cast<Interval>(y))); }
//!  \brief The binary multiplication function \c x*y.
inline Real mul(const Real& x, const Real& y) { return Real(mul(static_cast<Interval>(x),static_cast<Interval>(y))); }
//!  \brief The division function \c x/y.
inline Real div(const Real& x, const Real& y) { return Real(div(static_cast<Interval>(x),static_cast<Interval>(y))); }
//!  \brief The integer power function \c x^n.
inline Real pow(const Real& x, int n) { return Real(pow(static_cast<Interval>(x),n)); }
//!  \brief The square-root function.
inline Real sqrt(const Real& x) { return Real(sqrt(static_cast<Interval>(x))); }
//!  \brief The exponential function.
inline Real exp(const Real& x) { return Real(exp(static_cast<Interval>(x))); }
//!  \brief The natural logarithm function.
inline Real log(const Real& x) { return Real(log(static_cast<Interval>(x))); }
//!  \brief The constant \a pi.
template<> inline Real pi<>() { return Real(pi<Interval>()); }
//!  \brief The sine function.
inline Real sin(const Real& x) { return Real(sin(static_cast<Interval>(x))); }
//!  \brief The cosine function.
inline Real cos(const Real& x) { return Real(cos(static_cast<Interval>(x))); }
//!  \brief The tangent function.
inline Real tan(const Real& x) { return Real(tan(static_cast<Interval>(x))); }

//@}

//@{
//! \related Real \name Input/output operators

//! \related Real \brief Write to an output stream.
std::ostream& operator<<(std::ostream& os, const Real& x);
//@}

template<class R, class A> inline R numeric_cast(const A& a);
template<> inline double numeric_cast(const Real& a) { return a.get_d(); }
template<> inline Real numeric_cast(const Float& a) { return Real(a); }
template<> inline Real numeric_cast(const Interval& a) { return Real(a); }

//! \ingroup NumericModule \related Float \related Interval \related Real
//! \brief Cast one %Ariadne numerical type or builtin numerical type to another.
template<class R, class A> inline R numeric_cast(const A& a) { return R(a); }
template<> inline int numeric_cast(const Float& a) { return int(a.get_d()); }
template<> inline double numeric_cast(const Float& a) { return a.get_d(); }
template<> inline double numeric_cast(const Interval& a) { return a.get_d(); }
template<> inline Float numeric_cast(const Interval& a) { return a.midpoint(); }
template<> inline Interval numeric_cast(const Float& a) { return Interval(a); }

//! \ingroup NumericModule \related Float 
//! \brief Converts \a e to an object of type \a X, which may either be an
//! \c Float or \c Interval, with the semantics that \a e denotes and error bound.
//! Returns the float 0.0 (since floating-point computations do not keep track of errors)
//! and the interval [-e,+e].
template<class X> inline X convert_error(const Float& e);
template<> inline Float convert_error<Float>(const Float& e) { return 0.0; }
template<> inline Interval convert_error<Interval>(const Float& e) { return Interval(-e,+e); }

} // namespace Ariadne

#endif
