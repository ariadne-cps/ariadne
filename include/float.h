/***************************************************************************
 *            float.h
 *
 *  Copyright 2008-10  Pieter Collins
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

/*! \file float.h
 *  \brief Floating-point number class.
 */
#ifndef ARIADNE_FLOAT_H
#define ARIADNE_FLOAT_H

#include <iostream> // For std::floor std::ceil etc
#include <cmath> // For std::floor std::ceil etc
#include <algorithm> // For std::max, std::min
#include <limits> // For std::numeric_limits<double>

#include "rounding.h"
#include "rational.h"


// Simplifying typedef for unsigned integer type
typedef unsigned int uint;

namespace Ariadne {

class Real;

template<class X> X pi();
template<class X> X inf();

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

class Interval;
inline Interval sqr_ivl(Float x);
inline Interval rec_ivl(Float x);
inline Interval add_ivl(Float x, Float y);
inline Interval sub_ivl(Float x, Float y);
inline Interval mul_ivl(Float x, Float y);
inline Interval div_ivl(Float x, Float y);
inline Interval pow_ivl(Float x, int n);

inline Interval rad_ivl(Float x, Float y);
inline Interval med_ivl(Float x, Float y);


#ifdef HAVE_GMPXX_H
inline bool operator==(const Float& x, const Rational& q) { return x.get_d()==static_cast<const mpq_class&>(q); }
inline bool operator!=(const Float& x, const Rational& q) { return x.get_d()!=static_cast<const mpq_class&>(q); }
inline bool operator<=(const Float& x, const Rational& q) { return x.get_d()<=static_cast<const mpq_class&>(q); }
inline bool operator>=(const Float& x, const Rational& q) { return x.get_d()>=static_cast<const mpq_class&>(q); }
inline bool operator< (const Float& x, const Rational& q) { return x.get_d()< static_cast<const mpq_class&>(q); }
inline bool operator> (const Float& x, const Rational& q) { return x.get_d()> static_cast<const mpq_class&>(q); }

inline bool operator==(const Rational& q, const Float& x) { return static_cast<mpq_class>(q)==x.get_d(); }
inline bool operator!=(const Rational& q, const Float& x) { return static_cast<mpq_class>(q)!=x.get_d(); }
inline bool operator<=(const Rational& q, const Float& x) { return static_cast<mpq_class>(q)<=x.get_d(); }
inline bool operator>=(const Rational& q, const Float& x) { return static_cast<mpq_class>(q)>=x.get_d(); }
inline bool operator< (const Rational& q, const Float& x) { return static_cast<mpq_class>(q)< x.get_d(); }
inline bool operator> (const Rational& q, const Float& x) { return static_cast<mpq_class>(q)> x.get_d(); }
#endif // HAVE_GMPXX_H

} // namespace Ariadne

#endif
