/***************************************************************************
 *            float64.h
 *
 *  Copyright 2008-15  Pieter Collins
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

/*! \file float64.h
 *  \brief Raw floating-point number class based on double-precision floats.
 */

#ifndef ARIADNE_FLOAT64_H
#define ARIADNE_FLOAT64_H

#include <iostream> // For std::floor std::ceil etc
#include <iomanip> // For std::setprecision
#include <cmath> // For std::floor std::ceil etc
#include <algorithm> // For std::max, std::min
#include <limits> // For std::numeric_limits<double>

#include "utility/declarations.h"
#include "numeric/rounding.h"
#include "numeric/number.decl.h"
#include "numeric/float.decl.h"

namespace Ariadne {

class Float64;
typedef Float64 RawFloat64;

class Real;
class Dyadic;
class Decimal;

using std::min;
using std::max;

class Precision64 {
};

//! \ingroup NumericModule
//! \brief Floating point numbers (double precision) using approxiamate arithmetic.
//! \details
//! The \c %Float64 class represents floating-point numbers.
//! Unless otherwise mentioned, operations on floating-point numbers are performed approximately, with no guarantees
//! on the output.
//!
//! To implement <em>rounded arithmetic</em>, arithmetical operations of \c %Float64 can be performed with guaranteed rounding by
//! specifying \c _up and \c _down suffixes to arithmetical functions \c add, \c sub, \c mul and \c div.
//! Additionally, operations can be performed in the current <em>rounding mode</em> by using the \c _rnd suffix,
//! or with rounding reversed using the \c _opp suffix.
//! The \c _approx suffix is provided to specifically indicate that the operation is computed approximately.
//!
//! %Ariadne floating-point numbers can be constructed by conversion from built-in C++ types.
//! Note that the value of a built-in floating-point value may differ from the mathematical value of the literal.
//! For example, while <c>%Float64(3.25)</c> is represented exactly, <c>%Float64(3.3)</c> has a value of \f$3.2999999999999998224\ldots\f$.
//! \note In the future, the construction of a \c %Float64 from a string literal may be supported.
//! \sa Real, ExactFloat, ValidatedFloat, UpperFloat, LowerFloat, ApproximateFloat
class Float64 {
  public:
    volatile double dbl;
  public:
    typedef Raw Paradigm;
    typedef Float64 NumericType;
    typedef rounding_mode_t RoundingModeType;
  public:
    static const RoundingModeType downward;
    static const RoundingModeType upward;
    static const RoundingModeType to_nearest;
    static const RoundingModeType toward_zero;

    static RoundingModeType get_rounding_mode();
    static Void set_rounding_mode(RoundingModeType);
    static Void set_rounding_downward();
    static Void set_rounding_upward();
    static Void set_rounding_to_nearest();
    static Void set_rounding_toward_zero();
  public:
    //! \brief Default constructor creates an uninitialised number.
    Float64() : dbl() { }
    //! \brief Convert from a built-in double-precision floating-point number.
    Float64(double x) : dbl(x) { }
    //! \brief Copy constructor.
    Float64(const Float64& x) : dbl(x.dbl) { }
    //! \brief Construct from a rational number with given rounding
    explicit Float64(const Rational& q, RoundingModeType rnd);
    explicit operator volatile double& () { return const_cast<volatile double&>(dbl); }
    explicit operator const double& () const { return const_cast<const double&>(dbl); }
    Float64 const& raw() const { return *this; } // FIXME: Included for compatibility with user floats
    //! \brief An approximation by a built-in double-precision floating-point number.
    double get_d() const { return this->dbl; }
};

const Float64 inf = std::numeric_limits<double>::infinity();
const Float64 nan = (0.0/0.0);

template<class R, class A> inline R internal_cast(const A& a) { return static_cast<R>(a); }
template<> inline const double& internal_cast(const Float64& x) { return const_cast<const double&>(x.dbl); }
template<> inline double& internal_cast(const Float64& x) { return const_cast<double&>(x.dbl); }
template<> inline volatile double& internal_cast(const Float64& x) { return const_cast<volatile double&>(x.dbl); }
template<class R, class A> inline R internal_cast(A& a) { return static_cast<R>(a); }
template<> inline const double& internal_cast(Float64& x) { return const_cast<const double&>(x.dbl); }
template<> inline double& internal_cast(Float64& x) { return const_cast<double&>(x.dbl); }
template<> inline volatile double& internal_cast(Float64& x) { return const_cast<volatile double&>(x.dbl); }

template<class R, class A> inline R integer_cast(const A& a);
template<> inline Int integer_cast(const Float64& a) { return static_cast<Int>(a.dbl); }
template<> inline Nat integer_cast(const Float64& a) { return static_cast<Nat>(a.dbl); }

inline OutputStream& operator<<(OutputStream& os, const Float64& x) { return os << x.dbl; }
inline InputStream& operator>>(InputStream& is, Float64& x) { double dbl; is >> dbl; x=Float64(dbl); return is; }

inline Float64::RoundingModeType Float64::get_rounding_mode() { return Ariadne::get_rounding_mode(); }
inline Void Float64::set_rounding_mode(RoundingModeType rnd) { Ariadne::set_rounding_mode(rnd); }
inline Void Float64::set_rounding_upward() { Ariadne::set_rounding_upward(); }
inline Void Float64::set_rounding_downward() { Ariadne::set_rounding_downward(); }
inline Void Float64::set_rounding_to_nearest() { Ariadne::set_rounding_to_nearest(); }
inline Void Float64::set_rounding_toward_zero() { Ariadne::set_rounding_toward_zero(); }

// Exact raw data operations
inline Float64 operator+(Float64 x) { return +x.dbl; }
inline Float64 operator-(Float64 x) { return -x.dbl; }
inline Float64 pos(Float64 x) { return +x.dbl; }
inline Float64 neg(Float64 x) { return -x.dbl; }
inline Float64 abs(Float64 x) { return std::fabs(x.dbl); }
inline Float64 mag(Float64 x) { return std::fabs(x.dbl); }
inline Float64 half(Float64 x) { return x.dbl/2; }
inline Float64 max(Float64 x1, Float64 x2) { return std::max(x1,x2); }
inline Float64 min(Float64 x1, Float64 x2) { return std::min(x1,x2); }

// Raw data comparisons
inline Bool operator==(Float64 const& x1, Float64 const& x2) { return x1.dbl == x2.dbl; }
inline Bool operator!=(Float64 const& x1, Float64 const& x2) { return x1.dbl != x2.dbl; }
inline Bool operator<=(Float64 const& x1, Float64 const& x2) { return x1.dbl <= x2.dbl; }
inline Bool operator>=(Float64 const& x1, Float64 const& x2) { return x1.dbl >= x2.dbl; }
inline Bool operator< (Float64 const& x1, Float64 const& x2) { return x1.dbl <  x2.dbl; }
inline Bool operator> (Float64 const& x1, Float64 const& x2) { return x1.dbl >  x2.dbl; }

// Constants related to numerical limits
inline Float64 mx() { return std::numeric_limits<double>::max(); }
inline Float64 eps() { return std::numeric_limits<double>::epsilon(); }

// Checking whether a Float64 is not-a-number
inline Bool isnan(const Float64& x) { return std::isnan(x.dbl); }

// Operations for finding nearest representable values
inline Float64 down(Float64 x) { return x.dbl>0 ? x.dbl*(1-2e-16) : x.dbl*(1+2e-16); } // Deprecated
inline Float64 up(Float64 x) { return x.dbl>0 ? x.dbl*(1+2e-16) : x.dbl*(1-2e-16); } // Deprecated
//! \related Float64 \brief The next representable value above the given value.
inline Float64 above(Float64 x) { return x.dbl>0 ? x.dbl*(1-2e-16) : x.dbl*(1+2e-16); }
//! \related Float64 \brief The next representable value below the given value.
inline Float64 below(Float64 x) { return x.dbl>0 ? x.dbl*(1+2e-16) : x.dbl*(1-2e-16); }

// Discontinuous integer-valued functions
//! \related Float64 \brief The next lowest integer, represented as a floating-point type.
inline Float64 floor(Float64 x) { return std::floor(x.dbl); }
//! \related Float64 \brief The next highest integer, represented as a floating-point type.
inline Float64 ceil(Float64 x) { return std::ceil(x.dbl); }

inline Float64 pos_exact(Float64 x) { return +x.dbl; }
inline Float64 neg_exact(Float64 x) { return -x.dbl; }
inline Float64 half_exact(Float64 x) { return x.dbl/2; }

inline Float64 abs_exact(Float64 x) { return std::fabs(x.dbl); }
inline Float64 mag_exact(Float64 x) { return std::fabs(x.dbl); }

// Correctly rounded arithmetic
inline Float64 pos_rnd(const Float64& x) { return +x.dbl; }
inline Float64 neg_rnd(const Float64& x) { return -x.dbl; }
inline Float64 sqr_rnd(const Float64& x) { return x.dbl*x.dbl; }
inline Float64 rec_rnd(const Float64& x) { return 1.0/x.dbl; }
inline Float64 add_rnd(Float64 x, Float64 y) { return x.dbl+y.dbl; }
inline Float64 sub_rnd(Float64 x, Float64 y) { return x.dbl-y.dbl; }
inline Float64 mul_rnd(Float64 x, Float64 y) { return x.dbl*y.dbl; }
inline Float64 div_rnd(Float64 x, Float64 y) { return x.dbl/y.dbl; }

Float64 pow_rnd(Float64 x, Int n);

// Opposite rounded arithmetic
inline Float64 pos_opp(const Float64& x) { volatile double t=-x.dbl; return -t; }
inline Float64 neg_opp(const Float64& x) { volatile double t=x.dbl; return -t; }
inline Float64 sqr_opp(const Float64& x) { volatile double t=-x.dbl; t=t*x.dbl; return -t; }
inline Float64 rec_opp(const Float64& x) { volatile double t=-1.0/x.dbl; return -t; }
inline Float64 add_opp(Float64 x, Float64 y) { volatile double t=-x.dbl; t=t-y.dbl; return -t; }
inline Float64 sub_opp(Float64 x, Float64 y) { volatile double t=-x.dbl; t=t+y.dbl; return -t; }
inline Float64 mul_opp(Float64 x, Float64 y) { volatile double t=-x.dbl; t=t*y.dbl; return -t; }
inline Float64 div_opp(Float64 x, Float64 y) { volatile double t=x.dbl; t=t/y.dbl; return -t; }
Float64 pow_opp(Float64 x, Int n);

// Correctly rounded algebraic and transcendental functions
Float64 sqrt_rnd(Float64 x);
Float64 exp_rnd(Float64 x);
Float64 log_rnd(Float64 x);
Float64 sin_rnd(Float64 x);
Float64 cos_rnd(Float64 x);
Float64 tan_rnd(Float64 x);


inline Float64 next_down(Float64 x) {
    rounding_mode_t rounding_mode=Float64::get_rounding_mode(); Float64::set_rounding_mode(downward);
    Float64 r=add_rnd(x,Float64(-std::numeric_limits<double>::min())); Float64::set_rounding_mode(rounding_mode); return r; }
inline Float64 next_up(Float64 x) {
    rounding_mode_t rounding_mode=Float64::get_rounding_mode(); Float64::set_rounding_mode(upward);
    Float64 r=add_rnd(x,Float64(std::numeric_limits<double>::min())); Float64::set_rounding_mode(rounding_mode); return r; }

inline Float64 add_near(Float64 x, Float64 y) {
    rounding_mode_t rounding_mode=Float64::get_rounding_mode(); Float64::set_rounding_mode(to_nearest);
    Float64 r=add_rnd(x,y); Float64::set_rounding_mode(rounding_mode); return r; }
inline Float64 sub_near(Float64 x, Float64 y) {
    rounding_mode_t rounding_mode=Float64::get_rounding_mode(); Float64::set_rounding_mode(to_nearest);
    Float64 r=sub_rnd(x,y); Float64::set_rounding_mode(rounding_mode); return r; }
inline Float64 mul_near(Float64 x, Float64 y) {
    rounding_mode_t rounding_mode=Float64::get_rounding_mode(); Float64::set_rounding_mode(to_nearest);
    Float64 r=mul_rnd(x,y); Float64::set_rounding_mode(rounding_mode); return r; }
inline Float64 div_near(Float64 x, Float64 y) {
    rounding_mode_t rounding_mode=Float64::get_rounding_mode(); Float64::set_rounding_mode(to_nearest);
    Float64 r=div_rnd(x,y); Float64::set_rounding_mode(rounding_mode); return r; }


//! \related Float64 \brief The nearest floating-point approximation to the constant \a pi.
static Float64 pi_approx(Precision64) { return Float64(3.1415926535897931); }
static Float64 pi_near(Precision64) { return Float64(3.1415926535897931); }
static Float64 pi_down(Precision64) { return Float64(3.1415926535897931); }
static Float64 pi_up(Precision64) { return Float64(3.1415926535897936); }
static Float64 pi_approx() { return Float64(3.1415926535897931); }
static Float64 pi_near() { return Float64(3.1415926535897931); }
static Float64 pi_down() { return Float64(3.1415926535897931); }
static Float64 pi_up() { return Float64(3.1415926535897936); }

inline Float64 add_approx(Float64 x, Float64 y) {
    rounding_mode_t rounding_mode=Float64::get_rounding_mode(); Float64::set_rounding_mode(to_nearest);
    Float64 r=add_rnd(x,y); Float64::set_rounding_mode(rounding_mode); return r; }
inline Float64 sub_approx(Float64 x, Float64 y) {
    rounding_mode_t rounding_mode=Float64::get_rounding_mode(); Float64::set_rounding_mode(to_nearest);
    Float64 r=sub_rnd(x,y); Float64::set_rounding_mode(rounding_mode); return r; }
inline Float64 mul_approx(Float64 x, Float64 y) {
    rounding_mode_t rounding_mode=Float64::get_rounding_mode(); Float64::set_rounding_mode(to_nearest);
    Float64 r=mul_rnd(x,y); Float64::set_rounding_mode(rounding_mode); return r; }
inline Float64 div_approx(Float64 x, Float64 y) {
    rounding_mode_t rounding_mode=Float64::get_rounding_mode(); Float64::set_rounding_mode(to_nearest);
    Float64 r=div_rnd(x,y); Float64::set_rounding_mode(rounding_mode); return r; }
inline Float64 pow_approx(Float64 x, Int n) {
    rounding_mode_t rounding_mode=Float64::get_rounding_mode(); Float64::set_rounding_mode(to_nearest);
    Float64 r=pow_rnd(x,n); Float64::set_rounding_mode(rounding_mode); return r; }

inline Float64 add_up(Float64 x, Float64 y) {
    rounding_mode_t rounding_mode=Float64::get_rounding_mode(); Float64::set_rounding_mode(upward);
    Float64 r=add_rnd(x,y); Float64::set_rounding_mode(rounding_mode); return r; }
inline Float64 sub_up(Float64 x, Float64 y) {
    rounding_mode_t rounding_mode=Float64::get_rounding_mode(); Float64::set_rounding_mode(upward);
    Float64 r=sub_rnd(x,y); Float64::set_rounding_mode(rounding_mode); return r; }
inline Float64 mul_up(Float64 x, Float64 y) {
    rounding_mode_t rounding_mode=Float64::get_rounding_mode(); Float64::set_rounding_mode(upward);
    Float64 r=mul_rnd(x,y); Float64::set_rounding_mode(rounding_mode); return r; }
inline Float64 div_up(Float64 x, Float64 y) {
    rounding_mode_t rounding_mode=Float64::get_rounding_mode(); Float64::set_rounding_mode(upward);
    Float64 r=div_rnd(x,y); Float64::set_rounding_mode(rounding_mode); return r; }
inline Float64 pow_up(Float64 x, Int n) {
    rounding_mode_t rounding_mode=Float64::get_rounding_mode(); Float64::set_rounding_mode(upward);
    Float64 r=pow_rnd(x,n); Float64::set_rounding_mode(rounding_mode); return r; }

inline Float64 add_down(Float64 x, Float64 y) {
    rounding_mode_t rounding_mode=Float64::get_rounding_mode(); Float64::set_rounding_mode(downward);
    Float64 r=add_rnd(x,y); Float64::set_rounding_mode(rounding_mode); return r; }
inline Float64 sub_down(Float64 x, Float64 y) {
    rounding_mode_t rounding_mode=Float64::get_rounding_mode(); Float64::set_rounding_mode(downward);
    Float64 r=sub_rnd(x,y); Float64::set_rounding_mode(rounding_mode); return r; }
inline Float64 mul_down(Float64 x, Float64 y) {
    rounding_mode_t rounding_mode=Float64::get_rounding_mode(); Float64::set_rounding_mode(downward);
    Float64 r=mul_rnd(x,y); Float64::set_rounding_mode(rounding_mode); return r; }
inline Float64 div_down(Float64 x, Float64 y) {
    rounding_mode_t rounding_mode=Float64::get_rounding_mode(); Float64::set_rounding_mode(downward);
    Float64 r=div_rnd(x,y); Float64::set_rounding_mode(rounding_mode); return r; }
inline Float64 pow_down(Float64 x, Int n) {
    rounding_mode_t rounding_mode=Float64::get_rounding_mode(); Float64::set_rounding_mode(downward);
    Float64 r=pow_rnd(x,n); Float64::set_rounding_mode(rounding_mode); return r; }

//! \related Float64 \brief The average of two values, computed with nearest rounding. Also available with \c _ivl suffix.
inline Float64 med_approx(Float64 x, Float64 y) {
    rounding_mode_t rounding_mode=Float64::get_rounding_mode(); Float64::set_rounding_mode(to_nearest);
    Float64 r=half_exact(add_rnd(x,y)); Float64::set_rounding_mode(rounding_mode); return r; }
//! \related Float64 \brief Half of the difference of two values, computed with upward rounding. Also available with \c _ivl suffix.
inline Float64 rad_up(Float64 x, Float64 y) {
    rounding_mode_t rounding_mode=Float64::get_rounding_mode(); Float64::set_rounding_mode(upward);
    Float64 r=half_exact(sub_rnd(y,x)); Float64::set_rounding_mode(rounding_mode); return r; }

inline Float64 sqrt_approx(Float64 x) { return std::sqrt(x.dbl); }
inline Float64 exp_approx(Float64 x) { return std::exp(x.dbl); }
inline Float64 log_approx(Float64 x) { return std::log(x.dbl); }
inline Float64 sin_approx(Float64 x) { return std::sin(x.dbl); }
inline Float64 cos_approx(Float64 x) { return std::cos(x.dbl); }
inline Float64 tan_approx(Float64 x) { return std::tan(x.dbl); }
inline Float64 asin_approx(Float64 x) { return std::asin(x.dbl); }
inline Float64 acos_approx(Float64 x) { return std::acos(x.dbl); }
inline Float64 atan_approx(Float64 x) { return std::atan(x.dbl); }



// Deprecated approximate operations
inline Float64 operator+(Float64 x1, Float64 x2) { return x1.dbl+x2.dbl; }
inline Float64 operator-(Float64 x1, Float64 x2) { return x1.dbl-x2.dbl; }
inline Float64 operator*(Float64 x1, Float64 x2) { return x1.dbl*x2.dbl; }
inline Float64 operator/(Float64 x1, Float64 x2) { return x1.dbl/x2.dbl; }
inline Float64& operator+=(Float64& x1, Float64 x2) { x1.dbl+=x2.dbl; return x1; }
inline Float64& operator-=(Float64& x1, Float64 x2) { x1.dbl-=x2.dbl; return x1; }
inline Float64& operator*=(Float64& x1, Float64 x2) { x1.dbl*=x2.dbl; return x1; }
inline Float64& operator/=(Float64& x1, Float64 x2) { x1.dbl/=x2.dbl; return x1; }

inline Float64 add(Float64 x1, Float64 x2) { return x1.dbl+x2.dbl; }
inline Float64 sub(Float64 x1, Float64 x2) { return x1.dbl-x2.dbl; }
inline Float64 mul(Float64 x1, Float64 x2) { return x1.dbl*x2.dbl; }
inline Float64 div(Float64 x1, Float64 x2) { return x1.dbl/x2.dbl; }
inline Float64 pow(Float64 x, Nat n) { return std::pow(x.dbl,double(n)); }
inline Float64 pow(Float64 x, Int n) { return std::pow(x.dbl,double(n)); }
inline Float64 sqr(Float64 x) { return x.dbl * x.dbl; }
inline Float64 rec(Float64 x) { return 1.0/x.dbl; }
inline Float64 sqrt(Float64 x) { return std::sqrt(x.dbl); }
inline Float64 exp(Float64 x) { return std::exp(x.dbl); }
inline Float64 log(Float64 x) { return std::log(x.dbl); }
inline Float64 sin(Float64 x) { return std::sin(x.dbl); }
inline Float64 cos(Float64 x) { return std::cos(x.dbl); }
inline Float64 tan(Float64 x) { return std::tan(x.dbl); }
inline Float64 asin(Float64 x) { return std::asin(x.dbl); }
inline Float64 acos(Float64 x) { return std::acos(x.dbl); }
inline Float64 atan(Float64 x) { return std::atan(x.dbl); }

// Deprecated mixed operations
inline Float64 operator+(Float64 x1, double x2) { return x1.dbl+x2; }
inline Float64 operator-(Float64 x1, double x2) { return x1.dbl-x2; }
inline Float64 operator*(Float64 x1, double x2) { return x1.dbl*x2; }
inline Float64 operator/(Float64 x1, double x2) { return x1.dbl/x2; }
inline Float64 operator+(double x1, Float64 x2) { return x1+x2.dbl; }
inline Float64 operator-(double x1, Float64 x2) { return x1-x2.dbl; }
inline Float64 operator*(double x1, Float64 x2) { return x1*x2.dbl; }
inline Float64 operator/(double x1, Float64 x2) { return x1/x2.dbl; }
inline Float64& operator+=(Float64& x1, double x2) { x1.dbl+=x2; return x1; }
inline Float64& operator-=(Float64& x1, double x2) { x1.dbl-=x2; return x1; }
inline Float64& operator*=(Float64& x1, double x2) { x1.dbl*=x2; return x1; }
inline Float64& operator/=(Float64& x1, double x2) { x1.dbl/=x2; return x1; }

// Deprecated mixed comparisons
inline Bool operator==(Float64 x1, double x2) { return x1.dbl==x2; }
inline Bool operator!=(Float64 x1, double x2) { return x1.dbl!=x2; }
inline Bool operator<=(Float64 x1, double x2) { return x1.dbl<=x2; }
inline Bool operator>=(Float64 x1, double x2) { return x1.dbl>=x2; }
inline Bool operator< (Float64 x1, double x2) { return x1.dbl< x2; }
inline Bool operator> (Float64 x1, double x2) { return x1.dbl> x2; }
inline Bool operator==(double x1, Float64 x2) { return x1==x2.dbl; }
inline Bool operator!=(double x1, Float64 x2) { return x1!=x2.dbl; }
inline Bool operator<=(double x1, Float64 x2) { return x1<=x2.dbl; }
inline Bool operator>=(double x1, Float64 x2) { return x1>=x2.dbl; }
inline Bool operator< (double x1, Float64 x2) { return x1< x2.dbl; }
inline Bool operator> (double x1, Float64 x2) { return x1> x2.dbl; }

// Support for creating floating-point objects from builtin, concrete and generic numbers
ApproximateFloat create_float(Number<Approximate>);
LowerFloat create_float(Number<Lower>);
UpperFloat create_float(Number<Upper>);
ValidatedFloat create_float(Number<Validated>);
ValidatedFloat create_float(Number<Effective>);
ValidatedFloat create_float(Number<Exact>);
ValidatedFloat create_float(Real);
ValidatedFloat create_float(Rational);
ExactFloat create_float(Integer);
template<class N, EnableIf<IsIntegral<N>> =dummy> ExactFloat create_float(N n);
template<class F, EnableIf<IsFloatingPoint<F>> =dummy> ApproximateFloat create_float(F f);

template<class X> struct IsGenericNumber : And<IsNumber<X>,Not<IsFloat<X>>> { };

// Mixed Float64 - Generic operations
template<class F, class X, EnableIf<IsFloat<F>> =dummy, EnableIf<IsGenericNumber<X>> =dummy>
    decltype(declval<F>()+create_float(declval<X>())) operator+(F const& f, X const& x) { return f+create_float(x); }
template<class F, class X, EnableIf<IsFloat<F>> =dummy, EnableIf<IsGenericNumber<X>> =dummy>
    decltype(declval<F>()-create_float(declval<X>())) operator-(F const& f, X const& x) { return f-create_float(x); }
template<class F, class X, EnableIf<IsFloat<F>> =dummy, EnableIf<IsGenericNumber<X>> =dummy>
    decltype(declval<F>()*create_float(declval<X>())) operator*(F const& f, X const& x) { return f*create_float(x); }
template<class F, class X, EnableIf<IsFloat<F>> =dummy, EnableIf<IsGenericNumber<X>> =dummy>
    decltype(declval<F>()/create_float(declval<X>())) operator/(F const& f, X const& x) { return f/create_float(x); }
template<class F, class X, EnableIf<IsFloat<F>> =dummy, EnableIf<IsGenericNumber<X>> =dummy>
    decltype(create_float(declval<X>())+declval<F>()) operator+(X const& x, F const& f) { return create_float(x)+f; }
template<class F, class X, EnableIf<IsFloat<F>> =dummy, EnableIf<IsGenericNumber<X>> =dummy>
    decltype(create_float(declval<X>())-declval<F>()) operator-(X const& x, F const& f) { return create_float(x)-f; }
template<class F, class X, EnableIf<IsFloat<F>> =dummy, EnableIf<IsGenericNumber<X>> =dummy>
    decltype(create_float(declval<X>())*declval<F>()) operator*(X const& x, F const& f) { return create_float(x)*f; }
template<class F, class X, EnableIf<IsFloat<F>> =dummy, EnableIf<IsGenericNumber<X>> =dummy>
    decltype(create_float(declval<X>())/declval<F>()) operator/(X const& x, F const& f) { return create_float(x)/f; }

template<class F, class X, EnableIf<IsFloat<F>> =dummy, EnableIf<IsGenericNumber<X>> =dummy>
    decltype(declval<F>()==create_float(declval<X>())) operator==(F const& f, X const& x) { return f==create_float(x); }
template<class F, class X, EnableIf<IsFloat<F>> =dummy, EnableIf<IsGenericNumber<X>> =dummy>
    decltype(declval<F>()!=create_float(declval<X>())) operator!=(F const& f, X const& x) { return f!=create_float(x); }
template<class F, class X, EnableIf<IsFloat<F>> =dummy, EnableIf<IsGenericNumber<X>> =dummy>
    decltype(declval<F>()>=create_float(declval<X>())) operator>=(F const& f, X const& x) { return f>=create_float(x); }
template<class F, class X, EnableIf<IsFloat<F>> =dummy, EnableIf<IsGenericNumber<X>> =dummy>
    decltype(declval<F>()<=create_float(declval<X>())) operator<=(F const& f, X const& x) { return f<=create_float(x); }
template<class F, class X, EnableIf<IsFloat<F>> =dummy, EnableIf<IsGenericNumber<X>> =dummy>
    decltype(declval<F>()> create_float(declval<X>())) operator> (F const& f, X const& x) { return f> create_float(x); }
template<class F, class X, EnableIf<IsFloat<F>> =dummy, EnableIf<IsGenericNumber<X>> =dummy>
    decltype(declval<F>()> create_float(declval<X>())) operator< (F const& f, X const& x) { return f< create_float(x); }

} // namespace Ariadne

#endif
