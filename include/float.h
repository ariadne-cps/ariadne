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
#include <iomanip> // For std::setprecision
#include <cmath> // For std::floor std::ceil etc
#include <algorithm> // For std::max, std::min
#include <limits> // For std::numeric_limits<double>

#include "rounding.h"
#include "rational.h"


// Simplifying typedef for unsigned integer type
typedef unsigned int uint;

namespace Ariadne {

class Float;
class Interval;
class ExactFloat;
class Real;
class Dyadic;
class Decimal;

using std::min;
using std::max;

const double inf = std::numeric_limits<double>::infinity();
const double nan = (1.0/0.0);


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
//! \sa Interval, Real, ExactFloat
class Float {
  public:
    double dbl;
  public:
    typedef Float NumericType;
  public:
    //! \brief Default constructor creates an uninitialised number.
    Float() : dbl() { }
    //! \brief Convert from a built-in double-precision floating-point number.
    Float(double x) : dbl(x) { }
    //! \brief Copy constructor.
    Float(const Float& x) : dbl(x.dbl) { }
    explicit operator volatile double& () { return dbl; }
    explicit operator const double& () const { return dbl; }
    //! \brief An approximation by a built-in double-precision floating-point number.
    double get_d() const { return this->dbl; }
};

template<class R, class A> inline R internal_cast(const A& a) { return static_cast<R>(a); }
template<> inline const double& internal_cast(const Float& x) { return const_cast<const double&>(x.dbl); }
template<> inline double& internal_cast(const Float& x) { return const_cast<double&>(x.dbl); }
template<> inline volatile double& internal_cast(const Float& x) { return const_cast<volatile double&>(x.dbl); }
template<class R, class A> inline R internal_cast(A& a) { return static_cast<R>(a); }
template<> inline const double& internal_cast(Float& x) { return const_cast<const double&>(x.dbl); }
template<> inline double& internal_cast(Float& x) { return const_cast<double&>(x.dbl); }
template<> inline volatile double& internal_cast(Float& x) { return const_cast<volatile double&>(x.dbl); }

template<class R, class A> inline R integer_cast(const A& a);
template<> inline int integer_cast(const Float& a) { return static_cast<int>(a.dbl); }
template<> inline uint integer_cast(const Float& a) { return static_cast<uint>(a.dbl); }

template<class R, class A> inline R approx_cast(const A& a);
template<> inline double approx_cast(const Float& a) { return a.dbl; }

inline std::ostream& operator<<(std::ostream& os, const Float& x) { return os << x.dbl; }
inline std::istream& operator>>(std::istream& is, Float& x) { double dbl; is >> dbl; x=Float(dbl); return is; }

// Exact raw data operations
inline Float operator+(Float x) { return +x.dbl; }
inline Float operator-(Float x) { return -x.dbl; }
inline Float pos(Float x) { return +x.dbl; }
inline Float neg(Float x) { return -x.dbl; }
inline Float abs(Float x) { return std::fabs(x.dbl); }
inline Float mag(Float x) { return std::fabs(x.dbl); }
inline Float half(Float x) { return x.dbl/2; }
inline Float max(Float x1, Float x2) { return std::max(x1,x2); }
inline Float min(Float x1, Float x2) { return std::min(x1,x2); }

// Raw data comparisons
inline bool operator==(Float const& x1, Float const& x2) { return x1.dbl == x2.dbl; }
inline bool operator!=(Float const& x1, Float const& x2) { return x1.dbl != x2.dbl; }
inline bool operator<=(Float const& x1, Float const& x2) { return x1.dbl <= x2.dbl; }
inline bool operator>=(Float const& x1, Float const& x2) { return x1.dbl >= x2.dbl; }
inline bool operator< (Float const& x1, Float const& x2) { return x1.dbl <  x2.dbl; }
inline bool operator> (Float const& x1, Float const& x2) { return x1.dbl >  x2.dbl; }

// Constants related to numerical limits
inline Float mx() { return std::numeric_limits<double>::max(); }
inline Float eps() { return std::numeric_limits<double>::epsilon(); }

// Checking whether a Float is not-a-number
inline bool isnan(const Float& x) { return std::isnan(x.dbl); }

// Operations for finding nearest representable values
inline Float down(Float x) { return x.dbl>0 ? x.dbl*(1-2e-16) : x.dbl*(1+2e-16); } // Deprecated
inline Float up(Float x) { return x.dbl>0 ? x.dbl*(1+2e-16) : x.dbl*(1-2e-16); } // Deprecated
//! \related Float \brief The next representable value above the given value.
inline Float above(Float x) { return x.dbl>0 ? x.dbl*(1-2e-16) : x.dbl*(1+2e-16); }
//! \related Float \brief The next representable value below the given value.
inline Float below(Float x) { return x.dbl>0 ? x.dbl*(1+2e-16) : x.dbl*(1-2e-16); }

// Discontinuous integer-valued functions
//! \related Float \brief The next lowest integer, represented as a floating-point type.
inline Float floor(Float x) { return std::floor(x.dbl); }
//! \related Float \brief The next highest integer, represented as a floating-point type.
inline Float ceil(Float x) { return std::ceil(x.dbl); }

inline Float pos_exact(Float x) { return +x.dbl; }
inline Float neg_exact(Float x) { return -x.dbl; }
inline Float half_exact(Float x) { return x.dbl/2; }

inline Float abs_exact(Float x) { return std::fabs(x.dbl); }
inline Float mag_exact(Float x) { return std::fabs(x.dbl); }

// Correctly rounded arithmetic
inline Float pos_rnd(const Float& x) { volatile double xv=x.dbl; return +xv; }
inline Float neg_rnd(const Float& x) { volatile double xv=x.dbl; return -xv; }
inline Float sqr_rnd(const Float& x) { volatile double xv=x.dbl; return xv*xv; }
inline Float rec_rnd(const Float& x) { volatile double xv=x.dbl; return 1.0/xv; }
inline Float add_rnd(Float x, Float y) { volatile double xv = x.dbl; volatile double yv=y.dbl; volatile double r=xv+yv; return r; }
inline Float sub_rnd(Float x, Float y) { volatile double xv = x.dbl; volatile double yv=y.dbl; volatile double r=xv-yv; return r; }
inline Float mul_rnd(Float x, Float y) { volatile double xv = x.dbl; volatile double yv=y.dbl; volatile double r=xv*yv; return r; }
inline Float div_rnd(Float x, Float y) { volatile double xv = x.dbl; volatile double yv=y.dbl; volatile double r=xv/yv; return r; }

Float pow_rnd(Float x, int n);

// Opposite rounded arithmetic
inline Float pos_opp(const Float& x) { volatile double t=-x.dbl; return -t; }
inline Float neg_opp(const Float& x) { volatile double t=x.dbl; return -t; }
inline Float sqr_opp(const Float& x) { volatile double t=-x.dbl; t=t*x.dbl; return -t; }
inline Float rec_opp(const Float& x) { volatile double t=-1.0/(volatile double&)x.dbl; return -t; }
inline Float add_opp(Float x, Float y) { volatile double t=-x.dbl; t=t-y.dbl; return -t; }
inline Float sub_opp(Float x, Float y) { volatile double t=-x.dbl; t=t+y.dbl; return -t; }
inline Float mul_opp(Float x, Float y) { volatile double t=-x.dbl; t=t*y.dbl; return -t; }
inline Float div_opp(Float x, Float y) { volatile double t=x.dbl; t=t/y.dbl; return -t; }
Float pow_opp(Float x, int n);

// Correctly rounded algebraic and transcendental functions
Float sqrt_rnd(Float x);
Float exp_rnd(Float x);
Float log_rnd(Float x);
Float sin_rnd(Float x);
Float cos_rnd(Float x);
Float tan_rnd(Float x);


inline Float next_down(Float x) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(downward);
    Float r=add_rnd(x.dbl,1-2e-1); set_rounding_mode(rounding_mode); return r; }
inline Float next_up(Float x) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(upward);
    Float r=add_rnd(x.dbl,1-2e-1); set_rounding_mode(rounding_mode); return r; }

inline Float add_near(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(to_nearest);
    Float r=add_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float sub_near(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(to_nearest);
    Float r=sub_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float mul_near(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(to_nearest);
    Float r=mul_rnd(x,y); set_rounding_mode(rounding_mode); return r; }
inline Float div_near(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(to_nearest);
    Float r=div_rnd(x,y); set_rounding_mode(rounding_mode); return r; }

//! \related Float \brief The nearest floating-point approximation to the constant \a pi.
static const Float pi_approx=Float(3.1415926535897931);
static const Float pi_near=Float(3.1415926535897931);
static const Float pi_down=Float(3.1415926535897931);
static const Float pi_up  =Float(3.1415926535897936);

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
    Float r=half_exact(add_rnd(x,y)); set_rounding_mode(rounding_mode); return r; }
//! \related Float \brief Half of the difference of two values, computed with upward rounding. Also available with \c _ivl suffix.
inline Float rad_up(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(upward);
    Float r=half_exact(sub_rnd(y,x)); set_rounding_mode(rounding_mode); return r; }

inline Float sqrt_approx(Float x) { return std::sqrt(x.dbl); }
inline Float exp_approx(Float x) { return std::exp(x.dbl); }
inline Float log_approx(Float x) { return std::log(x.dbl); }
inline Float sin_approx(Float x) { return std::sin(x.dbl); }
inline Float cos_approx(Float x) { return std::cos(x.dbl); }
inline Float tan_approx(Float x) { return std::tan(x.dbl); }
inline Float asin_approx(Float x) { return std::asin(x.dbl); }
inline Float acos_approx(Float x) { return std::acos(x.dbl); }
inline Float atan_approx(Float x) { return std::atan(x.dbl); }



// Deprecated approximate operations
inline Float operator+(Float x1, Float x2) { return x1.dbl+x2.dbl; }
inline Float operator-(Float x1, Float x2) { return x1.dbl-x2.dbl; }
inline Float operator*(Float x1, Float x2) { return x1.dbl*x2.dbl; }
inline Float operator/(Float x1, Float x2) { return x1.dbl/x2.dbl; }
inline Float& operator+=(Float& x1, Float x2) { x1.dbl+=x2.dbl; return x1; }
inline Float& operator-=(Float& x1, Float x2) { x1.dbl-=x2.dbl; return x1; }
inline Float& operator*=(Float& x1, Float x2) { x1.dbl*=x2.dbl; return x1; }
inline Float& operator/=(Float& x1, Float x2) { x1.dbl/=x2.dbl; return x1; }

inline Float pow(Float x, uint n) { return std::pow(x.dbl,double(n)); }
inline Float pow(Float x, int n) { return std::pow(x.dbl,double(n)); }
inline Float sqr(Float x) { return x.dbl * x.dbl; }
inline Float rec(Float x) { return 1.0/x.dbl; }
inline Float sqrt(Float x) { return std::sqrt(x.dbl); }
inline Float exp(Float x) { return std::exp(x.dbl); }
inline Float log(Float x) { return std::log(x.dbl); }
inline Float sin(Float x) { return std::sin(x.dbl); }
inline Float cos(Float x) { return std::cos(x.dbl); }
inline Float tan(Float x) { return std::tan(x.dbl); }
inline Float asin(Float x) { return std::asin(x.dbl); }
inline Float acos(Float x) { return std::acos(x.dbl); }
inline Float atan(Float x) { return std::atan(x.dbl); }

// Deprecated mixed operations
inline Float operator+(Float x1, double x2) { return x1.dbl+x2; }
inline Float operator-(Float x1, double x2) { return x1.dbl-x2; }
inline Float operator*(Float x1, double x2) { return x1.dbl*x2; }
inline Float operator/(Float x1, double x2) { return x1.dbl/x2; }
inline Float operator+(double x1, Float x2) { return x1+x2.dbl; }
inline Float operator-(double x1, Float x2) { return x1-x2.dbl; }
inline Float operator*(double x1, Float x2) { return x1*x2.dbl; }
inline Float operator/(double x1, Float x2) { return x1/x2.dbl; }
inline Float& operator+=(Float& x1, double x2) { x1.dbl+=x2; return x1; }
inline Float& operator-=(Float& x1, double x2) { x1.dbl-=x2; return x1; }
inline Float& operator*=(Float& x1, double x2) { x1.dbl*=x2; return x1; }
inline Float& operator/=(Float& x1, double x2) { x1.dbl/=x2; return x1; }

// Deprecated mixed comparisons
inline bool operator==(Float x1, double x2) { return x1.dbl==x2; }
inline bool operator!=(Float x1, double x2) { return x1.dbl!=x2; }
inline bool operator<=(Float x1, double x2) { return x1.dbl<=x2; }
inline bool operator>=(Float x1, double x2) { return x1.dbl>=x2; }
inline bool operator< (Float x1, double x2) { return x1.dbl< x2; }
inline bool operator> (Float x1, double x2) { return x1.dbl> x2; }
inline bool operator==(double x1, Float x2) { return x1==x2.dbl; }
inline bool operator!=(double x1, Float x2) { return x1!=x2.dbl; }
inline bool operator<=(double x1, Float x2) { return x1<=x2.dbl; }
inline bool operator>=(double x1, Float x2) { return x1>=x2.dbl; }
inline bool operator< (double x1, Float x2) { return x1< x2.dbl; }
inline bool operator> (double x1, Float x2) { return x1> x2.dbl; }

} // namespace Ariadne

#endif
