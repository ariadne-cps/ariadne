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

#include <iosfwd> // For std::floor std::ceil etc
#include <cmath> // For std::floor std::ceil etc
#include <algorithm> // For std::max, std::min
#include <limits> // For std::numeric_limits<double>

#include "utility/declarations.h"
#include "expression/operators.h"
#include "numeric/rounding.h"
#include "numeric/number.decl.h"
#include "numeric/float.decl.h"

namespace Ariadne {

class Float64;
typedef Float64 RawFloat64;

class Rational;


class Precision64 {
};

using RoundingMode64 = RoundingModeType;

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
    typedef Precision64 PrecisionType;
    typedef RoundingMode64 RoundingModeType;
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
    static Precision64 get_default_precision();
    Precision64 precision() const;
    Void set_precision(Precision64);
  public:
    static Float64 inf();
    static Float64 max();
    static Float64 eps();
    static Float64 min();
  public:
    //! \brief Default constructor creates an uninitialised number.
    Float64() : dbl() { }
    //! \brief Convert from a built-in double-precision floating-point number.
    Float64(double x) : dbl(x) { }
    Float64(double x, Precision64) : dbl(x) { }
    //! \brief Copy constructor.
    Float64(const Float64& x) : dbl(x.dbl) { }

    //! \brief Construct from a rational number with given rounding
    explicit Float64(const Rational& q, RoundingModeType rnd);
    //! \brief Convert to a rational number.
    explicit operator Rational () const;
  public:
    Float64 const& raw() const { return *this; }
    //! \brief An approximation by a built-in double-precision floating-point number.
    double get_d() const { return this->dbl; }

  public:
    friend Float64 next_up(Float64 x);
    friend Float64 next_down(Float64 x);

    friend Float64 floor(Float64 x);
    friend Float64 ceil(Float64 x);

    friend Float64 nul(Float64 x);
    friend Float64 half(Float64 x);
    friend Float64 pos(Float64 x);
    friend Float64 neg(Float64 x);
    friend Float64 sqr(Float64 x);
    friend Float64 rec(Float64 x);
    friend Float64 add(Float64 x1, Float64 x2);
    friend Float64 sub(Float64 x1, Float64 x2);
    friend Float64 mul(Float64 x1, Float64 x2);
    friend Float64 div(Float64 x1, Float64 x2);
    friend Float64 fma(Float64 x1, Float64 x2, Float64 x3); // x1*x2+x3
    friend Float64 pow(Float64 x, Int n);
    friend Float64 sqrt(Float64 x);
    friend Float64 exp(Float64 x);
    friend Float64 log(Float64 x);
    friend Float64 sin(Float64 x);
    friend Float64 cos(Float64 x);
    friend Float64 tan(Float64 x);
    friend Float64 asin(Float64 x);
    friend Float64 acos(Float64 x);
    friend Float64 atan(Float64 x);

    static Float64 pi(PrecisionType pr=get_default_precision(), RoundingModeType rnd=get_rounding_mode());

    friend Float64 max(Float64 x1, Float64 x2);
    friend Float64 min(Float64 x1, Float64 x2);
    friend Float64 abs(Float64 x);
    friend Float64 mag(Float64 x);

    friend Bool is_nan(Float64 x);

    // Operators
    friend Float64 operator+(Float64 x);
    friend Float64 operator-(Float64 x);
    friend Float64 operator+(Float64 x1, Float64 x2);
    friend Float64 operator-(Float64 x1, Float64 x2);
    friend Float64 operator*(Float64 x1, Float64 x2);
    friend Float64 operator/(Float64 x1, Float64 x2);
    friend Float64& operator+=(Float64& x1, Float64 x2);
    friend Float64& operator-=(Float64& x1, Float64 x2);
    friend Float64& operator*=(Float64& x1, Float64 x2);
    friend Float64& operator/=(Float64& x1, Float64 x2);

    friend Bool operator==(Float64 x1, Float64 x2);
    friend Bool operator!=(Float64 x1, Float64 x2);
    friend Bool operator<=(Float64 x1, Float64 x2);
    friend Bool operator>=(Float64 x1, Float64 x2);
    friend Bool operator< (Float64 x1, Float64 x2);
    friend Bool operator> (Float64 x1, Float64 x2);

    friend OutputStream& operator<<(OutputStream& os, Float64 const&);
    friend InputStream& operator>>(InputStream& is, Float64&);
};

// Correctly rounded algebraic and transcendental functions
double pow_rnd(double x, int n);
double sqrt_rnd(double x);
double exp_rnd(double x);
double log_rnd(double x);
double sin_rnd(double x);
double cos_rnd(double x);
double tan_rnd(double x);
double atan_rnd(double x);


inline Float max(Float x1, Float x2) { return std::max(x1.dbl,x2.dbl); }
inline Float min(Float x1, Float x2) { return std::min(x1.dbl,x2.dbl); }
inline Float abs(Float x) { return std::fabs(x.dbl); }

// Correctly rounded arithmetic
inline Float nul(Float x) { return +0.0; }
inline Float pos(Float x) { volatile double xv=x.dbl; return +xv; }
inline Float neg(Float x) { volatile double xv=x.dbl; return -xv; }
inline Float half(Float x) { volatile double xv=x.dbl; return xv/2; }
inline Float sqr(Float x) { volatile double xv=x.dbl; return xv*xv; }
inline Float rec(Float x) { volatile double xv=x.dbl; return 1.0/xv; }
inline Float add(Float x1, Float x2) { volatile double xv = x1.dbl; volatile double yv=x2.dbl; volatile double r=xv+yv; return r; }
inline Float sub(Float x1, Float x2) { volatile double xv = x1.dbl; volatile double yv=x2.dbl; volatile double r=xv-yv; return r; }
inline Float mul(Float x1, Float x2) { volatile double xv = x1.dbl; volatile double yv=x2.dbl; volatile double r=xv*yv; return r; }
inline Float div(Float x1, Float x2) { volatile double xv = x1.dbl; volatile double yv=x2.dbl; volatile double r=xv/yv; return r; }
inline Float pow(Float x, Int n) { return pow_rnd(x.dbl,n); }
inline Float sqrt(Float x) { return sqrt_rnd(x.dbl); }
inline Float exp(Float x) { return exp_rnd(x.dbl); }
inline Float log(Float x) { return log_rnd(x.dbl); }
inline Float sin(Float x) { return sin_rnd(x.dbl); }
inline Float cos(Float x) { return cos_rnd(x.dbl); }
inline Float tan(Float x) { return tan_rnd(x.dbl); }
inline Float asin(Float x) { return std::asin(x.dbl); }
inline Float acos(Float x) { return std::acos(x.dbl); }
inline Float atan(Float x) { return atan_rnd(x.dbl); }

// Opposite rounded arithmetic
inline Float pos_opp(Float x) { volatile double t=-x.dbl; return -t; }
inline Float neg_opp(Float x) { volatile double t=x.dbl; return -t; }
inline Float sqr_opp(Float x) { volatile double t=-x.dbl; t=t*x.dbl; return -t; }
inline Float rec_opp(Float x) { volatile double t=-1.0/(volatile double&)x.dbl; return -t; }
inline Float add_opp(Float x, Float y) { volatile double t=-x.dbl; t=t-y.dbl; return -t; }
inline Float sub_opp(Float x, Float y) { volatile double t=-x.dbl; t=t+y.dbl; return -t; }
inline Float mul_opp(Float x, Float y) { volatile double t=-x.dbl; t=t*y.dbl; return -t; }
inline Float div_opp(Float x, Float y) { volatile double t=x.dbl; t=t/y.dbl; return -t; }
Float pow_opp(Float x, int n);

template<class OP> Float64 apply(OP op, Float64 const& x1, Float64 const& x2, RoundingMode64 rnd) {
    auto old_rnd=Float64::get_rounding_mode(); Float64::set_rounding_mode(rnd);
    Float64 r=op(x1,x2); Float64::set_rounding_mode(old_rnd); return std::move(r);
}

template<class OP> Float64 apply(OP op, Float64 const& x, RoundingMode64 rnd) {
    auto old_rnd=Float64::get_rounding_mode(); Float64::set_rounding_mode(rnd);
    Float64 r=op(x); Float64::set_rounding_mode(old_rnd); return std::move(r);
}

template<class OP> Float64 apply(OP op, Float64 const& x, Int n, RoundingMode64 rnd) {
    auto old_rnd=Float64::get_rounding_mode(); Float64::set_rounding_mode(rnd);
    Float64 r=op(x,n); Float64::set_rounding_mode(old_rnd); return std::move(r);
}

inline Float64 add(Float64 const& x1, Float64 const& x2, Float64::RoundingModeType rnd) { return apply(Add(),x1,x2,rnd); }
inline Float64 sub(Float64 const& x1, Float64 const& x2, Float64::RoundingModeType rnd) { return apply(Sub(),x1,x2,rnd); }
inline Float64 mul(Float64 const& x1, Float64 const& x2, Float64::RoundingModeType rnd) { return apply(Mul(),x1,x2,rnd); }
inline Float64 div(Float64 const& x1, Float64 const& x2, Float64::RoundingModeType rnd) { return apply(Div(),x1,x2,rnd); }
inline Float64 pow(Float64 const& x, Int n, Float64::RoundingModeType rnd) { return apply(Pow(),x,n,rnd); }
inline Float64 sqr(Float64 const& x, Float64::RoundingModeType rnd) { return apply(Sqr(),x,rnd); }
inline Float64 rec(Float64 const& x, Float64::RoundingModeType rnd) { return apply(Rec(),x,rnd); }
inline Float64 sqrt(Float64 const& x, Float64::RoundingModeType rnd) { return apply(Sqrt(),x,rnd); }
inline Float64 exp(Float64 const& x, Float64::RoundingModeType rnd) { return apply(Exp(),x,rnd); }
inline Float64 log(Float64 const& x, Float64::RoundingModeType rnd) { return apply(Log(),x,rnd); }
inline Float64 sin(Float64 const& x, Float64::RoundingModeType rnd) { return apply(Sin(),x,rnd); }
inline Float64 cos(Float64 const& x, Float64::RoundingModeType rnd) { return apply(Cos(),x,rnd); }
inline Float64 tan(Float64 const& x, Float64::RoundingModeType rnd) { return apply(Tan(),x,rnd); }
inline Float64 atan(Float64 const& x, Float64::RoundingModeType rnd) { return apply(Atan(),x,rnd); }


//! \related Float \brief The nearest floating-point approximation to the constant \a pi.
inline const Float pi_approx() { return Float(3.1415926535897931); }
inline const Float pi_near()   { return Float(3.1415926535897931); }
inline const Float pi_down()   { return Float(3.1415926535897931); }
inline const Float pi_up()     { return Float(3.1415926535897936); }

inline Float fma_approx(Float x, Float y, Float z) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(to_nearest);
    Float t=mul(x,y); Float r=add(t,z); set_rounding_mode(rounding_mode); return r; }

inline Float fma_up(Float x, Float y, Float z) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(upward);
    Float t=mul(x,y); Float r=add(t,z); set_rounding_mode(rounding_mode); return r; }

inline Float fma_down(Float x, Float y, Float z) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(downward);
    Float t=mul(x,y); Float r=add(t,z); set_rounding_mode(rounding_mode); return r; }

//! \related Float \brief The average of two values, computed with nearest rounding. Also available with \c _ivl suffix.
inline Float med_approx(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(to_nearest);
    Float r=half(add(x,y)); set_rounding_mode(rounding_mode); return r; }
inline Float med_near(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(to_nearest);
    Float r=half(add(x,y)); set_rounding_mode(rounding_mode); return r; }
//! \related Float \brief Half of the difference of two values, computed with upward rounding. Also available with \c _ivl suffix.
inline Float rad_up(Float x, Float y) {
    rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(upward);
    Float r=half(sub(y,x)); set_rounding_mode(rounding_mode); return r; }

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

// Checking whether a Float64 is not-a-number
inline Bool isnan(Float64 x) { return std::isnan(x.dbl); }

// Discontinuous integer-valued functions
inline Float64 floor(Float64 x) { return std::floor(x.dbl); }
inline Float64 ceil(Float64 x) { return std::ceil(x.dbl); }

} // namespace Ariadne

#endif
