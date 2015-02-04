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
#include "numeric/rounding.h"
#include "numeric/number.decl.h"
#include "numeric/float.decl.h"

namespace Ariadne {

class Float64;
typedef Float64 RawFloat64;

class Rational;


class Precision64 {
};

using RoundingMode64 = rounding_mode_t;

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
//    explicit operator volatile double& () { return const_cast<volatile double&>(dbl); }
//    explicit operator const double& () const { return const_cast<const double&>(dbl); }
    Float64 const& raw() const { return *this; }
    //! \brief An approximation by a built-in double-precision floating-point number.
    double get_d() const { return this->dbl; }

  public:
    friend Float64 next_up(Float64 const& x);
    friend Float64 next_down(Float64 const& x);

    friend Float64 floor(Float64 const& x);
    friend Float64 ceil(Float64 const& x);

    friend Float64 nul(Float64 const& x);
    friend Float64 half(Float64 const& x);
    friend Float64 pos(Float64 const& x);
    friend Float64 neg(Float64 const& x);
    friend Float64 sqr(Float64 const& x);
    friend Float64 rec(Float64 const& x);
    friend Float64 add(Float64 const& x1, Float64 const& x2);
    friend Float64 sub(Float64 const& x1, Float64 const& x2);
    friend Float64 mul(Float64 const& x1, Float64 const& x2);
    friend Float64 div(Float64 const& x1, Float64 const& x2);
    friend Float64 fma(Float64 const& x1, Float64 const& x2, Float64 const& x3); // x1*x2+x3
    friend Float64 pow(Float64 const& x, Int n);
    friend Float64 sqrt(Float64 const& x);
    friend Float64 exp(Float64 const& x);
    friend Float64 log(Float64 const& x);
    friend Float64 sin(Float64 const& x);
    friend Float64 cos(Float64 const& x);
    friend Float64 tan(Float64 const& x);
    friend Float64 asin(Float64 const& x);
    friend Float64 acos(Float64 const& x);
    friend Float64 atan(Float64 const& x);

    static Float64 pi(PrecisionType pr=get_default_precision(), RoundingModeType rnd=get_rounding_mode());

    friend Float64 max(Float64 const& x1, Float64 const& x2);
    friend Float64 min(Float64 const& x1, Float64 const& x2);
    friend Float64 abs(Float64 const& x);
    friend Float64 mag(Float64 const& x);

    friend Bool is_nan(Float64 const& x);

    // Operators
    friend Float64 operator+(Float64 const& x);
    friend Float64 operator-(Float64 const& x);
    friend Float64 operator+(Float64 const& x1, Float64 const& x2);
    friend Float64 operator-(Float64 const& x1, Float64 const& x2);
    friend Float64 operator*(Float64 const& x1, Float64 const& x2);
    friend Float64 operator/(Float64 const& x1, Float64 const& x2);
    friend Float64& operator+=(Float64& x1, Float64 const& x2);
    friend Float64& operator-=(Float64& x1, Float64 const& x2);
    friend Float64& operator*=(Float64& x1, Float64 const& x2);
    friend Float64& operator/=(Float64& x1, Float64 const& x2);

    friend Bool operator==(Float64 const& x1, Float64 const& x2);
    friend Bool operator!=(Float64 const& x1, Float64 const& x2);
    friend Bool operator<=(Float64 const& x1, Float64 const& x2);
    friend Bool operator>=(Float64 const& x1, Float64 const& x2);
    friend Bool operator< (Float64 const& x1, Float64 const& x2);
    friend Bool operator> (Float64 const& x1, Float64 const& x2);

    friend OutputStream& operator<<(OutputStream& os, Float64 const&);
    friend InputStream& operator>>(InputStream& is, Float64&);
};

template<class R, class A> R integer_cast(const A& a);
template<> inline Int integer_cast(const Float64& x) { return static_cast<Int>(x.get_d()); }
template<> inline Nat integer_cast(const Float64& x) { return static_cast<Nat>(x.get_d()); }

Float64 sqr_rnd(Float64 x);
Float64 add_rnd(Float64 x1, Float64 x2);
Float64 sub_rnd(Float64 x1, Float64 x2);
Float64 mul_rnd(Float64 x1, Float64 x2);
Float64 div_rnd(Float64 x1, Float64 x2);
Float64 pow_rnd(Float64 x, Int n);
Float64 sqrt_rnd(Float64 x);
Float64 exp_rnd(Float64 x);
Float64 log_rnd(Float64 x);
Float64 sin_rnd(Float64 x);
Float64 cos_rnd(Float64 x);
Float64 tan_rnd(Float64 x);
Float64 atan_rnd(Float64 x);



inline OutputStream& operator<<(OutputStream& os, const Float64& x) { return os << x.dbl; }
inline InputStream& operator>>(InputStream& is, Float64& x) { double dbl; is >> dbl; x=Float64(dbl); return is; }

inline Float64::RoundingModeType Float64::get_rounding_mode() { return Ariadne::get_rounding_mode(); }
inline Void Float64::set_rounding_mode(RoundingModeType rnd) { Ariadne::set_rounding_mode(rnd); }
inline Void Float64::set_rounding_upward() { Ariadne::set_rounding_upward(); }
inline Void Float64::set_rounding_downward() { Ariadne::set_rounding_downward(); }
inline Void Float64::set_rounding_to_nearest() { Ariadne::set_rounding_to_nearest(); }
inline Void Float64::set_rounding_toward_zero() { Ariadne::set_rounding_toward_zero(); }

inline Float64::PrecisionType Float64::get_default_precision() { return Float64::PrecisionType(); }
inline Float64::PrecisionType Float64::precision() const { return Float64::PrecisionType(); }
inline Void Float64::set_precision(Float64::PrecisionType) { }

// Exact raw data operations
inline Float64 nul(Float64 const& x) { return 0.0; }
inline Float64 pos(Float64 const& x) { return +x.dbl; }
inline Float64 neg(Float64 const& x) { return -x.dbl; }
inline Float64 half(Float64 const& x) { return x.dbl/2.0; }

inline Float64 max(Float64 const& x1, Float64 const& x2) { return std::max(x1,x2); }
inline Float64 min(Float64 const& x1, Float64 const& x2) { return std::min(x1,x2); }
inline Float64 abs(Float64 const& x) { return std::fabs(x.dbl); }
inline Float64 mag(Float64 const& x) { return std::fabs(x.dbl); }

inline Float64 add(Float64 const& x1, Float64 const& x2) { return x1.dbl+x2.dbl; }
inline Float64 sub(Float64 const& x1, Float64 const& x2) { return x1.dbl-x2.dbl; }
inline Float64 mul(Float64 const& x1, Float64 const& x2) { return x1.dbl*x2.dbl; }
inline Float64 div(Float64 const& x1, Float64 const& x2) { return x1.dbl/x2.dbl; }
inline Float64 fma(Float64 const& x1, Float64 const& x2, Float64 const& x3) { return std::fma(x1.dbl,x2.dbl,x3.dbl); }
inline Float64 pow(Float64 const& x, Int n) { return std::pow(x.dbl,double(n)); }
inline Float64 sqr(Float64 const& x) { return x.dbl * x.dbl; }
inline Float64 rec(Float64 const& x) { return 1.0/x.dbl; }
inline Float64 sqrt(Float64 const& x) { return std::sqrt(x.dbl); }
inline Float64 exp(Float64 const& x) { return std::exp(x.dbl); }
inline Float64 log(Float64 const& x) { return std::log(x.dbl); }
inline Float64 sin(Float64 const& x) { return std::sin(x.dbl); }
inline Float64 cos(Float64 const& x) { return cos_rnd(x); }
inline Float64 tan(Float64 const& x) { return std::tan(x.dbl); }
inline Float64 asin(Float64 const& x) { return std::asin(x.dbl); }
inline Float64 acos(Float64 const& x) { return std::acos(x.dbl); }
inline Float64 atan(Float64 const& x) { return std::atan(x.dbl); }

// Constants related to numerical limits
inline Float64 Float64::inf() { return std::numeric_limits<double>::infinity(); }
inline Float64 Float64::max() { return std::numeric_limits<double>::max(); }
inline Float64 Float64::min() { return std::numeric_limits<double>::min(); }
inline Float64 Float64::eps() { return std::numeric_limits<double>::epsilon(); }

// Arithmetic operators
inline Float64 operator+(Float64 const& x) { return +x.dbl; }
inline Float64 operator-(Float64 const& x) { return -x.dbl; }
inline Float64 operator+(Float64 const& x1, Float64 const& x2) { return x1.dbl+x2.dbl; }
inline Float64 operator-(Float64 const& x1, Float64 const& x2) { return x1.dbl-x2.dbl; }
inline Float64 operator*(Float64 const& x1, Float64 const& x2) { return x1.dbl*x2.dbl; }
inline Float64 operator/(Float64 const& x1, Float64 const& x2) { return x1.dbl/x2.dbl; }
inline Float64& operator+=(Float64& x1, Float64 const& x2) { x1.dbl+=x2.dbl; return x1; }
inline Float64& operator-=(Float64& x1, Float64 const& x2) { x1.dbl-=x2.dbl; return x1; }
inline Float64& operator*=(Float64& x1, Float64 const& x2) { x1.dbl*=x2.dbl; return x1; }
inline Float64& operator/=(Float64& x1, Float64 const& x2) { x1.dbl/=x2.dbl; return x1; }

// Comparison operators
inline Bool operator==(Float64 const& x1, Float64 const& x2) { return x1.dbl == x2.dbl; }
inline Bool operator!=(Float64 const& x1, Float64 const& x2) { return x1.dbl != x2.dbl; }
inline Bool operator<=(Float64 const& x1, Float64 const& x2) { return x1.dbl <= x2.dbl; }
inline Bool operator>=(Float64 const& x1, Float64 const& x2) { return x1.dbl >= x2.dbl; }
inline Bool operator< (Float64 const& x1, Float64 const& x2) { return x1.dbl <  x2.dbl; }
inline Bool operator> (Float64 const& x1, Float64 const& x2) { return x1.dbl >  x2.dbl; }

// Checking whether a Float64 is not-a-number
inline Bool isnan(Float64 const& x) { return std::isnan(x.dbl); }

// Discontinuous integer-valued functions
inline Float64 floor(Float64 const& x) { return std::floor(x.dbl); }
inline Float64 ceil(Float64 const& x) { return std::ceil(x.dbl); }

} // namespace Ariadne

#endif
