/***************************************************************************
 *            numeric/flt64.h
 *
 *  Copyright 2013-14  Pieter Collins
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

/*! \file numeric/flt64.h
 *  \brief
 */

#ifndef ARIADNE_FLT64_H
#define ARIADNE_FLT64_H

#include <cmath>
#include <iostream>

#include "utility/metaprogramming.h"
#include "utility/typedefs.h"
#include "numeric/paradigm.h"

#include "numeric/sign.h"
#include "numeric/is_number.h"

#include "numeric/logical.decl.h"
#include "numeric/float.decl.h"
#include "numeric/number.h"
#include "numeric/real.h"
#include "numeric/logical.h"

namespace Ariadne {

/************ Flt64 ********************************************************/

Void set_default_rounding();

Void Float64::set_rounding_to_nearest();
Void Float64::set_rounding_downward();
Void Float64::set_rounding_upward();
Void Float64::set_rounding_toward_zero();

typedef unsigned short RoundingModeType;
RoundingModeType Float64::get_rounding_mode();
Void Float64::set_rounding_mode(RoundingModeType);

extern const RoundingModeType ROUND_NEAR;
extern const RoundingModeType ROUND_DOWN;
extern const RoundingModeType ROUND_UP;
extern const RoundingModeType ROUND_ZERO;

// ----------------- RawTag Floatt64 class ---------------------------------------------------------

template<class Z, EnableIf<IsIntegral<Z>> =dummy> Z integer_cast(Flt64);

// ----------------- RawTag Float64 class ---------------------------------------------------------

//! \ingroup Float64SubModule
//! \brief Wrapper class for double-precision floating-point numbers.
//! \details
//! The main purpose of this class is as a base implementation for the \em safe floating-point number types
//! ApproximateFloat64, LowerFloat64, UpperFloat64,BoundedFloat64, MetricFloat64, ErrorFloat64 and ExactFloat64.
//! These classes
//!
//! Default arithmetic operations are approximate, and comparisons are exact, so this class is \em unsafe.
//! %Rounded operations on \b Float64 classes are provided using the \c _down
//! and \c _up suffixes. %ApproximateTag arithmetic without control of the
//! rounding mode can be specified by the \c _approx rounding suffix, or \c _near if round-to-nearest is available.
//! %ExactTag arithmetic can be specified explicitly using the \c _exact rounding suffix.
//!
//! Since the raw \b Float64 classes represent unsafe values, they cannot be converted to other numbers.
//! The exception is conversion to \b ApproximateFloat64 classes, since these represent approximate values with no control of the error.
//! Since the raw \b Float64 classes are used as the implementation of the safe \c Float64 classes,
//! they are valid arguments to constructors, but these constructors are all declared \c explicit.
//! Comparison of raw \b %Float64 data and other numbers is performed as if the %Float64 object were an \b %ApproximateFloat64.
//!
//! When testing, it is often useful to perform comparisons with \c %Flt64 or \c double values.
//! Although care must be taken, since compiler rounding may change the truth of certain comparisons.
//! comparisons with \c double are performed as if the value were \c ExactFloat64.


class Flt64 {
  public:
    volatile double d;
  public:
    typedef RawTag Paradigm;
    typedef Flt64 NumericType;
    Flt64() : d(0.0) { }
    Flt64(double dbl) : d(dbl) { }
    explicit operator double() const { return d; }
    Flt64 raw() const { return *this; }
    Flt64 get_flt() const { return *this; }
    double get_d() const { return d; }

    static Flt64 eps();
    static Flt64 inf();

    static Void Float64::set_rounding_mode(RoundingModeType);
    static RoundingModeType Float64::get_rounding_mode();

    static Void Float64::set_rounding_to_nearest();
    static Void Float64::set_rounding_downward();
    static Void Float64::set_rounding_upward();
    static Void Float64::set_rounding_toward_zero();

    friend inline Flt64 operator+(Flt64 x) { return Flt64{+x.d}; }
    friend inline Flt64 operator-(Flt64 x) { return Flt64{-x.d}; }
    friend inline Flt64 operator+(Flt64 x1, Flt64 x2) { return Flt64{x1.d+x2.d}; }
    friend inline Flt64 operator-(Flt64 x1, Flt64 x2) { return Flt64{x1.d-x2.d}; }
    friend inline Flt64 operator*(Flt64 x1, Flt64 x2) { return Flt64{x1.d*x2.d}; }
    friend inline Flt64 operator/(Flt64 x1, Flt64 x2) { return Flt64{x1.d/x2.d}; }
    friend inline Flt64 operator/(Flt64 x1, Int n2) { return Flt64{x1.d/n2}; }
    friend inline Flt64& operator+=(Flt64& x1, Flt64 x2) { x1.d+=x2.d; return x1; }
    friend inline Flt64& operator-=(Flt64& x1, Flt64 x2) { x1.d-=x2.d; return x1; }
    friend inline Flt64& operator*=(Flt64& x1, Flt64 x2) { x1.d*=x2.d; return x1; }
    friend inline Flt64& operator/=(Flt64& x1, Flt64 x2) { x1.d/=x2.d; return x1; }

    friend inline Flt64 min(Flt64 x1, Flt64 x2) { return Flt64{x2.d<x1.d?x2.d:x1.d}; }
    friend inline Flt64 max(Flt64 x1, Flt64 x2) { return Flt64{x2.d>x1.d?x2.d:x1.d}; }
    friend inline Flt64 abs(Flt64 x) { return Flt64{std::fabs(x.d)}; }
    friend inline Flt64 mig(Flt64 x) { return abs(x); }
    friend inline Flt64 mag(Flt64 x) { return abs(x); }

    friend inline Flt64 add(Flt64 x1, Flt64 x2) { return Flt64{x1.d+x2.d}; }
    friend inline Flt64 sub(Flt64 x1, Flt64 x2) { return Flt64{x1.d-x2.d}; }
    friend inline Flt64 mul(Flt64 x1, Flt64 x2) { return Flt64{x1.d*x2.d}; }
    friend inline Flt64 div(Flt64 x1, Flt64 x2) { return Flt64{x1.d/x2.d}; }

    friend inline Flt64 pos(Flt64 x) { return Flt64{x.d}; }
    friend inline Flt64 half(Flt64 x) { return Flt64{x.d/2}; }
    friend inline Flt64 sqr(Flt64 x) { return Flt64{x.d*x.d}; }

    friend inline Flt64 neg(Flt64 x) { return Flt64{-x.d}; }
    friend inline Flt64 rec(Flt64 x) { return Flt64{1.0/x.d}; }

    friend inline Flt64 pow(Flt64 x, Nat m) { return std::pow(x.d,m); }
    friend inline Flt64 pow(Flt64 x, Int n) { return std::pow(x.d,n); }
    friend inline Flt64 sqrt(Flt64 x) { return std::sqrt(x.d); }
    friend inline Flt64 exp(Flt64 x) { return std::exp(x.d); }
    friend inline Flt64 log(Flt64 x) { return std::log(x.d); }
    friend inline Flt64 sin(Flt64 x) { return std::sin(x.d); }
    friend inline Flt64 cos(Flt64 x) { return std::cos(x.d); }
    friend inline Flt64 tan(Flt64 x) { return std::tan(x.d); }
    friend inline Flt64 asin(Flt64 x) { return std::asin(x.d); }
    friend inline Flt64 acos(Flt64 x) { return std::acos(x.d); }
    friend inline Flt64 atan(Flt64 x) { return std::atan(x.d); }

    friend inline Int32Type floor(Flt64 x) { return int32_t(std::floor(x.d)); }
    friend inline Int32Type ceil(Flt64 x) { return int32_t(std::ceil(x.d)); }

    friend Bool isnan(Flt64 x);

    friend inline OutputStream& operator<<(OutputStream& os, Flt64 x) { return os << x.get_d(); }
    friend inline InputStream& operator>>(InputStream& is, Flt64& x) { return is >> reinterpret_cast<double&>(x); }

    friend inline Comparison cmp(Flt64 x1, Flt64 x2) {
        return (x1.d==x2.d)?Comparison::EQUAL:(x1.d>x2.d)?Comparison::GREATER:Comparison::LESS; }
    friend inline Bool operator==(Flt64 x1, Flt64 x2) { return x1.d==x2.d; }
    friend inline Bool operator!=(Flt64 x1, Flt64 x2) { return !(x1==x2); }
    friend inline Bool operator<=(Flt64 x1, Flt64 x2) { return x1.d<=x2.d; }
    friend inline Bool operator>=(Flt64 x1, Flt64 x2) { return (x2<=x1); }
    friend inline Bool operator< (Flt64 x1, Flt64 x2) { return !(x2<=x1); }
    friend inline Bool operator> (Flt64 x1, Flt64 x2) { return !(x1<=x2); }

    friend inline Bool operator==(Flt64 x1, double x2) { return x1.d==x2; }
    friend inline Bool operator!=(Flt64 x1, double x2) { return !(x1==x2); }
    friend inline Bool operator<=(Flt64 x1, double x2) { return x1.d<=x2; }
    friend inline Bool operator>=(Flt64 x1, double x2) { return (x2<=x1); }
    friend inline Bool operator< (Flt64 x1, double x2) { return !(x2<=x1); }
    friend inline Bool operator> (Flt64 x1, double x2) { return !(x1<=x2); }

};

inline Flt64 Flt64::inf() { return Flt64{1.0/0.0}; }

/************ Rounded arithmetic of Flt64 ***************************************************/

Flt64 half_exact(Flt64);

Flt64 next_down(Flt64);
Flt64 next_up(Flt64);

Flt64 add_near(Flt64, Flt64);
Flt64 add_down(Flt64, Flt64);
Flt64 add_up(Flt64, Flt64);

Flt64 sub_near(Flt64, Flt64);
Flt64 sub_down(Flt64, Flt64);
Flt64 sub_up(Flt64, Flt64);

Flt64 mul_near(Flt64, Flt64);
Flt64 mul_down(Flt64, Flt64);
Flt64 mul_up(Flt64, Flt64);

Flt64 div_near(Flt64, Flt64);
Flt64 div_down(Flt64, Flt64);
Flt64 div_up(Flt64, Flt64);

Flt64 rad_up(Flt64, Flt64);
Flt64 med_near(Flt64, Flt64);

Flt64 pow_approx(Flt64, Int);
Flt64 pow_down(Flt64, Int);
Flt64 pow_up(Flt64, Int);

Flt64 sqrt_approx(Flt64);
Flt64 sqrt_down(Flt64);
Flt64 sqrt_up(Flt64);
Flt64 exp_approx(Flt64);
Flt64 exp_down(Flt64);
Flt64 exp_up(Flt64);
Flt64 log_approx(Flt64);
Flt64 log_down(Flt64);
Flt64 log_up(Flt64);
Flt64 sin_approx(Flt64);
Flt64 sin_down(Flt64);
Flt64 sin_up(Flt64);
Flt64 cos_approx(Flt64);
Flt64 cos_down(Flt64);
Flt64 cos_up(Flt64);
Flt64 tan_approx(Flt64);
Flt64 tan_down(Flt64);
Flt64 tan_up(Flt64);
Flt64 atan_approx(Flt64);
Flt64 atan_down(Flt64);
Flt64 atan_up(Flt64);

// Arithmetic respecting the rounding mode
Flt64 add_rnd(Flt64 x, Flt64 y);
Flt64 sub_rnd(Flt64 x, Flt64 y);
Flt64 mul_rnd(Flt64 x, Flt64 y);
Flt64 div_rnd(Flt64 x, Flt64 y);
Flt64 add_opp(Flt64 x, Flt64 y);
Flt64 sub_opp(Flt64 x, Flt64 y);
Flt64 mul_opp(Flt64 x, Flt64 y);
Flt64 div_opp(Flt64 x, Flt64 y);

} // namespace Ariadne

#endif
