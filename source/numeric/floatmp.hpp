/***************************************************************************
 *            numeric/floatmp.hpp
 *
 *  Copyright 2013-17  Pieter Collins
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
 *  You should have received _a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file numeric/floatmp.hpp
 *  \brief RawTag floating-point numbers based on MPFT floats.
 */



#ifndef ARIADNE_FLOATMP_HPP
#define ARIADNE_FLOATMP_HPP

#include "paradigm.hpp"
#include "number.hpp"
#include <mpfr.h>

namespace Ariadne {

/************ FloatMP ********************************************************/

struct NoInit { };
struct RawPtr { };

//enum RoundingModeMP { NEAREST=MPFR_RNDN, UPWARD=MPFR_RNDU, DOWNWARD=MPFR_RNDD };
typedef mpfr_rnd_t RoundingModeMP;

class PrecisionMP {
    mpfr_prec_t prec;
  public:
    explicit PrecisionMP(mpfr_prec_t pr) : prec(pr) { }
    mpfr_prec_t bits() const { return prec; }
    operator mpfr_prec_t () const { return prec; }
    friend PrecisionMP max(PrecisionMP mp1, PrecisionMP mp2) { return PrecisionMP(std::max(mp1.bits(),mp2.bits())); }
    friend PrecisionMP min(PrecisionMP mp1, PrecisionMP mp2) { return PrecisionMP(std::min(mp1.bits(),mp2.bits())); }
    friend bool operator==(PrecisionMP mp1, PrecisionMP mp2) { return mp1.bits()==mp2.bits(); }
    friend bool operator<=(PrecisionMP mp1, PrecisionMP mp2) { return mp1.bits()<=mp2.bits(); }
    friend OutputStream& operator<<(OutputStream& os, PrecisionMP mp) { return os << "PrecisionMP("<<mp.bits()<<")"; }
};

//! \ingroup FltMPSubModule
//! \brief Multiple-precision floating-point numbers.
//! Currently defined as _a wrapper around \c mpfr_t from the MPFE library.
//! Default arithmetic operations are approximate, and comparisons are exact, so this class is \em unsafe.
class FloatMP {
  private:
    mpfr_t _mpfr;
    typedef decltype(_mpfr[0]) MpfrReference;
  public:
    typedef RawTag Paradigm;
    typedef FloatMP NumericType;
    typedef mpfr_exp_t ExponentType;
    typedef PrecisionMP PrecisionType;
    typedef RoundingModeMP RoundingModeType;
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
    static Void set_default_precision(PrecisionType prec);
    static PrecisionType get_default_precision();
  public:
    static FloatMP nan(PrecisionType);
    static FloatMP inf(PrecisionType);
    static FloatMP max(PrecisionType);
    static FloatMP eps(PrecisionType);
    static FloatMP min(PrecisionType);
  public:
    ~FloatMP();
    explicit FloatMP(NoInit);
    explicit FloatMP(PrecisionType, NoInit);

    FloatMP();
    FloatMP(double);

    explicit FloatMP(PrecisionType);
    explicit FloatMP(double, PrecisionType);
    explicit FloatMP(Float64, PrecisionType);
    explicit FloatMP(Dyadic const&, PrecisionType);

    FloatMP(const FloatMP&);
    FloatMP(FloatMP&&);

    FloatMP& operator=(const FloatMP&);
    FloatMP& operator=(FloatMP&&);

    FloatMP(Int32 n, PrecisionMP pr);
    FloatMP(double, RoundingModeType, PrecisionType);
    FloatMP(Float64, RoundingModeType, PrecisionType);
    FloatMP(Integer const&, RoundingModeType, PrecisionType);
    FloatMP(Dyadic const&, RoundingModeType, PrecisionType);
    FloatMP(Rational const&, RoundingModeType, PrecisionType);
    FloatMP(FloatMP const&, RoundingModeType, PrecisionType);
    explicit operator Dyadic() const;
    explicit operator Rational() const;

    ExponentType exponent() const;
    PrecisionMP precision() const;
    Void set_precision(PrecisionMP);
  public:
    FloatMP const& raw() const;
    MpfrReference get_mpfr();
    MpfrReference get_mpfr() const;
    double get_d() const;
  public:
    friend Bool is_nan(FloatMP const& x);
    friend Bool is_inf(FloatMP const& x);
    friend Bool is_finite(FloatMP const& x);

    friend FloatMP next_up(FloatMP const& x);
    friend FloatMP next_down(FloatMP const& x);

    friend FloatMP floor(FloatMP const& x);
    friend FloatMP ceil(FloatMP const& x);
    friend FloatMP round(FloatMP const& x);
  public:
    // Explcitly rounded operations
    friend FloatMP nul(FloatMP const& x, RoundingModeType rnd);
    friend FloatMP hlf(FloatMP const& x, RoundingModeType rnd);
    friend FloatMP pos(FloatMP const& x, RoundingModeType rnd);
    friend FloatMP neg(FloatMP const& x, RoundingModeType rnd);
    friend FloatMP sqr(FloatMP const& x, RoundingModeType rnd);
    friend FloatMP rec(FloatMP const& x, RoundingModeType rnd);
    friend FloatMP add(FloatMP const& x1, FloatMP const& x2, RoundingModeType rnd);
    friend FloatMP sub(FloatMP const& x1, FloatMP const& x2, RoundingModeType rnd);
    friend FloatMP mul(FloatMP const& x1, FloatMP const& x2, RoundingModeType rnd);
    friend FloatMP div(FloatMP const& x1, FloatMP const& x2, RoundingModeType rnd);
    friend FloatMP fma(FloatMP const& x1, FloatMP const& x2, FloatMP const& x3, RoundingModeType rnd);
    friend FloatMP pow(FloatMP const& x, Int n, RoundingModeType rnd);
    friend FloatMP sqrt(FloatMP const& x, RoundingModeType rnd);
    friend FloatMP exp(FloatMP const& x, RoundingModeType rnd);
    friend FloatMP log(FloatMP const& x, RoundingModeType rnd);
    friend FloatMP sin(FloatMP const& x, RoundingModeType rnd);
    friend FloatMP cos(FloatMP const& x, RoundingModeType rnd);
    friend FloatMP tan(FloatMP const& x, RoundingModeType rnd);
    friend FloatMP asin(FloatMP const& x, RoundingModeType rnd);
    friend FloatMP acos(FloatMP const& x, RoundingModeType rnd);
    friend FloatMP atan(FloatMP const& x, RoundingModeType rnd);
    static FloatMP pi(PrecisionMP pr, RoundingModeType rnd);

    friend FloatMP max(FloatMP const& x1, FloatMP const& x2, RoundingModeType rnd);
    friend FloatMP min(FloatMP const& x1, FloatMP const& x2, RoundingModeType rnd);
    friend FloatMP abs(FloatMP const& x, RoundingModeType rnd);
    friend FloatMP mag(FloatMP const& x, RoundingModeType rnd);

    // Mixed operations
    friend FloatMP add(FloatMP const& x1, Dbl x2, RoundingModeType rnd);
    friend FloatMP sub(FloatMP const& x1, Dbl x2, RoundingModeType rnd);
    friend FloatMP mul(FloatMP const& x1, Dbl x2, RoundingModeType rnd);
    friend FloatMP div(FloatMP const& x1, Dbl x2, RoundingModeType rnd);
    friend FloatMP add(Dbl x1, FloatMP const& x2, RoundingModeType rnd);
    friend FloatMP sub(Dbl x1, FloatMP const& x2, RoundingModeType rnd);
    friend FloatMP mul(Dbl x1, FloatMP const& x2, RoundingModeType rnd);
    friend FloatMP div(Dbl x1, FloatMP const& x2, RoundingModeType rnd);

    // Operations in current rounding mode
    friend FloatMP nul(FloatMP const& x);
    friend FloatMP hlf(FloatMP const& x);
    friend FloatMP pos(FloatMP const& x);
    friend FloatMP neg(FloatMP const& x);
    friend FloatMP sqr(FloatMP const& x);
    friend FloatMP rec(FloatMP const& x);
    friend FloatMP add(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP sub(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP mul(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP div(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP fma(FloatMP const& x1, FloatMP const& x2, FloatMP const& x3);
    friend FloatMP pow(FloatMP const& x, Int n);
    friend FloatMP sqrt(FloatMP const& x);
    friend FloatMP exp(FloatMP const& x);
    friend FloatMP log(FloatMP const& x);
    friend FloatMP sin(FloatMP const& x);
    friend FloatMP cos(FloatMP const& x);
    friend FloatMP tan(FloatMP const& x);
    friend FloatMP asin(FloatMP const& x);
    friend FloatMP acos(FloatMP const& x);
    friend FloatMP atan(FloatMP const& x);
    static FloatMP pi(PrecisionMP pr);

    friend FloatMP max(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP min(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP abs(FloatMP const& x);
    friend FloatMP mag(FloatMP const& x);

    // Mixed operations
    friend FloatMP add(FloatMP const& x1, Dbl x2);
    friend FloatMP sub(FloatMP const& x1, Dbl x2);
    friend FloatMP mul(FloatMP const& x1, Dbl x2);
    friend FloatMP div(FloatMP const& x1, Dbl x2);
    friend FloatMP add(Dbl x1, FloatMP const& x2);
    friend FloatMP sub(Dbl x1, FloatMP const& x2);
    friend FloatMP mul(Dbl x1, FloatMP const& x2);
    friend FloatMP div(Dbl x1, FloatMP const& x2);

    // Operators
    friend FloatMP operator+(FloatMP const& x);
    friend FloatMP operator-(FloatMP const& x);
    friend FloatMP operator+(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP operator-(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP operator*(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP operator/(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP& operator+=(FloatMP& x1, FloatMP const& x2);
    friend FloatMP& operator-=(FloatMP& x1, FloatMP const& x2);
    friend FloatMP& operator*=(FloatMP& x1, FloatMP const& x2);
    friend FloatMP& operator/=(FloatMP& x1, FloatMP const& x2);

    friend Comparison cmp(FloatMP const& x1, FloatMP const& x2);

    friend Bool operator==(FloatMP const& x1, FloatMP const& x2);
    friend Bool operator!=(FloatMP const& x1, FloatMP const& x2);
    friend Bool operator<=(FloatMP const& x1, FloatMP const& x2);
    friend Bool operator>=(FloatMP const& x1, FloatMP const& x2);
    friend Bool operator< (FloatMP const& x1, FloatMP const& x2);
    friend Bool operator> (FloatMP const& x1, FloatMP const& x2);

    friend OutputStream& operator<<(OutputStream& os, FloatMP const& x);
    friend InputStream& operator>>(InputStream& is, FloatMP& x);

    // Mixed operators
    friend FloatMP operator+(FloatMP const& x1, Dbl x2);
    friend FloatMP operator-(FloatMP const& x1, Dbl x2);
    friend FloatMP operator*(FloatMP const& x1, Dbl x2);
    friend FloatMP operator/(FloatMP const& x1, Dbl x2);
    friend FloatMP operator+(Dbl x1, FloatMP const& x2);
    friend FloatMP operator-(Dbl x1, FloatMP const& x2);
    friend FloatMP operator*(Dbl x1, FloatMP const& x2);
    friend FloatMP operator/(Dbl x1, FloatMP const& x2);

    friend Comparison cmp(FloatMP const& x1, Dbl x2);
    friend Bool operator==(FloatMP const& x1, Dbl x2);
    friend Bool operator!=(FloatMP const& x1, Dbl x2);
    friend Bool operator<=(FloatMP const& x1, Dbl x2);
    friend Bool operator>=(FloatMP const& x1, Dbl x2);
    friend Bool operator< (FloatMP const& x1, Dbl x2);
    friend Bool operator> (FloatMP const& x1, Dbl x2);

    friend Comparison cmp(FloatMP const& x1, Float64 const& x2);
    friend Bool operator==(FloatMP const& x1, Float64 const&  x2);
    friend Bool operator!=(FloatMP const& x1, Float64 const&  x2);
    friend Bool operator<=(FloatMP const& x1, Float64 const&  x2);
    friend Bool operator>=(FloatMP const& x1, Float64 const&  x2);
    friend Bool operator< (FloatMP const& x1, Float64 const&  x2);
    friend Bool operator> (FloatMP const& x1, Float64 const&  x2);

    friend Comparison cmp(FloatMP const& x, Rational const& q);
  private:
    friend OutputStream& write(OutputStream& os, FloatMP const& x, DecimalPlaces dgts, RoundingModeMP rnd);
    friend String print(FloatMP const& x, DecimalPlaces dgts, RoundingModeMP rnd);
    friend String print(FloatMP const& x, DecimalPrecision dgts, RoundingModeMP rnd);
};

template<class R, class A> R integer_cast(const A& a);
template<> inline Int integer_cast(const FloatMP& x) { return static_cast<Int>(x.get_d()); }
template<> inline Nat integer_cast(const FloatMP& x) { return static_cast<Nat>(x.get_d()); }


} // namespace Ariadne

#endif
