/***************************************************************************
 *            numeric/floatmp.h
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
 *  You should have received _a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file numeric/floatmp.h
 *  \brief Raw floating-point numbers based on MPFT floats.
 */



#ifndef ARIADNE_FLOATMP_H
#define ARIADNE_FLOATMP_H

#include "paradigm.h"
#include "number.h"
#include "mixins.h"
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
};
inline bool operator==(PrecisionMP mp1, PrecisionMP mp2) { return mp1.bits()==mp2.bits(); }

//! \ingroup FltMPSubModule
//! \brief Multiple-precision floating-point numbers.
//! Currently defined as _a wrapper around \c mpfr_t from the MPFE library.
//! Default arithmetic operations are approximate, and comparisons are exact, so this class is \em unsafe.
class FloatMP {
  private:
    mpfr_t _mpfr;
    typedef decltype(_mpfr[0]) MpfrReference;
  public:
    typedef Raw Paradigm;
    typedef FloatMP NumericType;
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
    static FloatMP nan(PrecisionType=get_default_precision());
    static FloatMP inf(PrecisionType=get_default_precision());
    static FloatMP max(PrecisionType=get_default_precision());
    static FloatMP eps(PrecisionType=get_default_precision());
    static FloatMP min(PrecisionType=get_default_precision());
  public:
    ~FloatMP();
    explicit FloatMP(NoInit);
    explicit FloatMP(PrecisionType, NoInit);

    FloatMP();
    FloatMP(double, PrecisionType=get_default_precision());
    explicit FloatMP(Float64, PrecisionType=get_default_precision());

    explicit FloatMP(PrecisionType);

    FloatMP(const FloatMP&);
    FloatMP(FloatMP&&);

    FloatMP& operator=(const FloatMP&);
    FloatMP& operator=(FloatMP&&);

    FloatMP(Int32 n, PrecisionMP pr);
    FloatMP(Integer const&, RoundingModeType=get_rounding_mode(), PrecisionType=get_default_precision());
    explicit FloatMP(Rational const&, RoundingModeType=get_rounding_mode(), PrecisionType=get_default_precision());
    explicit operator Rational() const;

    PrecisionMP precision() const;
    Void set_precision(PrecisionMP);
  public:
    FloatMP const& raw() const;
    MpfrReference get_mpfr();
    MpfrReference get_mpfr() const;
    double get_d() const;
  public:
    friend FloatMP next_up(FloatMP const& x);
    friend FloatMP next_down(FloatMP const& x);

    friend FloatMP floor(FloatMP const& x);
    friend FloatMP ceil(FloatMP const& x);

    friend FloatMP nul(FloatMP const& x, RoundingModeType=get_rounding_mode());
    friend FloatMP half(FloatMP const& x, RoundingModeType=get_rounding_mode());
    friend FloatMP pos(FloatMP const& x, RoundingModeType=get_rounding_mode());
    friend FloatMP neg(FloatMP const& x, RoundingModeType=get_rounding_mode());
    friend FloatMP sqr(FloatMP const& x, RoundingModeType=get_rounding_mode());
    friend FloatMP rec(FloatMP const& x, RoundingModeType=get_rounding_mode());
    friend FloatMP add(FloatMP const& x1, FloatMP const& x2, RoundingModeType=get_rounding_mode());
    friend FloatMP sub(FloatMP const& x1, FloatMP const& x2, RoundingModeType=get_rounding_mode());
    friend FloatMP mul(FloatMP const& x1, FloatMP const& x2, RoundingModeType=get_rounding_mode());
    friend FloatMP div(FloatMP const& x1, FloatMP const& x2, RoundingModeType=get_rounding_mode());
    friend FloatMP fma(FloatMP const& x1, FloatMP const& x2, FloatMP const& x3, RoundingModeType=get_rounding_mode());
    friend FloatMP pow(FloatMP const& x, Int n, RoundingModeType=get_rounding_mode());
    friend FloatMP sqrt(FloatMP const& x, RoundingModeType=get_rounding_mode());
    friend FloatMP exp(FloatMP const& x, RoundingModeType=get_rounding_mode());
    friend FloatMP log(FloatMP const& x, RoundingModeType=get_rounding_mode());
    friend FloatMP sin(FloatMP const& x, RoundingModeType=get_rounding_mode());
    friend FloatMP cos(FloatMP const& x, RoundingModeType=get_rounding_mode());
    friend FloatMP tan(FloatMP const& x, RoundingModeType=get_rounding_mode());
    friend FloatMP asin(FloatMP const& x, RoundingModeType=get_rounding_mode());
    friend FloatMP acos(FloatMP const& x, RoundingModeType=get_rounding_mode());
    friend FloatMP atan(FloatMP const& x, RoundingModeType=get_rounding_mode());
    static FloatMP pi(PrecisionMP pr=get_default_precision(), RoundingModeType rnd=get_rounding_mode());

    friend FloatMP max(FloatMP const& x1, FloatMP const& x2, RoundingModeType=get_rounding_mode());
    friend FloatMP min(FloatMP const& x1, FloatMP const& x2, RoundingModeType=get_rounding_mode());
    friend FloatMP abs(FloatMP const& x, RoundingModeType=get_rounding_mode());
    friend FloatMP mag(FloatMP const& x, RoundingModeType=get_rounding_mode());

    friend Bool is_nan(FloatMP const& x);
    friend Bool is_inf(FloatMP const& x);

    friend Comparison cmp(FloatMP const& x1, FloatMP const& x2);
    friend Comparison cmp(FloatMP const& x1, Float64 const& x2);

    // Mixed operations
    friend FloatMP add(FloatMP const& x1, Dbl x2, RoundingModeType=get_rounding_mode());
    friend FloatMP sub(FloatMP const& x1, Dbl x2, RoundingModeType=get_rounding_mode());
    friend FloatMP mul(FloatMP const& x1, Dbl x2, RoundingModeType=get_rounding_mode());
    friend FloatMP div(FloatMP const& x1, Dbl x2, RoundingModeType=get_rounding_mode());
    friend FloatMP add(Dbl x1, FloatMP const& x2, RoundingModeType=get_rounding_mode());
    friend FloatMP sub(Dbl x1, FloatMP const& x2, RoundingModeType=get_rounding_mode());
    friend FloatMP mul(Dbl x1, FloatMP const& x2, RoundingModeType=get_rounding_mode());
    friend FloatMP div(Dbl x1, FloatMP const& x2, RoundingModeType=get_rounding_mode());

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

    friend Bool operator==(FloatMP const& x1, FloatMP const& x2);
    friend Bool operator!=(FloatMP const& x1, FloatMP const& x2);
    friend Bool operator<=(FloatMP const& x1, FloatMP const& x2);
    friend Bool operator>=(FloatMP const& x1, FloatMP const& x2);
    friend Bool operator< (FloatMP const& x1, FloatMP const& x2);
    friend Bool operator> (FloatMP const& x1, FloatMP const& x2);

    friend OutputStream& operator<<(OutputStream& os, FloatMP const&);
    friend InputStream& operator>>(InputStream& is, FloatMP&);

    // Mixed operators
    friend FloatMP operator+(FloatMP const& x1, Dbl x2);
    friend FloatMP operator-(FloatMP const& x1, Dbl x2);
    friend FloatMP operator*(FloatMP const& x1, Dbl x2);
    friend FloatMP operator/(FloatMP const& x1, Dbl x2);
    friend FloatMP operator+(Dbl x1, FloatMP const& x2);
    friend FloatMP operator-(Dbl x1, FloatMP const& x2);
    friend FloatMP operator*(Dbl x1, FloatMP const& x2);
    friend FloatMP operator/(Dbl x1, FloatMP const& x2);

};

template<class R, class A> R integer_cast(const A& a);
template<> inline Int integer_cast(const FloatMP& x) { return static_cast<Int>(x.get_d()); }
template<> inline Nat integer_cast(const FloatMP& x) { return static_cast<Nat>(x.get_d()); }


} // namespace Ariadne

#endif
