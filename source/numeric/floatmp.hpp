/***************************************************************************
 *            numeric/floatmp.hpp
 *
 *  Copyright  2013-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file numeric/floatmp.hpp
 *  \brief RawTag floating-point numbers based on MPFT floats.
 */



#ifndef ARIADNE_FLOATMP_HPP
#define ARIADNE_FLOATMP_HPP

#include "paradigm.hpp"
#include "number.hpp"
#include "rounding.hpp"
#include "sign.hpp"
#include <mpfr.h>

namespace Ariadne {

/************ FloatMP ********************************************************/

struct NoInit { };
struct RawPtr { };

struct RoundUpward; struct RoundDownward;

typedef mpfr_rnd_t RoundingModeMP;

//! \ingroup NumericModule
//! \brief The precision of a FloatMP object.
//! \relates FloatMP
class MultiplePrecision {
    mpfr_prec_t prec;
    typedef std::make_unsigned<mpfr_prec_t>::type unsigned_mpfr_prec_t;
  public:
    typedef unsigned_mpfr_prec_t Type;
    //! \brief
    explicit MultiplePrecision(mpfr_prec_t pr) : prec(pr) { }
    //! \brief
    explicit MultiplePrecision(DoublePrecision const& pr) : prec(53u) { }
    //! \brief The number of binary digits of precision requested.
    unsigned_mpfr_prec_t bits() const { return static_cast<unsigned_mpfr_prec_t>(prec); }
    operator mpfr_prec_t () const { return prec; }
    //! \brief .
    friend MultiplePrecision max(MultiplePrecision mp1, MultiplePrecision mp2) { return MultiplePrecision(static_cast<mpfr_prec_t>(std::max(mp1.bits(),mp2.bits()))); }
    //! \brief .
    friend MultiplePrecision min(MultiplePrecision mp1, MultiplePrecision mp2) { return MultiplePrecision(static_cast<mpfr_prec_t>(std::min(mp1.bits(),mp2.bits()))); }
    //! \brief .
    friend Bool operator==(MultiplePrecision mp1, MultiplePrecision mp2) { return mp1.bits()==mp2.bits(); }
    //! \brief .
    friend Bool operator<=(MultiplePrecision mp1, MultiplePrecision mp2) { return mp1.bits()<=mp2.bits(); }
    //! \brief .
    friend OutputStream& operator<<(OutputStream& os, MultiplePrecision mp) { return os << "MultiplePrecision("<<mp.bits()<<")"; }
};
using MP = MultiplePrecision;
inline MultiplePrecision multiple_precision(mpfr_prec_t pr) { return MultiplePrecision(pr); }
inline MultiplePrecision precision(mpfr_prec_t pr) { return MultiplePrecision(pr); }
inline MP mp(mpfr_prec_t pr) { return MP(pr); }

//! \ingroup FltMPSubModule
//! \brief Multiple-precision floating-point numbers.
//! Currently defined as a wrapper around \c mpfr_t from the MPFR library.
//! Default arithmetic operations are approximate, and comparisons are exact, so this class is \em unsafe.
//! \see FloatDP, Real, FloatValue, FloatBounds, FloatUpperBound, FloatLowerBound, FloatApproximation, FloatError
class FloatMP {
  private:
    mpfr_t _mpfr;
  public:
    typedef RawTag Paradigm;
    typedef FloatMP NumericType;
    typedef mpfr_exp_t ExponentType;
    //! \brief The type representing the precision of the number.
    typedef MultiplePrecision PrecisionType;
    //! \brief The type representing the rounding modes available for arithmetic.
    typedef MPFRRoundingModeType RoundingModeType;
  public:
    static const RoundingModeType ROUND_TO_NEAREST;
    static const RoundingModeType ROUND_DOWNWARD;
    static const RoundingModeType ROUND_UPWARD;
    static const RoundingModeType ROUND_TOWARD_ZERO;

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
    static FloatMP inf(Sign, PrecisionType);
    static FloatMP inf(PrecisionType);
    static FloatMP nan(PrecisionType);

    static FloatMP max(PrecisionType);
    static FloatMP eps(PrecisionType);
    static FloatMP min(PrecisionType);
  public:
    ~FloatMP();
    explicit FloatMP(NoInit);
    explicit FloatMP(PrecisionType, NoInit);
    explicit FloatMP(const mpfr_t, RawPtr);

    FloatMP();
    FloatMP(double);

    explicit FloatMP(PrecisionType);
    explicit FloatMP(double, PrecisionType);
    explicit FloatMP(FloatDP const&, PrecisionType);
    explicit FloatMP(ExactDouble const& x, PrecisionType);
    explicit FloatMP(TwoExp const& x, PrecisionType);
    explicit FloatMP(Dyadic const&, PrecisionType);

    FloatMP(const FloatMP&);
    FloatMP(FloatMP&&);

    FloatMP& operator=(const FloatMP&);
    FloatMP& operator=(FloatMP&&);

    FloatMP(double, RoundingModeType, PrecisionType);
    FloatMP(FloatDP const&, RoundingModeType, PrecisionType);

    FloatMP(Integer const&, RoundingModeType, PrecisionType);
    FloatMP(Dyadic const&, RoundingModeType, PrecisionType);
    FloatMP(Decimal const&, RoundingModeType, PrecisionType);
    FloatMP(Rational const&, RoundingModeType, PrecisionType);
    FloatMP(FloatMP const&, RoundingModeType, PrecisionType);
    explicit operator Dyadic() const;
    explicit operator Rational() const;

    ExponentType exponent() const;
    MultiplePrecision precision() const;
    Void set_precision(MultiplePrecision);
  public:
    FloatMP const& raw() const;
    mpfr_t const& get_mpfr() const;
    mpfr_t& get_mpfr();
    double get_d() const;
  public:
    friend Bool is_nan(FloatMP const& x);
    friend Bool is_inf(FloatMP const& x);
    friend Bool is_finite(FloatMP const& x);
    friend Bool is_zero(FloatMP const& x);

    friend FloatMP next(RoundUpward rnd, FloatMP const& x);
    friend FloatMP next(RoundDownward rnd, FloatMP const& x);

    friend FloatMP floor(FloatMP const& x);
    friend FloatMP ceil(FloatMP const& x);
    friend FloatMP round(FloatMP const& x);
  public:
    // Exact lattice operations
    friend FloatMP max(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP min(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP abs(FloatMP const& x);
    friend FloatMP mag(FloatMP const& x);

    // Exact arithmetic
    friend FloatMP nul(FloatMP const& x);
    friend FloatMP pos(FloatMP const& x);
    friend FloatMP neg(FloatMP const& x);
    friend FloatMP hlf(FloatMP const& x);
    friend FloatMP shft(FloatMP const& x, Int n);

    // Exact arithmetic operators
    friend FloatMP operator+(FloatMP const& x);
    friend FloatMP operator-(FloatMP const& x);

    // Explcitly rounded lattice operations
    friend FloatMP max(RoundingModeType rnd, FloatMP const& x1, FloatMP const& x2);
    friend FloatMP min(RoundingModeType rnd, FloatMP const& x1, FloatMP const& x2);
    friend FloatMP abs(RoundingModeType rnd, FloatMP const& x);
    friend FloatMP mag(RoundingModeType rnd, FloatMP const& x);

    // Explcitly rounded inplace operations
    friend FloatMP& iadd(RoundingModeType rnd, FloatMP& x1, FloatMP const& x2);
    friend FloatMP& isub(RoundingModeType rnd, FloatMP& x1, FloatMP const& x2);
    friend FloatMP& imul(RoundingModeType rnd, FloatMP& x1, FloatMP const& x2);
    friend FloatMP& idiv(RoundingModeType rnd, FloatMP& x1, FloatMP const& x2);

    // Explcitly rounded operations
    friend FloatMP nul(RoundingModeType rnd, FloatMP const& x);
    friend FloatMP pos(RoundingModeType rnd, FloatMP const& x);
    friend FloatMP neg(RoundingModeType rnd, FloatMP const& x);
    friend FloatMP hlf(RoundingModeType rnd, FloatMP const& x);
    friend FloatMP shft(RoundingModeType rnd, FloatMP const& x, Int n);
    friend FloatMP add(RoundingModeType rnd, FloatMP const& x1, FloatMP const& x2);
    friend FloatMP sub(RoundingModeType rnd, FloatMP const& x1, FloatMP const& x2);
    friend FloatMP mul(RoundingModeType rnd, FloatMP const& x1, FloatMP const& x2);
    friend FloatMP div(RoundingModeType rnd, FloatMP const& x1, FloatMP const& x2);
    friend FloatMP fma(RoundingModeType rnd, FloatMP const& x1, FloatMP const& x2, FloatMP const& x3); // x1*x2+x3
    friend FloatMP pow(RoundingModeType rnd, FloatMP const& x, Int n);
    friend FloatMP sqr(RoundingModeType rnd, FloatMP const& x);
    friend FloatMP rec(RoundingModeType rnd, FloatMP const& x);
    friend FloatMP sqrt(RoundingModeType rnd, FloatMP const& x);
    friend FloatMP exp(RoundingModeType rnd, FloatMP const& x);
    friend FloatMP log(RoundingModeType rnd, FloatMP const& x);
    friend FloatMP sin(RoundingModeType rnd, FloatMP const& x);
    friend FloatMP cos(RoundingModeType rnd, FloatMP const& x);
    friend FloatMP tan(RoundingModeType rnd, FloatMP const& x);
    friend FloatMP asin(RoundingModeType rnd, FloatMP const& x);
    friend FloatMP acos(RoundingModeType rnd, FloatMP const& x);
    friend FloatMP atan(RoundingModeType rnd, FloatMP const& x);
    static FloatMP pi(RoundingModeType rnd, MultiplePrecision pr);

    friend FloatMP med(RoundingModeType rnd, FloatMP x1, FloatMP x2) { return hlf(add(rnd,x1,x2)); }
    friend FloatMP rad(RoundingModeType rnd, FloatMP x1, FloatMP x2) { return hlf(sub(rnd,x2,x1)); }

    // Mixed operations
    friend FloatMP add(RoundingModeType rnd, FloatMP const& x1, Dbl x2);
    friend FloatMP sub(RoundingModeType rnd, FloatMP const& x1, Dbl x2);
    friend FloatMP mul(RoundingModeType rnd, FloatMP const& x1, Dbl x2);
    friend FloatMP div(RoundingModeType rnd, FloatMP const& x1, Dbl x2);
    friend FloatMP add(RoundingModeType rnd, Dbl x1, FloatMP const& x2);
    friend FloatMP sub(RoundingModeType rnd, Dbl x1, FloatMP const& x2);
    friend FloatMP mul(RoundingModeType rnd, Dbl x1, FloatMP const& x2);
    friend FloatMP div(RoundingModeType rnd, Dbl x1, FloatMP const& x2);

    friend FloatMP add(RoundingModeType rnd, FloatMP const& x1, FloatDP const& x2);
    friend FloatMP sub(RoundingModeType rnd, FloatMP const& x1, FloatDP const& x2);
    friend FloatMP mul(RoundingModeType rnd, FloatMP const& x1, FloatDP const& x2);
    friend FloatMP div(RoundingModeType rnd, FloatMP const& x1, FloatDP const& x2);
    friend FloatMP add(RoundingModeType rnd, FloatDP const& x1, FloatMP const& x2);
    friend FloatMP sub(RoundingModeType rnd, FloatDP const& x1, FloatMP const& x2);
    friend FloatMP mul(RoundingModeType rnd, FloatDP const& x1, FloatMP const& x2);
    friend FloatMP div(RoundingModeType rnd, FloatDP const& x1, FloatMP const& x2);

    friend FloatMP add(RoundUpward rnd, FloatMP const& x1, FloatDP const& x2);
    friend FloatMP sub(RoundDownward rnd, FloatMP const& x1, FloatDP const& x2);

    // Correctly rounded arithmetic
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
    static FloatMP pi(MultiplePrecision pr);

    // Correctly rounded operators
    friend FloatMP operator+(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP operator-(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP operator*(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP operator/(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP& operator+=(FloatMP& x1, FloatMP const& x2);
    friend FloatMP& operator-=(FloatMP& x1, FloatMP const& x2);
    friend FloatMP& operator*=(FloatMP& x1, FloatMP const& x2);
    friend FloatMP& operator/=(FloatMP& x1, FloatMP const& x2);

    // Correctly rounded  mixed operators
    friend FloatMP operator+(FloatMP const& x1, Dbl x2);
    friend FloatMP operator-(FloatMP const& x1, Dbl x2);
    friend FloatMP operator*(FloatMP const& x1, Dbl x2);
    friend FloatMP operator/(FloatMP const& x1, Dbl x2);
    friend FloatMP operator+(Dbl x1, FloatMP const& x2);
    friend FloatMP operator-(Dbl x1, FloatMP const& x2);
    friend FloatMP operator*(Dbl x1, FloatMP const& x2);
    friend FloatMP operator/(Dbl x1, FloatMP const& x2);

    friend FloatMP operator+(FloatMP const& x1, FloatDP const& x2);
    friend FloatMP operator-(FloatMP const& x1, FloatDP const& x2);
    friend FloatMP operator*(FloatMP const& x1, FloatDP const& x2);
    friend FloatMP operator/(FloatMP const& x1, FloatDP const& x2);
    friend FloatMP operator+(FloatDP const& x1, FloatMP const& x2);
    friend FloatMP operator-(FloatDP const& x1, FloatMP const& x2);
    friend FloatMP operator*(FloatDP const& x1, FloatMP const& x2);
    friend FloatMP operator/(FloatDP const& x1, FloatMP const& x2);

    friend Comparison cmp(FloatMP const& x1, FloatMP const& x2);
    friend Bool operator==(FloatMP const& x1, FloatMP const& x2);
    friend Bool operator!=(FloatMP const& x1, FloatMP const& x2);
    friend Bool operator<=(FloatMP const& x1, FloatMP const& x2);
    friend Bool operator>=(FloatMP const& x1, FloatMP const& x2);
    friend Bool operator< (FloatMP const& x1, FloatMP const& x2);
    friend Bool operator> (FloatMP const& x1, FloatMP const& x2);

    friend Comparison cmp(FloatMP const& x1, Dbl x2);
    friend Bool operator==(FloatMP const& x1, Dbl x2) { return cmp(x1,x2)==Comparison::EQUAL; }
    friend Bool operator!=(FloatMP const& x1, Dbl x2) { return cmp(x1,x2)!=Comparison::EQUAL; }
    friend Bool operator<=(FloatMP const& x1, Dbl x2) { return cmp(x1,x2)<=Comparison::EQUAL; }
    friend Bool operator>=(FloatMP const& x1, Dbl x2) { return cmp(x1,x2)>=Comparison::EQUAL; }
    friend Bool operator< (FloatMP const& x1, Dbl x2) { return cmp(x1,x2)< Comparison::EQUAL; }
    friend Bool operator> (FloatMP const& x1, Dbl x2) { return cmp(x1,x2)> Comparison::EQUAL; }

    friend Comparison cmp(Dbl x1, FloatMP const& x2);
    friend Bool operator==(Dbl x1, FloatMP const& x2) { return cmp(x1,x2)==Comparison::EQUAL; }
    friend Bool operator!=(Dbl x1, FloatMP const& x2) { return cmp(x1,x2)!=Comparison::EQUAL; }
    friend Bool operator<=(Dbl x1, FloatMP const& x2) { return cmp(x1,x2)<=Comparison::EQUAL; }
    friend Bool operator>=(Dbl x1, FloatMP const& x2) { return cmp(x1,x2)>=Comparison::EQUAL; }
    friend Bool operator< (Dbl x1, FloatMP const& x2) { return cmp(x1,x2)< Comparison::EQUAL; }
    friend Bool operator> (Dbl x1, FloatMP const& x2) { return cmp(x1,x2)> Comparison::EQUAL; }

    friend Comparison cmp(FloatMP const& x1, Rational const& x2);
    friend Bool operator==(FloatMP const& x1, Rational const& x2) { return cmp(x1,x2)==Comparison::EQUAL; }
    friend Bool operator!=(FloatMP const& x1, Rational const& x2) { return cmp(x1,x2)!=Comparison::EQUAL; }
    friend Bool operator<=(FloatMP const& x1, Rational const& x2) { return cmp(x1,x2)<=Comparison::EQUAL; }
    friend Bool operator>=(FloatMP const& x1, Rational const& x2) { return cmp(x1,x2)>=Comparison::EQUAL; }
    friend Bool operator< (FloatMP const& x1, Rational const& x2) { return cmp(x1,x2)< Comparison::EQUAL; }
    friend Bool operator> (FloatMP const& x1, Rational const& x2) { return cmp(x1,x2)> Comparison::EQUAL; }
    friend Comparison cmp(Rational const& x1, FloatMP const& x2);
    friend Bool operator==(Rational const& x1, FloatMP const& x2) { return cmp(x1,x2)==Comparison::EQUAL; }
    friend Bool operator!=(Rational const& x1, FloatMP const& x2) { return cmp(x1,x2)!=Comparison::EQUAL; }
    friend Bool operator<=(Rational const& x1, FloatMP const& x2) { return cmp(x1,x2)<=Comparison::EQUAL; }
    friend Bool operator>=(Rational const& x1, FloatMP const& x2) { return cmp(x1,x2)>=Comparison::EQUAL; }
    friend Bool operator< (Rational const& x1, FloatMP const& x2) { return cmp(x1,x2)< Comparison::EQUAL; }
    friend Bool operator> (Rational const& x1, FloatMP const& x2) { return cmp(x1,x2)> Comparison::EQUAL; }

    friend Comparison cmp(FloatMP const& x1, FloatDP const& x2);
    friend Bool operator==(FloatMP const& x1, FloatDP const&  x2) { return cmp(x1,x2)==Comparison::EQUAL; }
    friend Bool operator!=(FloatMP const& x1, FloatDP const&  x2) { return cmp(x1,x2)!=Comparison::EQUAL; }
    friend Bool operator<=(FloatMP const& x1, FloatDP const&  x2) { return cmp(x1,x2)<=Comparison::EQUAL; }
    friend Bool operator>=(FloatMP const& x1, FloatDP const&  x2) { return cmp(x1,x2)>=Comparison::EQUAL; }
    friend Bool operator< (FloatMP const& x1, FloatDP const&  x2) { return cmp(x1,x2)< Comparison::EQUAL; }
    friend Bool operator> (FloatMP const& x1, FloatDP const&  x2) { return cmp(x1,x2)> Comparison::EQUAL; }
    friend Comparison cmp(FloatDP const& x1, FloatMP const& x2);
    friend Bool operator==(FloatDP const& x1, FloatMP const& x2) { return cmp(x1,x2)==Comparison::EQUAL; }
    friend Bool operator!=(FloatDP const& x1, FloatMP const& x2) { return cmp(x1,x2)!=Comparison::EQUAL; }
    friend Bool operator<=(FloatDP const& x1, FloatMP const& x2) { return cmp(x1,x2)<=Comparison::EQUAL; }
    friend Bool operator>=(FloatDP const& x1, FloatMP const& x2) { return cmp(x1,x2)>=Comparison::EQUAL; }
    friend Bool operator< (FloatDP const& x1, FloatMP const& x2) { return cmp(x1,x2)< Comparison::EQUAL; }
    friend Bool operator> (FloatDP const& x1, FloatMP const& x2) { return cmp(x1,x2)> Comparison::EQUAL; }

    friend OutputStream& operator<<(OutputStream& os, FloatMP const& x);
    friend InputStream& operator>>(InputStream& is, FloatMP& x);
  private:

    // Mixed operations
    friend FloatMP add(FloatMP const& x1, Dbl x2);
    friend FloatMP sub(FloatMP const& x1, Dbl x2);
    friend FloatMP mul(FloatMP const& x1, Dbl x2);
    friend FloatMP div(FloatMP const& x1, Dbl x2);
    friend FloatMP add(Dbl x1, FloatMP const& x2);
    friend FloatMP sub(Dbl x1, FloatMP const& x2);
    friend FloatMP mul(Dbl x1, FloatMP const& x2);
    friend FloatMP div(Dbl x1, FloatMP const& x2);

    // Correctly rounded operations
    friend FloatMP sqr_rnd(FloatMP const& x);
    friend FloatMP add_rnd(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP sub_rnd(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP mul_rnd(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP div_rnd(FloatMP const& x1, FloatMP const& x2);
    friend FloatMP pow_rnd(FloatMP const& x, Int n);
    friend FloatMP sqrt_rnd(FloatMP const& x);
    friend FloatMP exp_rnd(FloatMP const& x);
    friend FloatMP log_rnd(FloatMP const& x);
    friend FloatMP sin_rnd(FloatMP const& x);
    friend FloatMP cos_rnd(FloatMP const& x);
    friend FloatMP tan_rnd(FloatMP const& x);
    friend FloatMP atan_rnd(FloatMP const& x);

    friend OutputStream& write(OutputStream& os, FloatMP const& x, DecimalPlaces dgts, RoundingModeMP rnd);
    friend OutputStream& write(OutputStream& os, FloatMP const& x, DecimalPrecision dgts, RoundingModeMP rnd);
    friend String print(const mpfr_t x, int zdgts, int fdgts, mpfr_rnd_t rnd);
    friend String print(FloatMP const& x, DecimalPrecision figs, RoundingModeMP rnd);
    friend String print(FloatMP const& x, DecimalPlaces plcs, RoundingModeMP rnd);
};

template<class R, class A> R integer_cast(const A& a);
template<> inline Int integer_cast(const FloatMP& x) { return static_cast<Int>(x.get_d()); }
template<> inline Nat integer_cast(const FloatMP& x) { return static_cast<Nat>(x.get_d()); }


} // namespace Ariadne

#endif
