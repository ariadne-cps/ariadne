/***************************************************************************
 *            numeric/floatdp.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/*! \file numeric/floatdp.hpp
 *  \brief RawTag floating-point number class based on double-precision floats.
 */

#ifndef ARIADNE_FLOAT64_HPP
#define ARIADNE_FLOAT64_HPP

#include <iosfwd> // For std::floor std::ceil etc
#include <cmath> // For std::floor std::ceil etc
#include <algorithm> // For std::max, std::min
#include <limits> // For std::numeric_limits<double>

#include "../utility/declarations.hpp"
#include "../numeric/operators.hpp"
#include "../numeric/rounding.hpp"
#include "../numeric/sign.hpp"
#include "../numeric/number.decl.hpp"
#include "../numeric/float.decl.hpp"

namespace Ariadne {

class FloatDP;
typedef FloatDP RawFloatDP;

class Rational;
enum class Comparison : char;

//! \ingroup NumericModule
//! \brief The precision of a FloatDP object. Since this is fixed, the class is only a tag; all objects are equal.
//! \relates FloatDP
class DoublePrecision {
    //! \brief .
    friend constexpr DoublePrecision max(DoublePrecision, DoublePrecision) { return DoublePrecision(); }
    //! \brief .
    friend constexpr DoublePrecision min(DoublePrecision, DoublePrecision) { return DoublePrecision(); }
    //! \brief .
    friend constexpr Bool operator<=(DoublePrecision, DoublePrecision) { return true; }
    //! \brief Compare the two precisions for equality. Always \c true.
    friend constexpr Bool operator==(DoublePrecision, DoublePrecision) { return true; }
    //! \brief .
    friend OutputStream& operator<<(OutputStream& os, DoublePrecision) { return os << "DoublePrecision()"; }
};
using DP = DoublePrecision;
static const DoublePrecision double_precision = DoublePrecision();
static const DoublePrecision dp = DP();

// Correctly rounded functions
double sqr_rnd(double x);
double rec_rnd(double x);
double add_rnd(double x1, double x2);
double sub_rnd(double x1, double x2);
double mul_rnd(double x1, double x2);
double div_rnd(double x1, double x2);
double fma_rnd(double x1, double x2, double x3);
double pow_rnd(double x, Nat n);
double pow_rnd(double x, int n);
double sqrt_rnd(double x);
double exp_rnd(double x);
double log_rnd(double x);
double sin_rnd(double x);
double cos_rnd(double x);
double tan_rnd(double x);
double asin_rnd(double x);
double acos_rnd(double x);
double atan_rnd(double x);
double neg_rec_rnd(double x);
double atan_rnd_series(double x);
double pi_rnd();

double texp(double x);

double pi_opp();
double add_opp(double x, double y);
double sub_opp(double x, double y);
double mul_opp(double x, double y);
double div_opp(double x, double y);
double neg_rec_opp(double x);

//! \ingroup NumericModule
//! \brief Floating point numbers (double precision) using rounded arithmetic.
//! \details
//! The \c %FloatDP class represents floating-point numbers.
//! Unless otherwise mentioned, operations on floating-point numbers are performed approximately, with no guarantees
//! on the output.
//!
//! To implement <em>rounded arithmetic</em>, arithmetical operations of \c %FloatDP can be performed with guaranteed rounding by
//! specifying \c _up and \c _down suffixes to arithmetical functions \c add, \c sub, \c mul and \c div.
//! Additionally, operations can be performed in the current <em>rounding mode</em> by using the \c _rnd suffix,
//! or with rounding reversed using the \c _opp suffix.
//! The \c _approx suffix is provided to specifically indicate that the operation is computed approximately.
//!
//! %Ariadne floating-point numbers can be constructed by conversion from built-in C++ types.
//! Note that the value of a built-in floating-point value may differ from the mathematical value of the literal.
//! For example, while <c>%FloatDP(3.25)</c> is represented exactly, <c>%FloatDP(3.3)</c> has a value of \f$3.2999999999999998224\ldots\f$.
//! \note In the future, the construction of a \c %FloatDP from a string literal may be supported.
//! \sa FloatMP, Real, FloatDPValue, FloatDPBounds, FloatDPUpperBound, FloatDPLowerBound, FloatDPApproximation
class FloatDP {
  public:
    volatile double dbl;
  public:
    typedef RawTag Paradigm;
    typedef FloatDP NumericType;
    //! \brief The type representing the precision of the number.
    typedef DoublePrecision PrecisionType;
    //! \brief The type representing the rounding modes available for arithmetic.
    typedef BuiltinRoundingModeType RoundingModeType;
  public:
    static const RoundingModeType ROUND_TO_NEAREST = Ariadne::ROUND_TO_NEAREST;
    static const RoundingModeType ROUND_DOWNWARD = Ariadne::ROUND_DOWNWARD;
    static const RoundingModeType ROUND_UPWARD = Ariadne::ROUND_UPWARD;
    static const RoundingModeType ROUND_TOWARD_ZERO = Ariadne::ROUND_TOWARD_ZERO;

    static RoundingModeType get_rounding_mode();
    static Void set_rounding_mode(RoundingModeType);
    static Void set_rounding_downward();
    static Void set_rounding_upward();
    static Void set_rounding_to_nearest();
    static Void set_rounding_toward_zero();
  public:
    static DoublePrecision get_default_precision();
    DoublePrecision precision() const;
    Void set_precision(DoublePrecision);
  public:

    static FloatDP nan(DoublePrecision pr);
    static FloatDP inf(Sign sgn, DoublePrecision pr);
    static FloatDP inf(DoublePrecision pr);

    static FloatDP max(DoublePrecision pr);
    static FloatDP eps(DoublePrecision pr);
    static FloatDP min(DoublePrecision pr);
  public:
    //! \brief Default constructor creates an uninitialised number.
    FloatDP() : dbl() { }
    explicit FloatDP(DoublePrecision) : dbl() { }
    //! \brief Convert from a built-in double-precision floating-point number.
    FloatDP(double x) : dbl(x) { }
    explicit FloatDP(double x, DoublePrecision) : dbl(x) { }
    explicit FloatDP(ExactDouble const& x, DoublePrecision);
    explicit FloatDP(TwoExp const& x, DoublePrecision);
    explicit FloatDP(Dyadic const& x, DoublePrecision);
    //! \brief Copy constructor.
    FloatDP(const FloatDP& x) : dbl(x.dbl) { }
    //! \brief Copy assignment.
    FloatDP& operator=(const FloatDP& x) { this->dbl=x.dbl; return *this; }

    //! \brief Construct from a double number using given rounding
    explicit FloatDP(double d, RoundingModeType rnd, PrecisionType pr);
    //! \brief Construct from another FloatDP using given rounding
    explicit FloatDP(FloatDP const& d, RoundingModeType rnd, PrecisionType pr);
    //! \brief Construct from an integer number using given rounding
    explicit FloatDP(Integer const&, RoundingModeType rnd, PrecisionType pr);
    //! \brief Construct from a dyadic number with given rounding
    explicit FloatDP(Dyadic const& w, RoundingModeType rnd, PrecisionType pr);
    //! \brief Construct from a decimal number with given rounding
    explicit FloatDP(Decimal const& d, RoundingModeType rnd, PrecisionType pr);
    //! \brief Construct from a FloatMP using given rounding
    explicit FloatDP(FloatMP const& d, RoundingModeType rnd, PrecisionType pr);
    //! \brief Construct from a rational number with given rounding
    explicit FloatDP(Rational const& q, RoundingModeType rnd, PrecisionType pr);
    //! \brief Convert to a dyadic number.
    explicit operator Dyadic () const;
    //! \brief Convert to a rational number.
    explicit operator Rational () const;
  public:
    FloatDP const& raw() const { return *this; }
    //! \brief An approximation by a built-in double-precision floating-point number.
    double get_d() const { return this->dbl; }
  public:
    friend Bool is_nan(FloatDP x) { return std::isnan(x.dbl); }
    friend Bool is_inf(FloatDP x) { return std::isinf(x.dbl); }
    friend Bool is_finite(FloatDP x) { return std::isfinite(x.dbl); }
    friend Bool is_zero(FloatDP x) { return x.dbl==0.0; }

    friend FloatDP next(RoundUpward rnd, FloatDP x) { return add(rnd,x,FloatDP::min(x.precision())); }
    friend FloatDP next(RoundDownward rnd, FloatDP x) { return sub(rnd,x,FloatDP::min(x.precision())); }

    friend FloatDP floor(FloatDP x) { return std::floor(x.dbl); }
    friend FloatDP ceil(FloatDP x) { return std::ceil(x.dbl); }
    friend FloatDP round(FloatDP x) { return std::round(x.dbl); }
  public:
    // Exact lattic operations
    friend FloatDP max(FloatDP x1, FloatDP x2) { return std::max(x1.dbl,x2.dbl); }
    friend FloatDP min(FloatDP x1, FloatDP x2) { return std::min(x1.dbl,x2.dbl); }
    friend FloatDP abs(FloatDP x) { return std::fabs(x.dbl); }
    friend FloatDP mag(FloatDP x) { return std::fabs(x.dbl); }

    // Exact arithmetic
    friend FloatDP nul(FloatDP x) { return +0.0; }
    friend FloatDP pos(FloatDP x) { volatile double xv=x.dbl; return +xv; }
    friend FloatDP neg(FloatDP x) { volatile double xv=x.dbl; return -xv; }
    friend FloatDP hlf(FloatDP x) { volatile double xv=x.dbl; return xv/2; }
    friend FloatDP shft(FloatDP x, Int n) { volatile double xv=x.dbl; return ldexp(xv,n); }

    // Exact arithmetic operators
    friend FloatDP operator+(FloatDP x) { volatile double xv=x.dbl; return +xv; }
    friend FloatDP operator-(FloatDP x) { volatile double xv=x.dbl; return -xv; }

    // Explcitly rounded lattice operations
    friend FloatDP max(RoundingModeType rnd, FloatDP x1, FloatDP x2) { return std::max(x1.dbl,x2.dbl); }
    friend FloatDP min(RoundingModeType rnd, FloatDP x1, FloatDP x2) { return std::min(x1.dbl,x2.dbl); }
    friend FloatDP abs(RoundingModeType rnd, FloatDP x) { return std::fabs(x.dbl); }
    friend FloatDP mag(RoundingModeType rnd, FloatDP x) { return std::fabs(x.dbl); }

    // Explicitly rounded inplace operations
    friend FloatDP& iadd(RoundingModeType rnd, FloatDP& x1, FloatDP x2);
    friend FloatDP& isub(RoundingModeType rnd, FloatDP& x1, FloatDP x2);
    friend FloatDP& imul(RoundingModeType rnd, FloatDP& x1, FloatDP x2);
    friend FloatDP& idiv(RoundingModeType rnd, FloatDP& x1, FloatDP x2);

    static FloatDP pi(RoundingModeType rnd, PrecisionType pr);

    // Correctly rounded arithmetic
    friend FloatDP sqr(FloatDP x) { return sqr_rnd(x.dbl); }
    friend FloatDP rec(FloatDP x) { return rec_rnd(x.dbl); }
    friend FloatDP add(FloatDP x1, FloatDP x2) { return add_rnd(x1.dbl,x2.dbl); }
    friend FloatDP sub(FloatDP x1, FloatDP x2) { return sub_rnd(x1.dbl,x2.dbl); }
    friend FloatDP mul(FloatDP x1, FloatDP x2) { return mul_rnd(x1.dbl,x2.dbl); }
    friend FloatDP div(FloatDP x1, FloatDP x2) { return div_rnd(x1.dbl,x2.dbl); }
    friend FloatDP fma(FloatDP x1, FloatDP x2, FloatDP x3) { return fma_rnd(x1.dbl,x2.dbl,x3.dbl); }
    friend FloatDP pow(FloatDP x, Int n) { return pow_rnd(x.dbl,n); }
    friend FloatDP sqrt(FloatDP x) { return sqrt_rnd(x.dbl); }
    friend FloatDP exp(FloatDP x) { return exp_rnd(x.dbl); }
    friend FloatDP log(FloatDP x) { return log_rnd(x.dbl); }
    friend FloatDP sin(FloatDP x) { return sin_rnd(x.dbl); }
    friend FloatDP cos(FloatDP x) { return cos_rnd(x.dbl); }
    friend FloatDP tan(FloatDP x) { return tan_rnd(x.dbl); }
    friend FloatDP asin(FloatDP x) { return asin_rnd(x.dbl); }
    friend FloatDP acos(FloatDP x) { return acos_rnd(x.dbl); }
    friend FloatDP atan(FloatDP x) { return atan_rnd(x.dbl); }
    static FloatDP pi(PrecisionType pr) { return pi_rnd(); }

    // Correctly rounded operators
    friend FloatDP operator+(FloatDP x1, FloatDP x2) { volatile double x1v = x1.dbl; volatile double x2v=x2.dbl; volatile double r=x1v+x2v; return r; }
    friend FloatDP operator-(FloatDP x1, FloatDP x2) { volatile double x1v = x1.dbl; volatile double x2v=x2.dbl; volatile double r=x1v-x2v; return r; }
    friend FloatDP operator*(FloatDP x1, FloatDP x2) { volatile double x1v = x1.dbl; volatile double x2v=x2.dbl; volatile double r=x1v*x2v; return r; }
    friend FloatDP operator/(FloatDP x1, FloatDP x2) { volatile double x1v = x1.dbl; volatile double x2v=x2.dbl; volatile double r=x1v/x2v; return r; }
    friend FloatDP& operator+=(FloatDP& x1, FloatDP x2) { volatile double& x1v = x1.dbl; volatile double x2v=x2.dbl; x1v+=x2v; return x1; }
    friend FloatDP& operator-=(FloatDP& x1, FloatDP x2) { volatile double& x1v = x1.dbl; volatile double x2v=x2.dbl; x1v-=x2v; return x1; }
    friend FloatDP& operator*=(FloatDP& x1, FloatDP x2) { volatile double& x1v = x1.dbl; volatile double x2v=x2.dbl; x1v*=x2v; return x1; }
    friend FloatDP& operator/=(FloatDP& x1, FloatDP x2) { volatile double& x1v = x1.dbl; volatile double x2v=x2.dbl; x1v/=x2v; return x1; }

    template<class OP> friend FloatDP apply(OP op, RoundingModeType rnd, FloatDP x1, FloatDP x2, FloatDP x3) {
        auto old_rnd=FloatDP::get_rounding_mode(); FloatDP::set_rounding_mode(rnd);
        FloatDP r=op(x1,x2,x3); FloatDP::set_rounding_mode(old_rnd); return r;
    }

    template<class OP> friend FloatDP apply(OP op, RoundingModeType rnd, FloatDP x1, FloatDP x2) {
        auto old_rnd=FloatDP::get_rounding_mode(); FloatDP::set_rounding_mode(rnd);
        FloatDP r=op(x1,x2); FloatDP::set_rounding_mode(old_rnd); return r;
    }

    template<class OP> friend FloatDP apply(OP op, RoundingModeType rnd, FloatDP x) {
        auto old_rnd=FloatDP::get_rounding_mode(); FloatDP::set_rounding_mode(rnd);
        FloatDP r=op(x); FloatDP::set_rounding_mode(old_rnd); return r;
    }

    template<class OP> friend FloatDP apply(OP op, RoundingModeType rnd, FloatDP x, Int n) {
        auto old_rnd=FloatDP::get_rounding_mode(); FloatDP::set_rounding_mode(rnd);
        FloatDP r=op(x,n); FloatDP::set_rounding_mode(old_rnd); return r;
    }

    // Explicitly rounded operations
    friend FloatDP nul(RoundingModeType rnd, FloatDP x) { return apply(Nul(),rnd,x); }
    friend FloatDP pos(RoundingModeType rnd, FloatDP x) { return apply(Pos(),rnd,x); }
    friend FloatDP neg(RoundingModeType rnd, FloatDP x) { return apply(Neg(),rnd,x); }
    friend FloatDP hlf(RoundingModeType rnd, FloatDP x) { return apply(Hlf(),rnd,x); }
    friend FloatDP add(RoundingModeType rnd, FloatDP x1, FloatDP x2) { return apply(Add(),rnd,x1,x2); }
    friend FloatDP sub(RoundingModeType rnd, FloatDP x1, FloatDP x2) { return apply(Sub(),rnd,x1,x2); }
    friend FloatDP mul(RoundingModeType rnd, FloatDP x1, FloatDP x2) { return apply(Mul(),rnd,x1,x2); }
    friend FloatDP div(RoundingModeType rnd, FloatDP x1, FloatDP x2) { return apply(Div(),rnd,x1,x2); }
    friend FloatDP fma(RoundingModeType rnd, FloatDP x1, FloatDP x2, FloatDP x3) { return apply(Fma(),rnd,x1,x2,x3); }
    friend FloatDP pow(RoundingModeType rnd, FloatDP x, Int n) { return apply(Pow(),rnd,x,n); }
    friend FloatDP sqr(RoundingModeType rnd, FloatDP x) { return apply(Sqr(),rnd,x); }
    friend FloatDP rec(RoundingModeType rnd, FloatDP x) { return apply(Rec(),rnd,x); }
    friend FloatDP sqrt(RoundingModeType rnd, FloatDP x) { return apply(Sqrt(),rnd,x); }
    friend FloatDP exp(RoundingModeType rnd, FloatDP x) { return apply(Exp(),rnd,x); }
    friend FloatDP log(RoundingModeType rnd, FloatDP x) { return apply(Log(),rnd,x); }
    friend FloatDP sin(RoundingModeType rnd, FloatDP x) { return apply(Sin(),rnd,x); }
    friend FloatDP cos(RoundingModeType rnd, FloatDP x) { return apply(Cos(),rnd,x); }
    friend FloatDP tan(RoundingModeType rnd, FloatDP x) { return apply(Tan(),rnd,x); }
    friend FloatDP atan(RoundingModeType rnd, FloatDP x) { return apply(Atan(),rnd,x); }

    friend FloatDP med(RoundingModeType rnd, FloatDP x1, FloatDP x2) {
        rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(rnd);
        FloatDP r=hlf(add(rnd,x1,x2)); set_rounding_mode(rounding_mode); return r; }
    friend FloatDP rad(RoundingModeType rnd, FloatDP x1, FloatDP x2) {
        rounding_mode_t rounding_mode=get_rounding_mode(); set_rounding_mode(rnd);
        FloatDP r=hlf(sub(rnd,x2,x1)); set_rounding_mode(rounding_mode); return r; }

    friend FloatDP sqrt(RoundApprox, FloatDP x) { return std::sqrt(x.dbl); }
    friend FloatDP exp(RoundApprox, FloatDP x) { return std::exp(x.dbl); }
    friend FloatDP log(RoundApprox, FloatDP x) { return std::log(x.dbl); }
    friend FloatDP sin(RoundApprox, FloatDP x) { return std::sin(x.dbl); }
    friend FloatDP cos(RoundApprox, FloatDP x) { return std::cos(x.dbl); }
    friend FloatDP tan(RoundApprox, FloatDP x) { return std::tan(x.dbl); }
    friend FloatDP asin(RoundApprox, FloatDP x) { return std::asin(x.dbl); }
    friend FloatDP acos(RoundApprox, FloatDP x) { return std::acos(x.dbl); }
    friend FloatDP atan(RoundApprox, FloatDP x) { return std::atan(x.dbl); }

    friend Comparison cmp(FloatDP x1, FloatDP x2);
    friend Bool operator==(FloatDP x1, FloatDP x2) { return x1.dbl == x2.dbl; }
    friend Bool operator!=(FloatDP x1, FloatDP x2) { return x1.dbl != x2.dbl; }
    friend Bool operator<=(FloatDP x1, FloatDP x2) { return x1.dbl <= x2.dbl; }
    friend Bool operator>=(FloatDP x1, FloatDP x2) { return x1.dbl >= x2.dbl; }
    friend Bool operator< (FloatDP x1, FloatDP x2) { return x1.dbl <  x2.dbl; }
    friend Bool operator> (FloatDP x1, FloatDP x2) { return x1.dbl >  x2.dbl; }

    friend Comparison cmp(FloatDP x1, Dbl x2);
    friend Bool operator==(FloatDP x1, Dbl x2) { return x1.dbl == x2; }
    friend Bool operator!=(FloatDP x1, Dbl x2) { return x1.dbl != x2; }
    friend Bool operator<=(FloatDP x1, Dbl x2) { return x1.dbl <= x2; }
    friend Bool operator>=(FloatDP x1, Dbl x2) { return x1.dbl >= x2; }
    friend Bool operator< (FloatDP x1, Dbl x2) { return x1.dbl <  x2; }
    friend Bool operator> (FloatDP x1, Dbl x2) { return x1.dbl >  x2; }
    friend Comparison cmp(Dbl x1, FloatDP x2);
    friend Bool operator==(Dbl x1, FloatDP x2) { return x1 == x2.dbl; }
    friend Bool operator!=(Dbl x1, FloatDP x2) { return x1 != x2.dbl; }
    friend Bool operator<=(Dbl x1, FloatDP x2) { return x1 <= x2.dbl; }
    friend Bool operator>=(Dbl x1, FloatDP x2) { return x1 >= x2.dbl; }
    friend Bool operator< (Dbl x1, FloatDP x2) { return x1 <  x2.dbl; }
    friend Bool operator> (Dbl x1, FloatDP x2) { return x1 >  x2.dbl; }

    friend Comparison cmp(FloatDP x1, Rational const& x2);
    friend Bool operator==(FloatDP x1, Rational const& x2) { return cmp(x1,x2)==Comparison::EQUAL; }
    friend Bool operator!=(FloatDP x1, Rational const& x2) { return cmp(x1,x2)!=Comparison::EQUAL; }
    friend Bool operator<=(FloatDP x1, Rational const& x2) { return cmp(x1,x2)<=Comparison::EQUAL; }
    friend Bool operator>=(FloatDP x1, Rational const& x2) { return cmp(x1,x2)>=Comparison::EQUAL; }
    friend Bool operator< (FloatDP x1, Rational const& x2) { return cmp(x1,x2)< Comparison::EQUAL; }
    friend Bool operator> (FloatDP x1, Rational const& x2) { return cmp(x1,x2)> Comparison::EQUAL; }
    friend Comparison cmp(Rational const& x1, FloatDP x2);
    friend Bool operator==(Rational const& x1, FloatDP x2) { return cmp(x1,x2)==Comparison::EQUAL; }
    friend Bool operator!=(Rational const& x1, FloatDP x2) { return cmp(x1,x2)!=Comparison::EQUAL; }
    friend Bool operator<=(Rational const& x1, FloatDP x2) { return cmp(x1,x2)<=Comparison::EQUAL; }
    friend Bool operator>=(Rational const& x1, FloatDP x2) { return cmp(x1,x2)>=Comparison::EQUAL; }
    friend Bool operator< (Rational const& x1, FloatDP x2) { return cmp(x1,x2)< Comparison::EQUAL; }
    friend Bool operator> (Rational const& x1, FloatDP x2) { return cmp(x1,x2)> Comparison::EQUAL; }

    friend OutputStream& operator<<(OutputStream& os, FloatDP const&);
    friend InputStream& operator>>(InputStream& is, FloatDP&);
    friend OutputStream& write(OutputStream& os, FloatDP const& x, DecimalPlaces dgts, RoundingModeType rnd);
    friend OutputStream& write(OutputStream& os, FloatDP const& x, DecimalPrecision dgts, RoundingModeType rnd);
  private:
    // Rounded arithmetic
    friend FloatDP pow_rnd(FloatDP x, Int n);
    friend FloatDP sqrt_rnd(FloatDP x);
    friend FloatDP exp_rnd(FloatDP x);
    friend FloatDP log_rnd(FloatDP x);
    friend FloatDP sin_rnd(FloatDP x);
    friend FloatDP cos_rnd(FloatDP x);
    friend FloatDP tan_rnd(FloatDP x);
    friend FloatDP atan_rnd(FloatDP x);

    // Opposite rounded arithmetic
    friend FloatDP pos_opp(FloatDP x) { volatile double t=-x.dbl; return -t; }
    friend FloatDP neg_opp(FloatDP x) { volatile double t=x.dbl; return -t; }
    friend FloatDP sqr_opp(FloatDP x) { volatile double t=-x.dbl; t=t*x.dbl; return -t; }
    friend FloatDP rec_opp(FloatDP x) { volatile double t=-1.0/(volatile double&)x.dbl; return -t; }
    friend FloatDP add_opp(FloatDP x, FloatDP y) { volatile double t=-x.dbl; t=t-y.dbl; return -t; }
    friend FloatDP sub_opp(FloatDP x, FloatDP y) { volatile double t=-x.dbl; t=t+y.dbl; return -t; }
    friend FloatDP mul_opp(FloatDP x, FloatDP y) { volatile double t=-x.dbl; t=t*y.dbl; return -t; }
    friend FloatDP div_opp(FloatDP x, FloatDP y) { volatile double t=x.dbl; t=t/y.dbl; return -t; }
    friend FloatDP pow_opp(FloatDP x, int n);
};

template<class R, class A> R integer_cast(const A& a);
template<> Nat integer_cast<Nat,FloatDP>(FloatDP const&);
template<> Int integer_cast<Int,FloatDP>(FloatDP const&);

static const FloatDP inf = std::numeric_limits<double>::infinity();

struct Float32 {
    float flt;
  public:
    explicit Float32(FloatDP x, BuiltinRoundingModeType rnd) { set_rounding_mode(rnd); (volatile float&)flt = (volatile double&)x.dbl; }
    explicit operator FloatDP() const { return FloatDP((double)this->flt); }
};



} // namespace Ariadne

#endif
