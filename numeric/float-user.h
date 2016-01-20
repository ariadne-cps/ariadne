/***************************************************************************
 *            float-user.h
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

/*! \file float.h
 *  \brief Inclusion header for floating-point numbers.
 */

#ifndef ARIADNE_FLOAT_USER_H
#define ARIADNE_FLOAT_USER_H

#include "utility/macros.h"

#include "number.decl.h"
#include "float.decl.h"
#include "float64.h"
#include "floatmp.h"
#include "float-raw.h"
#include "twoexp.h"


namespace Ariadne {

template<class X> class Positive;

template<class X> struct DeclareNumericOperators {
    X nul(X const& x);
    X pos(X const& x);
    X neg(X const& x);
    X half(X const& x);
    X sqr(X const& x);
    X rec(X const& x);

    X add(X const& x1, X const& x2);
    X sub(X const& x1, X const& x2);
    X mul(X const& x1, X const& x2);
    X div(X const& x1, X const& x2);
    X fma(X const& x1, X const& x2, X const& y);
    X pow(X const& x, Nat m);
    X pow(X const& x, Int n);

    X sqrt(X const& x);
    X exp(X const& x);
    X log(X const& x);
    X sin(X const& x);
    X cos(X const& x);
    X tan(X const& x);
    X atan(X const& x);

    X max(X const& x1, X const& x2);
    X min(X const& x1, X const& x2);
    X abs(X const& x);
};


template<class PR, class P1, class P2> using FloatWeakerType = Float<Weaker<P1,P2>,PR>;

template<class PR, class P> using NegatedFloatType = Float<Negated<P>,PR>;
template<class PR, class P> using FloatNegateType = Float<Negated<P>,PR>;

template<class PR, class P1, class P2> using FloatSumType = Float<Widen<Weaker<P1,P2>>,PR>;
template<class PR, class P1, class P2> using FloatDifferenceType = Float<Widen<Weaker<P1,Negated<P2>>>,PR>;
template<class PR, class P1, class P2> using FloatProductType = Float<Widen<Weaker<P1,P2>>,PR>;
template<class PR, class P1, class P2> using FloatQuotientType = Float<Widen<Weaker<P1,Inverted<P2>>>,PR>;

template<class PR, class P1, class P2> using FloatEqualsType = Logical<Equality<Weaker<P1,Negated<P2>>>>;
template<class PR, class P1, class P2> using FloatLessType = Logical<Generic<Weaker<P1,Negated<P2>>>>;

//! \ingroup NumericModule
//! \brief Floating point numbers (double precision) using approxiamate arithmetic.
//! \details
//! The \c %Float64 class represents floating-point numbers.
//! Unless otherwise mentioned, operations on floating-point numbers are performed approximately, with no guarantees
//! on the output.
//!
//! To implement <em>interval arithmetic</em>, arithmetical operations of \c %Float64 can be performed with guaranteed rounding by
//! specifying \c _up and \c _down suffixes to arithmetical functions \c add, \c sub, \c mul and \c div.
//! Additionally, operations can be performed in the current <em>rounding mode</em> by using the \c _rnd suffix,
//! or with rounding reversed using the \c _opp suffix.
//! Operations can be specified to return an \c %ExactIntervalType answer by using the \c _ivl suffix.
//! The \c _approx suffix is provided to specifically indicate that the operation is computed approximately.
//!
//! %Ariadne floating-point numbers can be constructed by conversion from built-in C++ types.
//! Note that the value of a built-in floating-point value may differ from the mathematical value of the literal.
//! For example, while <c>%Float64(3.25)</c> is represented exactly, <c>%Float64(3.3)</c> has a value of \f$3.2999999999999998224\ldots\f$.
//! \note In the future, the construction of a \c %Float64 from a string literal may be supported.
//! \sa ExactIntervalType, Real, Float64Value
template<class PR> class Float<ApproximateTag,PR> {
    typedef ApproximateTag P; typedef RawFloat<PR> FLT;
  public:
    typedef ApproximateTag Paradigm;
    typedef Float<ApproximateTag,PR> NumericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
  public:
    Float<ApproximateTag,PR>() : _a(0.0) { }
    Float<ApproximateTag,PR>(PrecisionType pr) : _a(0.0,pr) { }
    explicit Float<ApproximateTag,PR>(RawFloatType const& a) : _a(a) { }
    template<class N, EnableIf<IsIntegral<N>> =dummy> Float<ApproximateTag,PR>(N n) : _a(n) { }
    template<class D, EnableIf<IsFloatingPoint<D>> =dummy> Float<ApproximateTag,PR>(D x) : _a(x) { }
    explicit Float<ApproximateTag,PR>(const Integer& z);
    explicit Float<ApproximateTag,PR>(const Dyadic& d);
    explicit Float<ApproximateTag,PR>(const Decimal& d);
    explicit Float<ApproximateTag,PR>(const Rational& q);
    explicit Float<ApproximateTag,PR>(const Real& r);
    explicit Float<ApproximateTag,PR>(const Number<ApproximateTag>& x);
    Float<ApproximateTag,PR>(const Rational& q, PR pr);
    Float<ApproximateTag,PR>(const Real& r, PR pr);
    Float<ApproximateTag,PR>(const Number<ApproximateTag>& x, PR pr);
    operator Number<ApproximateTag> () const;

    Float<ApproximateTag,PR>(Float<ExactTag,PR> const& x);
    Float<ApproximateTag,PR>(Float<MetricTag,PR> const& x);
    Float<ApproximateTag,PR>(Float<BoundedTag,PR> const& x);
    Float<ApproximateTag,PR>(Float<UpperTag,PR> const& x);
    Float<ApproximateTag,PR>(Float<LowerTag,PR> const& x);

    PrecisionType precision() const { return _a.precision(); }
    explicit operator RawFloatType () const { return this->_a; }
    RawFloatType const& raw() const { return this->_a; }
    RawFloatType& raw() { return this->_a; }
    double get_d() const { return this->_a.get_d(); }
  public:
    static Void set_output_precision(Nat p) { output_precision=p; }
    Float<ApproximateTag,PR> pm(Float<ApproximateTag,PR> _e) { return *this; }
  private: public:
    static Nat output_precision;
    RawFloatType _a;
};


//! \ingroup NumericModule
//! \brief Floating-point lower bounds for real numbers.
template<class PR> class Float<LowerTag,PR> {
    typedef LowerTag P; typedef RawFloat<PR> FLT;
  public:
    typedef LowerTag Paradigm;
    typedef Float<LowerTag,PR> NumericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
  public:
    Float<LowerTag,PR>() : _l(0.0) { }
    Float<LowerTag,PR>(PrecisionType pr) : _l(0.0,pr) { }
    explicit Float<LowerTag,PR>(RawFloatType const& l) : _l(l) { }

    template<class N, EnableIf<IsIntegral<N>> = dummy> Float<LowerTag,PR>(N n) : _l(n) { }
    template<class X, EnableIf<IsFloatingPoint<X>> = dummy> explicit Float<LowerTag,PR>(X x) : _l(x) { }

    Float<LowerTag,PR>(Float<BoundedTag,PR> const& x);
    Float<LowerTag,PR>(Float<MetricTag,PR> const& x);
    Float<LowerTag,PR>(Float<ExactTag,PR> const& x);

    Float<LowerTag,PR>(const Rational& q, PR pr);
    Float<LowerTag,PR>(const Real& r, PR pr);
    Float<LowerTag,PR>(const Number<LowerTag>& x, PR pr);
    operator Number<LowerTag> () const;

    explicit Float<LowerTag,PR>(const Integer& x);
    explicit Float<LowerTag,PR>(const Rational& x);
    explicit Float<LowerTag,PR>(const Real& x);
    explicit Float<LowerTag,PR>(const Number<LowerTag>& x);

    PrecisionType precision() const { return _l.precision(); }
    RawFloatType const& raw() const { return _l; }
    RawFloatType& raw() { return _l; }
    double get_d() const { return _l.get_d(); }
  private: public:
    static Nat output_precision;
    RawFloatType _l;
};


//! \ingroup NumericModule
//! \brief Floating-point upper bounds for real numbers.
template<class PR> class Float<UpperTag,PR> {
    typedef UpperTag P; typedef RawFloat<PR> FLT;
  public:
    typedef UpperTag Paradigm;
    typedef Float<UpperTag,PR> NumericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
  public:
    Float<UpperTag,PR>() : _u(0.0) { }
    Float<UpperTag,PR>(PrecisionType pr) : _u(0.0,pr) { }
    explicit Float<UpperTag,PR>(RawFloatType const& u) : _u(u) { }

    template<class N, EnableIf<IsIntegral<N>> = dummy> Float<UpperTag,PR>(N n) : _u(n) { }
    template<class X, EnableIf<IsFloatingPoint<X>> = dummy> explicit Float<UpperTag,PR>(X x) : _u(x) { }

    Float<UpperTag,PR>(Float<BoundedTag,PR> const& x);
    Float<UpperTag,PR>(Float<MetricTag,PR> const& x);
    Float<UpperTag,PR>(Float<ExactTag,PR> const& x);

    explicit Float<UpperTag,PR>(const Integer& x);
    explicit Float<UpperTag,PR>(const Rational& x);
    explicit Float<UpperTag,PR>(const Real& x);
    explicit Float<UpperTag,PR>(const Number<UpperTag>& x);

    Float<UpperTag,PR>(const Rational& q, PR pr);
    Float<UpperTag,PR>(const Real& r, PR pr);
    Float<UpperTag,PR>(const Number<UpperTag>& x, PR pr);
    operator Number<UpperTag> () const;

    PrecisionType precision() const { return _u.precision(); }
    RawFloatType const& raw() const { return _u; }
    RawFloatType& raw() { return _u; }
    double get_d() const { return _u.get_d(); }
  private: public:
    static Nat output_precision;
    RawFloatType _u;
};



//! \ingroup NumericModule
//! \brief ValidatedTag bounds on a number with floating-point endpoints supporting outwardly-rounded arithmetic.
//! \details
//! Note that direct construction from a floating-point number is prohibited, since <c>%Float64Bounds(3.3)</c> would the singleton interval \f$[3.2999999999999998224,3.2999999999999998224]\f$ (the constant is first interpreted by the C++ compiler to give a C++ \c double, whereas <c>%Float64Bounds(3.3_decimal)</c> yields the interval \f$[3.2999999999999998224,3.3000000000000002665]\f$ enclosing \f$3.3\f$.
//!
//! Comparison tests on \c FloatBounds use the idea that an interval represents a single number with an unknown value.
//! Hence the result is of type \c Kleenean, which can take values { \c True, \c False, \c Indeterminate }.
//! Hence a test \f$[\underline{x},\overline{x}]\leq [\underline{y},\overline{y}]\f$ returns \c True if \f$\overline{x}\leq \underline{y}\f$, since in this case \f$x\leq x\f$ whenever \f$x_1\in[\underline{x},\overline{x}]\f$ and \f$y\in[\underline{y},\overline{y}]\f$, \c False if \f$\underline{x}>\overline{y}\f$, since in this case we know \f$x>y\f$, and \c Indeterminate otherwise, since in this case we can find \f$x,y\f$ making the result either true or false.
//! In the case of equality, the comparison \f$[\underline{x},\overline{x}]\f$==\f$[\underline{y},\overline{y}]\f$ only returns \c True if both intervals are singletons, since otherwise we can find values making the result either true of false.
//! To test equality of representation, use \c same(x,y)
//!
//! To obtain the lower and upper bounds of the possible values, use \c x.lower() and \c x.upper().
//! To obtain a best estimate of the value, use \c x.value(), which has an error at most \a x.error().
//! If \f$v\f$ and \f$e\f$ are the returned value and error for the bounds \f$[l,u]\f$, then it is guaranteed that \f$v-e\leq l\f$ and \f$v+e\geq u\f$ in exact arithmetic.
//!
//! To test if the bounds contain a number , use \c models(FloatBounds,FloatValue), and to test if bounds are inconsistent use \c inconsistent(x,y), and to test if \c x provides a better approximation, use \c refines(x,y).
//! \sa Float64, FloatMP
//!
//! \par Python interface
//!
//! In the Python interface, %Ariadne validated bounds can be constructed from Python literals of the form \c {a:b} or (deprecated) \c [a,b] .
//! The former is preferred, as it cannot be confused with literals for other classes such as Vector and Array types.
//! Automatic conversion is used to convert FloatBounds literals of the form \c {a,b} to an FloatBounds in functions.
//!
//! Care must be taken when defining intervals using floating-point coefficients, since values are first converted to the nearest
//! representable value by the Python interpreter. <br><br>
//! \code
//!   Float64Bounds({1.1:2.3}) # Create the interval [1.1000000000000001, 2.2999999999999998]
//!   Float64Bounds({2.5:4.25}) # Create the interval [2.5, 4.25], which can be represented exactly
//!   Float64Bounds([2.5,4.25]) # Alternative syntax for creating the interval [2.5, 4.25]
//! \endcode
template<class PR> class Float<BoundedTag,PR> {
    typedef BoundedTag P; typedef RawFloat<PR> FLT;
  public:
    typedef BoundedTag Paradigm;
    typedef Float<BoundedTag,PR> NumericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
  public:
    Float<BoundedTag,PR>() : _l(0.0), _u(0.0) { }
    Float<BoundedTag,PR>(PrecisionType pr) : _l(0.0,pr), _u(0.0,pr) { }
    explicit Float<BoundedTag,PR>(RawFloatType const& v) : _l(v), _u(v) { }
    Float<BoundedTag,PR>(RawFloatType const& l, RawFloatType const& u) : _l(l), _u(u) { }
    Float<BoundedTag,PR>(Float<LowerTag,PR> const& lower, Float<UpperTag,PR> const& upper) : _l(lower.raw()), _u(upper.raw()) { }
    template<class N1, class N2, EnableIf<And<IsIntegral<N1>,IsIntegral<N2>>> = dummy> Float<BoundedTag,PR>(N1 n1, N2 n2) : _l(n1), _u(n2) { }
    Float<BoundedTag,PR>(Rational const& ql, Rational const& qu, PrecisionType pr);

    template<class N, EnableIf<IsIntegral<N>> = dummy> Float<BoundedTag,PR>(N n) : _l(n), _u(n) { }
    template<class X, EnableIf<IsFloatingPoint<X>> = dummy> explicit Float<BoundedTag,PR>(X x) : _l(x), _u(x) { }

    Float<BoundedTag,PR>(Float<MetricTag,PR> const& x);
    Float<BoundedTag,PR>(Float<ExactTag,PR> const& x);

    explicit Float<BoundedTag,PR>(const Dyadic& x);
    explicit Float<BoundedTag,PR>(const Decimal& x);
    explicit Float<BoundedTag,PR>(const Integer& z);
    explicit Float<BoundedTag,PR>(const Rational& q);
    explicit Float<BoundedTag,PR>(const Real& x);
    explicit Float<BoundedTag,PR>(const Number<ValidatedTag>& x);
    Float<BoundedTag,PR>(const Integer& z, PR pr);
    Float<BoundedTag,PR>(const Rational& q, PR pr);
    Float<BoundedTag,PR>(const Real& x, PR pr);
    Float<BoundedTag,PR>(const Number<ValidatedTag>& x, PR pr);
    operator Number<ValidatedTag> () const;

    Float<LowerTag,PR> const lower() const { return Float<LowerTag,PR>(lower_raw()); }
    Float<UpperTag,PR> const upper() const { return Float<UpperTag,PR>(upper_raw()); }
    Float<ExactTag,PR> const value() const;
    Float<ErrorTag,PR> const error() const;

    RawFloatType const& lower_raw() const { return _l; }
    RawFloatType const& upper_raw() const { return _u; }
    RawFloatType const value_raw() const { return half(add_near(_l,_u)); }
    RawFloatType const error_raw() const { RawFloatType v=value_raw(); return max(sub_up(_u,v),sub_up(v,_l)); }
    double get_d() const { return value_raw().get_d(); }

    PrecisionType precision() const { ARIADNE_DEBUG_ASSERT(_l.precision()==_u.precision()); return _u.precision(); }

    // DEPRECATED
    explicit operator RawFloatType () const { return value_raw(); }
    friend Float<ExactTag,PR> midpoint(Float<BoundedTag,PR> const& x);
  public:
    static Nat output_precision;
    static Void set_output_precision(Nat p) { output_precision=p; }
  private: public:
    RawFloatType _l, _u;
};


template<class PR> class Float<MetricTag,PR> {
    typedef MetricTag P; typedef RawFloat<PR> FLT;
  public:
    typedef MetricTag Paradigm;
    typedef Float<MetricTag,PR> NumericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
  public:
    Float<MetricTag,PR>() : _v(0.0), _e(0.0) { }
    Float<MetricTag,PR>(PrecisionType pr) : _v(0.0,pr), _e(0.0,pr) { }
    explicit Float<MetricTag,PR>(RawFloatType const& v) : _v(v), _e(0.0) { }
    Float<MetricTag,PR>(RawFloatType const& v, RawFloatType const& e) : _v(v), _e(e) { }
    Float<MetricTag,PR>(Float<ExactTag,PR> const& value, Float<ErrorTag,PR> const& error) : _v(value.raw()), _e(error.raw()) { }
    Float<MetricTag,PR>(Float<LowerTag,PR> const& lower, Float<UpperTag,PR> const& upper) =  delete;

    Float<MetricTag,PR>(Float<BoundedTag,PR> const& x);
    Float<MetricTag,PR>(Float<ExactTag,PR> const& x);

    template<class N, EnableIf<IsIntegral<N>> = dummy> Float<MetricTag,PR>(N n) : _v(n), _e(nul(_v)) { }
    template<class X, EnableIf<IsFloatingPoint<X>> = dummy> explicit Float<MetricTag,PR>(X x) : _v(x), _e(nul(_v)) { }
    explicit Float<MetricTag,PR>(const Integer& z);
    explicit Float<MetricTag,PR>(const Rational& q);
    explicit Float<MetricTag,PR>(const Real& x);
    explicit Float<MetricTag,PR>(const Number<ValidatedTag>& x);
    Float<MetricTag,PR>(const Rational& q, PR pr);
    Float<MetricTag,PR>(const Real& r, PR pr);
    Float<MetricTag,PR>(const Number<ValidatedTag>& x, PR pr);
    operator Number<ValidatedTag> () const;

    Float<LowerTag,PR> const lower() const { return Float<LowerTag,PR>(lower_raw()); }
    Float<UpperTag,PR> const upper() const { return Float<UpperTag,PR>(upper_raw()); }
    Float<ExactTag,PR> const value() const;
    Float<ErrorTag,PR> const error() const;

    RawFloatType const lower_raw() const { return sub_down(_v,_e); }
    RawFloatType const upper_raw() const { return add_up(_v,_e); }
    RawFloatType const& value_raw() const { return _v; }
    RawFloatType const& error_raw() const { return _e; }
    double get_d() const { return _v.get_d(); }

    PrecisionType precision() const { return _v.precision(); }
  private: public:
    RawFloatType _v, _e;
};

//! \ingroup NumericModule
//! \related Float64, Float64Bounds
//! \brief A floating-point number, which is taken to represent the \em exact value of a real quantity.
template<class PR> class Float<ExactTag,PR> {
    typedef ExactTag P; typedef RawFloat<PR> FLT;
  public:
    typedef ExactTag Paradigm;
    typedef Float<ExactTag,PR> NumericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
  public:
    Float<ExactTag,PR>() : _v(0.0) { }
    Float<ExactTag,PR>(PrecisionType pr) : _v(0.0,pr) { }
    explicit Float<ExactTag,PR>(RawFloatType const& v) : _v(v) { }
    template<class N, EnableIf<IsIntegral<N>> =dummy> Float<ExactTag,PR>(N n) : _v(n) { }
    template<class X, EnableIf<IsFloatingPoint<X>> =dummy> explicit Float<ExactTag,PR>(X x) : _v(x) { }


    explicit Float<ExactTag,PR>(const Integer& z);
    explicit Float<ExactTag,PR>(const Integer& z, PR pr);
    explicit Float<ExactTag,PR>(const TwoExp& ex, PR pr);
    explicit operator Rational () const;
    operator Number<ExactTag> () const;
    explicit operator RawFloatType () const { return _v; }

    PrecisionType precision() const { return _v.precision(); }
    RawFloatType const& raw() const { return _v; }
    RawFloatType& raw() { return _v; }
    double get_d() const { return _v.get_d(); }

    Float<MetricTag,PR> pm(Float<ErrorTag,PR> _e) const;
  public:
    static Nat output_precision;
    static Void set_output_precision(Nat p) { output_precision=p; }
  private: public:
    RawFloatType _v;
};

template<class PR> inline const Float<ExactTag,PR> Float<BoundedTag,PR>::value() const {
    return Float<ExactTag,PR>(med_near(this->_l,this->_u)); }

template<class PR> inline const Float<ErrorTag,PR> Float<BoundedTag,PR>::error() const {
    RawFloat<PR> _v=med_near(this->_l,this->_u); return Float<ErrorTag,PR>(max(sub_up(this->_u,_v),sub_up(_v,this->_l))); }

template<class PR> inline Float<ExactTag,PR> value(Float<BoundedTag,PR> const& x) {
    return x.value(); }

template<class PR> inline Float<ErrorTag,PR> error(Float<BoundedTag,PR> const& x) {
    return x.error(); }

template<class PR> inline const Float<ExactTag,PR> Float<MetricTag,PR>::value() const {
    return Float<ExactTag,PR>(this->_v); }

template<class PR> inline const Float<ErrorTag,PR> Float<MetricTag,PR>::error() const {
    return Float<ErrorTag,PR>(this->_e); }


template<class PR> class Float<PositiveExactTag,PR> : public Float<ExactTag,PR> {
  public:
    Float<PositiveExactTag,PR>() : Float<ExactTag,PR>() { }
    Float<PositiveExactTag,PR>(TwoExp ex) : Float<ExactTag,PR>(ex) { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy>
        Float<PositiveExactTag,PR>(M m) : Float<ExactTag,PR>(m) { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy>
        Float<PositiveExactTag,PR>(M m, PR pr) : Float<ExactTag,PR>(m,pr) { }
    explicit Float<PositiveExactTag,PR>(RawFloat<PR> const& x) : Float<ExactTag,PR>(x) { }
    explicit Float<PositiveExactTag,PR>(Float<ExactTag,PR> const& x) : Float<ExactTag,PR>(x) { }
};

template<class PR> class Float<PositiveBoundedTag,PR> : public Float<BoundedTag,PR> {
  public:
    Float<PositiveBoundedTag,PR>() : Float<BoundedTag,PR>() { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy>
        Float<PositiveBoundedTag,PR>(M m) : Float<BoundedTag,PR>(m) { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy>
        Float<PositiveBoundedTag,PR>(M m, PR pr) : Float<BoundedTag,PR>(m,pr) { }
    explicit Float<PositiveBoundedTag,PR>(RawFloat<PR> const& x) : Float<BoundedTag,PR>(x) { }
    explicit Float<PositiveBoundedTag,PR>(RawFloat<PR> const& l, RawFloat<PR> const& u) : Float<BoundedTag,PR>(l,u) { }
    explicit Float<PositiveBoundedTag,PR>(Float<BoundedTag,PR> const& x) : Float<BoundedTag,PR>(x) { }
};

template<class PR> class Float<PositiveUpperTag,PR> : public Float<UpperTag,PR> {
  public:
    Float<PositiveUpperTag,PR>() : Float<UpperTag,PR>() { }
    explicit Float<PositiveUpperTag,PR>(RawFloat<PR> const& x) : Float<UpperTag,PR>(x) {
        ARIADNE_PRECONDITION_MSG(!(x<0),"x="<<x); }
    explicit Float<PositiveUpperTag,PR>(Float<UpperTag,PR> const& x) : Float<UpperTag,PR>(x) { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy> Float<PositiveUpperTag,PR>(M m) : Float<UpperTag,PR>(m) { }
    template<class F, EnableIf<IsSame<F,Float<UpperTag,PR>>> =dummy>
        explicit Float<PositiveUpperTag,PR>(F const& x) : Float<UpperTag,PR>(x) { }
    Float<PositiveUpperTag,PR>(Float<PositiveExactTag,PR> const& x) : Float<UpperTag,PR>(x) { }
};
template<class PR> Float<PositiveUpperTag,PR> abs(Float<PositiveUpperTag,PR> const&);

template<class PR> class Float<PositiveLowerTag,PR> : public Float<LowerTag,PR> {
  public:
    Float<PositiveLowerTag,PR>() : Float<LowerTag,PR>() { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy>
        Float<PositiveLowerTag,PR>(M m) : Float<LowerTag,PR>(m) { }
    explicit Float<PositiveLowerTag,PR>(RawFloat<PR> const& x) : Float<LowerTag,PR>(x) { }
    explicit Float<PositiveLowerTag,PR>(Float<LowerTag,PR> const& x) : Float<LowerTag,PR>(x) { }
    Float<PositiveLowerTag,PR>(Float<PositiveExactTag,PR> const& x) : Float<LowerTag,PR>(x) { }
};

template<class PR> class Float<PositiveApproximateTag,PR> : public Float<ApproximateTag,PR> {
  public:
    Float<PositiveApproximateTag,PR>() : Float<ApproximateTag,PR>() { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy>
        Float<PositiveApproximateTag,PR>(M m) : Float<ApproximateTag,PR>(m) { }
    explicit Float<PositiveApproximateTag,PR>(RawFloat<PR> const& x) : Float<ApproximateTag,PR>(x) { }
    explicit Float<PositiveApproximateTag,PR>(Float<ApproximateTag,PR> const& x) : Float<ApproximateTag,PR>(x) { }
    Float<PositiveApproximateTag,PR>(Float<PositiveLowerTag,PR> const& x) : Float<ApproximateTag,PR>(x) { }
    Float<PositiveApproximateTag,PR>(Float<PositiveUpperTag,PR> const& x) : Float<ApproximateTag,PR>(x) { }
    Float<PositiveApproximateTag,PR>(Float<PositiveExactTag,PR> const& x) : Float<ApproximateTag,PR>(x) { }
};

template<class PR> inline Float<PositiveApproximateTag,PR> cast_positive(Float<ApproximateTag,PR> const& x) {
    return Float<PositiveApproximateTag,PR>(x); }

template<class PR> inline Float<PositiveLowerTag,PR> cast_positive(Float<LowerTag,PR> const& x) {
    return Float<PositiveLowerTag,PR>(x); }

template<class PR> inline Float<PositiveUpperTag,PR> cast_positive(Float<UpperTag,PR> const& x) {
    return Float<PositiveUpperTag,PR>(x); }

template<class PR> inline Float<PositiveExactTag,PR> cast_positive(Float<ExactTag,PR> const& x) {
    return Float<PositiveExactTag,PR>(x); }

template<class PR> inline OutputStream& operator<<(OutputStream& os, Float<PositiveApproximateTag,PR> const& x) {
    return os << static_cast<Float<ApproximateTag,PR>const&>(x); }

template<class R, class A> R integer_cast(const A& _a);

template<> Float<MetricTag,Precision64>::Float(Real const& x);
template<> Float<BoundedTag,Precision64>::Float(Real const& x);
template<> Float<UpperTag,Precision64>::Float(Real const& x);
template<> Float<LowerTag,Precision64>::Float(Real const& x);
template<> Float<ApproximateTag,Precision64>::Float(Real const& x);

template<class T, class F, EnableIf<Not<IsSame<T,F>>> =dummy> T convert(F const& x) { return T(x); }
template<class T> T const& convert(T const& x) { return x; }

template<class T> using NumericType = typename T::NumericType;

typedef Widen<ExactTag> Widened;

template<class PR, class P> auto
max(Float<P,PR> const& x1, Float<P,PR> const& x2) -> Float<P,PR>;
template<class PR, class P> auto
min(Float<P,PR> const& x1, Float<P,PR> const& x2) -> Float<P,PR>;
template<class PR, class P> auto
abs(Float<P,PR> const& x) -> Float<Weaker<P,Negated<P>>,PR>;

template<class PR, class P> auto
floor(Float<P,PR> const& x) -> Float<P,PR>;


template<class PR, class P> auto
nul(Float<P,PR> const&) -> Float<P,PR>;
template<class PR, class P> auto
pos(Float<P,PR> const&) -> Float<P,PR>;
template<class PR, class P> auto
neg(Float<P,PR> const&) -> Float<Negated<P>,PR>;
template<class PR, class P> auto
half(Float<P,PR> const&) -> Float<P,PR>;
template<class PR, class P> auto
sqr(Float<P,PR> const&) -> Float<Widen<P>,PR>;
template<class PR, class P> auto
rec(Float<P,PR> const&) -> Float<Widen<Inverted<P>>,PR>;

template<class PR, class P> auto
add(Float<P,PR> const&, Float<P,PR> const&) -> Float<Widen<P>,PR>;
template<class PR, class P> auto
sub(Float<P,PR> const&, Float<Negated<P>,PR> const&) -> Float<Widen<P>,PR>;
template<class PR, class P> auto
mul(Float<P,PR> const&, Float<P,PR> const&) -> Float<Widen<P>,PR>;
template<class PR, class P> auto
div(Float<P,PR> const&, Float<Inverted<P>,PR> const&) -> Float<Widen<P>,PR>;

template<class PR, class P> auto
pow(Float<P,PR> const&, Nat m) -> Float<Widen<P>,PR>;
template<class PR, class P> auto
pow(Float<P,PR> const&, Int n) -> Float<Widen<Undirect<P>>,PR>;

template<class PR, class P> auto
sqrt(Float<P,PR> const&) -> Float<Widen<P>,PR>;
template<class PR, class P> auto
exp(Float<P,PR> const&) -> Float<Widen<P>,PR>;
template<class PR, class P> auto
log(Float<P,PR> const&) -> Float<Widen<Signed<P>>,PR>;
template<class PR, class P> auto
sin(Float<P,PR> const&) -> Float<Widen<Unorder<P>>,PR>;
template<class PR, class P> auto
cos(Float<P,PR> const&) -> Float<Widen<Unorder<P>>,PR>;
template<class PR, class P> auto
tan(Float<P,PR> const&) -> Float<Widen<Unorder<P>>,PR>;
template<class PR, class P> auto
asin(Float<P,PR> const&) -> Float<Widen<Unorder<P>>,PR>;
template<class PR, class P> auto
acos(Float<P,PR> const&) -> Float<Widen<Unorder<P>>,PR>;
template<class PR, class P> auto
atan(Float<P,PR> const&) -> Float<Widen<P>,PR>;

template<class PR, class P> auto
floor(Float<P,PR> const& x) -> Float<P,PR>;

template<class PR, class P> auto
round(Float<P,PR> const& x) -> Float<P,PR>;

template<class PR, class P> auto mag(Float<P,PR> const& x) -> Float<Unsigned<Weaker<P,UpperTag>>,PR>;
template<class PR, class P> auto mig(Float<P,PR> const& x) -> Float<Unsigned<Weaker<P,LowerTag>>,PR>;

template<class PR, class P> auto is_zero(Float<P,PR> const&) -> Logical<Weaker<P,Opposite<P>>>;
template<class PR, class P> auto is_positive(Float<P,PR> const&) -> Logical<Opposite<P>>;
template<class PR, class P> auto same(Float<P,PR> const&, Float<P,PR> const&) -> Bool;

template<class PR, class P> auto eq(Float<P,PR> const& x1, Float<Negated<P>,PR> const& x2) -> FloatEqualsType<PR,P,Negated<P>>;
template<class PR, class P> auto leq(Float<P,PR> const& x1, Float<Negated<P>,PR> const& x2) -> FloatLessType<PR,P,Negated<P>>;


template<class PR, class P1, class P2> auto
max(Float<P1,PR> const& x1, Float<P2,PR> const& x2) -> Float<Weaker<P1,P2>,PR>;

template<class PR, class P1, class P2> auto
max(Float<P1,PR> const& x1, Float<P2,PR> const& x2) -> Float<Weaker<P1,P2>,PR> {
    typedef Float<Weaker<P1,P2>,PR> R; return max(convert<R>(x1),convert<R>(x2)); }

template<class PR, class P1, class P2> auto
min(Float<P1,PR> const& x1, Float<P2,PR> const& x2) -> Float<Weaker<P1,P2>,PR> {
    typedef Weaker<P1,P2> P0; typedef Float<P0,PR> X0;
    Float<P0,PR> rx1(x1); Float<P0,PR> rx2(x2); return min(rx1,rx2); }


template<class PR, class P> auto
operator+(Float<P,PR> const& x) -> Float<P,PR> {
    return pos(x); }

template<class PR, class P> auto
operator-(Float<P,PR> const& x) -> Float<Negated<P>,PR> {
    return neg(x); }

template<class PR, class P1, class P2> auto
operator+(Float<P1,PR> const& x1, Float<P2,PR> const& x2) -> FloatSumType<PR,P1,P2> {
    typedef Widen<Weaker<P1,P2>> P0; typedef Float<P0,PR> R;
    return add(convert<Float<P0,PR>>(x1),convert<Float<P0,PR>>(x2)); }

template<class PR, class P1, class P2> auto
operator-(Float<P1,PR> const& x1, Float<P2,PR> const& x2) -> Float<Widen<Weaker<P1,Negated<P2>>>,PR> {
    typedef Widen<Weaker<P1,Negated<P2>>> P0; typedef Float<P0,PR> R;
    return sub(convert<R>(x1),convert<Float<Negated<P0>,PR>>(x2)); }

template<class PR, class P1, class P2> auto
operator*(Float<P1,PR> const& x1, Float<P2,PR> const& x2) -> Float<Widen<Weaker<P1,P2>>,PR> {
    typedef Float<Widen<Weaker<P1,P2>>,PR> R;
    return mul(convert<R>(x1),convert<R>(x2)); }

template<class PR, class P1, class P2> auto
operator/(Float<P1,PR> const& x1, Float<P2,PR> const& x2) -> Float<Widen<Weaker<P1,Inverted<P2>>>,PR> {
    typedef Widen<Weaker<P1,Inverted<P2>>> P0; typedef Float<P0,PR> R;
    return div(convert<R>(x1),convert<Float<Inverted<P0>,PR>>(x2)); }


template<class PR, class P> auto
operator*(Float<P,PR> const& x1, TwoExp x2) -> Float<P,PR>;

template<class PR, class P> auto
operator/(Float<P,PR> const& x1, TwoExp x2) -> Float<P,PR>;

template<class P1, class Y2, class PR> auto
operator+=(Float<P1,PR>& x1, Y2 const& y2) -> decltype(x1=x1+y2) {
    return x1=x1+y2; }

template<class P1, class Y2, class PR> auto
operator-=(Float<P1,PR>& x1, Y2 const& y2) -> decltype(x1=x1-y2) {
    return x1=x1-y2; }

template<class P1, class Y2, class PR> auto
operator*=(Float<P1,PR>& x1, Y2 const& y2) -> decltype(x1=x1*y2) {
    return x1=x1*y2; }

template<class P1, class Y2, class PR> auto
operator/=(Float<P1,PR>& x1, Y2 const& y2) -> decltype(x1=x1/y2) {
    return x1=x1/y2; }

template<class PR, class P1, class P2> auto
operator==(Float<P1,PR> const& x1, Float<P2,PR> const& x2) -> FloatEqualsType<PR,P1,P2> {
    typedef Weaker<P1,Opposite<P2>> WP1; typedef Opposite<WP1> WP2; typedef Float<WP1,PR> X1; typedef Float<WP2,PR> X2;
    return eq(convert<X1>(x1),convert<X2>(x2)); }

template<class PR, class P1, class P2> auto
operator!=(Float<P1,PR> const& x1, Float<P2,PR> const& x2) -> decltype(not (x1==x2)) {
    return !(x1==x2); }

template<class PR, class P1, class P2> auto
operator<=(Float<P1,PR> const& x1, Float<P2,PR> const& x2) ->  FloatLessType<PR,P1,P2> {
    typedef Weaker<P1,Negated<P2>> WP1; typedef Negated<WP1> WP2; typedef Float<WP1,PR> X1; typedef Float<WP2,PR> X2;
    return leq(convert<X1>(x1),convert<X2>(x2)); }

template<class PR, class P1, class P2> auto
operator>=(Float<P1,PR> const& x1, Float<P2,PR> const& x2) -> decltype(x2<=x1) {
    return x2<=x1; }

template<class PR, class P1, class P2> auto
operator< (Float<P1,PR> const& x1, Float<P2,PR> const& x2) -> decltype(not (x2<=x1)) {
    return not (x2<=x1); }

template<class PR, class P1, class P2> auto
operator> (Float<P1,PR> const& x1, Float<P2,PR> const& x2) -> decltype(not (x1<=x2)){
    return not (x1<=x2); }

template<class PR, class P> auto
operator<<(OutputStream& os, Float<P,PR> const&) -> OutputStream&;

template<class PR, class P> auto
operator>>(InputStream& is, Float<P,PR>&) -> InputStream&;


extern const Float<ExactTag,Precision64> infty;

template<class P1, class P2, class PR> Float<Weaker<P1,P2>,PR> w(Float<P1,PR> x1, Float<P2,PR> x2) {
    return Float<Weaker<P1,P2>,PR>(x1); }


// Literals operations
Float64Value operator"" _exact(long double lx);
Float64Error operator"" _error(long double lx);
Float64Ball operator"" _near(long double lx);
Float64UpperBound operator"" _upper(long double lx);
Float64LowerBound operator"" _lower(long double lx);
Float64Approximation operator"" _approx(long double lx);


// ValidatedTag operations
template<class PR> Float<BoundedTag,PR> make_bounds(Float<ErrorTag,PR> const& e) {
    return Float<BoundedTag,PR>(-e.raw(),+e.raw()); }

//! \related Float, ValidatedTag \brief Tests if \_a x1 provides tighter bounds than \_a x2.
template<class PR> Bool refines(Float<MetricTag,PR> const& x1, Float<MetricTag,PR> const& x2);
template<class PR> Bool refines(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2);
template<class PR> Bool refines(Float<LowerTag,PR> const& x1, Float<LowerTag,PR> const& x2);
template<class PR> Bool refines(Float<UpperTag,PR> const& x1, Float<UpperTag,PR> const& x2);

//! \related Float, ValidatedTag \brief The common refinement of \_a x1 and \_a x2.
template<class PR> Float<BoundedTag,PR> refinement(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2);
template<class PR> Float<MetricTag,PR> refinement(Float<MetricTag,PR> const& x1, Float<MetricTag,PR> const& x2);

//! \related Float, ValidatedTag \brief Tests if \_a x1 and \_a x2 are consistent with representing the same number.
template<class PR> Bool consistent(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2);

//! \related Float, ValidatedTag \brief  Tests if \_a x1 and \_a x2 are inconsistent with representing the same number.
template<class PR> Bool inconsistent(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2);

//! \related Float, ValidatedTag \brief  Tests if \_a x1 is a model for the exact value \_a x2. number.
template<class PR> Bool models(Float<BoundedTag,PR> const& x1, Float<ExactTag,PR> const& x2);


template<class PR> inline Float<ApproximateTag,PR> make_float(Number<ApproximateTag> const& y, PR pr) { return Float<ApproximateTag,PR>(y,pr); }
template<class PR> Float<ApproximateTag,PR> make_float(Number<ApproximateTag> const& y, PR pr);
template<class PR> inline Float<LowerTag,PR> make_float(Number<LowerTag> const& y, PR pr) { return Float<LowerTag,PR>(y,pr); }
template<class PR> inline Float<UpperTag,PR> make_float(Number<UpperTag> const& y, PR pr) { return Float<UpperTag,PR>(y,pr); }
template<class PR> inline Float<BoundedTag,PR> make_float(Number<ValidatedTag> const& y, PR pr) { return Float<BoundedTag,PR>(y,pr); }
template<class PR> inline Float<BoundedTag,PR> make_float(Number<EffectiveTag> const& y, PR pr) { return Float<BoundedTag,PR>(y,pr); }
template<class PR> inline Float<BoundedTag,PR> make_float(Number<ExactTag> const& y, PR pr) { return Float<BoundedTag,PR>(y,pr); }
template<class PR> inline Float<BoundedTag,PR> make_float(Real const& y, PR pr) { return Float<BoundedTag,PR>(y,pr); }
template<class PR> inline Float<BoundedTag,PR> make_float(Rational const& y, PR pr) { return Float<BoundedTag,PR>(y,pr); }
template<class PR> inline Float<ExactTag,PR> make_float(Integer const& y, PR pr) { return Float<ExactTag,PR>(y,pr); }
template<class N, class PR, EnableIf<IsSignedIntegral<N>> =dummy> inline Float<ExactTag,PR> make_float(N const& y, PR pr) { return Float<ExactTag,PR>(y,pr); }
template<class M, class PR, EnableIf<IsUnsignedIntegral<M>> =dummy> inline Float<PositiveExactTag,PR> make_float(M const& y, PR pr) { return Float<PositiveExactTag,PR>(y,pr); }
template<class D, class PR, EnableIf<IsFloatingPoint<D>> =dummy> Float<ApproximateTag,PR> make_float(D const& y, PR pr){
    return Float<ApproximateTag,PR>(RawFloat<PR>(y,pr)); }

template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy> auto
operator+(X const& x, Y const& y) -> decltype(x+make_float(y,x.precision())) { return x+make_float(y,x.precision()); }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy> auto
operator-(X const& x, Y const& y) -> decltype(x-make_float(y,x.precision())) { return x-make_float(y,x.precision()); }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy> auto
operator*(X const& x, Y const& y) -> decltype(x*make_float(y,x.precision())) { return x*make_float(y,x.precision()); }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy> auto
operator/(X const& x, Y const& y) -> decltype(x/make_float(y,x.precision())) { return x/make_float(y,x.precision()); }

template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy> auto
operator+(Y const& y, X const& x) -> decltype(make_float(y,x.precision())+x) { return make_float(y,x.precision())+x; }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy> auto
operator-(Y const& y, X const& x) -> decltype(make_float(y,x.precision())-x) { return make_float(y,x.precision())-x; }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy> auto
operator*(Y const& y, X const& x) -> decltype(make_float(y,x.precision())*x) { return make_float(y,x.precision())*x; }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy> auto
operator/(Y const& y, X const& x) -> decltype(make_float(y,x.precision())/x) { return make_float(y,x.precision())/x; }

template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy> auto
operator==(X const& x, Y const& y) -> decltype(x==make_float(y,x.precision())) { return x==make_float(y,x.precision()); }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy> auto
operator!=(X const& x, Y const& y) -> decltype(x!=make_float(y,x.precision())) { return x!=make_float(y,x.precision()); }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy> auto
operator<=(X const& x, Y const& y) -> decltype(x<=make_float(y,x.precision())) { return x<=make_float(y,x.precision()); }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy> auto
operator>=(X const& x, Y const& y) -> decltype(x>=make_float(y,x.precision())) { return x>=make_float(y,x.precision()); }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy> auto
operator< (X const& x, Y const& y) -> decltype(x< make_float(y,x.precision())) { return x< make_float(y,x.precision()); }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy> auto
operator> (X const& x, Y const& y) -> decltype(x> make_float(y,x.precision())) { return x> make_float(y,x.precision()); }

template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy> auto
operator==(Y const& y, X const& x) -> decltype(make_float(y,x.precision())==x) { return make_float(y,x.precision())==x; }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy> auto
operator!=(Y const& y, X const& x) -> decltype(make_float(y,x.precision())!=x) { return make_float(y,x.precision())!=x; }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy> auto
operator<=(Y const& y, X const& x) -> decltype(make_float(y,x.precision())<=x) { return make_float(y,x.precision())<=x; }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy> auto
operator>=(Y const& y, X const& x) -> decltype(make_float(y,x.precision())>=x) { return make_float(y,x.precision())>=x; }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy> auto
operator< (Y const& y, X const& x) -> decltype(make_float(y,x.precision())< x) { return make_float(y,x.precision())< x; }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy> auto
operator> (Y const& y, X const& x) -> decltype(make_float(y,x.precision())> x) { return make_float(y,x.precision())> x; }

template<class PR> auto operator==(Float<ExactTag,PR> const& x, Rational const& q) -> decltype(Rational(x)==q) { return Rational(x)==q; }
template<class PR> auto operator!=(Float<ExactTag,PR> const& x, Rational const& q) -> decltype(Rational(x)!=q) { return Rational(x)!=q; }
template<class PR> auto operator<=(Float<ExactTag,PR> const& x, Rational const& q) -> decltype(Rational(x)<=q) { return Rational(x)<=q; }
template<class PR> auto operator>=(Float<ExactTag,PR> const& x, Rational const& q) -> decltype(Rational(x)>=q) { return Rational(x)>=q; }
template<class PR> auto operator< (Float<ExactTag,PR> const& x, Rational const& q) -> decltype(Rational(x)< q) { return Rational(x)< q; }
template<class PR> auto operator> (Float<ExactTag,PR> const& x, Rational const& q) -> decltype(Rational(x)> q) { return Rational(x)> q; }

template<class PR> auto operator==(Rational const& q, Float<ExactTag,PR> const& x) -> decltype(q==Rational(x)) { return q==Rational(x); }
template<class PR> auto operator!=(Rational const& q, Float<ExactTag,PR> const& x) -> decltype(q!=Rational(x)) { return q!=Rational(x); }
template<class PR> auto operator<=(Rational const& q, Float<ExactTag,PR> const& x) -> decltype(q<=Rational(x)) { return q<=Rational(x); }
template<class PR> auto operator>=(Rational const& q, Float<ExactTag,PR> const& x) -> decltype(q>=Rational(x)) { return q>=Rational(x); }
template<class PR> auto operator< (Rational const& q, Float<ExactTag,PR> const& x) -> decltype(q< Rational(x)) { return q< Rational(x); }
template<class PR> auto operator> (Rational const& q, Float<ExactTag,PR> const& x) -> decltype(q> Rational(x)) { return q> Rational(x); }

Float64Value cast_exact(const Real& x);

inline Float64Value const& cast_exact(RawFloat64 const& x) { return reinterpret_cast<Float64Value const&>(x); }
inline Float64Value const& cast_exact(Float64Approximation const& x) { return reinterpret_cast<Float64Value const&>(x); }
inline Float64Value const& cast_exact(Float64Value const& x) { return reinterpret_cast<Float64Value const&>(x); }

template<template<class>class T> inline const T<Float64Value>& cast_exact(const T<RawFloat64>& t) {
    return reinterpret_cast<const T<Float64Value>&>(t); }
template<template<class>class T> inline const T<Float64Value>& cast_exact(const T<Float64Approximation>& t) {
    return reinterpret_cast<const T<Float64Value>&>(t); }
template<template<class>class T> inline const T<Float64Value>& cast_exact(const T<Float64Value>& t) {
    return reinterpret_cast<const T<Float64Value>&>(t); }

inline RawFloat64 const& cast_raw(RawFloat64 const& x) { return reinterpret_cast<RawFloat64 const&>(x); }
inline RawFloat64 const& cast_raw(Float64Approximation const& x) { return reinterpret_cast<RawFloat64 const&>(x); }
inline RawFloat64 const& cast_raw(Float64Value const& x) { return reinterpret_cast<RawFloat64 const&>(x); }

template<template<class>class T> inline const T<RawFloat64>& cast_raw(const T<RawFloat64>& t) {
    return reinterpret_cast<const T<RawFloat64>&>(t); }
template<template<class>class T> inline const T<RawFloat64>& cast_raw(const T<Float64Approximation>& t) {
    return reinterpret_cast<const T<RawFloat64>&>(t); }
template<template<class>class T> inline const T<RawFloat64>& cast_raw(const T<Float64Value>& t) {
    return reinterpret_cast<const T<RawFloat64>&>(t); }

inline Float64Approximation const& cast_approximate(RawFloat64 const& x) { return reinterpret_cast<Float64Approximation const&>(x); }
inline Float64Approximation const& cast_approximate(Float64Approximation const& x) { return reinterpret_cast<Float64Approximation const&>(x); }
inline Float64Approximation const& cast_approximate(Float64Value const& x) { return reinterpret_cast<Float64Approximation const&>(x); }

template<template<class>class T> inline const T<Float64Approximation>& cast_approximate(const T<RawFloat64>& t) {
    return reinterpret_cast<const T<Float64Approximation>&>(t); }
template<template<class>class T> inline const T<Float64Approximation>& cast_approximate(const T<Float64Approximation>& t) {
    return reinterpret_cast<const T<Float64Approximation>&>(t); }
template<template<class>class T> inline const T<Float64Approximation>& cast_approximate(const T<Float64Value>& t) {
    return reinterpret_cast<const T<Float64Approximation>&>(t); }

} // namespace Ariadne

#endif
