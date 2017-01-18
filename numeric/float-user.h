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
#include "builtin.h"
#include "twoexp.h"


namespace Ariadne {


template<class X> class Positive;

template<class X> struct NumericTraits;

template<class PR> struct NumericTraits<FloatApproximation<PR>> {
    typedef ApproximateNumber GenericType;
    typedef PositiveFloatApproximation<PR> PositiveType;
    typedef Fuzzy LessType;
    typedef Fuzzy EqualsType;
};
template<class PR> struct NumericTraits<FloatLowerBound<PR>> {
    typedef ValidatedLowerNumber GenericType;
    typedef FloatUpperBound<PR> OppositeType;
    typedef PositiveFloatLowerBound<PR> PositiveType;
    typedef ValidatedLowerKleenean LessType;
    typedef ValidatedNegatedSierpinskian EqualsType;
};
template<class PR> struct NumericTraits<FloatUpperBound<PR>> {
    typedef ValidatedUpperNumber GenericType;
    typedef FloatLowerBound<PR> OppositeType;
    typedef PositiveFloatUpperBound<PR> PositiveType;
    typedef ValidatedUpperKleenean LessType;
    typedef ValidatedNegatedSierpinskian EqualsType;
};
template<class PR> struct NumericTraits<FloatBounds<PR>> {
    typedef ValidatedNumber GenericType;
    typedef PositiveFloatBounds<PR> PositiveType;
    typedef ValidatedKleenean LessType;
    typedef ValidatedKleenean EqualsType;
};
template<class PR> struct NumericTraits<FloatBall<PR>> {
    typedef ValidatedNumber GenericType;
    typedef PositiveFloatBall<PR> PositiveType;
    typedef ValidatedKleenean LessType;
    typedef ValidatedKleenean EqualsType;
};
template<class PR> struct NumericTraits<FloatValue<PR>> {
    typedef ExactNumber GenericType;
    typedef PositiveFloatValue<PR> PositiveType;
    typedef Boolean LessType;
    typedef Boolean EqualsType;
};

template<class X> using GenericTrait = typename NumericTraits<X>::GenericType;
template<class X> using PositiveTrait = typename NumericTraits<X>::PositiveType;
template<class X> using OppositeTrait = typename NumericTraits<X>::OppositeType;
template<class X> using LessTrait = typename NumericTraits<X>::LessType;
template<class X> using EqualsTrait = typename NumericTraits<X>::EqualsType;

template<class PR> inline FloatBall<PR> make_float_ball(Real const& y, PR pr) { return FloatBall<PR>(make_float(y,pr)); }

template<class X, class R=X> class DeclareFloatOperations
    : public DeclareRealOperations<X,R>
    , public DeclareComparisonOperations<X,LessTrait<X>,EqualsTrait<X>>
    , public DeclareMixedFieldOperators<X,GenericTrait<X>,R>
{
    friend OutputStream& operator<<(OutputStream&, X const&);
    friend InputStream& operator>>(InputStream&, X&);
};

template<class X, class NX=OppositeTrait<X>> class DeclareDirectedFloatOperations
    : DeclareDirectedNumericOperations<X,NX>
    , DeclareMixedDirectedGroupOperators<X,NX,GenericTrait<X>,GenericTrait<NX>>
{
    friend OutputStream& operator<<(OutputStream&, X const&);
    friend InputStream& operator>>(InputStream&, X&);
};

template<class PX, class QPX=OppositeTrait<PX>> class DeclarePositiveDirectedFloatOperations
    : DeclarePositiveDirectedNumericOperations<PX,QPX>
    , DeclareMixedDirectedSemifieldOperators<PX,QPX,GenericTrait<PX>,GenericTrait<QPX>>
{
    friend OutputStream& operator<<(OutputStream&, PX const&);
    friend InputStream& operator>>(InputStream&, PX&);
};

template<class PX> class DeclarePositiveFloatOperations
    : DeclarePositiveDirectedFloatOperations<PX,PX>
{
};

template<class X, class R=X> class DispatchFloatOperations
    : DispatchNumericOperations<X,R>
    , DispatchComparisonOperations<X,LessTrait<X>,EqualsTrait<X>>
    , ProvideConcreteGenericFieldOperations<X,GenericTrait<X>,R>
    , ProvideConcreteGenericComparisonOperations<X,GenericTrait<X>,LessTrait<X>,EqualsTrait<X>>
{
    friend OutputStream& operator<<(OutputStream& os, X const& x) { return Operations<X>::_write(os,x); }
    friend InputStream& operator>>(InputStream& is, X& x) { return Operations<X>::_read(is,x); }
};

template<class X> class DispatchDirectedFloatOperations
    : DispatchDirectedNumericOperations<X,OppositeTrait<X>>
    , DispatchDirectedNumericOperations<OppositeTrait<X>,X>
    , DispatchDirectedComparisonOperations<X,OppositeTrait<X>,LessTrait<X>,EqualsTrait<X>>
    , DispatchDirectedComparisonOperations<OppositeTrait<X>,X,LessTrait<OppositeTrait<X>>,EqualsTrait<OppositeTrait<X>>>
    , ProvideConcreteGenericDirectedGroupOperations<X,OppositeTrait<X>,GenericTrait<X>,GenericTrait<OppositeTrait<X>>>
    , ProvideConcreteGenericDirectedGroupOperations<OppositeTrait<X>,X,GenericTrait<OppositeTrait<X>>,GenericTrait<X>>
    , ProvideConcreteGenericDirectedComparisonOperations<X,GenericTrait<OppositeTrait<X>>,LessTrait<X>,EqualsTrait<X>>
    , ProvideConcreteGenericDirectedComparisonOperations<OppositeTrait<X>,GenericTrait<X>,LessTrait<OppositeTrait<X>>,EqualsTrait<OppositeTrait<X>>>
{
    friend OutputStream& operator<<(OutputStream& os, X const& x) { return Operations<X>::_write(os,x); }
    friend InputStream& operator>>(InputStream& is, X& x) { return Operations<X>::_read(is,x); }
};

template<class PX, class QPX=OppositeTrait<PX>> class DispatchPositiveDirectedFloatOperations
    : public DispatchPositiveDirectedNumericOperations<PX,QPX>
    , public DispatchPositiveDirectedNumericOperations<QPX,PX>
    , public ProvideConcreteGenericDirectedSemiFieldOperations<PX,QPX,Nat,Nat>
    , public ProvideConcreteGenericDirectedSemiFieldOperations<QPX,PX,Nat,Nat>
{
};

template<class PX> class DispatchPositiveFloatOperations
    : public DispatchPositiveDirectedNumericOperations<PX,PX>
    , public ProvideConcreteGenericDirectedSemiFieldOperations<PX,PX,Nat,Nat>
{
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
template<class PR> class FloatApproximation
    : public DispatchFloatOperations<FloatApproximation<PR>>
{
    typedef ApproximateTag P; typedef RawFloat<PR> FLT;
  public:
    typedef ApproximateTag Paradigm;
    typedef FloatApproximation<PR> NumericType;
    typedef ApproximateNumber GenericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
  public:
    FloatApproximation<PR>() : _a(0.0) { }
    explicit FloatApproximation<PR>(PrecisionType pr) : _a(0.0,pr) { }
    explicit FloatApproximation<PR>(RawFloatType const& a) : _a(a) { }

        FloatApproximation<PR>(double d, PR pr);
        FloatApproximation<PR>(ExactDouble d, PR pr);
        FloatApproximation<PR>(const Integer& z, PR pr);
        FloatApproximation<PR>(const Dyadic& w, PR pr);
        FloatApproximation<PR>(const Decimal& d, PR pr);
        FloatApproximation<PR>(const Rational& q, PR pr);
        FloatApproximation<PR>(const Real& r, PR pr);
        FloatApproximation<PR>(const FloatApproximation<PR>& r, PR pr);
    FloatApproximation<PR>(const ApproximateNumber& y, PR pr);

    FloatApproximation<PR>(FloatError<PR> const& x); // FIXME: Remove
    FloatApproximation<PR>(FloatValue<PR> const& x);
    FloatApproximation<PR>(FloatBall<PR> const& x);
    FloatApproximation<PR>(FloatBounds<PR> const& x);
    FloatApproximation<PR>(FloatUpperBound<PR> const& x);
    FloatApproximation<PR>(FloatLowerBound<PR> const& x);

    template<class N, EnableIf<IsIntegral<N>> =dummy> FloatApproximation<PR>& operator=(N n) { this->_a=n; return *this; }
    template<class D, EnableIf<IsFloatingPoint<D>> =dummy> FloatApproximation<PR>& operator=(D x) { this->_a=x; return *this; }
        FloatApproximation<PR>& operator=(const FloatLowerBound<PR>& x) { return *this=FloatApproximation<PR>(x); }
        FloatApproximation<PR>& operator=(const FloatUpperBound<PR>& x) { return *this=FloatApproximation<PR>(x); }
        FloatApproximation<PR>& operator=(const FloatBounds<PR>& x) { return *this=FloatApproximation<PR>(x); }
        FloatApproximation<PR>& operator=(const FloatValue<PR>& x) { return *this=FloatApproximation<PR>(x); }
    FloatApproximation<PR>& operator=(const ApproximateNumber& y);
    FloatApproximation<PR> create(const ApproximateNumber& y) const;

    operator ApproximateNumber () const;

    friend FloatApproximation<PR> round(FloatApproximation<PR> const& x);

    PrecisionType precision() const { return _a.precision(); }
    explicit operator RawFloatType () const { return this->_a; }
    RawFloatType const& raw() const { return this->_a; }
    RawFloatType& raw() { return this->_a; }
    double get_d() const { return this->_a.get_d(); }
  public:
    friend Bool same(FloatApproximation<PR> const&, FloatApproximation<PR> const&);
    friend PositiveFloatApproximation<PR> mag(FloatApproximation<PR> const&);
  public:
    static Void set_output_places(Nat p) { output_places=p; }
    FloatApproximation<PR> pm(FloatApproximation<PR> _e) { return *this; }
  private: public:
    static Nat output_places;
    RawFloatType _a;
};


//! \ingroup NumericModule
//! \brief Floating-point lower bounds for real numbers.
template<class PR> class FloatLowerBound
    : public DispatchDirectedFloatOperations<FloatLowerBound<PR>>
    , public DispatchFloatOperations<FloatApproximation<PR>>
{
    typedef LowerTag P; typedef RawFloat<PR> FLT;
  public:
    typedef LowerTag Paradigm;
    typedef FloatLowerBound<PR> NumericType;
    typedef ValidatedLowerNumber GenericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
  public:
    FloatLowerBound<PR>() : _l(0.0) { }
    explicit FloatLowerBound<PR>(PrecisionType pr) : _l(0.0,pr) { }
    explicit FloatLowerBound<PR>(RawFloatType const& l) : _l(l) { }

    template<class N, EnableIf<IsIntegral<N>> = dummy> FloatLowerBound<PR>(N n, PR pr) : FloatLowerBound<PR>(ExactDouble(n),pr) { }
    FloatLowerBound<PR>(ExactDouble d, PR pr);
        FloatLowerBound<PR>(const Integer& z, PR pr);
        FloatLowerBound<PR>(const Dyadic& w, PR pr);
        FloatLowerBound<PR>(const Decimal& d, PR pr);
        FloatLowerBound<PR>(const Rational& q, PR pr);
        FloatLowerBound<PR>(const Real& r, PR pr);
    FloatLowerBound<PR>(const FloatLowerBound<PR>& x, PR pr);
    FloatLowerBound<PR>(const ValidatedLowerNumber& y, PR pr);

    FloatLowerBound<PR>(FloatBounds<PR> const& x);
    FloatLowerBound<PR>(FloatBall<PR> const& x);
    FloatLowerBound<PR>(FloatValue<PR> const& x);

        FloatLowerBound<PR>& operator=(const FloatValue<PR>& x) { return *this=FloatLowerBound<PR>(x); }
    FloatLowerBound<PR>& operator=(const ValidatedLowerNumber&);
    FloatLowerBound<PR> create(const ValidatedLowerNumber& y) const;
    FloatUpperBound<PR> create(const ValidatedUpperNumber& y) const;

    operator ValidatedLowerNumber () const;

    PrecisionType precision() const { return _l.precision(); }
    RawFloatType const& raw() const { return _l; }
    RawFloatType& raw() { return _l; }
    double get_d() const { return _l.get_d(); }
  public: // To be removed
    friend Bool same(FloatLowerBound<PR> const&, FloatLowerBound<PR> const&);
    friend Bool refines(FloatLowerBound<PR> const&, FloatLowerBound<PR> const&);
    friend FloatLowerBound<PR> refinement(FloatLowerBound<PR> const&, FloatLowerBound<PR> const&);
  private: public:
    static Nat output_places;
    RawFloatType _l;
};


//! \ingroup NumericModule
//! \brief Floating-point upper bounds for real numbers.
template<class PR> class FloatUpperBound
    : public DispatchDirectedFloatOperations<FloatUpperBound<PR>>
    , public DispatchFloatOperations<FloatApproximation<PR>>
{
    typedef UpperTag P; typedef RawFloat<PR> FLT;
  public:
    typedef UpperTag Paradigm;
    typedef FloatUpperBound<PR> NumericType;
    typedef ValidatedUpperNumber GenericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
  public:
    FloatUpperBound<PR>() : _u(0.0) { }
    explicit FloatUpperBound<PR>(PrecisionType pr) : _u(0.0,pr) { }
    explicit FloatUpperBound<PR>(RawFloatType const& u) : _u(u) { }

    template<class N, EnableIf<IsIntegral<N>> = dummy> FloatUpperBound<PR>(N n, PR pr) : FloatUpperBound<PR>(ExactDouble(n),pr) { }
    FloatUpperBound<PR>(ExactDouble d, PR pr);
        FloatUpperBound<PR>(const Integer& z, PR pr);
        FloatUpperBound<PR>(const Dyadic& w, PR pr);
        FloatUpperBound<PR>(const Decimal& d, PR pr);
        FloatUpperBound<PR>(const Rational& q, PR pr);
        FloatUpperBound<PR>(const Real& r, PR pr);
    FloatUpperBound<PR>(const FloatUpperBound<PR>& x, PR pr);
    FloatUpperBound<PR>(const ValidatedUpperNumber& y, PR pr);

    FloatUpperBound<PR>(FloatBounds<PR> const& x);
    FloatUpperBound<PR>(FloatBall<PR> const& x);
    FloatUpperBound<PR>(FloatValue<PR> const& x);
    FloatUpperBound<PR>(FloatError<PR> const& x); // FIXME: Remove

        FloatUpperBound<PR>& operator=(const FloatValue<PR>& x) { return *this=FloatUpperBound<PR>(x); }
    FloatUpperBound<PR>& operator=(const ValidatedUpperNumber& y);
    FloatUpperBound<PR> create(const ValidatedUpperNumber& y) const;
    FloatLowerBound<PR> create(const ValidatedLowerNumber& y) const;

    operator ValidatedUpperNumber () const;

    PrecisionType precision() const { return _u.precision(); }
    RawFloatType const& raw() const { return _u; }
    RawFloatType& raw() { return _u; }
    double get_d() const { return _u.get_d(); }
  public: // To be removed
    friend Bool same(FloatUpperBound<PR> const&, FloatUpperBound<PR> const&);
    friend Bool refines(FloatUpperBound<PR> const&, FloatUpperBound<PR> const&);
    friend FloatUpperBound<PR> refinement(FloatUpperBound<PR> const&, FloatUpperBound<PR> const&);
  private: public:
    static Nat output_places;
    RawFloatType _u;
};



//! \ingroup NumericModule
//! \brief ValidatedTag bounds on a number with floating-point endpoints supporting outwardly-rounded arithmetic.
//! \details
//! Note that direct construction from a floating-point number is prohibited, since <c>%Float64Bounds(3.3)</c> would the singleton interval \f$[3.2999999999999998224,3.2999999999999998224]\f$ (the constant is first interpreted by the C++ compiler to give a C++ \c double, whereas <c>%Float64Bounds(3.3_decimal)</c> yields the interval \f$[3.2999999999999998224,3.3000000000000002665]\f$ enclosing \f$3.3\f$.
//!
//! Comparison tests on \c FloatBounds use the idea that an interval represents a single number with an unknown value.
//! Hence the result is of type \c ValidatedKleenean, which can take values { \c True, \c False, \c Indeterminate }.
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
template<class PR> class FloatBounds
    : public DispatchFloatOperations<FloatBounds<PR>>
//    , public ProvideConvertedFieldOperations<FloatBounds<PR>,FloatValue<PR>>
{
    typedef BoundedTag P; typedef RawFloat<PR> FLT;
  public:
    typedef BoundedTag Paradigm;
    typedef FloatBounds<PR> NumericType;
    typedef ValidatedNumber GenericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
  public:
    FloatBounds<PR>() : _l(0.0), _u(0.0) { }
    explicit FloatBounds<PR>(PrecisionType pr) : _l(0.0,pr), _u(0.0,pr) { }
    explicit FloatBounds<PR>(RawFloatType const& v) : _l(v), _u(v) { }
    explicit FloatBounds<PR>(RawFloatType const& l, RawFloatType const& u) : _l(l), _u(u) { }
    FloatBounds<PR>(FloatLowerBound<PR> const& lower, FloatUpperBound<PR> const& upper) : _l(lower.raw()), _u(upper.raw()) { }
    FloatBounds<PR>(FloatLowerBound<PR> const& lower, ValidatedUpperNumber const& upper) : FloatBounds<PR>(lower,lower.create(upper)) { }
    FloatBounds<PR>(ValidatedLowerNumber const& lower, FloatUpperBound<PR> const& upper) : FloatBounds<PR>(upper.create(lower),upper) { }
    template<class N1, class N2, EnableIf<And<IsIntegral<N1>,IsIntegral<N2>>> = dummy> FloatBounds<PR>(N1 n1, N2 n2, PR pr) : _l(n1,pr), _u(n2,pr) { }
    FloatBounds<PR>(ExactDouble const& dl, ExactDouble const& du, PrecisionType pr);
    FloatBounds<PR>(Dyadic const& wl, Dyadic const& wu, PrecisionType pr);
    FloatBounds<PR>(Rational const& ql, Rational const& qu, PrecisionType pr);
    FloatBounds<PR>(Pair<ExactDouble,ExactDouble> const& dlu, PrecisionType pr);

    template<class N, EnableIf<IsIntegral<N>> = dummy> FloatBounds<PR>(N n, PR pr) : FloatBounds<PR>(ExactDouble(n),pr) { }
    FloatBounds<PR>(ExactDouble d, PR pr);
        FloatBounds<PR>(const Integer& z, PR pr);
        FloatBounds<PR>(const Dyadic& w, PR pr);
        FloatBounds<PR>(const Decimal& w, PR pr);
        FloatBounds<PR>(const Rational& q, PR pr);
        FloatBounds<PR>(const Real& x, PR pr);
        FloatBounds<PR>(const FloatBounds<PR>& x, PR pr);
    FloatBounds<PR>(const ValidatedNumber& y, PR pr);

    FloatBounds<PR>(FloatBall<PR> const& x);
    FloatBounds<PR>(FloatValue<PR> const& x);

        FloatBounds<PR>& operator=(const FloatValue<PR>& x) { return *this=FloatBounds<PR>(x); }
    FloatBounds<PR>& operator=(const ValidatedNumber& y);
    template<class N, EnableIf<IsIntegral<N>> = dummy> FloatBounds<PR>& operator=(N n) { return *this=ValidatedNumber(n); }

    operator ValidatedNumber () const;

    FloatBounds<PR> create(const ValidatedNumber& y) const;

    FloatLowerBound<PR> const lower() const { return FloatLowerBound<PR>(lower_raw()); }
    FloatUpperBound<PR> const upper() const { return FloatUpperBound<PR>(upper_raw()); }
    FloatValue<PR> const value() const;
    FloatError<PR> const error() const;

    RawFloatType const& lower_raw() const { return _l; }
    RawFloatType const& upper_raw() const { return _u; }
    RawFloatType const value_raw() const { return hlf(add_near(_l,_u)); }
    RawFloatType const error_raw() const { RawFloatType v=value_raw(); return max(sub_up(_u,v),sub_up(v,_l)); }
    double get_d() const { return value_raw().get_d(); }

    PrecisionType precision() const { ARIADNE_DEBUG_ASSERT(_l.precision()==_u.precision()); return _u.precision(); }
    GenericType generic() const;

    FloatBounds<PR> pm(FloatError<PR> e) const;

    // DEPRECATED
    explicit operator RawFloatType () const { return value_raw(); }
    friend FloatApproximation<PR> round(FloatApproximation<PR> const& x);
    friend FloatValue<PR> midpoint(FloatBounds<PR> const& x);
  public:
    friend PositiveFloatUpperBound<PR> mag(FloatBounds<PR> const&);
    friend PositiveFloatLowerBound<PR> mig(FloatBounds<PR> const&);
    friend Bool same(FloatBounds<PR> const&, FloatBounds<PR> const&);
    friend Bool models(FloatBounds<PR> const&, FloatValue<PR> const&);
    friend Bool consistent(FloatBounds<PR> const&, FloatBounds<PR> const&);
    friend Bool inconsistent(FloatBounds<PR> const&, FloatBounds<PR> const&);
    friend Bool refines(FloatBounds<PR> const&, FloatBounds<PR> const&);
    friend FloatBounds<PR> refinement(FloatBounds<PR> const&, FloatBounds<PR> const&);
  public:
    static Nat output_places;
    static Void set_output_places(Nat p) { output_places=p; }
  private: public:
    RawFloatType _l, _u;
};


template<class PR> class FloatBall
    : public DispatchFloatOperations<FloatBall<PR>>
    , public ProvideConvertedFieldOperations<FloatBounds<PR>,FloatBall<PR>>
    , public ProvideConvertedFieldOperations<FloatBall<PR>,FloatValue<PR>>
{
    typedef MetricTag P; typedef RawFloat<PR> FLT;
  public:
    typedef MetricTag Paradigm;
    typedef FloatBall<PR> NumericType;
    typedef ValidatedNumber GenericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
  public:
    FloatBall<PR>() : _v(0.0), _e(0.0) { }
    explicit FloatBall<PR>(PrecisionType pr) : _v(0.0,pr), _e(0.0,pr) { }
    explicit FloatBall<PR>(RawFloatType const& v) : _v(v), _e(0.0) { }
    explicit FloatBall<PR>(RawFloatType const& v, RawFloatType const& e) : _v(v), _e(e) { }
    FloatBall<PR>(FloatValue<PR> const& value, FloatError<PR> const& error) : _v(value.raw()), _e(error.raw()) { }
    FloatBall<PR>(FloatLowerBound<PR> const& lower, FloatUpperBound<PR> const& upper) =  delete;

    FloatBall<PR>(ExactDouble d, PR pr);
        FloatBall<PR>(const Integer& z, PR pr);
        FloatBall<PR>(const Dyadic& w, PR pr);
        FloatBall<PR>(const Decimal& d, PR pr);
        FloatBall<PR>(const Rational& q, PR pr);
        FloatBall<PR>(const Real& r, PR pr);
        FloatBall<PR>(const FloatBall<PR>& x, PR pr);
    FloatBall<PR>(const ValidatedNumber& y, PR pr);

    explicit FloatBall<PR>(FloatBounds<PR> const& x);
    FloatBall<PR>(FloatValue<PR> const& x);

    FloatBall<PR>& operator=(const ValidatedNumber& y);

    operator ValidatedNumber () const;

    FloatBall<PR> create(const ValidatedNumber& y) const;

    FloatLowerBound<PR> const lower() const { return FloatLowerBound<PR>(lower_raw()); }
    FloatUpperBound<PR> const upper() const { return FloatUpperBound<PR>(upper_raw()); }
    FloatValue<PR> const value() const;
    FloatError<PR> const error() const;

    RawFloatType const lower_raw() const { return sub_down(_v,_e); }
    RawFloatType const upper_raw() const { return add_up(_v,_e); }
    RawFloatType const& value_raw() const { return _v; }
    RawFloatType const& error_raw() const { return _e; }
    double get_d() const { return _v.get_d(); }

    PrecisionType precision() const { return _v.precision(); }
    FloatBall<PR> pm(FloatError<PR> e) const;
  public:
    friend PositiveFloatUpperBound<PR> mag(FloatBall<PR> const&);
    friend PositiveFloatLowerBound<PR> mig(FloatBall<PR> const&);
    friend Bool same(FloatBall<PR> const&, FloatBall<PR> const&);
    friend Bool same(FloatBall<PR> const&, FloatBall<PR> const&);
    friend Bool models(FloatBall<PR> const&, FloatValue<PR> const&);
    friend Bool consistent(FloatBall<PR> const&, FloatBall<PR> const&);
    friend Bool refines(FloatBall<PR> const&, FloatBall<PR> const&);
    friend FloatBall<PR> refinement(FloatBall<PR> const&, FloatBall<PR> const&);
  private: public:
    RawFloatType _v, _e;
};

//! \ingroup NumericModule
//! \related Float64, Float64Bounds
//! \brief A floating-point number, which is taken to represent the \em exact value of a real quantity.
template<class PR> class FloatValue
    : DispatchNumericOperations<FloatValue<PR>,FloatBounds<PR>>
    , DispatchComparisonOperations<FloatValue<PR>,Boolean>
    , DefineMixedComparisonOperators<FloatValue<PR>,ExactNumber,Boolean>
{
    typedef ExactTag P; typedef RawFloat<PR> FLT;
  public:
    typedef ExactTag Paradigm;
    typedef FloatValue<PR> NumericType;
    typedef ExactNumber GenericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
  public:
    FloatValue<PR>() : _v(0.0) { }
    explicit FloatValue<PR>(PrecisionType pr) : _v(0.0,pr) { }
    explicit FloatValue<PR>(RawFloatType const& v) : _v(v) { }

    template<class N, EnableIf<IsIntegral<N>> = dummy> FloatValue<PR>(N n, PR pr) : FloatValue<PR>(ExactDouble(n),pr) { }
    FloatValue<PR>(ExactDouble d, PR pr);
    FloatValue<PR>(const Integer& z, PR pr);
    FloatValue<PR>(const TwoExp& t, PR pr);
    FloatValue<PR>(const Dyadic& w, PR pr);
    FloatValue<PR>(const FloatValue<PR>& x, PR pr);

    template<class N, EnableIf<IsIntegral<N>> = dummy> FloatValue<PR>& operator=(N n) { _v=n; return *this; }
    FloatValue<PR>& operator=(const Integer& z);
    FloatValue<PR>& operator=(const TwoExp& t);
    FloatValue<PR>& operator=(const Dyadic& w);

    operator ExactNumber () const;
    explicit operator Dyadic () const;
    explicit operator Rational () const;

    FloatBall<PR> create(ValidatedNumber const&) const;
//    explicit operator RawFloatType () const { return _v; }

    PrecisionType precision() const { return _v.precision(); }
    RawFloatType const& raw() const { return _v; }
    RawFloatType& raw() { return _v; }
    double get_d() const { return _v.get_d(); }

    FloatBall<PR> pm(FloatError<PR> _e) const;
  public:
    friend FloatValue<PR> operator*(FloatValue<PR> const&, TwoExp const&);
    friend FloatValue<PR> operator/(FloatValue<PR> const&, TwoExp const&);
    friend FloatValue<PR>& operator*=(FloatValue<PR>&, TwoExp const&);
    friend FloatValue<PR>& operator/=(FloatValue<PR>&, TwoExp const&);
    friend FloatError<PR> mag(FloatValue<PR> const&);
    friend FloatLowerBound<PR> mig(FloatValue<PR> const&);
    friend Bool same(FloatValue<PR> const&, FloatValue<PR> const&);
    friend OutputStream& operator<<(OutputStream&, FloatValue<PR> const&);
  public:
    static Nat output_places;
    static Void set_output_places(Nat p) { output_places=p; }
  private: public:
    RawFloatType _v;
  private:
    friend FloatValue<PR> operator*(FloatValue<PR> const& x, TwoExp const& y) {
        return FloatValue<PR>(x.raw()*RawFloat<PR>(y,x.precision())); }
    friend FloatValue<PR> operator/(FloatValue<PR> const& x, TwoExp const& y) {
        return FloatValue<PR>(x.raw()/RawFloat<PR>(y,x.precision())); }
    friend FloatValue<PR>& operator*=(FloatValue<PR>& x, TwoExp const& y) { return x=x*y; }
    friend FloatValue<PR>& operator/=(FloatValue<PR>& x, TwoExp const& y) { return x=x/y; }
    friend OutputStream& operator<<(OutputStream& os, FloatValue<PR> const& x) {
        return Operations<FloatValue<PR>>::_write(os,x); }
};


template<class PR> class PositiveFloatValue : public FloatValue<PR> {
  public:
    PositiveFloatValue<PR>() : FloatValue<PR>() { }
    explicit PositiveFloatValue<PR>(PR const& pr) : FloatValue<PR>(pr) { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy> PositiveFloatValue<PR>(M m, PR pr) : FloatValue<PR>(m,pr) { }
    PositiveFloatValue<PR>(TwoExp const& ex, PR pr) : FloatValue<PR>(ex,pr) { }
    explicit PositiveFloatValue<PR>(Dyadic const& w, PR pr) : FloatValue<PR>(w,pr) { }
    explicit PositiveFloatValue<PR>(RawFloat<PR> const& x) : FloatValue<PR>(x) { }
    explicit PositiveFloatValue<PR>(FloatValue<PR> const& x) : FloatValue<PR>(x) { }
  public:
    friend PositiveFloatValue<PR> hlf(PositiveFloatValue<PR> const&);
};

template<class PR> class PositiveFloatBall : public FloatBall<PR> {
  public:
    PositiveFloatBall<PR>() : FloatBounds<PR>() { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy>
        PositiveFloatBall<PR>(M m, PR pr) : FloatBall<PR>(m,pr) { }
    explicit PositiveFloatBall<PR>(FloatBall<PR> const& x) : FloatBall<PR>(x) { }
};

template<class PR> class PositiveFloatBounds : public FloatBounds<PR>
    , public DispatchPositiveFloatOperations<PositiveFloatBounds<PR>>
{
  public:
    PositiveFloatBounds<PR>() : FloatBounds<PR>() { }
    explicit PositiveFloatBounds<PR>(PR const& pr) : FloatBounds<PR>(pr) { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy> PositiveFloatBounds<PR>(M m, PR pr) : FloatBounds<PR>(m,pr) { }
    explicit PositiveFloatBounds<PR>(RawFloat<PR> const& x) : FloatBounds<PR>(x) { }
    explicit PositiveFloatBounds<PR>(RawFloat<PR> const& l, RawFloat<PR> const& u) : FloatBounds<PR>(l,u) { }
    explicit PositiveFloatBounds<PR>(FloatBounds<PR> const& x) : FloatBounds<PR>(x) { }
  public:
};

template<class PR> class PositiveFloatUpperBound : public FloatUpperBound<PR>
    , public DispatchPositiveDirectedFloatOperations<PositiveFloatUpperBound<PR>,PositiveFloatLowerBound<PR>>
{
  public:
    PositiveFloatUpperBound<PR>() : FloatUpperBound<PR>() { }
    explicit PositiveFloatUpperBound<PR>(PR const& pr) : FloatUpperBound<PR>(pr) { }
    explicit PositiveFloatUpperBound<PR>(RawFloat<PR> const& x) : FloatUpperBound<PR>(x) { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy> PositiveFloatUpperBound<PR>(M m, PR pr) : FloatUpperBound<PR>(m,pr) { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy> PositiveFloatValue<PR> create(M m) const { return PositiveFloatValue<PR>(m,this->precision()); }
    explicit PositiveFloatUpperBound<PR>(FloatUpperBound<PR> const& x) : FloatUpperBound<PR>(x) { ARIADNE_PRECONDITION_MSG(!(this->_u<0),"x="<<x); }
    explicit PositiveFloatUpperBound<PR>(ValidatedUpperNumber const& y, PR pr) : FloatUpperBound<PR>(y,pr) { ARIADNE_PRECONDITION_MSG(!(this->_u<0),"y="<<y); }
    PositiveFloatUpperBound<PR>(PositiveFloatValue<PR> const& x) : FloatUpperBound<PR>(x) { }
  public:
};

template<class PR> class PositiveFloatLowerBound : public FloatLowerBound<PR>
    , public DispatchPositiveDirectedFloatOperations<PositiveFloatLowerBound<PR>,PositiveFloatUpperBound<PR>>
{
  public:
    PositiveFloatLowerBound<PR>() : FloatLowerBound<PR>() { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy>
        PositiveFloatLowerBound<PR>(M m) : FloatLowerBound<PR>(m) { }
    explicit PositiveFloatLowerBound<PR>(RawFloat<PR> const& x) : FloatLowerBound<PR>(x) { }
    explicit PositiveFloatLowerBound<PR>(FloatLowerBound<PR> const& x) : FloatLowerBound<PR>(x) { }
    explicit PositiveFloatLowerBound<PR>(ValidatedLowerNumber const& y, PR pr) : FloatLowerBound<PR>(y,pr) { }
    PositiveFloatLowerBound<PR>(PositiveFloatValue<PR> const& x) : FloatLowerBound<PR>(x) { }
  public:
};

template<class PR> class PositiveFloatApproximation : public FloatApproximation<PR>
    , public DispatchPositiveFloatOperations<PositiveFloatApproximation<PR>>
{
  public:
    PositiveFloatApproximation<PR>() : FloatApproximation<PR>() { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy>
        PositiveFloatApproximation<PR>(M m) : FloatApproximation<PR>(m) { }
    explicit PositiveFloatApproximation<PR>(RawFloat<PR> const& x) : FloatApproximation<PR>(x) { }
    explicit PositiveFloatApproximation<PR>(FloatApproximation<PR> const& x) : FloatApproximation<PR>(x) { }
    explicit PositiveFloatApproximation<PR>(ApproximateNumber const& y, PR pr) : FloatApproximation<PR>(y,pr) { }
    PositiveFloatApproximation<PR>(PositiveFloatLowerBound<PR> const& x) : FloatApproximation<PR>(x) { }
    PositiveFloatApproximation<PR>(PositiveFloatUpperBound<PR> const& x) : FloatApproximation<PR>(x) { }
    PositiveFloatApproximation<PR>(PositiveFloatValue<PR> const& x) : FloatApproximation<PR>(x) { }
    PositiveFloatApproximation<PR>(FloatError<PR> const& x) : FloatApproximation<PR>(x) { }
  public:
};

template<class PR> inline PositiveFloatApproximation<PR> cast_positive(FloatApproximation<PR> const& x) {
    return PositiveFloatApproximation<PR>(x); }

template<class PR> inline PositiveFloatLowerBound<PR> cast_positive(FloatLowerBound<PR> const& x) {
    return PositiveFloatLowerBound<PR>(x); }

template<class PR> inline PositiveFloatUpperBound<PR> cast_positive(FloatUpperBound<PR> const& x) {
    return PositiveFloatUpperBound<PR>(x); }

template<class PR> inline PositiveFloatBounds<PR> cast_positive(FloatBounds<PR> const& x) {
    return PositiveFloatBounds<PR>(x); }

template<class PR> inline PositiveFloatBall<PR> cast_positive(FloatBall<PR> const& x) {
    return PositiveFloatBall<PR>(x); }

template<class PR> inline PositiveFloatValue<PR> cast_positive(FloatValue<PR> const& x) {
    return PositiveFloatValue<PR>(x); }


template<class PR> class FloatError
    : public DispatchDirectedFloatOperations<FloatUpperBound<PR>>
    , public DispatchPositiveDirectedNumericOperations<PositiveFloatUpperBound<PR>,PositiveFloatLowerBound<PR>>
    , public ProvideConcreteGenericDirectedSemiFieldOperations<PositiveFloatUpperBound<PR>,PositiveFloatLowerBound<PR>,Nat,Nat>
{
  private: public:
    RawFloat<PR> _e;
  public:
    typedef PR PrecisionType;
  public:
    FloatError<PR>(PositiveFloatUpperBound<PR> const& x) : _e(x._u) { }
    operator PositiveFloatUpperBound<PR> const& () const { return reinterpret_cast<PositiveFloatUpperBound<PR>const&>(*this); }
    operator PositiveFloatUpperBound<PR>& () { return reinterpret_cast<PositiveFloatUpperBound<PR>&>(*this); }
  public:
    FloatError<PR>() : _e() { }
    explicit FloatError<PR>(PR const& pr) : _e(pr) { }
    explicit FloatError<PR>(RawFloat<PR> const& x) : _e(x) { ARIADNE_PRECONDITION_MSG((this->_e>=0),"e="<<*this); }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy> FloatError<PR>(M m, PR pr) : _e(m,pr) { }
    explicit FloatError<PR>(FloatUpperBound<PR> const& x) : FloatError<PR>(x._u) { }
    explicit FloatError<PR>(ValidatedUpperNumber const& y, PR pr) : FloatError(FloatUpperBound<PR>(y,pr)) { }
    FloatError<PR>(PositiveFloatValue<PR> const& x) : _e(x._v) { }
    FloatError<PR>& operator=(Nat m) { reinterpret_cast<FloatUpperBound<PR>&>(*this)=m; return *this; }
  public:
    PrecisionType precision() const { return _e.precision(); }
    RawFloat<PR> const& raw() const { return _e; }
    RawFloat<PR>& raw() { return _e; }
  public:
    friend FloatError<PR> mag(FloatError<PR> const& x) { return x; }
    friend FloatUpperBound<PR> operator+(FloatError<PR> const& x) { return FloatUpperBound<PR>(+x._e); }
    friend FloatLowerBound<PR> operator-(FloatError<PR> const& x) { return FloatLowerBound<PR>(-x._e); }
    friend FloatUpperBound<PR> operator+(FloatValue<PR> const& x1, FloatError<PR> const& x2) { return FloatUpperBound<PR>(add_up(x1._v,x2._e)); }
    friend FloatLowerBound<PR> operator-(FloatValue<PR> const& x1, FloatError<PR> const& x2) { return FloatLowerBound<PR>(sub_down(x1._v,x2._e)); }

    friend Bool same(FloatError<PR> const& x1, FloatError<PR> const& x2) { return x1._e==x2._e; }
    friend Bool refines(FloatError<PR> const& x1, FloatError<PR> const& x2) { return x1._e<=x2._e; }
    friend FloatError<PR> refinement(FloatError<PR> const& x1, FloatError<PR> const& x2) { return FloatError<PR>(min(x1._e,x2._e)); }
    friend OutputStream& operator<<(OutputStream& os, FloatError<PR> const& x) { return Operations<FloatError<PR>>::_write(os,x); }
  public:
    static Nat output_places;
    static Void set_output_places(Nat p) { output_places=p; }
};


template<class PR> inline const FloatValue<PR> FloatBounds<PR>::value() const {
    return FloatValue<PR>(med_near(this->_l,this->_u)); }

template<class PR> inline const FloatError<PR> FloatBounds<PR>::error() const {
    RawFloat<PR> _v=med_near(this->_l,this->_u); return FloatError<PR>(max(sub_up(this->_u,_v),sub_up(_v,this->_l))); }

template<class PR> inline FloatValue<PR> value(FloatBounds<PR> const& x) {
    return x.value(); }

template<class PR> inline FloatError<PR> error(FloatBounds<PR> const& x) {
    return x.error(); }

template<class PR> inline const FloatValue<PR> FloatBall<PR>::value() const {
    return FloatValue<PR>(this->_v); }

template<class PR> inline const FloatError<PR> FloatBall<PR>::error() const {
    return FloatError<PR>(this->_e); }


template<class R, class A> R integer_cast(const A& _a);

template<class T, class F, EnableIf<Not<IsSame<T,F>>> =dummy> T convert(F const& x) { return T(x); }
template<class T> T const& convert(T const& x) { return x; }

template<class T> using NumericType = typename T::NumericType;



extern const FloatValue<Precision64> infty;



// Literals operations
Float64Value operator"" _exact(long double lx);
Float64Error operator"" _error(long double lx);
Float64Ball operator"" _near(long double lx);
Float64UpperBound operator"" _upper(long double lx);
Float64LowerBound operator"" _lower(long double lx);
Float64Approximation operator"" _approx(long double lx);


// ValidatedTag operations
template<class PR> FloatBounds<PR> make_bounds(FloatError<PR> const& e) {
    return FloatBounds<PR>(-e.raw(),+e.raw()); }

template<class PR> class FloatFactory {
    PR _pr;
  public:
    typedef PR PrecisionType;
    FloatFactory(PR const& pr) : _pr(pr) { }
    PR precision() const { return this->_pr; }
  public:
    FloatApproximation<PR> create(Number<ApproximateTag> const& y) { return FloatApproximation<PR>(y,_pr); }
    FloatLowerBound<PR> create(Number<ValidatedLowerTag> const& y) { return FloatLowerBound<PR>(y,_pr); }
    FloatUpperBound<PR> create(Number<ValidatedUpperTag> const& y) { return FloatUpperBound<PR>(y,_pr); }
    FloatBounds<PR> create(Number<ValidatedTag> const& y) { return FloatBounds<PR>(y,_pr); }
    FloatBounds<PR> create(Number<EffectiveTag> const& y) { return FloatBounds<PR>(y,_pr); }
    FloatBounds<PR> create(Number<ExactTag> const& y) { return FloatBounds<PR>(y,_pr); }
    FloatBounds<PR> create(Real const& y) { return FloatBounds<PR>(y,_pr); }
    FloatBounds<PR> create(Rational const& y) { return FloatBounds<PR>(y,_pr); }
    FloatBounds<PR> create(Decimal const& y) { return FloatBounds<PR>(y,_pr); }
    FloatValue<PR> create(Dyadic const& y) { return FloatValue<PR>(y,_pr); }
    FloatValue<PR> create(Integer const& y) { return FloatValue<PR>(y,_pr); }
    template<class N, EnableIf<IsSignedIntegral<N>> =dummy> FloatValue<PR> create(N const& y) { return FloatValue<PR>(y,_pr); }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy> PositiveFloatValue<PR> create(M const& y) { return PositiveFloatValue<PR>(y,_pr); }
    template<class D, EnableIf<IsFloatingPoint<D>> =dummy> FloatApproximation<PR> create(D const& y) { return FloatApproximation<PR>(RawFloat<PR>(y,_pr)); }
};
template<class Y, class PR> using ConcreteType = decltype(declval<FloatFactory<PR>>().create(declval<Y>()));

template<class PR> inline FloatFactory<PR> float_factory(PR pr) { return FloatFactory<PR>(pr); }
template<template<class>class FLT, class PR> inline FloatFactory<PR> factory(FLT<PR> flt) { return FloatFactory<PR>(flt.precision()); }

template<class Y, class PR> inline decltype(auto) make_float(Y const& y, PR pr) { return float_factory(pr).create(y); }

template<class X, class Y> using AreConcreteGenericNumbers = And<IsFloat<X>,IsGenericNumericType<Y>>;


template<class X, class Y, EnableIf<AreConcreteGenericNumbers<X,Y>> =dummy>
decltype(auto) operator+(X const& x, Y const& y) { return x+factory(x).create(y); }
template<class X, class Y, EnableIf<AreConcreteGenericNumbers<X,Y>> =dummy>
decltype(auto) operator-(X const& x, Y const& y) { return x-factory(x).create(y); }
template<class X, class Y, EnableIf<AreConcreteGenericNumbers<X,Y>> =dummy>
decltype(auto) operator*(X const& x, Y const& y) { return x*factory(x).create(y); }
template<class X, class Y, EnableIf<AreConcreteGenericNumbers<X,Y>> =dummy>
decltype(auto) operator/(X const& x, Y const& y) { return x/factory(x).create(y); }

template<class X, class Y, EnableIf<AreConcreteGenericNumbers<X,Y>> =dummy>
decltype(auto) operator+(Y const& y, X const& x) { return factory(x).create(y)+x; }
template<class X, class Y, EnableIf<AreConcreteGenericNumbers<X,Y>> =dummy>
decltype(auto) operator-(Y const& y, X const& x) { return factory(x).create(y)-x; }
template<class X, class Y, EnableIf<AreConcreteGenericNumbers<X,Y>> =dummy>
decltype(auto) operator*(Y const& y, X const& x) { return factory(x).create(y)*x; }
template<class X, class Y, EnableIf<AreConcreteGenericNumbers<X,Y>> =dummy>
decltype(auto) operator/(Y const& y, X const& x) { return factory(x).create(y)/x; }

template<class X, class Y, EnableIf<AreConcreteGenericNumbers<X,Y>> =dummy>
decltype(auto) operator==(X const& x, Y const& y) { return x==factory(x).create(y); }
template<class X, class Y, EnableIf<AreConcreteGenericNumbers<X,Y>> =dummy>
decltype(auto) operator!=(X const& x, Y const& y) { return x!=factory(x).create(y); }
template<class X, class Y, EnableIf<AreConcreteGenericNumbers<X,Y>> =dummy>
decltype(auto) operator<=(X const& x, Y const& y) { return x<=factory(x).create(y); }
template<class X, class Y, EnableIf<AreConcreteGenericNumbers<X,Y>> =dummy>
decltype(auto) operator>=(X const& x, Y const& y) { return x>=factory(x).create(y); }
template<class X, class Y, EnableIf<AreConcreteGenericNumbers<X,Y>> =dummy>
decltype(auto) operator< (X const& x, Y const& y) { return x< factory(x).create(y); }
template<class X, class Y, EnableIf<AreConcreteGenericNumbers<X,Y>> =dummy>
decltype(auto) operator> (X const& x, Y const& y) { return x> factory(x).create(y); }

template<class X, class Y, EnableIf<AreConcreteGenericNumbers<X,Y>> =dummy>
decltype(auto) operator==(Y const& y, X const& x) { return factory(x).create(y)==x; }
template<class X, class Y, EnableIf<AreConcreteGenericNumbers<X,Y>> =dummy>
decltype(auto) operator!=(Y const& y, X const& x) { return factory(x).create(y)!=x; }
template<class X, class Y, EnableIf<AreConcreteGenericNumbers<X,Y>> =dummy>
decltype(auto) operator<=(Y const& y, X const& x) { return factory(x).create(y)<=x; }
template<class X, class Y, EnableIf<AreConcreteGenericNumbers<X,Y>> =dummy>
decltype(auto) operator>=(Y const& y, X const& x) { return factory(x).create(y)>=x; }
template<class X, class Y, EnableIf<AreConcreteGenericNumbers<X,Y>> =dummy>
decltype(auto) operator< (Y const& y, X const& x) { return factory(x).create(y)< x; }
template<class X, class Y, EnableIf<AreConcreteGenericNumbers<X,Y>> =dummy>
decltype(auto) operator> (Y const& y, X const& x) { return factory(x).create(y)> x; }

template<class PR> Bool operator==(FloatValue<PR> const& x, Rational const& q) { return cmp(x.raw(),q)==Comparison::EQUAL; }
template<class PR> Bool operator!=(FloatValue<PR> const& x, Rational const& q) { return cmp(x.raw(),q)!=Comparison::EQUAL; }
template<class PR> Bool operator<=(FloatValue<PR> const& x, Rational const& q) { return cmp(x.raw(),q)!=Comparison::GREATER; }
template<class PR> Bool operator>=(FloatValue<PR> const& x, Rational const& q) { return cmp(x.raw(),q)!=Comparison::LESS; }
template<class PR> Bool operator< (FloatValue<PR> const& x, Rational const& q) { return cmp(x.raw(),q)==Comparison::LESS; }
template<class PR> Bool operator> (FloatValue<PR> const& x, Rational const& q) { return cmp(x.raw(),q)==Comparison::GREATER; }
template<class PR> Bool operator>=(Rational const& q, FloatValue<PR> const& x);


// FIXME: Should be able to use cmp directly in >=
template<class PR> Bool operator==(Rational const& q, FloatValue<PR> const& x ) { return cmp(x.raw(),q)==Comparison::EQUAL; }
template<class PR> Bool operator!=(Rational const& q, FloatValue<PR> const& x ) { return cmp(x.raw(),q)!=Comparison::EQUAL; }
template<class PR> Bool operator>=(Rational const& q, FloatValue<PR> const& x ) { return cmp(x.raw(),q)!=Comparison::GREATER; return q>=Rational(x); }
template<class PR> Bool operator<=(Rational const& q, FloatValue<PR> const& x ) { return cmp(x.raw(),q)!=Comparison::LESS; }
template<class PR> Bool operator< (Rational const& q, FloatValue<PR> const& x ) { return cmp(x.raw(),q)==Comparison::GREATER; }
template<class PR> Bool operator> (Rational const& q, FloatValue<PR> const& x ) { return cmp(x.raw(),q)==Comparison::LESS; }



Float64Value cast_exact(const Real& x);

inline Float64Value const& cast_exact(RawFloat64 const& x) { return reinterpret_cast<Float64Value const&>(x); }
inline Float64Value const& cast_exact(Float64Approximation const& x) { return reinterpret_cast<Float64Value const&>(x); }
inline Float64Value const& cast_exact(Float64Value const& x) { return reinterpret_cast<Float64Value const&>(x); }
inline Float64Value const& cast_exact(Float64Error const& x) { return reinterpret_cast<Float64Value const&>(x); }

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
