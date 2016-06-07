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

template<class X, class NX=X> struct DeclareNumericOperators {
    friend X operator+(X const& x);
    friend X operator-(NX const& x);
    friend NX operator-(X const& x);
    friend X operator+(X const& x1, X const& x2);
    friend X operator-(X const& x1, NX const& x2);
    friend NX operator-(NX const& x1, X const& x2);
    friend X& operator+=(X& x1, X const& x2);
    friend X& operator-=(X& x1, NX const& x2);

    friend X nul(X const& x);
    friend X pos(X const& x);
    friend X neg(NX const& x);
    friend NX neg(X const& x);
    friend X half(X const& x);
    friend X add(X const& x1, X const& x2);
    friend X sub(X const& x1, NX const& x2);
    friend NX sub(NX const& x1, X const& x2);

    friend X max(X const& x1, X const& x2);
    friend X min(X const& x1, X const& x2);
};

template<class X> struct DeclareNumericOperators<X> {
    friend X operator+(X const& x);
    friend X operator-(X const& x);
    friend X operator+(X const& x1, X const& x2);
    friend X operator-(X const& x1, X const& x2);
    friend X operator*(X const& x1, X const& x2);
    friend X operator/(X const& x1, X const& x2);
    friend X& operator+=(X& x1, X const& x2);
    friend X& operator-=(X& x1, X const& x2);
    friend X& operator*=(X& x1, X const& x2);
    friend X& operator/=(X& x1, X const& x2);

    friend X nul(X const& x);
    friend X pos(X const& x);
    friend X neg(X const& x);
    friend X half(X const& x);
    friend X sqr(X const& x);
    friend X rec(X const& x);

    friend X add(X const& x1, X const& x2);
    friend X sub(X const& x1, X const& x2);
    friend X mul(X const& x1, X const& x2);
    friend X div(X const& x1, X const& x2);
    friend X fma(X const& x1, X const& x2, X const& y);
    friend X pow(X const& x, Nat m);
    friend X pow(X const& x, Int n);

    friend X sqrt(X const& x);
    friend X exp(X const& x);
    friend X log(X const& x);
    friend X sin(X const& x);
    friend X cos(X const& x);
    friend X tan(X const& x);
    friend X asin(X const& x);
    friend X acos(X const& x);
    friend X atan(X const& x);

    friend X max(X const& x1, X const& x2);
    friend X min(X const& x1, X const& x2);
    friend X abs(X const& x);

    friend X round(X const&);
};

template<class X, class Y, class R=X, class NR=R> class DeclareMixedOperators {
    friend R operator+(Y const& y1, X const& x2);
    friend R operator+(X const& x1, Y const& y2);
    friend NR operator-(Y const& y1, X const& x2);
    friend R operator-(X const& x1, Y const& y2);
    friend R operator*(Y const& y1, X const& x2);
    friend R operator*(X const& x1, Y const& y2);
    friend NR operator/(Y const& y1, X const& x2);
    friend R operator/(X const& x1, Y const& y2);
/*
    friend R add(Y const& y1, X const& x2);
    friend R add(X const& x1, Y const& y2);
    friend NR sub(Y const& y1, X const& x2);
    friend R sub(X const& x1, Y const& y2);
    friend R mul(Y const& y1, X const& x2);
    friend R mul(X const& x1, Y const& y2);
    friend NR div(Y const& y1, X const& x2);
    friend R div(X const& x1, Y const& y2);
*/
};

template<class X, class NX=X> class DeclareFloatOperators : DeclareNumericOperators<X,NX>, DeclareMixedOperators<X,Real,X,NX> {
    friend OutputStream& operator<<(OutputStream&, X const&);
    friend InputStream& operator>>(InputStream&, X&);
};

template<class X,class R> class ProvideInplaceOperators
{
};
template<class X> class ProvideInplaceOperators<X,X>
{
    friend X& operator+=(X& x1, X const& x2) { return x1=add(x1,x2); }
    friend X& operator-=(X& x1, X const& x2) { return x1=sub(x1,x2); }
    friend X& operator*=(X& x1, X const& x2) { return x1=mul(x1,x2); }
    friend X& operator/=(X& x1, X const& x2) { return x1=div(x1,x2); }
};
template<class X, class R=X> class ProvideOperators
    : ProvideInplaceOperators<X,R>
{
    friend X operator+(X const& x) { return pos(x); }
    friend X operator-(X const& x) { return neg(x); }
    friend R operator+(X const& x1, X const& x2) { return add(x1,x2); }
    friend R operator-(X const& x1, X const& x2) { return sub(x1,x2); }
    friend R operator*(X const& x1, X const& x2) { return mul(x1,x2); }
    friend R operator/(X const& x1, X const& x2) { return div(x1,x2); }
};

template<class X, class NX> class ProvideDirectedOperators
{
    friend NX operator-(X const&);
    friend X operator+(X const& x) { return pos(x); }
    friend X operator-(NX const& x) { return neg(x); }
    friend X operator+(X const& x1, X const& x2) { return add(x1,x2); }
    friend X operator-(X const& x1, NX const& x2) { return sub(x1,x2); }
    friend X& operator+=(X& x1, X const& x2) { return x1=add(x1,x2); }
    friend X& operator-=(X& x1, NX const& x2) { return x1=sub(x1,x2); }
};

template<class X, class Y, class R=X, class NR=R> class ProvideMixedOperators {
    friend R operator+(Y const& y1, X const& x2) { return make_float(y1,x2.precision())+x2; }
    friend R operator+(X const& x1, Y const& y2) { return x1+make_float(y2,x1.precision()); }
    friend NR operator-(Y const& y1, X const& x2) { return make_float(y1,x2.precision())-x2; }
    friend R operator-(X const& x1, Y const& y2) { return x1-make_float(y2,x1.precision()); }
    friend R operator*(Y const& y1, X const& x2) { return make_float(y1,x2.precision())*x2; }
    friend R operator*(X const& x1, Y const& y2) { return x1*make_float(y2,x1.precision()); }
    friend NR operator/(Y const& y1, X const& x2) { return make_float(y1,x2.precision())/x2; }
    friend R operator/(X const& x1, Y const& y2) { return x1/make_float(y2,x1.precision()); }
};

template<class PR> inline FloatBall<PR> make_float_ball(Real const& y, PR pr) { return FloatBall<PR>(make_float(y,pr)); }
template<class PR> class ProvideMixedOperators<FloatBall<PR>,Real> {
    typedef FloatBall<PR> X; typedef Real Y; typedef X R; typedef R NR;
//friend FloatBall<PR> make_float_ball(Real const& y, PR pr) { return FloatBall<PR>(make_float(y,pr)); }
    friend R operator+(Y const& y1, X const& x2) { return make_float_ball(y1,x2.precision())+x2; }
    friend R operator+(X const& x1, Y const& y2) { return x1+make_float_ball(y2,x1.precision()); }
    friend NR operator-(Y const& y1, X const& x2) { return make_float_ball(y1,x2.precision())-x2; }
    friend R operator-(X const& x1, Y const& y2) { return x1-make_float_ball(y2,x1.precision()); }
    friend R operator*(Y const& y1, X const& x2) { return make_float_ball(y1,x2.precision())*x2; }
    friend R operator*(X const& x1, Y const& y2) { return x1*make_float_ball(y2,x1.precision()); }
    friend NR operator/(Y const& y1, X const& x2) { return make_float_ball(y1,x2.precision())/x2; }
    friend R operator/(X const& x1, Y const& y2) { return x1/make_float_ball(y2,x1.precision()); }
};

template<class X, class Y, class NX, class NY=Y> class ProvideMixedDirectedOperators
{
    friend X operator+(X const& x1, Y const& y2) { return x1+make_float(y2,x1.precision()); }
    friend X operator+(Y const& y1, X const& x2) { return make_float(y1,x2.precision())+x2; }
    friend X operator-(X const& x1, NY const& y2) { return x1-make_float(y2,x1.precision()); }
    friend X operator-(Y const& y1, NX const& x2) { return make_float(y1,x2.precision())-x2; }
    friend X& operator+=(X& x1, Y const& y2) { return x1=x1+y2; }
    friend X& operator-=(X& x1, NY const& y2) { return x1=x1-y2; }
};

template<class X, class C, class A=C> class DeclareFloatComparisons {
    typedef decltype(not declval<A>()) E;
    friend E operator==(X const&, X const&);
    friend A operator!=(X const&, X const&);
    friend C operator<=(X const&, X const&);
    friend C operator>=(X const&, X const&);
    friend C operator< (X const&, X const&);
    friend C operator> (X const&, X const&);
};

template<class X, class NX, class C> class DeclareDirectedFloatComparisons {
    typedef decltype(not declval<C>()) NC;
    friend C operator<=(X const&, NX const&);
    friend NC operator<=(NX const&, X const&);
    friend NC operator>=(X const&, NX const&);
    friend C operator>=(NX const&, X const&);
    friend C operator< (X const&, NX const&);
    friend NC operator< (NX const&, X const&);
    friend NC operator> (X const&, NX const&);
    friend C operator> (NX const&, X const&);
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
    : public DeclareFloatOperators<FloatApproximation<PR>>
    , public DeclareFloatComparisons<FloatApproximation<PR>,ApproximateKleenean>
    , public ProvideOperators<FloatApproximation<PR>>
    , public ProvideMixedOperators<FloatApproximation<PR>,Real>
{
    typedef ApproximateTag P; typedef RawFloat<PR> FLT;
  public:
    typedef ApproximateTag Paradigm;
    typedef FloatApproximation<PR> NumericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
  public:
    FloatApproximation<PR>() : _a(0.0) { }
    FloatApproximation<PR>(PrecisionType pr) : _a(0.0,pr) { }
    explicit FloatApproximation<PR>(RawFloatType const& a) : _a(a) { }
    template<class N, EnableIf<IsIntegral<N>> =dummy> FloatApproximation<PR>(N n) : _a(n) { }
    template<class D, EnableIf<IsFloatingPoint<D>> =dummy> FloatApproximation<PR>(D x) : _a(x) { }
    explicit FloatApproximation<PR>(const Integer& z);
    explicit FloatApproximation<PR>(const Dyadic& d);
    explicit FloatApproximation<PR>(const Decimal& d);
    explicit FloatApproximation<PR>(const Rational& q);
    explicit FloatApproximation<PR>(const Real& r);
    explicit FloatApproximation<PR>(const Number<ApproximateTag>& x);
    FloatApproximation<PR>(const Rational& q, PR pr);
    FloatApproximation<PR>(const Real& r, PR pr);
    FloatApproximation<PR>(const Number<ApproximateTag>& x, PR pr);
    operator Number<ApproximateTag> () const;

    FloatApproximation<PR>(FloatError<PR> const& x); // FIXME: Remove
    FloatApproximation<PR>(FloatValue<PR> const& x);
    FloatApproximation<PR>(FloatBall<PR> const& x);
    FloatApproximation<PR>(FloatBounds<PR> const& x);
    FloatApproximation<PR>(FloatUpperBound<PR> const& x);
    FloatApproximation<PR>(FloatLowerBound<PR> const& x);

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
    : public DeclareFloatOperators<FloatLowerBound<PR>,FloatUpperBound<PR>>
    , public DeclareFloatOperators<FloatApproximation<PR>>
    , public DeclareFloatComparisons<FloatApproximation<PR>,ApproximateKleenean>
    , public DeclareDirectedFloatComparisons<FloatUpperBound<PR>,FloatLowerBound<PR>,ValidatedSierpinskian>
    , public ProvideDirectedOperators<FloatLowerBound<PR>,FloatUpperBound<PR>>
    , public ProvideMixedDirectedOperators<FloatLowerBound<PR>,Real,FloatUpperBound<PR>>
{
    typedef LowerTag P; typedef RawFloat<PR> FLT;
  public:
    typedef LowerTag Paradigm;
    typedef FloatLowerBound<PR> NumericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
  public:
    FloatLowerBound<PR>() : _l(0.0) { }
    FloatLowerBound<PR>(PrecisionType pr) : _l(0.0,pr) { }
    explicit FloatLowerBound<PR>(RawFloatType const& l) : _l(l) { }

    template<class N, EnableIf<IsIntegral<N>> = dummy> FloatLowerBound<PR>(N n) : _l(n) { }
    template<class X, EnableIf<IsFloatingPoint<X>> = dummy> explicit FloatLowerBound<PR>(X x) : _l(x) { }

    FloatLowerBound<PR>(FloatBounds<PR> const& x);
    FloatLowerBound<PR>(FloatBall<PR> const& x);
    FloatLowerBound<PR>(FloatValue<PR> const& x);

    FloatLowerBound<PR>(const Rational& q, PR pr);
    FloatLowerBound<PR>(const Real& r, PR pr);
    FloatLowerBound<PR>(const Number<ValidatedLowerTag>& x, PR pr);
    operator Number<ValidatedLowerTag> () const;

    explicit FloatLowerBound<PR>(const Integer& x);
    explicit FloatLowerBound<PR>(const Rational& x);
    explicit FloatLowerBound<PR>(const Real& x);
    explicit FloatLowerBound<PR>(const Number<ValidatedLowerTag>& x);

    PrecisionType precision() const { return _l.precision(); }
    RawFloatType const& raw() const { return _l; }
    RawFloatType& raw() { return _l; }
    double get_d() const { return _l.get_d(); }
  public: // To be removed
    friend FloatUpperBound<PR> rec(FloatLowerBound<PR> const&);
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
    : public DeclareFloatOperators<FloatUpperBound<PR>,FloatLowerBound<PR>>
    , public DeclareFloatOperators<FloatApproximation<PR>>
    , public DeclareFloatComparisons<FloatApproximation<PR>,ApproximateKleenean>
    , public DeclareDirectedFloatComparisons<FloatUpperBound<PR>,FloatLowerBound<PR>,ValidatedSierpinskian>
    , public ProvideDirectedOperators<FloatUpperBound<PR>,FloatLowerBound<PR>>
    , public ProvideMixedDirectedOperators<FloatUpperBound<PR>,Real,FloatLowerBound<PR>>
{
    typedef UpperTag P; typedef RawFloat<PR> FLT;
  public:
    typedef UpperTag Paradigm;
    typedef FloatUpperBound<PR> NumericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
  public:
    FloatUpperBound<PR>() : _u(0.0) { }
    FloatUpperBound<PR>(PrecisionType pr) : _u(0.0,pr) { }
    explicit FloatUpperBound<PR>(RawFloatType const& u) : _u(u) { }

    template<class N, EnableIf<IsIntegral<N>> = dummy> FloatUpperBound<PR>(N n) : _u(n) { }
    template<class X, EnableIf<IsFloatingPoint<X>> = dummy> explicit FloatUpperBound<PR>(X x) : _u(x) { }

    FloatUpperBound<PR>(FloatBounds<PR> const& x);
    FloatUpperBound<PR>(FloatBall<PR> const& x);
    FloatUpperBound<PR>(FloatValue<PR> const& x);
    FloatUpperBound<PR>(FloatError<PR> const& x); // FIXME: Remove

    explicit FloatUpperBound<PR>(const Integer& x);
    explicit FloatUpperBound<PR>(const Rational& x);
    explicit FloatUpperBound<PR>(const Real& x);
    explicit FloatUpperBound<PR>(const Number<ValidatedUpperTag>& x);

    FloatUpperBound<PR>(const Rational& q, PR pr);
    FloatUpperBound<PR>(const Real& r, PR pr);
    FloatUpperBound<PR>(const Number<ValidatedUpperTag>& x, PR pr);
    operator Number<ValidatedUpperTag> () const;

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
    : public DeclareFloatOperators<FloatBounds<PR>>
    , public DeclareFloatComparisons<FloatBounds<PR>,ValidatedKleenean>
    , public ProvideOperators<FloatBounds<PR>>
    , public ProvideMixedOperators<FloatBounds<PR>,Real>
{
    typedef BoundedTag P; typedef RawFloat<PR> FLT;
  public:
    typedef BoundedTag Paradigm;
    typedef FloatBounds<PR> NumericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
  public:
    FloatBounds<PR>() : _l(0.0), _u(0.0) { }
    FloatBounds<PR>(PrecisionType pr) : _l(0.0,pr), _u(0.0,pr) { }
    explicit FloatBounds<PR>(RawFloatType const& v) : _l(v), _u(v) { }
    FloatBounds<PR>(RawFloatType const& l, RawFloatType const& u) : _l(l), _u(u) { }
    FloatBounds<PR>(FloatLowerBound<PR> const& lower, FloatUpperBound<PR> const& upper) : _l(lower.raw()), _u(upper.raw()) { }
    template<class N1, class N2, EnableIf<And<IsIntegral<N1>,IsIntegral<N2>>> = dummy> FloatBounds<PR>(N1 n1, N2 n2) : _l(n1), _u(n2) { }
    FloatBounds<PR>(Rational const& ql, Rational const& qu, PrecisionType pr);

    template<class N, EnableIf<IsIntegral<N>> = dummy> FloatBounds<PR>(N n) : _l(n), _u(n) { }
    template<class X, EnableIf<IsFloatingPoint<X>> = dummy> explicit FloatBounds<PR>(X x) : _l(x), _u(x) { }

    FloatBounds<PR>(FloatBall<PR> const& x);
    FloatBounds<PR>(FloatValue<PR> const& x);

    explicit FloatBounds<PR>(const Dyadic& x);
    explicit FloatBounds<PR>(const Decimal& x);
    explicit FloatBounds<PR>(const Integer& z);
    explicit FloatBounds<PR>(const Rational& q);
    explicit FloatBounds<PR>(const Real& x);
    explicit FloatBounds<PR>(const Number<ValidatedTag>& x);
    FloatBounds<PR>(const Integer& z, PR pr);
    FloatBounds<PR>(const Rational& q, PR pr);
    FloatBounds<PR>(const Real& x, PR pr);
    FloatBounds<PR>(const Number<ValidatedTag>& x, PR pr);
    operator Number<ValidatedTag> () const;

    FloatLowerBound<PR> const lower() const { return FloatLowerBound<PR>(lower_raw()); }
    FloatUpperBound<PR> const upper() const { return FloatUpperBound<PR>(upper_raw()); }
    FloatValue<PR> const value() const;
    FloatError<PR> const error() const;

    RawFloatType const& lower_raw() const { return _l; }
    RawFloatType const& upper_raw() const { return _u; }
    RawFloatType const value_raw() const { return half(add_near(_l,_u)); }
    RawFloatType const error_raw() const { RawFloatType v=value_raw(); return max(sub_up(_u,v),sub_up(v,_l)); }
    double get_d() const { return value_raw().get_d(); }

    PrecisionType precision() const { ARIADNE_DEBUG_ASSERT(_l.precision()==_u.precision()); return _u.precision(); }

    // DEPRECATED
    explicit operator RawFloatType () const { return value_raw(); }
    friend FloatValue<PR> midpoint(FloatBounds<PR> const& x);
  public:
    friend FloatError<PR> mag(FloatBounds<PR> const&);
    friend FloatLowerBound<PR> mig(FloatBounds<PR> const&);
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
    : public DeclareFloatOperators<FloatBall<PR>>
    , public ProvideOperators<FloatBall<PR>>
    , public ProvideMixedOperators<FloatBall<PR>,Real>
{
    typedef MetricTag P; typedef RawFloat<PR> FLT;
  public:
    typedef MetricTag Paradigm;
    typedef FloatBall<PR> NumericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
  public:
    FloatBall<PR>() : _v(0.0), _e(0.0) { }
    FloatBall<PR>(PrecisionType pr) : _v(0.0,pr), _e(0.0,pr) { }
    explicit FloatBall<PR>(RawFloatType const& v) : _v(v), _e(0.0) { }
    FloatBall<PR>(RawFloatType const& v, RawFloatType const& e) : _v(v), _e(e) { }
    FloatBall<PR>(FloatValue<PR> const& value, FloatError<PR> const& error) : _v(value.raw()), _e(error.raw()) { }
    FloatBall<PR>(FloatLowerBound<PR> const& lower, FloatUpperBound<PR> const& upper) =  delete;

    explicit FloatBall<PR>(FloatBounds<PR> const& x);
    FloatBall<PR>(FloatValue<PR> const& x);

    template<class N, EnableIf<IsIntegral<N>> = dummy> FloatBall<PR>(N n) : _v(n), _e(nul(_v)) { }
    template<class X, EnableIf<IsFloatingPoint<X>> = dummy> explicit FloatBall<PR>(X x) : _v(x), _e(nul(_v)) { }
    explicit FloatBall<PR>(const Integer& z);
    explicit FloatBall<PR>(const Rational& q);
    explicit FloatBall<PR>(const Real& x);
    explicit FloatBall<PR>(const Number<ValidatedTag>& x);
    FloatBall<PR>(const Rational& q, PR pr);
    FloatBall<PR>(const Real& r, PR pr);
    FloatBall<PR>(const Number<ValidatedTag>& x, PR pr);
    operator Number<ValidatedTag> () const;

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
  public:
    friend FloatError<PR> mag(FloatBall<PR> const&);
    friend FloatLowerBound<PR> mig(FloatBall<PR> const&);
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
    : public DeclareFloatOperators<FloatBounds<PR>>
    , public DeclareMixedOperators<FloatValue<PR>, Real, FloatBounds<PR>>
    , public DeclareFloatComparisons<FloatValue<PR>,Boolean>
    , public ProvideOperators<FloatValue<PR>,FloatBounds<PR>>
    , public ProvideMixedOperators<FloatValue<PR>,Real,FloatBounds<PR>>
{
    typedef ExactTag P; typedef RawFloat<PR> FLT;
  public:
    typedef ExactTag Paradigm;
    typedef FloatValue<PR> NumericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
  public:
    FloatValue<PR>() : _v(0.0) { }
    FloatValue<PR>(PrecisionType pr) : _v(0.0,pr) { }
    explicit FloatValue<PR>(RawFloatType const& v) : _v(v) { }
    template<class N, EnableIf<IsIntegral<N>> =dummy> FloatValue<PR>(N n) : _v(n) { }
    template<class X, EnableIf<IsFloatingPoint<X>> =dummy> explicit FloatValue<PR>(X x) : _v(x) { }


    explicit FloatValue<PR>(const Integer& z);
    explicit FloatValue<PR>(const Integer& z, PR pr);
    explicit FloatValue<PR>(const TwoExp& ex, PR pr);
    explicit operator Rational () const;
    operator Number<ExactTag> () const;
    explicit operator RawFloatType () const { return _v; }

    PrecisionType precision() const { return _v.precision(); }
    RawFloatType const& raw() const { return _v; }
    RawFloatType& raw() { return _v; }
    double get_d() const { return _v.get_d(); }

    FloatBall<PR> pm(FloatError<PR> _e) const;
  public:
    friend FloatValue<PR> pos(FloatValue<PR> const&);
    friend FloatValue<PR> neg(FloatValue<PR> const&);
    friend FloatValue<PR> half(FloatValue<PR> const&);
    friend FloatValue<PR> operator+(FloatValue<PR> const&);
    friend FloatValue<PR> operator-(FloatValue<PR> const&);
    friend FloatValue<PR>& operator*=(FloatValue<PR>&, TwoExp);
    friend FloatValue<PR>& operator/=(FloatValue<PR>&, TwoExp);
    friend FloatValue max(FloatValue<PR> const&, FloatValue<PR> const&);
    friend FloatValue min(FloatValue<PR> const&, FloatValue<PR> const&);
    friend FloatValue abs(FloatValue<PR> const&);
    friend FloatError<PR> mag(FloatValue<PR> const&);
    friend FloatLowerBound<PR> mig(FloatValue<PR> const&);
    friend Bool same(FloatValue<PR> const&, FloatValue<PR> const&);
    friend OutputStream& operator<<(OutputStream&, FloatValue<PR> const&);
  public:
    static Nat output_places;
    static Void set_output_places(Nat p) { output_places=p; }
  private: public:
    RawFloatType _v;
};

template<class PR> class FloatError
    : public DeclareDirectedFloatComparisons<FloatUpperBound<PR>,FloatLowerBound<PR>,ValidatedSierpinskian>
{
  private: public:
    RawFloat<PR> _e;
  public:
    typedef PR PrecisionType;
  public:
    FloatError<PR>() : _e() { }
    explicit FloatError<PR>(RawFloat<PR> const& x) : _e(x) {
        ARIADNE_PRECONDITION_MSG(!(x<0),"x="<<x); }
    explicit FloatError<PR>(FloatUpperBound<PR> const& x) : FloatError<PR>(x.raw()) { }
    FloatError<PR>(PositiveFloatUpperBound<PR> const& x) : _e(x.raw()) { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy> FloatError<PR>(M m) : _e(m) { }
    template<class F, EnableIf<IsSame<F,FloatUpperBound<PR>>> =dummy>
        explicit FloatError<PR>(F const& x) : _e(x) { }
//    operator PositiveFloatUpperBound<PR>() const;
    PrecisionType precision() const { return _e.precision(); }
    RawFloat<PR> const& raw() const { return _e; }
    RawFloat<PR>& raw() { return _e; }
  public:
    friend FloatError<PR> mag(FloatError<PR> const&);
    friend FloatError<PR> max(FloatError<PR> const&, FloatError<PR> const&);
    friend FloatError<PR> min(FloatError<PR> const&, FloatError<PR> const&);
    friend FloatUpperBound<PR> operator+(FloatError<PR> const&);
    friend FloatLowerBound<PR> operator-(FloatError<PR> const&);
    friend FloatError<PR> operator+(FloatError<PR> const&, FloatError<PR> const&);
    friend FloatError<PR> operator*(FloatError<PR> const&, FloatError<PR> const&);
    friend FloatError<PR> operator*(Nat, FloatError<PR> const&);
    friend FloatError<PR> operator/(FloatError<PR> const&, PositiveFloatLowerBound<PR> const&);
    friend FloatError<PR>& operator+=(FloatError<PR>&, FloatError<PR> const&);
    friend FloatError<PR>& operator*=(FloatError<PR>&, FloatError<PR> const&);
    friend FloatError<PR> sqr(FloatError<PR> const&);
    friend FloatError<PR> pow(FloatError<PR> const&, Nat);
    friend FloatError<PR> abs(FloatError<PR> const&);
    friend OutputStream& operator<<(OutputStream&, FloatError<PR> const&);

    friend Bool same(FloatError<PR> const&, FloatError<PR> const&);

    friend FloatUpperBound<PR> operator+(FloatUpperBound<PR> const&, FloatUpperBound<PR> const&);
    friend FloatLowerBound<PR> operator-(FloatLowerBound<PR> const&, FloatUpperBound<PR> const&);
    friend FloatUpperBound<PR> operator+(Real const&, FloatUpperBound<PR> const&);
    friend FloatLowerBound<PR> operator-(Real const&, FloatUpperBound<PR> const&);
    friend FloatLowerBound<PR> operator-(Int, FloatUpperBound<PR> const&);
    friend FloatUpperBound<PR> log(FloatError<PR> const&);
    friend FloatLowerBound<PR> rec(FloatError<PR> const&);

    friend ValidatedSierpinskian operator<=(FloatUpperBound<PR> const&, FloatLowerBound<PR> const&);
    friend ValidatedNegatedSierpinskian operator> (FloatUpperBound<PR> const&, FloatLowerBound<PR> const&);
    friend ValidatedNegatedSierpinskian operator> (FloatUpperBound<PR> const&, double); // FIXME: Remove
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



template<class PR> class PositiveFloatValue : public FloatValue<PR> {
  public:
    PositiveFloatValue<PR>() : FloatValue<PR>() { }
    PositiveFloatValue<PR>(TwoExp ex) : FloatValue<PR>(ex) { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy>
        PositiveFloatValue<PR>(M m) : FloatValue<PR>(m) { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy>
        PositiveFloatValue<PR>(M m, PR pr) : FloatValue<PR>(m,pr) { }
    explicit PositiveFloatValue<PR>(RawFloat<PR> const& x) : FloatValue<PR>(x) { }
    explicit PositiveFloatValue<PR>(FloatValue<PR> const& x) : FloatValue<PR>(x) { }
  public:
    friend PositiveFloatValue<PR> max(PositiveFloatValue<PR> const&, PositiveFloatValue<PR> const&);
    friend PositiveFloatBounds<PR> operator*(PositiveFloatBounds<PR> const&, PositiveFloatBounds<PR> const&);
};


template<class PR> class PositiveFloatBounds : public FloatBounds<PR> {
  public:
    PositiveFloatBounds<PR>() : FloatBounds<PR>() { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy>
        PositiveFloatBounds<PR>(M m) : FloatBounds<PR>(m) { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy>
        PositiveFloatBounds<PR>(M m, PR pr) : FloatBounds<PR>(m,pr) { }
    explicit PositiveFloatBounds<PR>(RawFloat<PR> const& x) : FloatBounds<PR>(x) { }
    explicit PositiveFloatBounds<PR>(RawFloat<PR> const& l, RawFloat<PR> const& u) : FloatBounds<PR>(l,u) { }
    explicit PositiveFloatBounds<PR>(FloatBounds<PR> const& x) : FloatBounds<PR>(x) { }
  public:
    friend PositiveFloatBounds<PR> max(PositiveFloatBounds<PR> const&, PositiveFloatBounds<PR> const&);
    friend PositiveFloatBounds<PR> operator*(PositiveFloatBounds<PR> const&, PositiveFloatBounds<PR> const&);
};

template<class PR> class PositiveFloatUpperBound : public FloatUpperBound<PR> {
  public:
    PositiveFloatUpperBound<PR>() : FloatUpperBound<PR>() { }
    explicit PositiveFloatUpperBound<PR>(RawFloat<PR> const& x) : FloatUpperBound<PR>(x) {
        ARIADNE_PRECONDITION_MSG(!(x<0),"x="<<x); }
    explicit PositiveFloatUpperBound<PR>(FloatUpperBound<PR> const& x) : FloatUpperBound<PR>(x) { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy> PositiveFloatUpperBound<PR>(M m) : FloatUpperBound<PR>(m) { }
    template<class F, EnableIf<IsSame<F,FloatUpperBound<PR>>> =dummy>
        explicit PositiveFloatUpperBound<PR>(F const& x) : FloatUpperBound<PR>(x) { }
    PositiveFloatUpperBound<PR>(PositiveFloatValue<PR> const& x) : FloatUpperBound<PR>(x) { }
    PositiveFloatUpperBound<PR>(FloatError<PR> const& x) : FloatUpperBound<PR>(x.raw()) { }
  public:
    friend PositiveFloatUpperBound<PR> max(PositiveFloatUpperBound<PR> const&, PositiveFloatUpperBound<PR> const&);
    friend PositiveFloatUpperBound<PR> operator*(PositiveFloatUpperBound<PR> const&, PositiveFloatUpperBound<PR> const&);
    friend PositiveFloatUpperBound<PR> operator/(PositiveFloatUpperBound<PR> const&, PositiveFloatLowerBound<PR> const&);
    friend PositiveFloatUpperBound<PR>& operator*=(PositiveFloatUpperBound<PR>&, PositiveFloatUpperBound<PR> const&);
    friend PositiveFloatUpperBound<PR>& operator*=(PositiveFloatUpperBound<PR>&, PositiveFloatLowerBound<PR> const&);
  public:
    static Nat output_places;
    static Void set_output_places(Nat p) { output_places=p; }
};
template<class PR> PositiveFloatUpperBound<PR> abs(PositiveFloatUpperBound<PR> const&);

template<class PR> class PositiveFloatLowerBound : public FloatLowerBound<PR> {
  public:
    PositiveFloatLowerBound<PR>() : FloatLowerBound<PR>() { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy>
        PositiveFloatLowerBound<PR>(M m) : FloatLowerBound<PR>(m) { }
    explicit PositiveFloatLowerBound<PR>(RawFloat<PR> const& x) : FloatLowerBound<PR>(x) { }
    explicit PositiveFloatLowerBound<PR>(FloatLowerBound<PR> const& x) : FloatLowerBound<PR>(x) { }
    PositiveFloatLowerBound<PR>(PositiveFloatValue<PR> const& x) : FloatLowerBound<PR>(x) { }
  public:
    friend PositiveFloatLowerBound<PR> max(PositiveFloatLowerBound<PR> const&, PositiveFloatLowerBound<PR> const&);
    friend PositiveFloatLowerBound<PR> operator*(PositiveFloatLowerBound<PR> const&, PositiveFloatLowerBound<PR> const&);
};

template<class PR> class PositiveFloatApproximation : public FloatApproximation<PR> {
  public:
    PositiveFloatApproximation<PR>() : FloatApproximation<PR>() { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy>
        PositiveFloatApproximation<PR>(M m) : FloatApproximation<PR>(m) { }
    explicit PositiveFloatApproximation<PR>(RawFloat<PR> const& x) : FloatApproximation<PR>(x) { }
    explicit PositiveFloatApproximation<PR>(FloatApproximation<PR> const& x) : FloatApproximation<PR>(x) { }
    PositiveFloatApproximation<PR>(PositiveFloatLowerBound<PR> const& x) : FloatApproximation<PR>(x) { }
    PositiveFloatApproximation<PR>(PositiveFloatUpperBound<PR> const& x) : FloatApproximation<PR>(x) { }
    PositiveFloatApproximation<PR>(PositiveFloatValue<PR> const& x) : FloatApproximation<PR>(x) { }
    PositiveFloatApproximation<PR>(FloatError<PR> const& x) : FloatApproximation<PR>(x) { }
  public:
    friend PositiveFloatApproximation<PR> max(PositiveFloatApproximation<PR> const&, PositiveFloatApproximation<PR> const&);
    friend PositiveFloatApproximation<PR> operator*(PositiveFloatApproximation<PR> const&, PositiveFloatApproximation<PR> const&);
};

template<class PR> inline PositiveFloatApproximation<PR> cast_positive(FloatApproximation<PR> const& x) {
    return PositiveFloatApproximation<PR>(x); }

template<class PR> inline PositiveFloatLowerBound<PR> cast_positive(FloatLowerBound<PR> const& x) {
    return PositiveFloatLowerBound<PR>(x); }

template<class PR> inline PositiveFloatUpperBound<PR> cast_positive(FloatUpperBound<PR> const& x) {
    return PositiveFloatUpperBound<PR>(x); }

template<class PR> inline PositiveFloatBounds<PR> cast_positive(FloatBounds<PR> const& x) {
    return PositiveFloatBounds<PR>(x); }

template<class PR> inline PositiveFloatValue<PR> cast_positive(FloatValue<PR> const& x) {
    return PositiveFloatValue<PR>(x); }

template<class PR> inline OutputStream& operator<<(OutputStream& os, PositiveFloatApproximation<PR> const& x) {
    return os << static_cast<FloatApproximation<PR>const&>(x); }




template<class R, class A> R integer_cast(const A& _a);

template<> FloatBall<Precision64>::FloatBall(Real const& x);
template<> FloatBounds<Precision64>::FloatBounds(Real const& x);
template<> FloatUpperBound<Precision64>::FloatUpperBound(Real const& x);
template<> FloatLowerBound<Precision64>::FloatLowerBound(Real const& x);
template<> FloatApproximation<Precision64>::FloatApproximation(Real const& x);

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


extern const FloatValue<Precision64> infty;

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
template<class PR> FloatBounds<PR> make_bounds(FloatError<PR> const& e) {
    return FloatBounds<PR>(-e.raw(),+e.raw()); }


template<class PR> inline FloatApproximation<PR> make_float(Number<ApproximateTag> const& y, PR pr) { return FloatApproximation<PR>(y,pr); }
template<class PR> FloatApproximation<PR> make_float(Number<ApproximateTag> const& y, PR pr);
template<class PR> inline FloatLowerBound<PR> make_float(Number<ValidatedLowerTag> const& y, PR pr) { return FloatLowerBound<PR>(y,pr); }
template<class PR> inline FloatUpperBound<PR> make_float(Number<ValidatedUpperTag> const& y, PR pr) { return FloatUpperBound<PR>(y,pr); }
template<class PR> inline FloatBounds<PR> make_float(Number<ValidatedTag> const& y, PR pr) { return FloatBounds<PR>(y,pr); }
template<class PR> inline FloatBounds<PR> make_float(Number<EffectiveTag> const& y, PR pr) { return FloatBounds<PR>(y,pr); }
template<class PR> inline FloatBounds<PR> make_float(Number<ExactTag> const& y, PR pr) { return FloatBounds<PR>(y,pr); }
template<class PR> inline FloatBounds<PR> make_float(Real const& y, PR pr) { return FloatBounds<PR>(y,pr); }
template<class PR> inline FloatBounds<PR> make_float(Rational const& y, PR pr) { return FloatBounds<PR>(y,pr); }
template<class PR> inline FloatValue<PR> make_float(Integer const& y, PR pr) { return FloatValue<PR>(y,pr); }
template<class N, class PR, EnableIf<IsSignedIntegral<N>> =dummy> inline FloatValue<PR> make_float(N const& y, PR pr) { return FloatValue<PR>(y,pr); }
template<class M, class PR, EnableIf<IsUnsignedIntegral<M>> =dummy> inline PositiveFloatValue<PR> make_float(M const& y, PR pr) { return PositiveFloatValue<PR>(y,pr); }
template<class D, class PR, EnableIf<IsFloatingPoint<D>> =dummy> FloatApproximation<PR> make_float(D const& y, PR pr){
    return FloatApproximation<PR>(RawFloat<PR>(y,pr)); }

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

template<class PR> auto operator==(FloatValue<PR> const& x, Rational const& q) -> decltype(Rational(x)==q) { return Rational(x)==q; }
template<class PR> auto operator!=(FloatValue<PR> const& x, Rational const& q) -> decltype(Rational(x)!=q) { return Rational(x)!=q; }
template<class PR> auto operator<=(FloatValue<PR> const& x, Rational const& q) -> decltype(Rational(x)<=q) { return Rational(x)<=q; }
template<class PR> auto operator>=(FloatValue<PR> const& x, Rational const& q) -> decltype(Rational(x)>=q) { return Rational(x)>=q; }
template<class PR> auto operator< (FloatValue<PR> const& x, Rational const& q) -> decltype(Rational(x)< q) { return Rational(x)< q; }
template<class PR> auto operator> (FloatValue<PR> const& x, Rational const& q) -> decltype(Rational(x)> q) { return Rational(x)> q; }

template<class PR> auto operator==(Rational const& q, FloatValue<PR> const& x) -> decltype(q==Rational(x)) { return q==Rational(x); }
template<class PR> auto operator!=(Rational const& q, FloatValue<PR> const& x) -> decltype(q!=Rational(x)) { return q!=Rational(x); }
template<class PR> auto operator<=(Rational const& q, FloatValue<PR> const& x) -> decltype(q<=Rational(x)) { return q<=Rational(x); }
template<class PR> auto operator>=(Rational const& q, FloatValue<PR> const& x) -> decltype(q>=Rational(x)) { return q>=Rational(x); }
template<class PR> auto operator< (Rational const& q, FloatValue<PR> const& x) -> decltype(q< Rational(x)) { return q< Rational(x); }
template<class PR> auto operator> (Rational const& q, FloatValue<PR> const& x) -> decltype(q> Rational(x)) { return q> Rational(x); }

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
