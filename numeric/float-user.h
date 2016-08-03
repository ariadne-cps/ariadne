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

template<class X, class NX=X> struct DeclareDirectedNumericOperators {
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

template<class X, class NX=X> class DeclareFloatOperations : DeclareNumericOperations<X,NX>, DeclareMixedOperators<X,Real,X,NX> {
    friend OutputStream& operator<<(OutputStream&, X const&);
    friend InputStream& operator>>(InputStream&, X&);
};

template<class X, class NX=X> class DeclareDirectedFloatOperations : DeclareDirectedNumericOperations<X,NX>, DeclareMixedOperators<X,Real,X,NX> {
    friend OutputStream& operator<<(OutputStream&, X const&);
    friend InputStream& operator>>(InputStream&, X&);
};


template<class X, class Y, class R=X> class ProvideMixedOperators : ProvideInplaceOperators<X,Y,R> {
    friend R operator+(X const& x1, Y const& y2) { return x1+make_float(y2,x1.precision()); }
    friend R operator-(X const& x1, Y const& y2) { return x1-make_float(y2,x1.precision()); }
    friend R operator*(X const& x1, Y const& y2) { return x1*make_float(y2,x1.precision()); }
    friend R operator/(X const& x1, Y const& y2) { return x1/make_float(y2,x1.precision()); }
    friend R operator+(Y const& y1, X const& x2) { return make_float(y1,x2.precision())+x2; }
    friend R operator-(Y const& y1, X const& x2) { return make_float(y1,x2.precision())-x2; }
    friend R operator*(Y const& y1, X const& x2) { return make_float(y1,x2.precision())*x2; }
    friend R operator/(Y const& y1, X const& x2) { return make_float(y1,x2.precision())/x2; }
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

template<class X, class Y, class NX, class NY=Y> struct ProvideMixedDirectedOperators
{
    friend X operator+(X const& x1, Y const& y2) { return x1+make_float(y2,x1.precision()); }
    friend X operator+(Y const& y1, X const& x2) { return make_float(y1,x2.precision())+x2; }
    friend X operator-(X const& x1, NY const& y2) { return x1-make_float(y2,x1.precision()); }
    friend X operator-(Y const& y1, NX const& x2) { return make_float(y1,x2.precision())-x2; }
    friend X& operator+=(X& x1, Y const& y2) { return x1=x1+y2; }
    friend X& operator-=(X& x1, NY const& y2) { return x1=x1-y2; }
};

template<class X, class C, class E=C> class DeclareComparisons {
    typedef decltype(not declval<E>()) NE;
    typedef decltype(not declval<C>()) NC;
    friend E operator==(X const&, X const&);
    friend NE operator!=(X const&, X const&);
    friend C operator<=(X const&, X const&);
    friend C operator>=(X const&, X const&);
    friend NC operator< (X const&, X const&);
    friend NC operator> (X const&, X const&);
};

template<class X, class NX, class C> class DeclareDirectedComparisons {
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

template<class X, class C, class E=C> struct ProvideComparisons {
    typedef decltype(not declval<E>()) NE;
    typedef decltype(not declval<C>()) NC;
    friend E eq(X const&, X const&);
    friend C leq(X const&, X const&);
    friend C geq(X const& x1, X const& x2) { return leq(x2,x1); }
    friend E operator==(X const& x1, X const& x2){ return eq(x1,x2); }
    friend NE operator!=(X const& x1, X const& x2){ return not eq(x1,x2); }
    friend C operator<=(X const& x1, X const& x2){ return leq(x1,x2); }
    friend C operator>=(X const& x1, X const& x2){ return geq(x1,x2); }
    friend NC operator< (X const& x1, X const& x2){ return not geq(x1,x2); }
    friend NC operator> (X const& x1, X const& x2){ return not leq(x1,x2); }
};

template<class X1, class X2, class C1, class C2=C1, class E=C1> struct ProvideMixedComparisons {
    typedef decltype(not declval<E>()) NE;
    typedef decltype(not declval<C1>()) NC1;
    typedef decltype(not declval<C2>()) NC2;
    friend E eq(X1 const&, X2 const&);
    friend C1 leq(X1 const&, X2 const&);
    friend C2 geq(X1 const& x1, X2 const& x2) { return leq(x2,x1); }
    friend E operator==(X1 const& x1, X2 const& x2){ return eq(x1,x2); }
    friend NE operator!=(X1 const& x1, X2 const& x2){ return not eq(x1,x2); }
    friend C1 operator<=(X1 const& x1, X2 const& x2){ return leq(x1,x2); }
    friend C2 operator>=(X1 const& x1, X2 const& x2){ return geq(x1,x2); }
    friend NC2 operator< (X1 const& x1, X2 const& x2){ return not geq(x1,x2); }
    friend NC1 operator> (X1 const& x1, X2 const& x2){ return not leq(x1,x2); }
    friend E operator==(X2 const& x2, X1 const& x1){ return eq(x1,x2); }
    friend NE operator!=(X2 const& x2, X1 const& x1){ return not eq(x1,x2); }
    friend C2 operator<=(X2 const& x2, X1 const& x1){ return geq(x1,x2); }
    friend C1 operator>=(X2 const& x2, X1 const& x1){ return leq(x1,x2); }
    friend NC1 operator< (X2 const& x2, X1 const& x1){ return not leq(x1,x2); }
    friend NC2 operator> (X2 const& x2, X1 const& x1){ return not geq(x1,x2); }
};

template<class X, class NX, class C> struct ProvideDirectedComparisons {
    typedef decltype(not declval<C>()) NC;
    friend C leq(X const&, NX const&);
    friend C operator<=(X const& x1, NX const& x2) { return leq(x1,x2); }
    friend C operator>=(NX const& x1, X const& x2) { return leq(x2,x1); }
    friend NC operator< (NX const& x1, X const& x2) { return not leq(x2,x1); }
    friend NC operator> (X const& x1, NX const& x2) { return not leq(x1,x2); }
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
    : public DeclareFloatOperations<FloatApproximation<PR>>
    , public DeclareComparisons<FloatApproximation<PR>,ApproximateKleenean>
    , public DispatchNumericOperations<FloatApproximation<PR>>
    , public ProvideMixedOperators<FloatApproximation<PR>,Real>
    , public ProvideMixedOperators<FloatApproximation<PR>,double>
    , public ProvideComparisons<FloatApproximation<PR>,ApproximateKleenean>
    , public DeclareFieldComparisons<FloatApproximation<PR>,double,ApproximateKleenean>
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
        FloatApproximation<PR>(const Integer& z, PR pr);
        FloatApproximation<PR>(const Dyadic& w, PR pr);
        FloatApproximation<PR>(const Rational& q, PR pr);
        FloatApproximation<PR>(const Real& r, PR pr);
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
    : public DeclareDirectedFloatOperations<FloatLowerBound<PR>,FloatUpperBound<PR>>
    , public DeclareDirectedComparisons<FloatUpperBound<PR>,FloatLowerBound<PR>,ValidatedSierpinskian>
    , public DispatchDirectedNumericOperations<FloatLowerBound<PR>,FloatUpperBound<PR>>
    , public ProvideMixedDirectedOperators<FloatLowerBound<PR>,Real,FloatUpperBound<PR>>
    , public ProvideDirectedComparisons<FloatUpperBound<PR>,FloatLowerBound<PR>,ValidatedSierpinskian>
    , public ProvideDirectedComparisons<FloatLowerBound<PR>,ValidatedUpperNumber,ValidatedNegatedSierpinskian>
    , public DeclareFloatOperations<FloatApproximation<PR>>
    , public DeclareComparisons<FloatApproximation<PR>,ApproximateKleenean>
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
    FloatLowerBound<PR>(const Rational& q, PR pr);
    FloatLowerBound<PR>(const Real& r, PR pr);
    FloatLowerBound<PR>(const ValidatedLowerNumber& y, PR pr);

    FloatLowerBound<PR>(FloatBounds<PR> const& x);
    FloatLowerBound<PR>(FloatBall<PR> const& x);
    FloatLowerBound<PR>(FloatValue<PR> const& x);

        FloatLowerBound<PR>& operator=(const FloatValue<PR>& x);
    FloatLowerBound<PR>& operator=(const ValidatedLowerNumber&);
    FloatLowerBound<PR> create(const ValidatedLowerNumber& y) const;
    FloatUpperBound<PR> create(const ValidatedUpperNumber& y) const;

    operator ValidatedLowerNumber () const;

    PrecisionType precision() const { return _l.precision(); }
    RawFloatType const& raw() const { return _l; }
    RawFloatType& raw() { return _l; }
    double get_d() const { return _l.get_d(); }
  public: // To be removed
    friend FloatError<PR> rec(FloatLowerBound<PR> const&);
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
    : public DeclareDirectedFloatOperations<FloatUpperBound<PR>,FloatLowerBound<PR>>
    , public DeclareDirectedComparisons<FloatUpperBound<PR>,FloatLowerBound<PR>,ValidatedSierpinskian>
    , public DispatchDirectedNumericOperations<FloatUpperBound<PR>,FloatLowerBound<PR>>
    , public ProvideMixedDirectedOperators<FloatUpperBound<PR>,Real,FloatLowerBound<PR>>
    , public ProvideDirectedComparisons<FloatLowerBound<PR>,FloatUpperBound<PR>,ValidatedNegatedSierpinskian>
    , public ProvideDirectedComparisons<FloatUpperBound<PR>,ValidatedLowerNumber,ValidatedSierpinskian>
    , public DeclareFloatOperations<FloatApproximation<PR>>
    , public DeclareComparisons<FloatApproximation<PR>,ApproximateKleenean>
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
    FloatUpperBound<PR>(const Rational& q, PR pr);
    FloatUpperBound<PR>(const Real& r, PR pr);
    FloatUpperBound<PR>(const ValidatedUpperNumber& y, PR pr);

    FloatUpperBound<PR>(FloatBounds<PR> const& x);
    FloatUpperBound<PR>(FloatBall<PR> const& x);
    FloatUpperBound<PR>(FloatValue<PR> const& x);
    FloatUpperBound<PR>(FloatError<PR> const& x); // FIXME: Remove

        FloatUpperBound<PR>& operator=(const FloatValue<PR>& x);
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
    : public DeclareFloatOperations<FloatBounds<PR>>
    , public DispatchNumericOperations<FloatBounds<PR>>
    , public ProvideMixedOperators<FloatBounds<PR>,Real>
    , public ProvideComparisons<FloatBounds<PR>,ValidatedKleenean>
    , public ProvideMixedComparisons<FloatBounds<PR>,ValidatedNumber,ValidatedKleenean>
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
    FloatBounds<PR>(ValidatedLowerNumber const& lower, FloatLowerBound<PR> const& upper) : FloatBounds<PR>(upper.create(lower),upper) { }
    template<class N1, class N2, EnableIf<And<IsIntegral<N1>,IsIntegral<N2>>> = dummy> FloatBounds<PR>(N1 n1, N2 n2, PR pr) : _l(n1,pr), _u(n2,pr) { }
    FloatBounds<PR>(ExactDouble const& dl, ExactDouble const& du, PrecisionType pr);
    FloatBounds<PR>(Rational const& ql, Rational const& qu, PrecisionType pr);

    template<class N, EnableIf<IsIntegral<N>> = dummy> FloatBounds<PR>(N n, PR pr) : FloatBounds<PR>(ExactDouble(n),pr) { }
    FloatBounds<PR>(ExactDouble d, PR pr);
        FloatBounds<PR>(const Integer& z, PR pr);
        FloatBounds<PR>(const Dyadic& w, PR pr);
        FloatBounds<PR>(const Rational& q, PR pr);
        FloatBounds<PR>(const Real& x, PR pr);
    FloatBounds<PR>(const ValidatedNumber& y, PR pr);

    FloatBounds<PR>(FloatBall<PR> const& x);
    FloatBounds<PR>(FloatValue<PR> const& x);

        FloatBounds<PR>& operator=(const FloatValue<PR>& x) { return *this=FloatBounds<PR>(x); }
    FloatBounds<PR>& operator=(const ValidatedNumber& y);
    template<class N, EnableIf<IsIntegral<N>> = dummy> FloatBounds<PR>& operator=(N n) { return *this=ValidatedNumber(n); }

    operator ValidatedNumber () const;


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
    friend FloatApproximation<PR> round(FloatApproximation<PR> const& x);
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
    : public DeclareFloatOperations<FloatBall<PR>>
    , public DispatchNumericOperations<FloatBall<PR>>
    , public ProvideMixedOperators<FloatBall<PR>,Real>
    , public ProvideComparisons<FloatBall<PR>,ValidatedKleenean>
    , public ProvideMixedComparisons<FloatBall<PR>,ValidatedNumber,ValidatedKleenean>
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
        FloatBall<PR>(const Rational& q, PR pr);
        FloatBall<PR>(const Real& r, PR pr);
    FloatBall<PR>(const ValidatedNumber& y, PR pr);

    explicit FloatBall<PR>(FloatBounds<PR> const& x);
    FloatBall<PR>(FloatValue<PR> const& x);

    FloatBall<PR>& operator=(const ValidatedNumber& y);

    operator ValidatedNumber () const;

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
    : public DeclareFloatOperations<FloatBounds<PR>>
    , public DeclareMixedOperators<FloatValue<PR>, Real, FloatBounds<PR>>
    , public DeclareComparisons<FloatValue<PR>,Boolean>
    , public DispatchNumericOperations<FloatValue<PR>,FloatBounds<PR>>
    , public ProvideMixedOperators<FloatValue<PR>,Real,FloatBounds<PR>>
    , public ProvideComparisons<FloatValue<PR>,Boolean>
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

    template<class N, EnableIf<IsIntegral<N>> = dummy> FloatValue<PR>& operator=(N n) { _v=n; }
    FloatValue<PR>& operator=(const Integer& z);
    FloatValue<PR>& operator=(const TwoExp& t);
    FloatValue<PR>& operator=(const Dyadic& w);

    operator ExactNumber () const;
    explicit operator Dyadic () const;
    explicit operator Rational () const;

//    explicit operator RawFloatType () const { return _v; }

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
    friend FloatValue<PR>& operator*=(FloatValue<PR>&, TwoExp const&);
    friend FloatValue<PR>& operator/=(FloatValue<PR>&, TwoExp const&);
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
    : public DispatchPositiveDirectedNumericOperations<FloatError<PR>,FloatLowerBound<PR>>
    , public DeclareDirectedNumericOperations<FloatUpperBound<PR>,FloatLowerBound<PR>>
    , public DeclareDirectedComparisons<FloatUpperBound<PR>,FloatLowerBound<PR>,ValidatedSierpinskian>
{
  private: public:
    RawFloat<PR> _e;
  public:
    typedef PR PrecisionType;
  public:
    FloatError<PR>() : _e() { }
    explicit FloatError<PR>(PR pr) : _e(pr) { }
    explicit FloatError<PR>(RawFloat<PR> const& x) : _e(x) {
        ARIADNE_PRECONDITION_MSG(!(x<0),"x="<<x); }
    explicit FloatError<PR>(FloatUpperBound<PR> const& x) : FloatError<PR>(x.raw()) { }
    FloatError<PR>(PositiveFloatUpperBound<PR> const& x) : _e(x.raw()) { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy> FloatError<PR>(M m, PR pr) : _e(m,pr) { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy> FloatError<PR>& operator=(M m) { _e=m; return *this; }
    PrecisionType precision() const { return _e.precision(); }
    RawFloat<PR> const& raw() const { return _e; }
    RawFloat<PR>& raw() { return _e; }
  public:
    friend FloatError<PR> mag(FloatError<PR> const& x);
    friend FloatUpperBound<PR> operator+(FloatError<PR> const& x) { return FloatUpperBound<PR>(+x._e); }
    friend FloatLowerBound<PR> operator-(FloatError<PR> const& x) { return FloatLowerBound<PR>(-x._e); }
    friend FloatUpperBound<PR> log(FloatError<PR> const& x);
    friend FloatError<PR> exp(FloatError<PR> const& x);

    friend FloatError<PR> operator*(Nat m1, FloatError<PR> const& x2);
    friend FloatError<PR> operator*(FloatError<PR> const& x1, Nat m2);
    friend FloatError<PR>& operator*=(FloatError<PR>& x1, Nat m2);

    friend Bool same(FloatError<PR> const&, FloatError<PR> const&);
    friend Bool refines(FloatError<PR> const&, FloatError<PR> const&);
    friend OutputStream& operator<<(OutputStream&, FloatError<PR> const&);
    friend OutputStream& operator<<(OutputStream&, FloatError<PR> const&);
  public:
    friend FloatUpperBound<PR> operator+(Real const&, FloatUpperBound<PR> const&);
    friend FloatLowerBound<PR> operator-(Real const&, FloatUpperBound<PR> const&);
    friend FloatLowerBound<PR> operator-(Int, FloatUpperBound<PR> const&);
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
    PositiveFloatValue<PR>(Nat m, PR pr) : FloatValue<PR>(m,pr) { }
    template<class N, EnableIf<IsIntegral<N>> =dummy> PositiveFloatValue<PR>(N n, PR pr) : FloatValue<PR>(n,pr) { }
    PositiveFloatValue<PR>(TwoExp const& ex, PR pr) : FloatValue<PR>(ex,pr) { }
    explicit PositiveFloatValue<PR>(Dyadic const& w, PR pr) : FloatValue<PR>(w,pr) { }
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
    explicit PositiveFloatUpperBound<PR>(ValidatedUpperNumber const& y, PR pr) : FloatUpperBound<PR>(y,pr) { }
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
    explicit PositiveFloatLowerBound<PR>(ValidatedLowerNumber const& y, PR pr) : FloatLowerBound<PR>(y,pr) { }
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
    explicit PositiveFloatApproximation<PR>(ApproximateNumber const& y, PR pr) : FloatApproximation<PR>(y,pr) { }
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


template<class PR> inline FloatApproximation<PR> make_float(Number<ApproximateTag> const& y, PR pr) { return FloatApproximation<PR>(y,pr); }
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
