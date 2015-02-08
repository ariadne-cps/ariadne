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
//! Operations can be specified to return an \c %ExactInterval answer by using the \c _ivl suffix.
//! The \c _approx suffix is provided to specifically indicate that the operation is computed approximately.
//!
//! %Ariadne floating-point numbers can be constructed by conversion from built-in C++ types.
//! Note that the value of a built-in floating-point value may differ from the mathematical value of the literal.
//! For example, while <c>%Float64(3.25)</c> is represented exactly, <c>%Float64(3.3)</c> has a value of \f$3.2999999999999998224\ldots\f$.
//! \note In the future, the construction of a \c %Float64 from a string literal may be supported.
//! \sa ExactInterval, Real, ExactFloat64
template<class PR> class Float<Approximate,PR> {
    typedef Approximate P; typedef RawFloatType<PR> FLT;
  public:
    typedef Approximate Paradigm;
    typedef ApproximateFloat64 NumericType;
    typedef PR PrecisionType;
  public:
    Float<Approximate,PR>() : a() { }
    template<class N, EnableIf<IsIntegral<N>> =dummy> Float<Approximate,PR>(N n) : a(n) { }
    template<class D, EnableIf<IsFloatingPoint<D>> =dummy> Float<Approximate,PR>(D x) : a(x) { }
    explicit Float<Approximate,PR>(Float64 const& x) : a(x) { }
    explicit Float<Approximate,PR>(const Dyadic& d);
    explicit Float<Approximate,PR>(const Decimal& d);

    explicit Float<Approximate,PR>(const Integer& z);
    explicit Float<Approximate,PR>(const Rational& q);
    explicit Float<Approximate,PR>(const Real& r);
    explicit Float<Approximate,PR>(const Number<Approximate>& x);
    Float<Approximate,PR>(const Number<Approximate>& x, PR pr);
    operator Number<Approximate> () const;

    Float<Approximate,PR>(Float<Exact,PR> const& x);
    Float<Approximate,PR>(Float<Validated,PR> const& x);
    Float<Approximate,PR>(Float<Upper,PR> const& x);
    Float<Approximate,PR>(Float<Lower,PR> const& x);

    PrecisionType precision() const { return a.precision(); }
    explicit operator Float64 () const { return this->a; }
    Float64 const& raw() const { return this->a; }
    Float64& raw() { return this->a; }
    double get_d() const { return this->a.get_d(); }
  public:
    static Void set_output_precision(Nat p) { output_precision=p; }
    Float<Approximate,PR> pm(Float<Approximate,PR> e) { return *this; }
  private: public:
    static Nat output_precision;
    Float64 a;
};


//! \ingroup NumericModule
//! \brief Floating-point lower bounds for real numbers.
template<class PR> class Float<Lower,PR> {
    typedef Lower P; typedef RawFloatType<PR> FLT;
  public:
    typedef Lower Paradigm;
    typedef Float<Lower,PR> NumericType;
    typedef PR PrecisionType;
  public:
    Float<Lower,PR>() : l(0.0) { }
    template<class N, EnableIf<IsIntegral<N>> = dummy> Float<Lower,PR>(N n) : l(n) { }
    template<class X, EnableIf<IsFloatingPoint<X>> = dummy> explicit Float<Lower,PR>(X x) : l(x) { }
    explicit Float<Lower,PR>(Float64 const& x) : l(x) { }

    Float<Lower,PR>(Float<Validated,PR> const& x);
    Float<Lower,PR>(Float<Exact,PR> const& x);

    explicit Float<Lower,PR>(const Number<Lower>& x);
    Float<Lower,PR>(const Number<Lower>& x, PR pr);
    operator Number<Lower> () const;

    explicit Float<Lower,PR>(const Real& x);
    explicit Float<Lower,PR>(const Rational& x);
    explicit Float<Lower,PR>(const Integer& x);

    PrecisionType precision() const { return l.precision(); }
    Float64 const& raw() const { return l; }
    Float64& raw() { return l; }
    double get_d() const { return l.get_d(); }
  private: public:
    static Nat output_precision;
    Float64 l;
};


//! \ingroup NumericModule
//! \brief Floating-point upper bounds for real numbers.
template<class PR> class Float<Upper,PR> {
    typedef Upper P; typedef RawFloatType<PR> FLT;
  public:
    typedef Upper Paradigm;
    typedef Float<Upper,PR> NumericType;
    typedef PR PrecisionType;
  public:
    Float<Upper,PR>() : u(0.0) { }
    template<class N, EnableIf<IsIntegral<N>> = dummy> Float<Upper,PR>(N n) : u(n) { }
    template<class X, EnableIf<IsFloatingPoint<X>> = dummy> explicit Float<Upper,PR>(X x) : u(x) { }

    explicit Float<Upper,PR>(Float64 const& x) : u(x) { }

    Float<Upper,PR>(Float<Validated,PR> const& x);
    Float<Upper,PR>(Float<Exact,PR> const& x);

    explicit Float<Upper,PR>(const Real& x);
    explicit Float<Upper,PR>(const Rational& x);
    explicit Float<Upper,PR>(const Integer& x);
    explicit Float<Upper,PR>(const Number<Upper>& x);

    Float<Upper,PR>(const Number<Upper>& x, PR pr);
    operator Number<Upper> () const;

    PrecisionType precision() const { return u.precision(); }
    Float64 const& raw() const { return u; }
    Float64& raw() { return u; }
    double get_d() const { return u.get_d(); }
  private: public:
    static Nat output_precision;
    Float64 u;
};



//! \ingroup NumericModule
//! \brief Validated bounds on a number with floating-point endpoints supporting outwardly-rounded arithmetic.
//! \details
//! Note that <c>%ValidatedFloat64(3.3)</c> yields the singleton interval \f$[3.2999999999999998224,3.2999999999999998224]\f$ (the constant is first interpreted by the C++ compiler to give a C++ \c double, whereas <c>%ValidatedFloat64("3.3")</c> yields the interval \f$[3.2999999999999998224,3.3000000000000002665]\f$ enclosing \f$3.3\f$.
//!
//! Comparison tests on \c ValidatedFloat64 use the idea that an interval represents a single number with an unknown value.
//! Hence the result is of type \c Tribool, which can take values { \c True, \c False, \c Indeterminate }.
//! Hence a test \f$[l_1,u_1]\leq [l_2,u_2]\f$ returns \c True if \f$u_1\leq u_2\f$, since in this case \f$x_1\leq x_2\f$ whenever \f$x_1\in[l_1,u_2]\f$ and \f$x_2\in[l_2,u_2]\f$, \c False if \f$l_1>u_2\f$, since in this case we know \f$x_1>x_2\f$, and \c Indeterminate otherwise, since in this case we can find \f$x_1,x_2\f$ making the result either true or false.
//! In the case of equality, the comparison \f$[l_1,u_1]\f$==\f$[l_2,u_2]\f$ only returns \c True if both intervals are singletons, since otherwise we can find values making the result either true of false.
//!
//! To obtain the lower and upper bounds of an interval, use \c ivl.lower() and \c ivl.upper().
//! To obtain the midpoint and radius, use \c ivl.midpoint() and \c ivl.radius().
//! Alternatives \c midpoint(ivl) and \c radius(ivl) are also provided.
//! Note that \c midpoint and \c radius return approximations to the true midpoint and radius of the interval. If \f$m\f$ and \f$r\f$ are the returned midpoint and radius of the interval \f$[l,u]\f$, the using exact arithmetic, we guarentee \f$m-r\leq l\f$ and \f$m+r\geq u\f$
//!
//! To test if an interval contains a point or another interval, use \c encloses(ValidatedFloat64,Float64) or \c encloses(ValidatedFloat64,ValidatedFloat64).
//! The test \c refines(ValidatedFloat64,ValidatedFloat64) can also be used.
//! \sa Float64
//!
//! \par Python interface
//!
//! In the Python interface, %Ariadne intervals can be constructed from Python literals of the form \c {a:b} or (deprecated) \c [a,b] .
//! The former is preferred, as it cannot be confused with literals for other classes such as Vector and Array types.
//! Automatic conversion is used to convert ValidatedFloat64 literals of the form \c {a,b} to an ValidatedFloat64 in functions.
//!
//! Care must be taken when defining intervals using floating-point coefficients, since values are first converted to the nearest
//! representable value by the Python interpreter. <br><br>
//! \code
//!   ValidatedFloat64({1.1:2.3}) # Create the interval [1.1000000000000001, 2.2999999999999998]
//!   ValidatedFloat64({2.5:4.25}) # Create the interval [2.5, 4.25], which can be represented exactly
//!   ValidatedFloat64([2.5,4.25]) # Alternative syntax for creating the interval [2.5, 4.25]
//! \endcode
template<class PR> class Float<Validated,PR> {
    typedef Validated P; typedef RawFloatType<PR> FLT;
  public:
    typedef Validated Paradigm;
    typedef Float<Validated,PR> NumericType;
    typedef PR PrecisionType;
  public:
    Float<Validated,PR>() : l(0.0), u(0.0) { }
    template<class N, EnableIf<IsIntegral<N>> = dummy> Float<Validated,PR>(N n) : l(n), u(n) { }
    template<class X, EnableIf<IsFloatingPoint<X>> = dummy> explicit Float<Validated,PR>(X x) : l(x), u(x) { }
    explicit Float<Validated,PR>(Float64 const& x) : l(x), u(x) { }

    Float<Validated,PR>(Float<Exact,PR> const& x);

    explicit Float<Validated,PR>(const Dyadic& x);
    explicit Float<Validated,PR>(const Decimal& x);
    explicit Float<Validated,PR>(const Integer& z);
    explicit Float<Validated,PR>(const Rational& q);
    explicit Float<Validated,PR>(const Real& x);
    explicit Float<Validated,PR>(const Number<Validated>& x);
    Float<Validated,PR>(const Number<Validated>& x, PR pr);
    operator Number<Validated> () const;

    template<class N1, class N2, EnableIf<And<IsIntegral<N1>,IsIntegral<N2>>> =dummy>
        Float<Validated,PR>(N1 lower, N2 upper) : l(lower), u(upper) { }
    Float<Validated,PR>(Float64 const& lower, Float64 const& upper) : l(lower), u(upper) { }
    Float<Validated,PR>(Float<Lower,PR> const& lower, Float<Upper,PR> const& upper) : l(lower.raw()), u(upper.raw()) { }
    Float<Validated,PR>(const Rational& lower, const Rational& upper);

    Float64 const& lower_raw() const { return l; }
    Float64 const& upper_raw() const { return u; }
    double get_d() const { return (l.get_d()+u.get_d())/2; }

    Float<Lower,PR> lower() const { return Float<Lower,PR>(l); }
    Float<Upper,PR> upper() const { return Float<Upper,PR>(u); }
    const Float<Exact,PR> value() const;
    const Float<PositiveUpper,PR> error() const;

    PrecisionType precision() const { ARIADNE_DEBUG_ASSERT(l.precision()==u.precision()); return u.precision(); }

    // DEPRECATED
    explicit operator Float64 () const { return (l+u)/2; }
    friend Float<Exact,PR> midpoint(Float<Validated,PR> const& x);
  public:
    static Nat output_precision;
    static Void set_output_precision(Nat p) { output_precision=p; }
  private: public:
    Float64 l, u;
};



//! \ingroup NumericModule
//! \related Float64, ValidatedFloat64
//! \brief A floating-point number, which is taken to represent the \em exact value of a real quantity.
template<class PR> class Float<Exact,PR> {
    typedef Exact P; typedef RawFloatType<PR> FLT;
  public:
    typedef Exact Paradigm;
    typedef Float<Exact,PR> NumericType;
    typedef PR PrecisionType;
  public:
    Float<Exact,PR>() : v(0) { }
    template<class N, EnableIf<IsIntegral<N>> =dummy> Float<Exact,PR>(N n) : v(n) { }
    template<class X, EnableIf<IsFloatingPoint<X>> =dummy> explicit Float<Exact,PR>(X x) : v(x) { }

    explicit Float<Exact,PR>(Float64 const& x) : v(x) { }

    explicit Float<Exact,PR>(const Integer& z);
    explicit Float<Exact,PR>(const Integer& z, PR pr);
    explicit operator Rational () const;
    operator Number<Exact> () const;
    explicit operator Float64 () const { return v; }

    PrecisionType precision() const { return v.precision(); }
    Float64 const& raw() const { return v; }
    Float64& raw() { return v; }
    double get_d() const { return v.get_d(); }

    Float<Validated,PR> pm(Float<Error,PR> e) const;
  public:
    static Nat output_precision;
    static Void set_output_precision(Nat p) { output_precision=p; }
  private: public:
    Float64 v;
};

template<class PR> inline const Float<Exact,PR> Float<Validated,PR>::value() const {
    return Float<Exact,PR>(med_near(this->l,this->u)); }

template<class PR> inline const Float<Error,PR> Float<Validated,PR>::error() const {
    Float64 v=med_near(this->l,this->u); return Float<Error,PR>(max(sub_up(this->u,v),sub_up(v,this->l))); }


template<class PR> class Float<PositiveExact,PR> : public Float<Exact,PR> {
  public:
    Float<PositiveExact,PR>() : Float<Exact,PR>() { }
    template<class M, EnableIf<IsIntegral<M>>, EnableIf<IsUnsigned<M>> =dummy>
        Float<PositiveExact,PR>(M m) : Float<Exact,PR>(m) { }
    explicit Float<PositiveExact,PR>(RawFloatType<PR> const& x) : Float<Exact,PR>(x) { }
    explicit Float<PositiveExact,PR>(Float<Exact,PR> const& x) : Float<Exact,PR>(x) { }
};

template<class PR> class Float<PositiveUpper,PR> : public Float<Upper,PR> {
  public:
    Float<PositiveUpper,PR>() : Float<Upper,PR>() { }
    explicit Float<PositiveUpper,PR>(RawFloatType<PR> const& x) : Float<Upper,PR>(x) { ARIADNE_PRECONDITION(x>=0); }
    explicit Float<PositiveUpper,PR>(Float<Upper,PR> const& x) : Float<Upper,PR>(x) { }
    template<class M, EnableIf<IsUnsigned<M>> =dummy> Float<PositiveUpper,PR>(M m) : Float<Upper,PR>(m) { }
    template<class F, EnableIf<IsSame<F,Float<Upper,PR>>> =dummy>
        explicit Float<PositiveUpper,PR>(F const& x) : Float<Upper,PR>(x) { }
    Float<PositiveUpper,PR>(Float<PositiveExact,PR> const& x) : Float<Upper,PR>(x) { }
};

template<class PR> class Float<PositiveLower,PR> : public Float<Lower,PR> {
  public:
    Float<PositiveLower,PR>() : Float<Lower,PR>() { }
    template<class M, EnableIf<IsSigned<M>> =dummy>
        Float<PositiveLower,PR>(M m) : Float<Lower,PR>(m) { }
    explicit Float<PositiveLower,PR>(RawFloatType<PR> const& x) : Float<Lower,PR>(x) { }
    explicit Float<PositiveLower,PR>(Float<Lower,PR> const& x) : Float<Lower,PR>(x) { }
    Float<PositiveLower,PR>(Float<PositiveExact,PR> const& x) : Float<Lower,PR>(x) { }
};

template<class PR> class Float<PositiveApproximate,PR> : public Float<Approximate,PR> {
  public:
    Float<PositiveApproximate,PR>() : Float<Approximate,PR>() { }
    template<class M, EnableIf<IsSigned<M>> =dummy>
        Float<PositiveApproximate,PR>(M m) : Float<Approximate,PR>(m) { }
    explicit Float<PositiveApproximate,PR>(RawFloatType<PR> const& x) : Float<Approximate,PR>(x) { }
    explicit Float<PositiveApproximate,PR>(Float<Approximate,PR> const& x) : Float<Approximate,PR>(x) { }
    Float<PositiveApproximate,PR>(Float<PositiveLower,PR> const& x) : Float<Approximate,PR>(x) { }
    Float<PositiveApproximate,PR>(Float<PositiveUpper,PR> const& x) : Float<Approximate,PR>(x) { }
    Float<PositiveApproximate,PR>(Float<PositiveExact,PR> const& x) : Float<Approximate,PR>(x) { }
};

template<class PR> inline Float<PositiveUpper,PR> cast_positive(Float<Upper,PR> const& x) {
    return Float<PositiveUpper,PR>(x); }

template<class PR> inline Float<PositiveExact,PR> cast_positive(Float<Exact,PR> const& x) {
    return Float<PositiveExact,PR>(x); }

template<class R, class A> R integer_cast(const A& a);

template<> Float<Validated,Precision64>::Float(Real const& x);
template<> Float<Upper,Precision64>::Float(Real const& x);
template<> Float<Lower,Precision64>::Float(Real const& x);
template<> Float<Approximate,Precision64>::Float(Real const& x);

inline ApproximateFloat64 floor(ApproximateFloat64 const& x) { return ApproximateFloat64(floor(x.a)); }
inline ApproximateFloat64 ceil(ApproximateFloat64 const& x) { return ApproximateFloat64(ceil(x.a)); }

inline ApproximateFloat64 abs(ApproximateFloat64 const& x) { return ApproximateFloat64(abs_exact(x.a)); }
inline ApproximateFloat64 max(ApproximateFloat64 const& x, ApproximateFloat64 y) { return ApproximateFloat64(max_exact(x.a,y.a)); }
inline ApproximateFloat64 min(ApproximateFloat64 const& x, ApproximateFloat64 y) { return ApproximateFloat64(min_exact(x.a,y.a)); }

inline ApproximateFloat64 nul(ApproximateFloat64 const& x) { return ApproximateFloat64(nul_exact(x.a)); }
inline ApproximateFloat64 pos(ApproximateFloat64 const& x) { return ApproximateFloat64(pos_exact(x.a)); }
inline ApproximateFloat64 neg(ApproximateFloat64 const& x) { return ApproximateFloat64(neg_exact(x.a)); }
inline ApproximateFloat64 half(ApproximateFloat64 const& x) { return ApproximateFloat64(half_exact(x.a)); }
inline ApproximateFloat64 sqr(ApproximateFloat64 const& x) { return ApproximateFloat64(mul_near(x.a,x.a)); }
inline ApproximateFloat64 rec(ApproximateFloat64 const& x) { return ApproximateFloat64(div_near(1.0,x.a)); }

inline ApproximateFloat64 add(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return ApproximateFloat64(add_near(x1.a,x2.a)); }
inline ApproximateFloat64 sub(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return ApproximateFloat64(sub_near(x1.a,x2.a)); }
inline ApproximateFloat64 mul(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return ApproximateFloat64(mul_near(x1.a,x2.a)); }
inline ApproximateFloat64 div(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return ApproximateFloat64(div_near(x1.a,x2.a)); }

inline ApproximateFloat64 pow(ApproximateFloat64 const& x, Nat m) { return ApproximateFloat64(pow_approx(x.a,m)); }
inline ApproximateFloat64 pow(ApproximateFloat64 const& x, Int n) { return ApproximateFloat64(pow_approx(x.a,n)); }

inline ApproximateFloat64 sqrt(ApproximateFloat64 const& x) { return ApproximateFloat64(sqrt_approx(x.a)); }
inline ApproximateFloat64 exp(ApproximateFloat64 const& x) { return ApproximateFloat64(exp_approx(x.a)); }
inline ApproximateFloat64 log(ApproximateFloat64 const& x) { return ApproximateFloat64(log_approx(x.a)); }
inline ApproximateFloat64 sin(ApproximateFloat64 const& x) { return ApproximateFloat64(sin_approx(x.a)); }
inline ApproximateFloat64 cos(ApproximateFloat64 const& x) { return ApproximateFloat64(cos_approx(x.a)); }
inline ApproximateFloat64 tan(ApproximateFloat64 const& x) { return ApproximateFloat64(tan_approx(x.a)); }
inline ApproximateFloat64 asin(ApproximateFloat64 const& x) { return ApproximateFloat64(asin_approx(x.a)); }
inline ApproximateFloat64 acos(ApproximateFloat64 const& x) { return ApproximateFloat64(acos_approx(x.a)); }
inline ApproximateFloat64 atan(ApproximateFloat64 const& x) { return ApproximateFloat64(atan_approx(x.a)); }

inline ApproximateFloat64 operator+(ApproximateFloat64 const& x) { return pos(x); }
inline ApproximateFloat64 operator-(ApproximateFloat64 const& x) { return neg(x); }
inline ApproximateFloat64 operator+(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return add(x1,x2); }
inline ApproximateFloat64 operator-(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return sub(x1,x2); }
inline ApproximateFloat64 operator*(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return mul(x1,x2); }
inline ApproximateFloat64 operator/(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return div(x1,x2); }
inline ApproximateFloat64& operator+=(ApproximateFloat64& x1, ApproximateFloat64 const& x2) { x1.a+=x2.a; return x1; }
inline ApproximateFloat64& operator-=(ApproximateFloat64& x1, ApproximateFloat64 const& x2) { x1.a-=x2.a; return x1; }
inline ApproximateFloat64& operator*=(ApproximateFloat64& x1, ApproximateFloat64 const& x2) { x1.a*=x2.a; return x1; }
inline ApproximateFloat64& operator/=(ApproximateFloat64& x1, ApproximateFloat64 const& x2) { x1.a/=x2.a; return x1; }

inline Bool operator==(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return x1.a==x2.a; }
inline Bool operator!=(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return x1.a!=x2.a; }
inline Bool operator<=(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return x1.a<=x2.a; }
inline Bool operator>=(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return x1.a>=x2.a; }
inline Bool operator< (ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return x1.a< x2.a; }
inline Bool operator> (ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return x1.a> x2.a; }

OutputStream& operator<<(OutputStream& os, ApproximateFloat64 const& x);
InputStream& operator>>(InputStream& is, ApproximateFloat64& x);



inline LowerFloat64 max(LowerFloat64 const& x1, LowerFloat64 const& x2) { return LowerFloat64(max_exact(x1.l,x2.l)); }
inline LowerFloat64 min(LowerFloat64 const& x1, LowerFloat64 const& x2) { return LowerFloat64(min_exact(x1.l,x2.l)); }

inline LowerFloat64 nul(LowerFloat64 const& x) { return LowerFloat64(pos_exact(x.l)); }
inline LowerFloat64 pos(LowerFloat64 const& x) { return LowerFloat64(pos_exact(x.l)); }
inline LowerFloat64 neg(UpperFloat64 const& x) { return LowerFloat64(neg_exact(x.u)); }
inline LowerFloat64 half(LowerFloat64 const& x) { return LowerFloat64(half_exact(x.l)); }
LowerFloat64 sqr(LowerFloat64 const& x);
LowerFloat64 rec(UpperFloat64 const& x);

inline LowerFloat64 add(LowerFloat64 const& x1, LowerFloat64 const& x2) { return LowerFloat64(add_down(x1.l,x2.l)); }
inline LowerFloat64 sub(LowerFloat64 const& x1, UpperFloat64 const& x2) { return LowerFloat64(sub_down(x1.l,x2.u)); }
LowerFloat64 mul(LowerFloat64 const& x1, LowerFloat64 const& x2);
LowerFloat64 div(LowerFloat64 const& x1, UpperFloat64 const& x2);
LowerFloat64 pow(LowerFloat64 const& x, Nat m);

LowerFloat64 sqrt(LowerFloat64 const& x);
LowerFloat64 exp(LowerFloat64 const& x);
LowerFloat64 log(LowerFloat64 const& x);
LowerFloat64 atan(LowerFloat64 const& x);

inline LowerFloat64 operator+(LowerFloat64 const& x) { return pos(x); }
inline LowerFloat64 operator-(UpperFloat64 const& x) { return neg(x); }
inline LowerFloat64 operator+(LowerFloat64 const& x1, LowerFloat64 const& x2) { return add(x1,x2); }
inline LowerFloat64 operator-(LowerFloat64 const& x1, UpperFloat64 const& x2) { return sub(x1,x2); }
inline LowerFloat64 operator*(LowerFloat64 const& x1, LowerFloat64 const& x2) { return mul(x1,x2); }
inline LowerFloat64 operator/(LowerFloat64 const& x1, UpperFloat64 const& x2) { return div(x1,x2); }
inline LowerFloat64& operator+=(LowerFloat64& x1, LowerFloat64 const& x2) { return x1=x1+x2; }
inline LowerFloat64& operator-=(LowerFloat64& x1, UpperFloat64 const& x2) { return x1=x1-x2; }
inline LowerFloat64& operator*=(LowerFloat64& x1, LowerFloat64 const& x2) { return x1=x1*x2; }
inline LowerFloat64& operator/=(LowerFloat64& x1, UpperFloat64 const& x2) { return x1=x1/x2; }

OutputStream& operator<<(OutputStream& os, LowerFloat64 const& x);
InputStream& operator>>(InputStream& is, LowerFloat64& x);

inline UpperFloat64 max(UpperFloat64 const& x1, UpperFloat64 const& x2) { return UpperFloat64(max_exact(x1.u,x2.u)); }
inline UpperFloat64 min(UpperFloat64 const& x1, UpperFloat64 const& x2) { return UpperFloat64(min_exact(x1.u,x2.u)); }

inline UpperFloat64 nul(UpperFloat64 const& x) { return UpperFloat64(pos_exact(x.u)); }
inline UpperFloat64 pos(UpperFloat64 const& x) { return UpperFloat64(pos_exact(x.u)); }
inline UpperFloat64 neg(LowerFloat64 const& x) { return UpperFloat64(neg_exact(x.l)); }
inline UpperFloat64 half(UpperFloat64 const& x) { return UpperFloat64(half_exact(x.u)); }
UpperFloat64 sqr(UpperFloat64 const& x);
UpperFloat64 rec(LowerFloat64 const& x);

inline UpperFloat64 add(UpperFloat64 const& x1, UpperFloat64 const& x2) { return UpperFloat64(add_up(x1.u,x2.u)); }
inline UpperFloat64 sub(UpperFloat64 const& x1, LowerFloat64 const& x2) { return UpperFloat64(sub_up(x1.u,x2.l)); }
UpperFloat64 mul(UpperFloat64 const& x1, UpperFloat64 const& x2);
UpperFloat64 div(UpperFloat64 const& x1, LowerFloat64 const& x2);
UpperFloat64 pow(UpperFloat64 const& x, Nat m);

UpperFloat64 sqrt(UpperFloat64 const& x);
UpperFloat64 exp(UpperFloat64 const& x);
UpperFloat64 log(UpperFloat64 const& x);
UpperFloat64 atan(UpperFloat64 const& x);

inline UpperFloat64 operator+(UpperFloat64 const& x) { return pos(x); }
inline UpperFloat64 operator-(LowerFloat64 const& x) { return neg(x); }
inline UpperFloat64 operator+(UpperFloat64 const& x1, UpperFloat64 const& x2) { return add(x1,x2); }
inline UpperFloat64 operator-(UpperFloat64 const& x1, LowerFloat64 const& x2) { return sub(x1,x2); }
inline UpperFloat64 operator*(UpperFloat64 const& x1, UpperFloat64 const& x2) { return mul(x1,x2); }
inline UpperFloat64 operator/(UpperFloat64 const& x1, LowerFloat64 const& x2) { return div(x1,x2); }
inline UpperFloat64& operator+=(UpperFloat64& x1, UpperFloat64 const& x2) { return x1=x1+x2; }
inline UpperFloat64& operator-=(UpperFloat64& x1, LowerFloat64 const& x2) { return x1=x1-x2; }
inline UpperFloat64& operator*=(UpperFloat64& x1, UpperFloat64 const& x2) { return x1=x1*x2; }
inline UpperFloat64& operator/=(UpperFloat64& x1, LowerFloat64 const& x2) { return x1=x1/x2; }

OutputStream& operator<<(OutputStream& os, UpperFloat64 const& x);
InputStream& operator>>(InputStream& is, UpperFloat64& x);





ValidatedFloat64 max(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2);
ValidatedFloat64 min(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2);
ValidatedFloat64 abs(ValidatedFloat64 const& x);

ValidatedFloat64 nul(ValidatedFloat64 const& x);
ValidatedFloat64 pos(ValidatedFloat64 const& x);
ValidatedFloat64 neg(ValidatedFloat64 const& x);
ValidatedFloat64 half(ValidatedFloat64 const& x);
ValidatedFloat64 sqr(ValidatedFloat64 const& x);
ValidatedFloat64 rec(ValidatedFloat64 const& x);

ValidatedFloat64 add(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2);
ValidatedFloat64 sub(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2);
ValidatedFloat64 mul(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2);
ValidatedFloat64 div(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2);
ValidatedFloat64 pow(ValidatedFloat64 const& x, Nat m);
ValidatedFloat64 pow(ValidatedFloat64 const& x, Int m);

ValidatedFloat64 sqrt(ValidatedFloat64 const& x);
ValidatedFloat64 exp(ValidatedFloat64 const& x);
ValidatedFloat64 log(ValidatedFloat64 const& x);
ValidatedFloat64 sin(ValidatedFloat64 const& x);
ValidatedFloat64 cos(ValidatedFloat64 const& x);
ValidatedFloat64 tan(ValidatedFloat64 const& x);
ValidatedFloat64 asin(ValidatedFloat64 const& x);
ValidatedFloat64 acos(ValidatedFloat64 const& x);
ValidatedFloat64 atan(ValidatedFloat64 const& x);

inline ValidatedFloat64 operator+(ValidatedFloat64 const& x) { return pos(x); }
inline ValidatedFloat64 operator-(ValidatedFloat64 const& x) { return neg(x); }
inline ValidatedFloat64 operator+(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) { return add(x1,x2); }
inline ValidatedFloat64 operator-(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) { return sub(x1,x2); }
inline ValidatedFloat64 operator*(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) { return mul(x1,x2); }
inline ValidatedFloat64 operator/(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) { return div(x1,x2); }
inline ValidatedFloat64& operator+=(ValidatedFloat64& x1, ValidatedFloat64 const& x2) { return x1=x1+x2; }
inline ValidatedFloat64& operator-=(ValidatedFloat64& x1, ValidatedFloat64 const& x2) { return x1=x1-x2; }
inline ValidatedFloat64& operator*=(ValidatedFloat64& x1, ValidatedFloat64 const& x2) { return x1=x1*x2; }
inline ValidatedFloat64& operator/=(ValidatedFloat64& x1, ValidatedFloat64 const& x2) { return x1=x1/x2; }

OutputStream& operator<<(OutputStream& os, ValidatedFloat64 const& x);
InputStream& operator>>(InputStream& is, ValidatedFloat64& x);



extern const ExactFloat64 infty;
inline ExactFloat64 operator"" _exact(long double lx) { double x=lx; assert(x==lx); return ExactFloat64(x); }
inline TwoExp::operator ExactFloat64 () const { return ExactFloat64(this->get_d()); }

inline ExactFloat64 max(ExactFloat64 const& x1,  ExactFloat64 const& x2) { return ExactFloat64(max(x1.v,x2.v)); }
inline ExactFloat64 min(ExactFloat64 const& x1,  ExactFloat64 const& x2) { return ExactFloat64(min(x1.v,x2.v)); }
inline ExactFloat64 abs(ExactFloat64 const& x) { return ExactFloat64(abs(x.v)); }

inline ExactFloat64 nul(ExactFloat64 const& x) { return ExactFloat64(nul(x.v)); }
inline ExactFloat64 pos(ExactFloat64 const& x) { return ExactFloat64(pos(x.v)); }
inline ExactFloat64 neg(ExactFloat64 const& x) { return ExactFloat64(neg(x.v)); }
inline ExactFloat64 half(ExactFloat64 const& x) { return ExactFloat64(half(x.v)); }

inline ValidatedFloat64 sqr(ExactFloat64 const& x);
inline ValidatedFloat64 rec(ExactFloat64 const& x);
inline ValidatedFloat64 add(ExactFloat64 const& x1,  ExactFloat64 const& x2);
inline ValidatedFloat64 sub(ExactFloat64 const& x1,  ExactFloat64 const& x2);
inline ValidatedFloat64 mul(ExactFloat64 const& x1,  ExactFloat64 const& x2);
inline ValidatedFloat64 div(ExactFloat64 const& x1,  ExactFloat64 const& x2);
inline ValidatedFloat64 pow(ExactFloat64 const& x, Int n);

inline ValidatedFloat64 sqrt(ExactFloat64 const& x);
inline ValidatedFloat64 exp(ExactFloat64 const& x);
inline ValidatedFloat64 log(ExactFloat64 const& x);
inline ValidatedFloat64 sin(ExactFloat64 const& x);
inline ValidatedFloat64 cos(ExactFloat64 const& x);
inline ValidatedFloat64 tan(ExactFloat64 const& x);
inline ValidatedFloat64 atan(ExactFloat64 const& x);

inline ValidatedFloat64 rad(ExactFloat64 const& x1, ExactFloat64 const& x2);
inline ValidatedFloat64 med(ExactFloat64 const& x1, ExactFloat64 const& x2);

inline ExactFloat64 operator+(ExactFloat64 const& x) { return pos(x); }
inline ExactFloat64 operator-(ExactFloat64 const& x) { return neg(x); }
inline ValidatedFloat64 operator+(ExactFloat64 const& x1,  ExactFloat64 const& x2) { return add(x1,x2); }
inline ValidatedFloat64 operator-(ExactFloat64 const& x1,  ExactFloat64 const& x2) { return sub(x1,x2); }
inline ValidatedFloat64 operator*(ExactFloat64 const& x1,  ExactFloat64 const& x2) { return mul(x1,x2); }
inline ValidatedFloat64 operator/(ExactFloat64 const& x1,  ExactFloat64 const& x2) { return div(x1,x2); }

/*
inline ValidatedFloat64 operator+(ValidatedFloat64 const& x1,  ExactFloat64 const& x2);
inline ValidatedFloat64 operator-(ValidatedFloat64 const& x1,  ExactFloat64 const& x2);
inline ValidatedFloat64 operator*(ValidatedFloat64 const& x1,  ExactFloat64 const& x2);
inline ValidatedFloat64 operator/(ValidatedFloat64 const& x1,  ExactFloat64 const& x2);
inline ValidatedFloat64 operator+(ExactFloat64 const& x1,  ValidatedFloat64 const& x2);
inline ValidatedFloat64 operator-(ExactFloat64 const& x1,  ValidatedFloat64 const& x2);
inline ValidatedFloat64 operator*(ExactFloat64 const& x1,  ValidatedFloat64 const& x2);
inline ValidatedFloat64 operator/(ExactFloat64 const& x1,  ValidatedFloat64 const& x2);
*/

inline ExactFloat64 operator*(ExactFloat64 const& x, TwoExp y) { ExactFloat64 yv=y; return ExactFloat64(x.raw()*yv.raw()); }
inline ExactFloat64 operator/(ExactFloat64 const& x, TwoExp y) { ExactFloat64 yv=y; return ExactFloat64(x.raw()/yv.raw()); }
inline ExactFloat64& operator*=(ExactFloat64& x, TwoExp y) { ExactFloat64 yv=y; return x=ExactFloat64(x.raw()*yv.raw()); }
inline ExactFloat64& operator/=(ExactFloat64& x, TwoExp y) { ExactFloat64 yv=y; return x=ExactFloat64(x.raw()/yv.raw()); }

inline Boolean operator==(ExactFloat64 const& x1, ExactFloat64 const& x2) { return x1.raw()==x2.raw(); }
inline Boolean operator!=(ExactFloat64 const& x1, ExactFloat64 const& x2) { return x1.raw()!=x2.raw(); }
inline Boolean operator<=(ExactFloat64 const& x1, ExactFloat64 const& x2) { return x1.raw()<=x2.raw(); }
inline Boolean operator>=(ExactFloat64 const& x1, ExactFloat64 const& x2) { return x1.raw()>=x2.raw(); }
inline Boolean operator< (ExactFloat64 const& x1, ExactFloat64 const& x2) { return x1.raw()< x2.raw(); }
inline Boolean operator> (ExactFloat64 const& x1, ExactFloat64 const& x2) { return x1.raw()> x2.raw(); }

template<class N, EnableIf<IsIntegral<N>> =dummy> inline Bool operator==(ExactFloat64 const& x1, N n2) { return x1.raw()==n2; }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Bool operator!=(ExactFloat64 const& x1, N n2) { return x1.raw()!=n2; }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Bool operator<=(ExactFloat64 const& x1, N n2) { return x1.raw()<=n2; }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Bool operator>=(ExactFloat64 const& x1, N n2) { return x1.raw()>=n2; }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Bool operator< (ExactFloat64 const& x1, N n2) { return x1.raw()< n2; }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Bool operator> (ExactFloat64 const& x1, N n2) { return x1.raw()> n2; }

OutputStream& operator<<(OutputStream& os, ExactFloat64 const& x);


inline ExactFloat64 const& make_exact(RawFloat64 const& x) { return reinterpret_cast<ExactFloat64 const&>(x); }
inline ExactFloat64 const& make_exact(ApproximateFloat64 const& x) { return reinterpret_cast<ExactFloat64 const&>(x); }
inline ExactFloat64 const& make_exact(ExactFloat64 const& x) { return reinterpret_cast<ExactFloat64 const&>(x); }
ExactFloat64 make_exact(const Real& x);

template<template<class>class T> inline const T<ExactFloat64>& make_exact(const T<RawFloat64>& t) {
    return reinterpret_cast<const T<ExactFloat64>&>(t); }
template<template<class>class T> inline const T<ExactFloat64>& make_exact(const T<ApproximateFloat64>& t) {
    return reinterpret_cast<const T<ExactFloat64>&>(t); }
template<template<class>class T> inline const T<ExactFloat64>& make_exact(const T<ExactFloat64>& t) {
    return reinterpret_cast<const T<ExactFloat64>&>(t); }

inline RawFloat64 const& make_raw(RawFloat64 const& x) { return reinterpret_cast<RawFloat64 const&>(x); }
inline RawFloat64 const& make_raw(ApproximateFloat64 const& x) { return reinterpret_cast<RawFloat64 const&>(x); }
inline RawFloat64 const& make_raw(ExactFloat64 const& x) { return reinterpret_cast<RawFloat64 const&>(x); }

template<template<class>class T> inline const T<RawFloat64>& make_raw(const T<RawFloat64>& t) {
    return reinterpret_cast<const T<RawFloat64>&>(t); }
template<template<class>class T> inline const T<RawFloat64>& make_raw(const T<ApproximateFloat64>& t) {
    return reinterpret_cast<const T<RawFloat64>&>(t); }
template<template<class>class T> inline const T<RawFloat64>& make_raw(const T<ExactFloat64>& t) {
    return reinterpret_cast<const T<RawFloat64>&>(t); }

inline ApproximateFloat64 const& make_approximate(RawFloat64 const& x) { return reinterpret_cast<ApproximateFloat64 const&>(x); }
inline ApproximateFloat64 const& make_approximate(ApproximateFloat64 const& x) { return reinterpret_cast<ApproximateFloat64 const&>(x); }
inline ApproximateFloat64 const& make_approximate(ExactFloat64 const& x) { return reinterpret_cast<ApproximateFloat64 const&>(x); }

template<template<class>class T> inline const T<ApproximateFloat64>& make_approximate(const T<RawFloat64>& t) {
    return reinterpret_cast<const T<ApproximateFloat64>&>(t); }
template<template<class>class T> inline const T<ApproximateFloat64>& make_approximate(const T<ApproximateFloat64>& t) {
    return reinterpret_cast<const T<ApproximateFloat64>&>(t); }
template<template<class>class T> inline const T<ApproximateFloat64>& make_approximate(const T<ExactFloat64>& t) {
    return reinterpret_cast<const T<ApproximateFloat64>&>(t); }


inline Bool operator==(ExactFloat64 const& x, const Rational& q) { return Rational(x)==q; }
inline Bool operator!=(ExactFloat64 const& x, const Rational& q) { return Rational(x)!=q; }
inline Bool operator<=(ExactFloat64 const& x, const Rational& q) { return Rational(x)<=q; }
inline Bool operator>=(ExactFloat64 const& x, const Rational& q) { return Rational(x)>=q; }
inline Bool operator< (ExactFloat64 const& x, const Rational& q) { return Rational(x)< q; }
inline Bool operator> (ExactFloat64 const& x, const Rational& q) { return Rational(x)> q; }

inline Bool operator==(const Rational& q, ExactFloat64 const& x) { return q==Rational(x); }
inline Bool operator!=(const Rational& q, ExactFloat64 const& x) { return q!=Rational(x); }
inline Bool operator<=(const Rational& q, ExactFloat64 const& x) { return q<=Rational(x); }
inline Bool operator>=(const Rational& q, ExactFloat64 const& x) { return q>=Rational(x); }
inline Bool operator< (const Rational& q, ExactFloat64 const& x) { return q< Rational(x); }
inline Bool operator> (const Rational& q, ExactFloat64 const& x) { return q> Rational(x); }


inline PositiveExactFloat64 mag(ExactFloat64 const& x) {
    return PositiveExactFloat64(abs(x.raw())); }
// FIXME: Unsafe since x may be negative
inline PositiveUpperFloat64 mag(UpperFloat64 const& x) {
    return PositiveUpperFloat64(abs(x.raw())); }
inline PositiveUpperFloat64 mag(ValidatedFloat64 const& x) {
    return PositiveUpperFloat64(max(neg(x.lower_raw()),x.upper_raw())); }
inline PositiveLowerFloat64 mig(ValidatedFloat64 const& x) {
    Float64 r=max(x.lower_raw(),neg(x.upper_raw()));
    return PositiveLowerFloat64(max(r,nul(r))); }
inline PositiveApproximateFloat64 mag(ApproximateFloat64 const& x) {
    return PositiveApproximateFloat64(abs(x.raw())); }


inline ValidatedFloat64 make_bounds(PositiveUpperFloat64 const& e) {
    return ValidatedFloat64(-e.raw(),+e.raw());
}

inline ExactFloat64 value(ValidatedFloat64 const& x) {
    return ExactFloat64(half_exact(add_near(x.lower_raw(),x.upper_raw())));
}

inline PositiveUpperFloat64 error(ValidatedFloat64 const& x) {
    return PositiveUpperFloat64(half_exact(sub_up(x.upper_raw(),x.lower_raw())));
}



template<class PR> inline Float<Approximate,PR>::Float(Float<Lower,PR> const& x) : a(x.raw()) {
}

template<class PR> inline Float<Approximate,PR>::Float(Float<Upper,PR> const& x) : a(x.raw()) {
}

template<class PR> inline Float<Approximate,PR>::Float(Float<Validated,PR> const& x) : a(half_exact(add_near(x.lower_raw(),x.upper_raw()))) {
}

template<class PR> inline Float<Approximate,PR>::Float(Float<Exact,PR> const& x) : a(x.raw()) {
}

template<class PR> inline Float<Lower,PR>::Float(Float<Validated,PR> const& x) : l(x.lower_raw()) {
}

template<class PR> inline Float<Lower,PR>::Float(Float<Exact,PR> const& x) : l(x.raw()) {
}

template<class PR> inline Float<Upper,PR>::Float(Float<Validated,PR> const& x) : u(x.upper_raw()) {
}

template<class PR> inline Float<Upper,PR>::Float(Float<Exact,PR> const& x) : u(x.raw()) {
}

template<class PR> inline Float<Validated,PR>::Float(Float<Exact,PR> const& x) : l(x.raw()), u(x.raw()) {
}


inline Bool same(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) {
    return x1.raw()==x2.raw(); }

inline Bool same(LowerFloat64 const& x1, LowerFloat64 const& x2) {
    return x1.raw()==x2.raw(); }

inline Bool same(UpperFloat64 const& x1, UpperFloat64 const& x2) {
    return x1.raw()==x2.raw(); }

inline Bool same(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) {
    return x1.lower_raw()==x2.lower_raw() && x1.upper_raw()==x2.upper_raw(); }

inline Bool same(ExactFloat64 const& x1, ExactFloat64 const& x2) {
    return x1.raw()==x2.raw(); }





inline ExactFloat64 midpoint(ValidatedFloat64 const& x) { return x.value(); }






inline ValidatedFloat64 max(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2)
{
    return ValidatedFloat64(max(x1.lower_raw(),x2.lower_raw()),max(x1.upper_raw(),x2.upper_raw()));
}

inline ValidatedFloat64 min(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2)
{
    return ValidatedFloat64(min(x1.lower_raw(),x2.lower_raw()),min(x1.upper_raw(),x2.upper_raw()));
}


inline ValidatedFloat64 abs(ValidatedFloat64 const& x)
{
    if(x.lower_raw()>=0) {
        return ValidatedFloat64(x.lower_raw(),x.upper_raw());
    } else if(x.upper_raw()<=0) {
        return ValidatedFloat64(neg(x.upper_raw()),neg(x.lower_raw()));
    } else {
        return ValidatedFloat64(static_cast<Float64>(0.0),max(neg(x.lower_raw()),x.upper_raw()));
    }
}

inline ValidatedFloat64 pos(ValidatedFloat64 const& x)
{
    return ValidatedFloat64(pos(x.lower_raw()),pos(x.upper_raw()));
}

inline ValidatedFloat64 neg(ValidatedFloat64 const& x)
{
    return ValidatedFloat64(neg(x.upper_raw()),neg(x.lower_raw()));
}

inline ValidatedFloat64 half(ValidatedFloat64 const& x) {
    return ValidatedFloat64(half(x.lower_raw()),half(x.upper_raw()));
}

inline ValidatedFloat64 sqr(ExactFloat64 const& x)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64 const& xv=x.raw();
    Float64::set_rounding_downward();
    Float64 rl=mul(xv,xv);
    Float64::set_rounding_upward();
    Float64 ru=mul(xv,xv);
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}

inline ValidatedFloat64 rec(ExactFloat64 const& x)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64 const& xv=x.raw();
    Float64::set_rounding_downward();
    Float64 rl=1.0/xv;
    Float64::set_rounding_upward();
    Float64 ru=1.0/xv;
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}



inline ValidatedFloat64 add(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64 const& x1l=x1.lower_raw();
    Float64 const& x1u=x1.upper_raw();
    Float64 const& x2l=x2.lower_raw();
    Float64 const& x2u=x2.upper_raw();
    Float64::set_rounding_downward();
    Float64 rl=add(x1l,x2l);
    Float64::set_rounding_upward();
    Float64 ru=add(x1u,x2u);
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}

inline ValidatedFloat64 add(ValidatedFloat64 const& x1, ExactFloat64 const& x2)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64 const& x1l=x1.lower_raw();
    Float64 const& x1u=x1.upper_raw();
    Float64 const& x2v=x2.raw();
    Float64::set_rounding_downward();
    Float64 rl=add(x1l,x2v);
    Float64::set_rounding_upward();
    Float64 ru=add(x1u,x2v);
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}

inline ValidatedFloat64 add(ExactFloat64 const& x1, ValidatedFloat64 const& x2)
{
    return add(x2,x1);
}

inline ValidatedFloat64 add(ExactFloat64 const& x1, ExactFloat64 const& x2)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64 const& x1v=x1.raw();
    Float64 const& x2v=x2.raw();
    Float64::set_rounding_downward();
    Float64 rl=add(x1v,x2v);
    Float64::set_rounding_upward();
    Float64 ru=add(x1v,x2v);
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}

inline ValidatedFloat64 sub(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64 const& x1l=x1.lower_raw();
    Float64 const& x1u=x1.upper_raw();
    Float64 const& x2l=x2.lower_raw();
    Float64 const& x2u=x2.upper_raw();
    Float64::set_rounding_downward();
    Float64 rl=sub(x1l,x2u);
    Float64::set_rounding_upward();
    Float64 ru=sub(x1u,x2l);
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}

inline ValidatedFloat64 sub(ValidatedFloat64 const& x1, ExactFloat64 const& x2)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64 const& x1l=x1.lower_raw();
    Float64 const& x1u=x1.upper_raw();
    Float64 const& x2v=x2.raw();
    Float64::set_rounding_downward();
    Float64 rl=sub(x1l,x2v);
    Float64::set_rounding_upward();
    Float64 ru=sub(x1u,x2v);
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}

inline ValidatedFloat64 sub(ExactFloat64 const& x1, ValidatedFloat64 const& x2)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64 const& x1v=x1.raw();
    Float64 const& x2l=x2.lower_raw();
    Float64 const& x2u=x2.upper_raw();
    Float64::set_rounding_downward();
    Float64 rl=sub(x1v,x2u);
    Float64::set_rounding_upward();
    Float64 ru=sub(x1v,x2l);
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}

inline ValidatedFloat64 sub(ExactFloat64 const& x1, ExactFloat64 const& x2)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64 const& x1v=x1.raw();
    Float64 const& x2v=x2.raw();
    Float64::set_rounding_downward();
    Float64 rl=sub(x1v,x2v);
    Float64::set_rounding_upward();
    Float64 ru=sub(x1v,x2v);
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}

inline ValidatedFloat64 mul(ExactFloat64 const& x1, ExactFloat64 const& x2)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64 const& x1v=x1.raw();
    Float64 const& x2v=x2.raw();
    Float64::set_rounding_downward();
    Float64 rl=mul(x1v,x2v);
    Float64::set_rounding_upward();
    Float64 ru=mul(x1v,x2v);
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}

inline ValidatedFloat64 div(ExactFloat64 const& x1, ExactFloat64 const& x2)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64 const& x1v=x1.raw();
    Float64 const& x2v=x2.raw();
    Float64::set_rounding_downward();
    Float64 rl=div(x1v,x2v);
    Float64::set_rounding_upward();
    Float64 ru=div(x1v,x2v);
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}

inline ValidatedFloat64 pow(ExactFloat64 const& x1, Int n2) {
    return pow(ValidatedFloat64(x1),n2);
}

inline ValidatedFloat64 med(ExactFloat64 const& x1, ExactFloat64 const& x2) {
    return add(half(x1),half(x2));
}

inline ValidatedFloat64 rad(ExactFloat64 const& x1, ExactFloat64 const& x2) {
    return sub(half(x2),half(x1));
}

inline ValidatedFloat64 sqrt(ExactFloat64 const& x) {
    return sqrt(ValidatedFloat64(x));
}

inline ValidatedFloat64 exp(ExactFloat64 const& x) {
    return exp(ValidatedFloat64(x));
}

inline ValidatedFloat64 log(ExactFloat64 const& x) {
    return log(ValidatedFloat64(x));
}

inline ValidatedFloat64 sin(ExactFloat64 const& x) {
    return sin(ValidatedFloat64(x));
}

inline ValidatedFloat64 cos(ExactFloat64 const& x) {
    return cos(ValidatedFloat64(x));
}



inline ValidatedFloat64 med(ValidatedFloat64 const& x);

inline ValidatedFloat64 rad(ValidatedFloat64 const& x);


/*
inline ValidatedFloat64 operator+(ValidatedFloat64 const& x1, ExactFloat64 const& x2) { return add(x1,x2); }
inline ValidatedFloat64 operator-(ValidatedFloat64 const& x1, ExactFloat64 const& x2) { return sub(x1,x2); }
inline ValidatedFloat64 operator*(ValidatedFloat64 const& x1, ExactFloat64 const& x2) { return mul(x1,x2); }
inline ValidatedFloat64 operator/(ValidatedFloat64 const& x1, ExactFloat64 const& x2) { return div(x1,x2); }
inline ValidatedFloat64 operator+(ExactFloat64 const& x1, ValidatedFloat64 const& x2) { return add(x2,x1); }
inline ValidatedFloat64 operator-(ExactFloat64 const& x1, ValidatedFloat64 const& x2) { return sub(x1,x2); }
inline ValidatedFloat64 operator*(ExactFloat64 const& x1, ValidatedFloat64 const& x2) { return mul(x2,x1); }
inline ValidatedFloat64 operator/(ExactFloat64 const& x1, ValidatedFloat64 const& x2) { return div(x1,x2); }
*/

// Standard equality operators
//! \related ValidatedFloat64 \brief Tests if \a x1 provides tighter bounds than \a x2.
inline Bool refines(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) {
    return x1.lower_raw()>=x2.lower_raw() && x1.upper_raw()<=x2.upper_raw(); }

//! \related ValidatedFloat64 \brief The common refinement of \a x1 and \a x2.
inline ValidatedFloat64 refinement(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) {
    return ValidatedFloat64(max(x1.lower_raw(),x2.lower_raw()),min(x1.upper_raw(),x2.upper_raw())); }

//! \related ValidatedFloat64 \brief Tests if \a x1 and \a x2 are consistent with representing the same number.
inline Bool consistent(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) {
    return x1.lower_raw()<=x2.upper_raw() && x1.upper_raw()>=x2.lower_raw(); }

//! \related ValidatedFloat64 \brief  Tests if \a x1 and \a x2 are inconsistent with representing the same number.
inline Bool inconsistent(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) {
    return not consistent(x1,x2); }

//! \related ValidatedFloat64 \brief  Tests if \a x1 is a model for the exact value \a x2. number.
inline Bool models(ValidatedFloat64 const& x1, ExactFloat64 const& x2) {
    return x1.lower_raw()<=x2.raw() && x1.upper_raw()>=x2.raw(); }

// Standard equality operators
//! \related ValidatedFloat64 \brief Equality operator. Tests equality of intervals as geometric objects, so \c [0,1]==[0,1] returns \c true.
inline Bool operator==(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) { return x1.lower_raw()==x2.lower_raw() && x1.upper_raw()==x2.upper_raw(); }
//! \related ValidatedFloat64 \brief Inequality operator. Tests equality of intervals as geometric objects, so \c [0,2]!=[1,3] returns \c true,
//! even though the intervals possibly represent the same exact real value.
inline Bool operator!=(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) { return x1.lower_raw()!=x2.lower_raw() || x1.upper_raw()!=x2.upper_raw(); }



//! \related ValidatedFloat64 \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
//! Hence \c [1.0,3.0]>[0.0,2.0] yields \c indeterminate since the first interval could represent the number 1.25 and the second 1.75.
inline Tribool operator> (ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) {
    if(x1.lower_raw()> x2.upper_raw()) { return true; }
    else if(x1.upper_raw()<=x2.lower_raw()) { return false; }
    else { return indeterminate; }
}

//! \related ValidatedFloat64 \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
inline Tribool operator< (ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) {
    if(x1.upper_raw()< x2.lower_raw()) { return true; }
    else if(x1.lower_raw()>=x2.upper_raw()) { return false; }
    else { return indeterminate; }
}

//! \related ValidatedFloat64 \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
inline Tribool operator>=(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) {
    if(x1.lower_raw()>=x2.upper_raw()) { return true; }
    else if(x1.upper_raw()< x2.lower_raw()) { return false; }
    else { return indeterminate; }
}

//! \related ValidatedFloat64 \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
inline Tribool operator<=(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) {
    if(x1.upper_raw()<=x2.lower_raw()) { return true; }
    else if(x1.lower_raw()> x2.upper_raw()) { return false; }
    else { return indeterminate; }
}

template<class N, EnableIf<IsIntegral<N>> =dummy> inline Tribool operator==(ValidatedFloat64 const& x1, N n2) { return x1==ValidatedFloat64(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Tribool operator!=(ValidatedFloat64 const& x1, N n2) { return x1!=ValidatedFloat64(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Tribool operator<=(ValidatedFloat64 const& x1, N n2) { return x1<=ValidatedFloat64(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Tribool operator>=(ValidatedFloat64 const& x1, N n2) { return x1>=ValidatedFloat64(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Tribool operator< (ValidatedFloat64 const& x1, N n2) { return x1< ValidatedFloat64(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Tribool operator> (ValidatedFloat64 const& x1, N n2) { return x1> ValidatedFloat64(n2); }

#ifdef ARIADNE_ENABLE_SERIALIZATION
  template<class A> Void serialize(A& a, ValidatedFloat64& ivl, const Nat version) {
    a & ivl.lower_raw() & ivl.upper_raw(); }
#endif

OutputStream& operator<<(OutputStream&, ValidatedFloat64 const&);
InputStream& operator>>(InputStream&, ValidatedFloat64&);


PositiveUpperFloat64 operator+(PositiveUpperFloat64 const& x1, PositiveUpperFloat64 const& x2);
PositiveUpperFloat64 operator-(PositiveUpperFloat64 const& x1, LowerFloat64 const& x2);
PositiveUpperFloat64 operator*(PositiveUpperFloat64 const& x1, PositiveUpperFloat64 const& x2);
PositiveUpperFloat64 operator/(PositiveUpperFloat64 const& x1, LowerFloat64 const& x2);
PositiveUpperFloat64 pow(PositiveUpperFloat64 const& x, Nat m);
PositiveUpperFloat64 half(PositiveUpperFloat64 const& x);


inline ErrorFloat64 operator"" _error(long double lx) { double x=lx; assert(x==lx); return ErrorFloat64(Float64(x)); }

template<class PR> inline  Float<Approximate,PR> make_float(Number<Approximate> const& y, PR pr) { return Float<Approximate,PR>(y,pr); }
template<class PR> inline  Float<Lower,PR> make_float(Float<Lower,PR> const& y, PR pr) { return Float<Lower,PR>(y,pr); }
template<class PR> inline  Float<Upper,PR> make_float(Number<Upper> const& y, PR pr) { return Float<Upper,PR>(y,pr); }
template<class PR> inline  Float<Validated,PR> make_float(Number<Validated> const& y, PR pr) { return Float<Validated,PR>(y,pr); }
template<class PR> inline  Float<Validated,PR> make_float(Number<Effective> const& y, PR pr) { return Float<Validated,PR>(y,pr); }
template<class PR> inline  Float<Validated,PR> make_float(Number<Exact> const& y, PR pr) { return Float<Validated,PR>(y,pr); }
template<class PR> inline  Float<Validated,PR> make_float(Real const& y, PR pr) { return Float<Validated,PR>(y,pr); }

template<class X> struct IsGenericNumber : IsConvertible<X,Real> { };
template<> struct IsGenericNumber<Real> : True { };
template<class P> struct IsGenericNumber<Number<P>> : True { };

template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumber<Y>> =dummy> auto
operator+(X const& x, Y const& y) -> decltype(x+make_float(y,x.precision())) { return x+make_float(y,x.precision()); }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumber<Y>> =dummy> auto
operator-(X const& x, Y const& y) -> decltype(x-make_float(y,x.precision())) { return x-make_float(y,x.precision()); }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumber<Y>> =dummy> auto
operator*(X const& x, Y const& y) -> decltype(x*make_float(y,x.precision())) { return x*make_float(y,x.precision()); }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumber<Y>> =dummy> auto
operator/(X const& x, Y const& y) -> decltype(x/make_float(y,x.precision())) { return x/make_float(y,x.precision()); }

template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumber<Y>> =dummy> auto
operator+(Y const& y, X const& x) -> decltype(make_float(y,x.raw().precision())+x) { return make_float(y,x.precision())+x; }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumber<Y>> =dummy> auto
operator-(Y const& y, X const& x) -> decltype(make_float(y,x.precision())-x) { return make_float(y,x.precision())-x; }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumber<Y>> =dummy> auto
operator*(Y const& y, X const& x) -> decltype(make_float(y,x.precision())*x) { return make_float(y,x.precision())*x; }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumber<Y>> =dummy> auto
operator/(Y const& y, X const& x) -> decltype(make_float(y,x.precision())/x) { return make_float(y,x.precision())/x; }


} // namespace Ariadne

#endif
