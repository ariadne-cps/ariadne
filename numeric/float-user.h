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
//! The \c %Float class represents floating-point numbers.
//! Unless otherwise mentioned, operations on floating-point numbers are performed approximately, with no guarantees
//! on the output.
//!
//! To implement <em>interval arithmetic</em>, arithmetical operations of \c %Float can be performed with guaranteed rounding by
//! specifying \c _up and \c _down suffixes to arithmetical functions \c add, \c sub, \c mul and \c div.
//! Additionally, operations can be performed in the current <em>rounding mode</em> by using the \c _rnd suffix,
//! or with rounding reversed using the \c _opp suffix.
//! Operations can be specified to return an \c %ExactInterval answer by using the \c _ivl suffix.
//! The \c _approx suffix is provided to specifically indicate that the operation is computed approximately.
//!
//! %Ariadne floating-point numbers can be constructed by conversion from built-in C++ types.
//! Note that the value of a built-in floating-point value may differ from the mathematical value of the literal.
//! For example, while <c>%Float(3.25)</c> is represented exactly, <c>%Float(3.3)</c> has a value of \f$3.2999999999999998224\ldots\f$.
//! \note In the future, the construction of a \c %Float from a string literal may be supported.
//! \sa ExactInterval, Real, ExactFloat
class ApproximateFloat {
  public:
    typedef Approximate Paradigm;
    typedef ApproximateFloat NumericType;
  public:
    ApproximateFloat() : a() { }
    template<class N, EnableIf<IsIntegral<N>> =dummy> ApproximateFloat(N n) : a(n) { }
    template<class D, EnableIf<IsFloatingPoint<D>> =dummy> ApproximateFloat(D x) : a(x) { }
    explicit ApproximateFloat(Float const& x) : a(x) { }
    explicit ApproximateFloat(const Dyadic& d);
    explicit ApproximateFloat(const Decimal& d);

    explicit ApproximateFloat(const Integer& z);
    explicit ApproximateFloat(const Rational& q);
    explicit ApproximateFloat(const Real& r);
    explicit ApproximateFloat(const Number<Approximate>& x);
    operator Number<Approximate> () const;

    ApproximateFloat(ExactFloat const& x);
    ApproximateFloat(ValidatedFloat const& x);
    ApproximateFloat(UpperFloat const& x);
    ApproximateFloat(LowerFloat const& x);

    explicit operator Float () const { return this->a; }
    Float const& raw() const { return this->a; }
    Float& raw() { return this->a; }
    double get_d() const { return this->a.get_d(); }
  public:
    static Void set_output_precision(Nat p) { output_precision=p; }
    ApproximateFloat pm(ApproximateFloat e) { return *this; }
  private: public:
    static Nat output_precision;
    Float a;
};

template<class R, class A> R integer_cast(const A& a);


inline ApproximateFloat floor(ApproximateFloat const& x) { return ApproximateFloat(floor(x.a)); }
inline ApproximateFloat ceil(ApproximateFloat const& x) { return ApproximateFloat(ceil(x.a)); }

inline ApproximateFloat abs(ApproximateFloat const& x) { return ApproximateFloat(abs_exact(x.a)); }
inline ApproximateFloat max(ApproximateFloat const& x, ApproximateFloat y) { return ApproximateFloat(max_exact(x.a,y.a)); }
inline ApproximateFloat min(ApproximateFloat const& x, ApproximateFloat y) { return ApproximateFloat(min_exact(x.a,y.a)); }

inline ApproximateFloat nul(ApproximateFloat const& x) { return ApproximateFloat(nul_exact(x.a)); }
inline ApproximateFloat pos(ApproximateFloat const& x) { return ApproximateFloat(pos_exact(x.a)); }
inline ApproximateFloat neg(ApproximateFloat const& x) { return ApproximateFloat(neg_exact(x.a)); }
inline ApproximateFloat half(ApproximateFloat const& x) { return ApproximateFloat(half_exact(x.a)); }
inline ApproximateFloat sqr(ApproximateFloat const& x) { return ApproximateFloat(mul_near(x.a,x.a)); }
inline ApproximateFloat rec(ApproximateFloat const& x) { return ApproximateFloat(div_near(1.0,x.a)); }

inline ApproximateFloat add(ApproximateFloat const& x1, ApproximateFloat const& x2) { return ApproximateFloat(add_near(x1.a,x2.a)); }
inline ApproximateFloat sub(ApproximateFloat const& x1, ApproximateFloat const& x2) { return ApproximateFloat(sub_near(x1.a,x2.a)); }
inline ApproximateFloat mul(ApproximateFloat const& x1, ApproximateFloat const& x2) { return ApproximateFloat(mul_near(x1.a,x2.a)); }
inline ApproximateFloat div(ApproximateFloat const& x1, ApproximateFloat const& x2) { return ApproximateFloat(div_near(x1.a,x2.a)); }

inline ApproximateFloat pow(ApproximateFloat const& x, Nat m) { return ApproximateFloat(pow_approx(x.a,m)); }
inline ApproximateFloat pow(ApproximateFloat const& x, Int n) { return ApproximateFloat(pow_approx(x.a,n)); }

inline ApproximateFloat sqrt(ApproximateFloat const& x) { return ApproximateFloat(sqrt_approx(x.a)); }
inline ApproximateFloat exp(ApproximateFloat const& x) { return ApproximateFloat(exp_approx(x.a)); }
inline ApproximateFloat log(ApproximateFloat const& x) { return ApproximateFloat(log_approx(x.a)); }
inline ApproximateFloat sin(ApproximateFloat const& x) { return ApproximateFloat(sin_approx(x.a)); }
inline ApproximateFloat cos(ApproximateFloat const& x) { return ApproximateFloat(cos_approx(x.a)); }
inline ApproximateFloat tan(ApproximateFloat const& x) { return ApproximateFloat(tan_approx(x.a)); }
inline ApproximateFloat asin(ApproximateFloat const& x) { return ApproximateFloat(asin_approx(x.a)); }
inline ApproximateFloat acos(ApproximateFloat const& x) { return ApproximateFloat(acos_approx(x.a)); }
inline ApproximateFloat atan(ApproximateFloat const& x) { return ApproximateFloat(atan_approx(x.a)); }

inline ApproximateFloat operator+(ApproximateFloat const& x) { return pos(x); }
inline ApproximateFloat operator-(ApproximateFloat const& x) { return neg(x); }
inline ApproximateFloat operator+(ApproximateFloat const& x1, ApproximateFloat const& x2) { return add(x1,x2); }
inline ApproximateFloat operator-(ApproximateFloat const& x1, ApproximateFloat const& x2) { return sub(x1,x2); }
inline ApproximateFloat operator*(ApproximateFloat const& x1, ApproximateFloat const& x2) { return mul(x1,x2); }
inline ApproximateFloat operator/(ApproximateFloat const& x1, ApproximateFloat const& x2) { return div(x1,x2); }
inline ApproximateFloat& operator+=(ApproximateFloat& x1, ApproximateFloat const& x2) { x1.a+=x2.a; return x1; }
inline ApproximateFloat& operator-=(ApproximateFloat& x1, ApproximateFloat const& x2) { x1.a-=x2.a; return x1; }
inline ApproximateFloat& operator*=(ApproximateFloat& x1, ApproximateFloat const& x2) { x1.a*=x2.a; return x1; }
inline ApproximateFloat& operator/=(ApproximateFloat& x1, ApproximateFloat const& x2) { x1.a/=x2.a; return x1; }

inline Bool operator==(ApproximateFloat const& x1, ApproximateFloat const& x2) { return x1.a==x2.a; }
inline Bool operator!=(ApproximateFloat const& x1, ApproximateFloat const& x2) { return x1.a!=x2.a; }
inline Bool operator<=(ApproximateFloat const& x1, ApproximateFloat const& x2) { return x1.a<=x2.a; }
inline Bool operator>=(ApproximateFloat const& x1, ApproximateFloat const& x2) { return x1.a>=x2.a; }
inline Bool operator< (ApproximateFloat const& x1, ApproximateFloat const& x2) { return x1.a< x2.a; }
inline Bool operator> (ApproximateFloat const& x1, ApproximateFloat const& x2) { return x1.a> x2.a; }

OutputStream& operator<<(OutputStream& os, ApproximateFloat const& x);
InputStream& operator>>(InputStream& is, ApproximateFloat& x);



//! \ingroup NumericModule
//! \brief Floating-point lower bounds for real numbers.
class LowerFloat {
  public:
    typedef Lower Paradigm;
    typedef LowerFloat NumericType;
  public:
    LowerFloat() : l(0.0) { }
    template<class N, EnableIf<IsIntegral<N>> = dummy> LowerFloat(N n) : l(n) { }
    template<class X, EnableIf<IsFloatingPoint<X>> = dummy> explicit LowerFloat(X x) : l(x) { }
    explicit LowerFloat(Float const& x) : l(x) { }

    LowerFloat(ValidatedFloat const& x);
    LowerFloat(ExactFloat const& x);

    explicit LowerFloat(const Number<Lower>& x);
    operator Number<Lower> () const;

    explicit LowerFloat(const Real& x);
    explicit LowerFloat(const Rational& x);
    explicit LowerFloat(const Integer& x);

    Float const& raw() const { return l; }
    Float& raw() { return l; }
    double get_d() const { return l.get_d(); }
  private: public:
    static Nat output_precision;
    Float l;
};


//! \ingroup NumericModule
//! \brief Floating-point upper bounds for real numbers.
class UpperFloat {
  public:
    typedef Upper Paradigm;
    typedef UpperFloat NumericType;
  public:
    UpperFloat() : u(0.0) { }
    template<class N, EnableIf<IsIntegral<N>> = dummy> UpperFloat(N n) : u(n) { }
    template<class X, EnableIf<IsFloatingPoint<X>> = dummy> explicit UpperFloat(X x) : u(x) { }

    explicit UpperFloat(Float const& x) : u(x) { }

    UpperFloat(ValidatedFloat const& x);
    UpperFloat(ExactFloat const& x);

    explicit UpperFloat(const Real& x);
    explicit UpperFloat(const Rational& x);
    explicit UpperFloat(const Integer& x);
    explicit UpperFloat(const Number<Upper>& x);

    operator Number<Upper> () const;

    Float const& raw() const { return u; }
    Float& raw() { return u; }
    double get_d() const { return u.get_d(); }
  private: public:
    static Nat output_precision;
    Float u;
};


inline LowerFloat max(LowerFloat const& x1, LowerFloat const& x2) { return LowerFloat(max_exact(x1.l,x2.l)); }
inline LowerFloat min(LowerFloat const& x1, LowerFloat const& x2) { return LowerFloat(min_exact(x1.l,x2.l)); }

inline LowerFloat nul(LowerFloat const& x) { return LowerFloat(pos_exact(x.l)); }
inline LowerFloat pos(LowerFloat const& x) { return LowerFloat(pos_exact(x.l)); }
inline LowerFloat neg(UpperFloat const& x) { return LowerFloat(neg_exact(x.u)); }
inline LowerFloat half(LowerFloat const& x) { return LowerFloat(half_exact(x.l)); }
LowerFloat sqr(LowerFloat const& x);
LowerFloat rec(UpperFloat const& x);

inline LowerFloat add(LowerFloat const& x1, LowerFloat const& x2) { return LowerFloat(add_down(x1.l,x2.l)); }
inline LowerFloat sub(LowerFloat const& x1, UpperFloat const& x2) { return LowerFloat(sub_down(x1.l,x2.u)); }
LowerFloat mul(LowerFloat const& x1, LowerFloat const& x2);
LowerFloat div(LowerFloat const& x1, UpperFloat const& x2);
LowerFloat pow(LowerFloat const& x, Nat m);

LowerFloat sqrt(LowerFloat const& x);
LowerFloat exp(LowerFloat const& x);
LowerFloat log(LowerFloat const& x);
LowerFloat atan(LowerFloat const& x);

inline LowerFloat operator+(LowerFloat const& x) { return pos(x); }
inline LowerFloat operator-(UpperFloat const& x) { return neg(x); }
inline LowerFloat operator+(LowerFloat const& x1, LowerFloat const& x2) { return add(x1,x2); }
inline LowerFloat operator-(LowerFloat const& x1, UpperFloat const& x2) { return sub(x1,x2); }
inline LowerFloat operator*(LowerFloat const& x1, LowerFloat const& x2) { return mul(x1,x2); }
inline LowerFloat operator/(LowerFloat const& x1, UpperFloat const& x2) { return div(x1,x2); }
inline LowerFloat& operator+=(LowerFloat& x1, LowerFloat const& x2) { return x1=x1+x2; }
inline LowerFloat& operator-=(LowerFloat& x1, UpperFloat const& x2) { return x1=x1-x2; }
inline LowerFloat& operator*=(LowerFloat& x1, LowerFloat const& x2) { return x1=x1*x2; }
inline LowerFloat& operator/=(LowerFloat& x1, UpperFloat const& x2) { return x1=x1/x2; }

OutputStream& operator<<(OutputStream& os, LowerFloat const& x);
InputStream& operator>>(InputStream& is, LowerFloat& x);

inline UpperFloat max(UpperFloat const& x1, UpperFloat const& x2) { return UpperFloat(max_exact(x1.u,x2.u)); }
inline UpperFloat min(UpperFloat const& x1, UpperFloat const& x2) { return UpperFloat(min_exact(x1.u,x2.u)); }

inline UpperFloat nul(UpperFloat const& x) { return UpperFloat(pos_exact(x.u)); }
inline UpperFloat pos(UpperFloat const& x) { return UpperFloat(pos_exact(x.u)); }
inline UpperFloat neg(LowerFloat const& x) { return UpperFloat(neg_exact(x.l)); }
inline UpperFloat half(UpperFloat const& x) { return UpperFloat(half_exact(x.u)); }
UpperFloat sqr(UpperFloat const& x);
UpperFloat rec(LowerFloat const& x);

inline UpperFloat add(UpperFloat const& x1, UpperFloat const& x2) { return UpperFloat(add_up(x1.u,x2.u)); }
inline UpperFloat sub(UpperFloat const& x1, LowerFloat const& x2) { return UpperFloat(sub_up(x1.u,x2.l)); }
UpperFloat mul(UpperFloat const& x1, UpperFloat const& x2);
UpperFloat div(UpperFloat const& x1, LowerFloat const& x2);
UpperFloat pow(UpperFloat const& x, Nat m);

UpperFloat sqrt(UpperFloat const& x);
UpperFloat exp(UpperFloat const& x);
UpperFloat log(UpperFloat const& x);
UpperFloat atan(UpperFloat const& x);

inline UpperFloat operator+(UpperFloat const& x) { return pos(x); }
inline UpperFloat operator-(LowerFloat const& x) { return neg(x); }
inline UpperFloat operator+(UpperFloat const& x1, UpperFloat const& x2) { return add(x1,x2); }
inline UpperFloat operator-(UpperFloat const& x1, LowerFloat const& x2) { return sub(x1,x2); }
inline UpperFloat operator*(UpperFloat const& x1, UpperFloat const& x2) { return mul(x1,x2); }
inline UpperFloat operator/(UpperFloat const& x1, LowerFloat const& x2) { return div(x1,x2); }
inline UpperFloat& operator+=(UpperFloat& x1, UpperFloat const& x2) { return x1=x1+x2; }
inline UpperFloat& operator-=(UpperFloat& x1, LowerFloat const& x2) { return x1=x1-x2; }
inline UpperFloat& operator*=(UpperFloat& x1, UpperFloat const& x2) { return x1=x1*x2; }
inline UpperFloat& operator/=(UpperFloat& x1, LowerFloat const& x2) { return x1=x1/x2; }

OutputStream& operator<<(OutputStream& os, UpperFloat const& x);
InputStream& operator>>(InputStream& is, UpperFloat& x);





//! \ingroup NumericModule
//! \brief Validated bounds on a number with floating-point endpoints supporting outwardly-rounded arithmetic.
//! \details
//! Note that <c>%ValidatedFloat(3.3)</c> yields the singleton interval \f$[3.2999999999999998224,3.2999999999999998224]\f$ (the constant is first interpreted by the C++ compiler to give a C++ \c double, whereas <c>%ValidatedFloat("3.3")</c> yields the interval \f$[3.2999999999999998224,3.3000000000000002665]\f$ enclosing \f$3.3\f$.
//!
//! Comparison tests on \c ValidatedFloat use the idea that an interval represents a single number with an unknown value.
//! Hence the result is of type \c Tribool, which can take values { \c True, \c False, \c Indeterminate }.
//! Hence a test \f$[l_1,u_1]\leq [l_2,u_2]\f$ returns \c True if \f$u_1\leq u_2\f$, since in this case \f$x_1\leq x_2\f$ whenever \f$x_1\in[l_1,u_2]\f$ and \f$x_2\in[l_2,u_2]\f$, \c False if \f$l_1>u_2\f$, since in this case we know \f$x_1>x_2\f$, and \c Indeterminate otherwise, since in this case we can find \f$x_1,x_2\f$ making the result either true or false.
//! In the case of equality, the comparison \f$[l_1,u_1]\f$==\f$[l_2,u_2]\f$ only returns \c True if both intervals are singletons, since otherwise we can find values making the result either true of false.
//!
//! To obtain the lower and upper bounds of an interval, use \c ivl.lower() and \c ivl.upper().
//! To obtain the midpoint and radius, use \c ivl.midpoint() and \c ivl.radius().
//! Alternatives \c midpoint(ivl) and \c radius(ivl) are also provided.
//! Note that \c midpoint and \c radius return approximations to the true midpoint and radius of the interval. If \f$m\f$ and \f$r\f$ are the returned midpoint and radius of the interval \f$[l,u]\f$, the using exact arithmetic, we guarentee \f$m-r\leq l\f$ and \f$m+r\geq u\f$
//!
//! To test if an interval contains a point or another interval, use \c encloses(ValidatedFloat,Float) or \c encloses(ValidatedFloat,ValidatedFloat).
//! The test \c refines(ValidatedFloat,ValidatedFloat) can also be used.
//! \sa Float
//!
//! \par Python interface
//!
//! In the Python interface, %Ariadne intervals can be constructed from Python literals of the form \c {a:b} or (deprecated) \c [a,b] .
//! The former is preferred, as it cannot be confused with literals for other classes such as Vector and Array types.
//! Automatic conversion is used to convert ValidatedFloat literals of the form \c {a,b} to an ValidatedFloat in functions.
//!
//! Care must be taken when defining intervals using floating-point coefficients, since values are first converted to the nearest
//! representable value by the Python interpreter. <br><br>
//! \code
//!   ValidatedFloat({1.1:2.3}) # Create the interval [1.1000000000000001, 2.2999999999999998]
//!   ValidatedFloat({2.5:4.25}) # Create the interval [2.5, 4.25], which can be represented exactly
//!   ValidatedFloat([2.5,4.25]) # Alternative syntax for creating the interval [2.5, 4.25]
//! \endcode
class ValidatedFloat {
  public:
    typedef Validated Paradigm;
    typedef ValidatedFloat NumericType;
  public:
    ValidatedFloat() : l(0.0), u(0.0) { }
    template<class N, EnableIf<IsIntegral<N>> = dummy> ValidatedFloat(N n) : l(n), u(n) { }
    template<class X, EnableIf<IsFloatingPoint<X>> = dummy> explicit ValidatedFloat(X x) : l(x), u(x) { }
    explicit ValidatedFloat(Float const& x) : l(x), u(x) { }

    ValidatedFloat(ExactFloat const& x);

    explicit ValidatedFloat(const Dyadic& x);
    explicit ValidatedFloat(const Decimal& x);
    explicit ValidatedFloat(const Integer& z);
    explicit ValidatedFloat(const Rational& q);
    explicit ValidatedFloat(const Real& x);
    explicit ValidatedFloat(const Number<Validated>& x);
    operator Number<Validated> () const;

    template<class N1, class N2, EnableIf<And<IsIntegral<N1>,IsIntegral<N2>>> =dummy>
        ValidatedFloat(N1 lower, N2 upper) : l(lower), u(upper) { }
    ValidatedFloat(Float const& lower, Float const& upper) : l(lower), u(upper) { }
    ValidatedFloat(LowerFloat const& lower, UpperFloat const& upper) : l(lower.raw()), u(upper.raw()) { }
    ValidatedFloat(const Rational& lower, const Rational& upper);

    Float const& lower_raw() const { return l; }
    Float const& upper_raw() const { return u; }
    double get_d() const { return (l.get_d()+u.get_d())/2; }

    LowerFloat lower() const { return LowerFloat(l); }
    UpperFloat upper() const { return UpperFloat(u); }
    const ExactFloat value() const;
    const PositiveUpperFloat error() const;

    // DEPRECATED
    explicit operator Float () const { return (l+u)/2; }
    friend ExactFloat midpoint(ValidatedFloat const& x);
  public:
    static Nat output_precision;
    static Void set_output_precision(Nat p) { output_precision=p; }
  private: public:
    Float l, u;
};

ValidatedFloat max(ValidatedFloat const& x1, ValidatedFloat const& x2);
ValidatedFloat min(ValidatedFloat const& x1, ValidatedFloat const& x2);
ValidatedFloat abs(ValidatedFloat const& x);

ValidatedFloat nul(ValidatedFloat const& x);
ValidatedFloat pos(ValidatedFloat const& x);
ValidatedFloat neg(ValidatedFloat const& x);
ValidatedFloat half(ValidatedFloat const& x);
ValidatedFloat sqr(ValidatedFloat const& x);
ValidatedFloat rec(ValidatedFloat const& x);

ValidatedFloat add(ValidatedFloat const& x1, ValidatedFloat const& x2);
ValidatedFloat sub(ValidatedFloat const& x1, ValidatedFloat const& x2);
ValidatedFloat mul(ValidatedFloat const& x1, ValidatedFloat const& x2);
ValidatedFloat div(ValidatedFloat const& x1, ValidatedFloat const& x2);
ValidatedFloat pow(ValidatedFloat const& x, Nat m);
ValidatedFloat pow(ValidatedFloat const& x, Int m);

ValidatedFloat sqrt(ValidatedFloat const& x);
ValidatedFloat exp(ValidatedFloat const& x);
ValidatedFloat log(ValidatedFloat const& x);
ValidatedFloat sin(ValidatedFloat const& x);
ValidatedFloat cos(ValidatedFloat const& x);
ValidatedFloat tan(ValidatedFloat const& x);
ValidatedFloat asin(ValidatedFloat const& x);
ValidatedFloat acos(ValidatedFloat const& x);
ValidatedFloat atan(ValidatedFloat const& x);

inline ValidatedFloat operator+(ValidatedFloat const& x) { return pos(x); }
inline ValidatedFloat operator-(ValidatedFloat const& x) { return neg(x); }
inline ValidatedFloat operator+(ValidatedFloat const& x1, ValidatedFloat const& x2) { return add(x1,x2); }
inline ValidatedFloat operator-(ValidatedFloat const& x1, ValidatedFloat const& x2) { return sub(x1,x2); }
inline ValidatedFloat operator*(ValidatedFloat const& x1, ValidatedFloat const& x2) { return mul(x1,x2); }
inline ValidatedFloat operator/(ValidatedFloat const& x1, ValidatedFloat const& x2) { return div(x1,x2); }
inline ValidatedFloat& operator+=(ValidatedFloat& x1, ValidatedFloat const& x2) { return x1=x1+x2; }
inline ValidatedFloat& operator-=(ValidatedFloat& x1, ValidatedFloat const& x2) { return x1=x1-x2; }
inline ValidatedFloat& operator*=(ValidatedFloat& x1, ValidatedFloat const& x2) { return x1=x1*x2; }
inline ValidatedFloat& operator/=(ValidatedFloat& x1, ValidatedFloat const& x2) { return x1=x1/x2; }

OutputStream& operator<<(OutputStream& os, ValidatedFloat const& x);
InputStream& operator>>(InputStream& is, ValidatedFloat& x);


//! \ingroup NumericModule
//! \related Float, ValidatedFloat
//! \brief A floating-point number, which is taken to represent the \em exact value of a real quantity.
class ExactFloat {
  public:
    typedef Exact Paradigm;
    typedef ExactFloat NumericType;

    ExactFloat() : v(0) { }
    template<class N, EnableIf<IsIntegral<N>> =dummy> ExactFloat(N n) : v(n) { }
    template<class X, EnableIf<IsFloatingPoint<X>> =dummy> explicit ExactFloat(X x) : v(x) { }

    explicit ExactFloat(Float const& x) : v(x) { }

    explicit ExactFloat(const Integer& z);
    explicit operator Rational () const;
    operator Number<Exact> () const;
    explicit operator Float () const { return v; }

    Float const& raw() const { return v; }
    Float& raw() { return v; }
    double get_d() const { return v.get_d(); }

    ValidatedFloat pm(ErrorFloat e) const;
  public:
    static Nat output_precision;
    static Void set_output_precision(Nat p) { output_precision=p; }
  private: public:
    Float v;
};

extern const ExactFloat infty;
inline ExactFloat operator"" _exact(long double lx) { double x=lx; assert(x==lx); return ExactFloat(x); }
inline TwoExp::operator ExactFloat () const { return ExactFloat(this->get_d()); }

inline ExactFloat max(ExactFloat const& x1,  ExactFloat const& x2) { return ExactFloat(max(x1.v,x2.v)); }
inline ExactFloat min(ExactFloat const& x1,  ExactFloat const& x2) { return ExactFloat(min(x1.v,x2.v)); }
inline ExactFloat abs(ExactFloat const& x) { return ExactFloat(abs(x.v)); }

inline ExactFloat nul(ExactFloat const& x) { return ExactFloat(nul(x.v)); }
inline ExactFloat pos(ExactFloat const& x) { return ExactFloat(pos(x.v)); }
inline ExactFloat neg(ExactFloat const& x) { return ExactFloat(neg(x.v)); }
inline ExactFloat half(ExactFloat const& x) { return ExactFloat(half(x.v)); }

inline ValidatedFloat sqr(ExactFloat const& x);
inline ValidatedFloat rec(ExactFloat const& x);
inline ValidatedFloat add(ExactFloat const& x1,  ExactFloat const& x2);
inline ValidatedFloat sub(ExactFloat const& x1,  ExactFloat const& x2);
inline ValidatedFloat mul(ExactFloat const& x1,  ExactFloat const& x2);
inline ValidatedFloat div(ExactFloat const& x1,  ExactFloat const& x2);
inline ValidatedFloat pow(ExactFloat const& x, Int n);

inline ValidatedFloat sqrt(ExactFloat const& x);
inline ValidatedFloat exp(ExactFloat const& x);
inline ValidatedFloat log(ExactFloat const& x);
inline ValidatedFloat sin(ExactFloat const& x);
inline ValidatedFloat cos(ExactFloat const& x);
inline ValidatedFloat tan(ExactFloat const& x);
inline ValidatedFloat atan(ExactFloat const& x);

inline ValidatedFloat rad(ExactFloat const& x1, ExactFloat const& x2);
inline ValidatedFloat med(ExactFloat const& x1, ExactFloat const& x2);

inline ExactFloat operator+(ExactFloat const& x) { return pos(x); }
inline ExactFloat operator-(ExactFloat const& x) { return neg(x); }
inline ValidatedFloat operator+(ExactFloat const& x1,  ExactFloat const& x2) { return add(x1,x2); }
inline ValidatedFloat operator-(ExactFloat const& x1,  ExactFloat const& x2) { return sub(x1,x2); }
inline ValidatedFloat operator*(ExactFloat const& x1,  ExactFloat const& x2) { return mul(x1,x2); }
inline ValidatedFloat operator/(ExactFloat const& x1,  ExactFloat const& x2) { return div(x1,x2); }

/*
inline ValidatedFloat operator+(ValidatedFloat const& x1,  ExactFloat const& x2);
inline ValidatedFloat operator-(ValidatedFloat const& x1,  ExactFloat const& x2);
inline ValidatedFloat operator*(ValidatedFloat const& x1,  ExactFloat const& x2);
inline ValidatedFloat operator/(ValidatedFloat const& x1,  ExactFloat const& x2);
inline ValidatedFloat operator+(ExactFloat const& x1,  ValidatedFloat const& x2);
inline ValidatedFloat operator-(ExactFloat const& x1,  ValidatedFloat const& x2);
inline ValidatedFloat operator*(ExactFloat const& x1,  ValidatedFloat const& x2);
inline ValidatedFloat operator/(ExactFloat const& x1,  ValidatedFloat const& x2);
*/

inline ExactFloat operator*(ExactFloat const& x, TwoExp y) { ExactFloat yv=y; return ExactFloat(x.raw()*yv.raw()); }
inline ExactFloat operator/(ExactFloat const& x, TwoExp y) { ExactFloat yv=y; return ExactFloat(x.raw()/yv.raw()); }
inline ExactFloat& operator*=(ExactFloat& x, TwoExp y) { ExactFloat yv=y; return x=ExactFloat(x.raw()*yv.raw()); }
inline ExactFloat& operator/=(ExactFloat& x, TwoExp y) { ExactFloat yv=y; return x=ExactFloat(x.raw()/yv.raw()); }

inline Boolean operator==(ExactFloat const& x1, ExactFloat const& x2) { return x1.raw()==x2.raw(); }
inline Boolean operator!=(ExactFloat const& x1, ExactFloat const& x2) { return x1.raw()!=x2.raw(); }
inline Boolean operator<=(ExactFloat const& x1, ExactFloat const& x2) { return x1.raw()<=x2.raw(); }
inline Boolean operator>=(ExactFloat const& x1, ExactFloat const& x2) { return x1.raw()>=x2.raw(); }
inline Boolean operator< (ExactFloat const& x1, ExactFloat const& x2) { return x1.raw()< x2.raw(); }
inline Boolean operator> (ExactFloat const& x1, ExactFloat const& x2) { return x1.raw()> x2.raw(); }

template<class N, EnableIf<IsIntegral<N>> =dummy> inline Bool operator==(ExactFloat const& x1, N n2) { return x1.raw()==n2; }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Bool operator!=(ExactFloat const& x1, N n2) { return x1.raw()!=n2; }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Bool operator<=(ExactFloat const& x1, N n2) { return x1.raw()<=n2; }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Bool operator>=(ExactFloat const& x1, N n2) { return x1.raw()>=n2; }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Bool operator< (ExactFloat const& x1, N n2) { return x1.raw()< n2; }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Bool operator> (ExactFloat const& x1, N n2) { return x1.raw()> n2; }

OutputStream& operator<<(OutputStream& os, ExactFloat const& x);


inline ExactFloat const& make_exact(RawFloat const& x) { return reinterpret_cast<ExactFloat const&>(x); }
inline ExactFloat const& make_exact(ApproximateFloat const& x) { return reinterpret_cast<ExactFloat const&>(x); }
inline ExactFloat const& make_exact(ExactFloat const& x) { return reinterpret_cast<ExactFloat const&>(x); }
ExactFloat make_exact(const Real& x);

template<template<class>class T> inline const T<ExactFloat>& make_exact(const T<RawFloat>& t) {
    return reinterpret_cast<const T<ExactFloat>&>(t); }
template<template<class>class T> inline const T<ExactFloat>& make_exact(const T<ApproximateFloat>& t) {
    return reinterpret_cast<const T<ExactFloat>&>(t); }
template<template<class>class T> inline const T<ExactFloat>& make_exact(const T<ExactFloat>& t) {
    return reinterpret_cast<const T<ExactFloat>&>(t); }

inline RawFloat const& make_raw(RawFloat const& x) { return reinterpret_cast<RawFloat const&>(x); }
inline RawFloat const& make_raw(ApproximateFloat const& x) { return reinterpret_cast<RawFloat const&>(x); }
inline RawFloat const& make_raw(ExactFloat const& x) { return reinterpret_cast<RawFloat const&>(x); }

template<template<class>class T> inline const T<RawFloat>& make_raw(const T<RawFloat>& t) {
    return reinterpret_cast<const T<RawFloat>&>(t); }
template<template<class>class T> inline const T<RawFloat>& make_raw(const T<ApproximateFloat>& t) {
    return reinterpret_cast<const T<RawFloat>&>(t); }
template<template<class>class T> inline const T<RawFloat>& make_raw(const T<ExactFloat>& t) {
    return reinterpret_cast<const T<RawFloat>&>(t); }

inline ApproximateFloat const& make_approximate(RawFloat const& x) { return reinterpret_cast<ApproximateFloat const&>(x); }
inline ApproximateFloat const& make_approximate(ApproximateFloat const& x) { return reinterpret_cast<ApproximateFloat const&>(x); }
inline ApproximateFloat const& make_approximate(ExactFloat const& x) { return reinterpret_cast<ApproximateFloat const&>(x); }

template<template<class>class T> inline const T<ApproximateFloat>& make_approximate(const T<RawFloat>& t) {
    return reinterpret_cast<const T<ApproximateFloat>&>(t); }
template<template<class>class T> inline const T<ApproximateFloat>& make_approximate(const T<ApproximateFloat>& t) {
    return reinterpret_cast<const T<ApproximateFloat>&>(t); }
template<template<class>class T> inline const T<ApproximateFloat>& make_approximate(const T<ExactFloat>& t) {
    return reinterpret_cast<const T<ApproximateFloat>&>(t); }


inline Bool operator==(ExactFloat const& x, const Rational& q) { return Rational(x)==q; }
inline Bool operator!=(ExactFloat const& x, const Rational& q) { return Rational(x)!=q; }
inline Bool operator<=(ExactFloat const& x, const Rational& q) { return Rational(x)<=q; }
inline Bool operator>=(ExactFloat const& x, const Rational& q) { return Rational(x)>=q; }
inline Bool operator< (ExactFloat const& x, const Rational& q) { return Rational(x)< q; }
inline Bool operator> (ExactFloat const& x, const Rational& q) { return Rational(x)> q; }

inline Bool operator==(const Rational& q, ExactFloat const& x) { return q==Rational(x); }
inline Bool operator!=(const Rational& q, ExactFloat const& x) { return q!=Rational(x); }
inline Bool operator<=(const Rational& q, ExactFloat const& x) { return q<=Rational(x); }
inline Bool operator>=(const Rational& q, ExactFloat const& x) { return q>=Rational(x); }
inline Bool operator< (const Rational& q, ExactFloat const& x) { return q< Rational(x); }
inline Bool operator> (const Rational& q, ExactFloat const& x) { return q> Rational(x); }


class PositiveExactFloat : public ExactFloat {
  public:
    PositiveExactFloat() : ExactFloat() { }
    template<class M, EnableIf<IsIntegral<M>>, EnableIf<IsUnsigned<M>> =dummy>
        PositiveExactFloat(M m) : ExactFloat(m) { }
    explicit PositiveExactFloat(Float const& x) : ExactFloat(x) { }
};

class PositiveUpperFloat : public UpperFloat {
  public:
    PositiveUpperFloat() : UpperFloat() { }
    explicit PositiveUpperFloat(Float const& x) : UpperFloat(x) { ARIADNE_PRECONDITION(x>=0); }
    template<class M, EnableIf<IsUnsigned<M>> =dummy> PositiveUpperFloat(M m) : UpperFloat(m) { }
    template<class F, EnableIf<IsSame<F,UpperFloat>> =dummy>
        explicit PositiveUpperFloat(F const& x) : UpperFloat(x) { }
    PositiveUpperFloat(PositiveExactFloat const& x) : UpperFloat(x) { }
};

class PositiveLowerFloat : public LowerFloat {
  public:
    PositiveLowerFloat() : LowerFloat() { }
    template<class M, EnableIf<IsSigned<M>> =dummy>
        PositiveLowerFloat(M m) : LowerFloat(m) { }
    explicit PositiveLowerFloat(Float const& x) : LowerFloat(x) { }
    PositiveLowerFloat(PositiveExactFloat const& x) : LowerFloat(x) { }
};

class PositiveApproximateFloat : public ApproximateFloat {
  public:
    PositiveApproximateFloat() : ApproximateFloat() { }
    template<class M, EnableIf<IsSigned<M>> =dummy>
        PositiveApproximateFloat(M m) : ApproximateFloat(m) { }
    explicit PositiveApproximateFloat(Float const& x) : ApproximateFloat(x) { }
    PositiveApproximateFloat(PositiveLowerFloat const& x) : ApproximateFloat(x) { }
    PositiveApproximateFloat(PositiveUpperFloat const& x) : ApproximateFloat(x) { }
    PositiveApproximateFloat(PositiveExactFloat const& x) : ApproximateFloat(x) { }
};

inline PositiveExactFloat mag(ExactFloat const& x) {
    return PositiveExactFloat(abs(x.raw())); }
// FIXME: Unsafe since x may be negative
inline PositiveUpperFloat mag(UpperFloat const& x) {
    return PositiveUpperFloat(abs(x.raw())); }
inline PositiveUpperFloat mag(ValidatedFloat const& x) {
    return PositiveUpperFloat(max(neg(x.lower_raw()),x.upper_raw())); }
inline PositiveLowerFloat mig(ValidatedFloat const& x) {
    Float r=max(x.lower_raw(),neg(x.upper_raw()));
    return PositiveLowerFloat(max(r,nul(r))); }
inline PositiveApproximateFloat mag(ApproximateFloat const& x) {
    return PositiveApproximateFloat(abs(x.raw())); }


inline ValidatedFloat make_bounds(PositiveUpperFloat const& e) {
    return ValidatedFloat(-e.raw(),+e.raw());
}

inline ExactFloat value(ValidatedFloat const& x) {
    return ExactFloat(half_exact(add_near(x.lower_raw(),x.upper_raw())));
}

inline PositiveUpperFloat error(ValidatedFloat const& x) {
    return PositiveUpperFloat(half_exact(sub_up(x.upper_raw(),x.lower_raw())));
}



inline ApproximateFloat::ApproximateFloat(LowerFloat const& x) : a(x.raw()) {
}

inline ApproximateFloat::ApproximateFloat(UpperFloat const& x) : a(x.raw()) {
}

inline ApproximateFloat::ApproximateFloat(ValidatedFloat const& x) : a(half_exact(add_near(x.lower_raw(),x.upper_raw()))) {
}

inline ApproximateFloat::ApproximateFloat(ExactFloat const& x) : a(x.raw()) {
}

inline LowerFloat::LowerFloat(ValidatedFloat const& x) : l(x.lower_raw()) {
}

inline LowerFloat::LowerFloat(ExactFloat const& x) : l(x.raw()) {
}

inline UpperFloat::UpperFloat(ValidatedFloat const& x) : u(x.upper_raw()) {
}

inline UpperFloat::UpperFloat(ExactFloat const& x) : u(x.raw()) {
}

inline ValidatedFloat::ValidatedFloat(ExactFloat const& x) : l(x.raw()), u(x.raw()) {
}


inline Bool same(ApproximateFloat const& x1, ApproximateFloat const& x2) {
    return x1.raw()==x2.raw(); }

inline Bool same(LowerFloat const& x1, LowerFloat const& x2) {
    return x1.raw()==x2.raw(); }

inline Bool same(UpperFloat const& x1, UpperFloat const& x2) {
    return x1.raw()==x2.raw(); }

inline Bool same(ValidatedFloat const& x1, ValidatedFloat const& x2) {
    return x1.lower_raw()==x2.lower_raw() && x1.upper_raw()==x2.upper_raw(); }

inline Bool same(ExactFloat const& x1, ExactFloat const& x2) {
    return x1.raw()==x2.raw(); }





inline const ExactFloat ValidatedFloat::value() const {
    return ExactFloat(med_near(this->l,this->u)); }

inline const ErrorFloat ValidatedFloat::error() const {
    Float v=med_near(this->l,this->u); return ErrorFloat(max(sub_up(this->u,v),sub_up(v,this->l))); }

inline ExactFloat midpoint(ValidatedFloat const& x) { return x.value(); }






inline ValidatedFloat max(ValidatedFloat const& x1, ValidatedFloat const& x2)
{
    return ValidatedFloat(max(x1.lower_raw(),x2.lower_raw()),max(x1.upper_raw(),x2.upper_raw()));
}

inline ValidatedFloat min(ValidatedFloat const& x1, ValidatedFloat const& x2)
{
    return ValidatedFloat(min(x1.lower_raw(),x2.lower_raw()),min(x1.upper_raw(),x2.upper_raw()));
}


inline ValidatedFloat abs(ValidatedFloat const& x)
{
    if(x.lower_raw()>=0) {
        return ValidatedFloat(x.lower_raw(),x.upper_raw());
    } else if(x.upper_raw()<=0) {
        return ValidatedFloat(neg(x.upper_raw()),neg(x.lower_raw()));
    } else {
        return ValidatedFloat(static_cast<Float>(0.0),max(neg(x.lower_raw()),x.upper_raw()));
    }
}

inline ValidatedFloat pos(ValidatedFloat const& x)
{
    return ValidatedFloat(pos(x.lower_raw()),pos(x.upper_raw()));
}

inline ValidatedFloat neg(ValidatedFloat const& x)
{
    return ValidatedFloat(neg(x.upper_raw()),neg(x.lower_raw()));
}

inline ValidatedFloat half(ValidatedFloat const& x) {
    return ValidatedFloat(half(x.lower_raw()),half(x.upper_raw()));
}

inline ValidatedFloat sqr(ExactFloat const& x)
{
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    Float const& xv=x.raw();
    Float::set_rounding_downward();
    Float rl=mul(xv,xv);
    Float::set_rounding_upward();
    Float ru=mul(xv,xv);
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

inline ValidatedFloat rec(ExactFloat const& x)
{
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    Float const& xv=x.raw();
    Float::set_rounding_downward();
    Float rl=1.0/xv;
    Float::set_rounding_upward();
    Float ru=1.0/xv;
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}



inline ValidatedFloat add(ValidatedFloat const& x1, ValidatedFloat const& x2)
{
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    Float const& x1l=x1.lower_raw();
    Float const& x1u=x1.upper_raw();
    Float const& x2l=x2.lower_raw();
    Float const& x2u=x2.upper_raw();
    Float::set_rounding_downward();
    Float rl=add(x1l,x2l);
    Float::set_rounding_upward();
    Float ru=add(x1u,x2u);
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

inline ValidatedFloat add(ValidatedFloat const& x1, ExactFloat const& x2)
{
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    Float const& x1l=x1.lower_raw();
    Float const& x1u=x1.upper_raw();
    Float const& x2v=x2.raw();
    Float::set_rounding_downward();
    Float rl=add(x1l,x2v);
    Float::set_rounding_upward();
    Float ru=add(x1u,x2v);
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

inline ValidatedFloat add(ExactFloat const& x1, ValidatedFloat const& x2)
{
    return add(x2,x1);
}

inline ValidatedFloat add(ExactFloat const& x1, ExactFloat const& x2)
{
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    Float const& x1v=x1.raw();
    Float const& x2v=x2.raw();
    Float::set_rounding_downward();
    Float rl=add(x1v,x2v);
    Float::set_rounding_upward();
    Float ru=add(x1v,x2v);
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

inline ValidatedFloat sub(ValidatedFloat const& x1, ValidatedFloat const& x2)
{
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    Float const& x1l=x1.lower_raw();
    Float const& x1u=x1.upper_raw();
    Float const& x2l=x2.lower_raw();
    Float const& x2u=x2.upper_raw();
    Float::set_rounding_downward();
    Float rl=sub(x1l,x2u);
    Float::set_rounding_upward();
    Float ru=sub(x1u,x2l);
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

inline ValidatedFloat sub(ValidatedFloat const& x1, ExactFloat const& x2)
{
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    Float const& x1l=x1.lower_raw();
    Float const& x1u=x1.upper_raw();
    Float const& x2v=x2.raw();
    Float::set_rounding_downward();
    Float rl=sub(x1l,x2v);
    Float::set_rounding_upward();
    Float ru=sub(x1u,x2v);
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

inline ValidatedFloat sub(ExactFloat const& x1, ValidatedFloat const& x2)
{
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    Float const& x1v=x1.raw();
    Float const& x2l=x2.lower_raw();
    Float const& x2u=x2.upper_raw();
    Float::set_rounding_downward();
    Float rl=sub(x1v,x2u);
    Float::set_rounding_upward();
    Float ru=sub(x1v,x2l);
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

inline ValidatedFloat sub(ExactFloat const& x1, ExactFloat const& x2)
{
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    Float const& x1v=x1.raw();
    Float const& x2v=x2.raw();
    Float::set_rounding_downward();
    Float rl=sub(x1v,x2v);
    Float::set_rounding_upward();
    Float ru=sub(x1v,x2v);
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

inline ValidatedFloat mul(ExactFloat const& x1, ExactFloat const& x2)
{
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    Float const& x1v=x1.raw();
    Float const& x2v=x2.raw();
    Float::set_rounding_downward();
    Float rl=mul(x1v,x2v);
    Float::set_rounding_upward();
    Float ru=mul(x1v,x2v);
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

inline ValidatedFloat div(ExactFloat const& x1, ExactFloat const& x2)
{
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    Float const& x1v=x1.raw();
    Float const& x2v=x2.raw();
    Float::set_rounding_downward();
    Float rl=div(x1v,x2v);
    Float::set_rounding_upward();
    Float ru=div(x1v,x2v);
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

inline ValidatedFloat pow(ExactFloat const& x1, Int n2) {
    return pow(ValidatedFloat(x1),n2);
}

inline ValidatedFloat med(ExactFloat const& x1, ExactFloat const& x2) {
    return add(half(x1),half(x2));
}

inline ValidatedFloat rad(ExactFloat const& x1, ExactFloat const& x2) {
    return sub(half(x2),half(x1));
}

inline ValidatedFloat sqrt(ExactFloat const& x) {
    return sqrt(ValidatedFloat(x));
}

inline ValidatedFloat exp(ExactFloat const& x) {
    return exp(ValidatedFloat(x));
}

inline ValidatedFloat log(ExactFloat const& x) {
    return log(ValidatedFloat(x));
}

inline ValidatedFloat sin(ExactFloat const& x) {
    return sin(ValidatedFloat(x));
}

inline ValidatedFloat cos(ExactFloat const& x) {
    return cos(ValidatedFloat(x));
}



inline ValidatedFloat med(ValidatedFloat const& x);

inline ValidatedFloat rad(ValidatedFloat const& x);


/*
inline ValidatedFloat operator+(ValidatedFloat const& x1, ExactFloat const& x2) { return add(x1,x2); }
inline ValidatedFloat operator-(ValidatedFloat const& x1, ExactFloat const& x2) { return sub(x1,x2); }
inline ValidatedFloat operator*(ValidatedFloat const& x1, ExactFloat const& x2) { return mul(x1,x2); }
inline ValidatedFloat operator/(ValidatedFloat const& x1, ExactFloat const& x2) { return div(x1,x2); }
inline ValidatedFloat operator+(ExactFloat const& x1, ValidatedFloat const& x2) { return add(x2,x1); }
inline ValidatedFloat operator-(ExactFloat const& x1, ValidatedFloat const& x2) { return sub(x1,x2); }
inline ValidatedFloat operator*(ExactFloat const& x1, ValidatedFloat const& x2) { return mul(x2,x1); }
inline ValidatedFloat operator/(ExactFloat const& x1, ValidatedFloat const& x2) { return div(x1,x2); }
*/

// Standard equality operators
//! \related ValidatedFloat \brief Tests if \a x1 provides tighter bounds than \a x2.
inline Bool refines(ValidatedFloat const& x1, ValidatedFloat const& x2) {
    return x1.lower_raw()>=x2.lower_raw() && x1.upper_raw()<=x2.upper_raw(); }

//! \related ValidatedFloat \brief The common refinement of \a x1 and \a x2.
inline ValidatedFloat refinement(ValidatedFloat const& x1, ValidatedFloat const& x2) {
    return ValidatedFloat(max(x1.lower_raw(),x2.lower_raw()),min(x1.upper_raw(),x2.upper_raw())); }

//! \related ValidatedFloat \brief Tests if \a x1 and \a x2 are consistent with representing the same number.
inline Bool consistent(ValidatedFloat const& x1, ValidatedFloat const& x2) {
    return x1.lower_raw()<=x2.upper_raw() && x1.upper_raw()>=x2.lower_raw(); }

//! \related ValidatedFloat \brief  Tests if \a x1 and \a x2 are inconsistent with representing the same number.
inline Bool inconsistent(ValidatedFloat const& x1, ValidatedFloat const& x2) {
    return not consistent(x1,x2); }

//! \related ValidatedFloat \brief  Tests if \a x1 is a model for the exact value \a x2. number.
inline Bool models(ValidatedFloat const& x1, ExactFloat const& x2) {
    return x1.lower_raw()<=x2.raw() && x1.upper_raw()>=x2.raw(); }

// Standard equality operators
//! \related ValidatedFloat \brief Equality operator. Tests equality of intervals as geometric objects, so \c [0,1]==[0,1] returns \c true.
inline Bool operator==(ValidatedFloat const& x1, ValidatedFloat const& x2) { return x1.lower_raw()==x2.lower_raw() && x1.upper_raw()==x2.upper_raw(); }
//! \related ValidatedFloat \brief Inequality operator. Tests equality of intervals as geometric objects, so \c [0,2]!=[1,3] returns \c true,
//! even though the intervals possibly represent the same exact real value.
inline Bool operator!=(ValidatedFloat const& x1, ValidatedFloat const& x2) { return x1.lower_raw()!=x2.lower_raw() || x1.upper_raw()!=x2.upper_raw(); }



//! \related ValidatedFloat \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
//! Hence \c [1.0,3.0]>[0.0,2.0] yields \c indeterminate since the first interval could represent the number 1.25 and the second 1.75.
inline Tribool operator> (ValidatedFloat const& x1, ValidatedFloat const& x2) {
    if(x1.lower_raw()> x2.upper_raw()) { return true; }
    else if(x1.upper_raw()<=x2.lower_raw()) { return false; }
    else { return indeterminate; }
}

//! \related ValidatedFloat \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
inline Tribool operator< (ValidatedFloat const& x1, ValidatedFloat const& x2) {
    if(x1.upper_raw()< x2.lower_raw()) { return true; }
    else if(x1.lower_raw()>=x2.upper_raw()) { return false; }
    else { return indeterminate; }
}

//! \related ValidatedFloat \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
inline Tribool operator>=(ValidatedFloat const& x1, ValidatedFloat const& x2) {
    if(x1.lower_raw()>=x2.upper_raw()) { return true; }
    else if(x1.upper_raw()< x2.lower_raw()) { return false; }
    else { return indeterminate; }
}

//! \related ValidatedFloat \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
inline Tribool operator<=(ValidatedFloat const& x1, ValidatedFloat const& x2) {
    if(x1.upper_raw()<=x2.lower_raw()) { return true; }
    else if(x1.lower_raw()> x2.upper_raw()) { return false; }
    else { return indeterminate; }
}

template<class N, EnableIf<IsIntegral<N>> =dummy> inline Tribool operator==(ValidatedFloat const& x1, N n2) { return x1==ValidatedFloat(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Tribool operator!=(ValidatedFloat const& x1, N n2) { return x1!=ValidatedFloat(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Tribool operator<=(ValidatedFloat const& x1, N n2) { return x1<=ValidatedFloat(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Tribool operator>=(ValidatedFloat const& x1, N n2) { return x1>=ValidatedFloat(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Tribool operator< (ValidatedFloat const& x1, N n2) { return x1< ValidatedFloat(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Tribool operator> (ValidatedFloat const& x1, N n2) { return x1> ValidatedFloat(n2); }

#ifdef ARIADNE_ENABLE_SERIALIZATION
  template<class A> Void serialize(A& a, ValidatedFloat& ivl, const Nat version) {
    a & ivl.lower_raw() & ivl.upper_raw(); }
#endif

OutputStream& operator<<(OutputStream&, ValidatedFloat const&);
InputStream& operator>>(InputStream&, ValidatedFloat&);


PositiveUpperFloat operator+(PositiveUpperFloat const& x1, PositiveUpperFloat const& x2);
PositiveUpperFloat operator-(PositiveUpperFloat const& x1, LowerFloat const& x2);
PositiveUpperFloat operator*(PositiveUpperFloat const& x1, PositiveUpperFloat const& x2);
PositiveUpperFloat operator/(PositiveUpperFloat const& x1, LowerFloat const& x2);
PositiveUpperFloat pow(PositiveUpperFloat const& x, Nat m);
PositiveUpperFloat half(PositiveUpperFloat const& x);


inline ErrorFloat operator"" _error(long double lx) { double x=lx; assert(x==lx); return ErrorFloat(Float(x)); }



} // namespace Ariadne

#endif
