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
//! \note In the future, the construction of _a \c %Float from a string literal may be supported.
//! \sa ExactInterval, Real, ExactFloat
template<class PR> class FloatTemplate<Approximate,PR> {
    typedef Approximate P; typedef FloatType<PR> FLT;
  public:
    typedef Approximate Paradigm;
    typedef ApproximateFloat NumericType;
  public:
    FloatTemplate<Approximate,PR>() : _a() { }
    template<class N, EnableIf<IsIntegral<N>> =dummy> FloatTemplate<Approximate,PR>(N n) : _a(n) { }
    template<class D, EnableIf<IsFloatingPoint<D>> =dummy> FloatTemplate<Approximate,PR>(D x) : _a(x) { }
    explicit FloatTemplate<Approximate,PR>(Float const& x) : _a(x) { }
    explicit FloatTemplate<Approximate,PR>(const Dyadic& d);
    explicit FloatTemplate<Approximate,PR>(const Decimal& d);

    explicit FloatTemplate<Approximate,PR>(const Integer& z);
    explicit FloatTemplate<Approximate,PR>(const Rational& q);
    explicit FloatTemplate<Approximate,PR>(const Real& r);
    explicit FloatTemplate<Approximate,PR>(const Number<Approximate>& x);
    operator Number<Approximate> () const;

    FloatTemplate<Approximate,PR>(FloatTemplate<Exact,PR> const& x);
    FloatTemplate<Approximate,PR>(FloatTemplate<Validated,PR> const& x);
    FloatTemplate<Approximate,PR>(FloatTemplate<Upper,PR> const& x);
    FloatTemplate<Approximate,PR>(FloatTemplate<Lower,PR> const& x);

    FloatTemplate<Approximate,PR> pm(FloatTemplate<Approximate,PR> e) { return *this; }

    explicit operator Float () const { return this->_a; }
    Float const& raw() const { return this->_a; }
    Float& raw() { return this->_a; }
    double get_d() const { return this->_a.get_d(); }
  public:
    static Void set_output_precision(Nat p) { output_precision=p; }
    static Nat get_output_precision() { return output_precision; }
    static Nat output_precision;
  private: public:
    Float _a;
};




//! \ingroup NumericModule
//! \brief Floating-point lower bounds for real numbers.
template<class PR> class FloatTemplate<Lower,PR> {
    typedef Lower P; typedef FloatType<PR> FLT;
  public:
    typedef Lower Paradigm;
    typedef FloatTemplate<Lower,PR> NumericType;
  public:
    FloatTemplate<Lower,PR>() : _l(0.0) { }
    template<class N, EnableIf<IsIntegral<N>> = dummy> FloatTemplate<Lower,PR>(N n) : _l(n) { }
    template<class X, EnableIf<IsFloatingPoint<X>> = dummy> explicit FloatTemplate<Lower,PR>(X x) : _l(x) { }
    explicit FloatTemplate<Lower,PR>(Float const& x) : _l(x) { }

    FloatTemplate<Lower,PR>(FloatTemplate<Validated,PR> const& x);
    FloatTemplate<Lower,PR>(FloatTemplate<Exact,PR> const& x);

    explicit FloatTemplate<Lower,PR>(const Number<Lower>& x);
    operator Number<Lower> () const;

    explicit FloatTemplate<Lower,PR>(const Real& x);
    explicit FloatTemplate<Lower,PR>(const Rational& x);
    explicit FloatTemplate<Lower,PR>(const Integer& x);

    Float const& raw() const { return _l; }
    Float& raw() { return _l; }
    double get_d() const { return _l.get_d(); }
  private: public:
    static Nat output_precision;
    Float _l;
};


//! \ingroup NumericModule
//! \brief Floating-point upper bounds for real numbers.
template<class PR> class FloatTemplate<Upper,PR> {
  public:
    typedef Upper Paradigm;
    typedef FloatTemplate<Upper,PR> NumericType;
  public:
    FloatTemplate<Upper,PR>() : _u(0.0) { }
    template<class N, EnableIf<IsIntegral<N>> = dummy> FloatTemplate<Upper,PR>(N n) : _u(n) { }
    template<class X, EnableIf<IsFloatingPoint<X>> = dummy> explicit FloatTemplate<Upper,PR>(X x) : _u(x) { }

    explicit FloatTemplate<Upper,PR>(Float const& x) : _u(x) { }

    FloatTemplate<Upper,PR>(FloatTemplate<Validated,PR> const& x);
    FloatTemplate<Upper,PR>(FloatTemplate<Exact,PR> const& x);

    explicit FloatTemplate<Upper,PR>(const Real& x);
    explicit FloatTemplate<Upper,PR>(const Rational& x);
    explicit FloatTemplate<Upper,PR>(const Integer& x);
    explicit FloatTemplate<Upper,PR>(const Number<Upper>& x);

    operator Number<Upper> () const;

    Float const& raw() const { return _u; }
    Float& raw() { return _u; }
    double get_d() const { return _u.get_d(); }
  private: public:
    static Nat output_precision;
    Float _u;
};


template<class PR> class FloatTemplate<Validated,PR> {
  public:
    typedef Validated Paradigm;
    typedef FloatTemplate<Validated,PR> NumericType;
  public:
    FloatTemplate<Validated,PR>() : _l(0.0), _u(0.0) { }
    template<class N, EnableIf<IsIntegral<N>> = dummy> FloatTemplate<Validated,PR>(N n) : _l(n), _u(n) { }
    template<class X, EnableIf<IsFloatingPoint<X>> = dummy> explicit FloatTemplate<Validated,PR>(X x) : _l(x), _u(x) { }
    explicit FloatTemplate<Validated,PR>(Float const& x) : _l(x), _u(x) { }

    FloatTemplate<Validated,PR>(FloatTemplate<Exact,PR> const& x);

    explicit FloatTemplate<Validated,PR>(const Dyadic& x);
    explicit FloatTemplate<Validated,PR>(const Decimal& x);
    explicit FloatTemplate<Validated,PR>(const Integer& z);
    explicit FloatTemplate<Validated,PR>(const Rational& q);
    explicit FloatTemplate<Validated,PR>(const Real& x);
    explicit FloatTemplate<Validated,PR>(const Number<Validated>& x);
    operator Number<Validated> () const;

    template<class N1, class N2, EnableIf<And<IsIntegral<N1>,IsIntegral<N2>>> =dummy>
        FloatTemplate<Validated,PR>(N1 lower, N2 upper) : _l(lower), _u(upper) { }
    FloatTemplate<Validated,PR>(Float const& lower, Float const& upper) : _l(lower), _u(upper) { }
    FloatTemplate<Validated,PR>(LowerFloat const& lower, UpperFloat const& upper) : _l(lower.raw()), _u(upper.raw()) { }
    FloatTemplate<Validated,PR>(const Rational& lower, const Rational& upper);

    Float const& lower_raw() const { return _l; }
    Float const& upper_raw() const { return _u; }
    Float const value_raw() const { return med_near(_l,_u); }
    Float const error_raw() const { return rad_up(_l,_u); }
    double get_d() const { return (_l.get_d()+_u.get_d())/2; }

    FloatTemplate<Lower,PR> lower() const { return FloatTemplate<Lower,PR>(_l); }
    FloatTemplate<Upper,PR> upper() const { return FloatTemplate<Upper,PR>(_u); }
    const FloatTemplate<Exact,PR> value() const;
    const FloatTemplate<PositiveUpper,PR> error() const;

    // DEPRECATED
    explicit operator Float () const { return (_l+_u)/2; }
    friend ExactFloat midpoint(ValidatedFloat const& x);
  public:
    static Nat output_precision;
    static Void set_output_precision(Nat p) { output_precision=p; }
    static Nat get_output_precision() { return output_precision; }
  private: public:
    Float _l, _u;
};


//! \ingroup NumericModule
//! \related Float, ValidatedFloat
//! \brief A floating-point number, which is taken to represent the \em exact value of a real quantity.
template<class PR> class FloatTemplate<Exact,PR> {
  public:
    typedef Exact Paradigm;
    typedef FloatTemplate<Exact,PR> NumericType;

    FloatTemplate<Exact,PR>() : _v(0) { }
    template<class N, EnableIf<IsIntegral<N>> =dummy> FloatTemplate<Exact,PR>(N n) : _v(n) { }
    template<class X, EnableIf<IsFloatingPoint<X>> =dummy> explicit FloatTemplate<Exact,PR>(X x) : _v(x) { }

    explicit FloatTemplate<Exact,PR>(Float const& x) : _v(x) { }

    explicit FloatTemplate<Exact,PR>(const Integer& z);
    explicit operator Rational () const;
    operator Number<Exact> () const;
    explicit operator Float () const { return _v; }

    Float const& raw() const { return _v; }
    Float& raw() { return _v; }
    double get_d() const { return _v.get_d(); }

    ValidatedFloat pm(ErrorFloat e) const;
  public:
    static Nat output_precision;
    static Void set_output_precision(Nat p) { output_precision=p; }
    static Nat get_output_precision() { return output_precision; }
  private: public:
    Float _v;
};


class PositiveExactFloat : public ExactFloat {
  public:
    PositiveExactFloat() : ExactFloat() { }
    template<class M, EnableIf<IsIntegral<M>>, EnableIf<IsUnsigned<M>> =dummy>
        PositiveExactFloat(M m) : ExactFloat(m) { }
    explicit PositiveExactFloat(Float const& x) : ExactFloat(x) { }
};

template<class PR> class FloatTemplate<PositiveUpper,PR> : public FloatTemplate<Upper,PR> {
    typedef FloatTemplate<Upper,PR> UpperFloat;
  public:
    FloatTemplate<PositiveUpper,PR>() : UpperFloat() { }
    explicit FloatTemplate<PositiveUpper,PR>(Float const& x) : UpperFloat(x) { ARIADNE_PRECONDITION(x>=0); }
    template<class M, EnableIf<IsUnsigned<M>> =dummy> FloatTemplate<PositiveUpper,PR>(M m) : UpperFloat(m) { }
    template<class F, EnableIf<IsSame<F,UpperFloat>> =dummy>
        explicit FloatTemplate<PositiveUpper,PR>(F const& x) : UpperFloat(x) { }
    FloatTemplate<PositiveUpper,PR>(PositiveExactFloat const& x) : UpperFloat(x) { }
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



template<class R, class A> R integer_cast(const A& _a);


inline ApproximateFloat floor(ApproximateFloat const& x) { return ApproximateFloat(floor(x._a)); }
inline ApproximateFloat ceil(ApproximateFloat const& x) { return ApproximateFloat(ceil(x._a)); }

inline ApproximateFloat abs(ApproximateFloat const& x) { return ApproximateFloat(abs_exact(x._a)); }
inline ApproximateFloat max(ApproximateFloat const& x, ApproximateFloat y) { return ApproximateFloat(max_exact(x._a,y._a)); }
inline ApproximateFloat min(ApproximateFloat const& x, ApproximateFloat y) { return ApproximateFloat(min_exact(x._a,y._a)); }

inline ApproximateFloat nul(ApproximateFloat const& x) { return ApproximateFloat(nul_exact(x._a)); }
inline ApproximateFloat pos(ApproximateFloat const& x) { return ApproximateFloat(pos_exact(x._a)); }
inline ApproximateFloat neg(ApproximateFloat const& x) { return ApproximateFloat(neg_exact(x._a)); }
inline ApproximateFloat half(ApproximateFloat const& x) { return ApproximateFloat(half_exact(x._a)); }
inline ApproximateFloat sqr(ApproximateFloat const& x) { return ApproximateFloat(mul_near(x._a,x._a)); }
inline ApproximateFloat rec(ApproximateFloat const& x) { return ApproximateFloat(div_near(1.0,x._a)); }

inline ApproximateFloat add(ApproximateFloat const& x1, ApproximateFloat const& x2) { return ApproximateFloat(add_near(x1._a,x2._a)); }
inline ApproximateFloat sub(ApproximateFloat const& x1, ApproximateFloat const& x2) { return ApproximateFloat(sub_near(x1._a,x2._a)); }
inline ApproximateFloat mul(ApproximateFloat const& x1, ApproximateFloat const& x2) { return ApproximateFloat(mul_near(x1._a,x2._a)); }
inline ApproximateFloat div(ApproximateFloat const& x1, ApproximateFloat const& x2) { return ApproximateFloat(div_near(x1._a,x2._a)); }

inline ApproximateFloat pow(ApproximateFloat const& x, Nat m) { return ApproximateFloat(pow_approx(x._a,m)); }
inline ApproximateFloat pow(ApproximateFloat const& x, Int n) { return ApproximateFloat(pow_approx(x._a,n)); }

inline ApproximateFloat sqrt(ApproximateFloat const& x) { return ApproximateFloat(sqrt_approx(x._a)); }
inline ApproximateFloat exp(ApproximateFloat const& x) { return ApproximateFloat(exp_approx(x._a)); }
inline ApproximateFloat log(ApproximateFloat const& x) { return ApproximateFloat(log_approx(x._a)); }
inline ApproximateFloat sin(ApproximateFloat const& x) { return ApproximateFloat(sin_approx(x._a)); }
inline ApproximateFloat cos(ApproximateFloat const& x) { return ApproximateFloat(cos_approx(x._a)); }
inline ApproximateFloat tan(ApproximateFloat const& x) { return ApproximateFloat(tan_approx(x._a)); }
inline ApproximateFloat asin(ApproximateFloat const& x) { return ApproximateFloat(asin_approx(x._a)); }
inline ApproximateFloat acos(ApproximateFloat const& x) { return ApproximateFloat(acos_approx(x._a)); }
inline ApproximateFloat atan(ApproximateFloat const& x) { return ApproximateFloat(atan_approx(x._a)); }

inline ApproximateFloat operator+(ApproximateFloat const& x) { return pos(x); }
inline ApproximateFloat operator-(ApproximateFloat const& x) { return neg(x); }
inline ApproximateFloat operator+(ApproximateFloat const& x1, ApproximateFloat const& x2) { return add(x1,x2); }
inline ApproximateFloat operator-(ApproximateFloat const& x1, ApproximateFloat const& x2) { return sub(x1,x2); }
inline ApproximateFloat operator*(ApproximateFloat const& x1, ApproximateFloat const& x2) { return mul(x1,x2); }
inline ApproximateFloat operator/(ApproximateFloat const& x1, ApproximateFloat const& x2) { return div(x1,x2); }
inline ApproximateFloat& operator+=(ApproximateFloat& x1, ApproximateFloat const& x2) { x1._a+=x2._a; return x1; }
inline ApproximateFloat& operator-=(ApproximateFloat& x1, ApproximateFloat const& x2) { x1._a-=x2._a; return x1; }
inline ApproximateFloat& operator*=(ApproximateFloat& x1, ApproximateFloat const& x2) { x1._a*=x2._a; return x1; }
inline ApproximateFloat& operator/=(ApproximateFloat& x1, ApproximateFloat const& x2) { x1._a/=x2._a; return x1; }

inline Bool operator==(ApproximateFloat const& x1, ApproximateFloat const& x2) { return x1._a==x2._a; }
inline Bool operator!=(ApproximateFloat const& x1, ApproximateFloat const& x2) { return x1._a!=x2._a; }
inline Bool operator<=(ApproximateFloat const& x1, ApproximateFloat const& x2) { return x1._a<=x2._a; }
inline Bool operator>=(ApproximateFloat const& x1, ApproximateFloat const& x2) { return x1._a>=x2._a; }
inline Bool operator< (ApproximateFloat const& x1, ApproximateFloat const& x2) { return x1._a< x2._a; }
inline Bool operator> (ApproximateFloat const& x1, ApproximateFloat const& x2) { return x1._a> x2._a; }

OutputStream& operator<<(OutputStream& os, ApproximateFloat const& x);
InputStream& operator>>(InputStream& is, ApproximateFloat& x);



inline LowerFloat max(LowerFloat const& x1, LowerFloat const& x2) { return LowerFloat(max_exact(x1._l,x2._l)); }
inline LowerFloat min(LowerFloat const& x1, LowerFloat const& x2) { return LowerFloat(min_exact(x1._l,x2._l)); }

inline LowerFloat nul(LowerFloat const& x) { return LowerFloat(pos_exact(x._l)); }
inline LowerFloat pos(LowerFloat const& x) { return LowerFloat(pos_exact(x._l)); }
inline LowerFloat neg(UpperFloat const& x) { return LowerFloat(neg_exact(x._u)); }
inline LowerFloat half(LowerFloat const& x) { return LowerFloat(half_exact(x._l)); }
LowerFloat sqr(LowerFloat const& x);
LowerFloat rec(UpperFloat const& x);

inline LowerFloat add(LowerFloat const& x1, LowerFloat const& x2) { return LowerFloat(add_down(x1._l,x2._l)); }
inline LowerFloat sub(LowerFloat const& x1, UpperFloat const& x2) { return LowerFloat(sub_down(x1._l,x2._u)); }
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

inline UpperFloat max(UpperFloat const& x1, UpperFloat const& x2) { return UpperFloat(max_exact(x1._u,x2._u)); }
inline UpperFloat min(UpperFloat const& x1, UpperFloat const& x2) { return UpperFloat(min_exact(x1._u,x2._u)); }

inline UpperFloat nul(UpperFloat const& x) { return UpperFloat(pos_exact(x._u)); }
inline UpperFloat pos(UpperFloat const& x) { return UpperFloat(pos_exact(x._u)); }
inline UpperFloat neg(LowerFloat const& x) { return UpperFloat(neg_exact(x._l)); }
inline UpperFloat half(UpperFloat const& x) { return UpperFloat(half_exact(x._u)); }
UpperFloat sqr(UpperFloat const& x);
UpperFloat rec(LowerFloat const& x);

inline UpperFloat add(UpperFloat const& x1, UpperFloat const& x2) { return UpperFloat(add_up(x1._u,x2._u)); }
inline UpperFloat sub(UpperFloat const& x1, LowerFloat const& x2) { return UpperFloat(sub_up(x1._u,x2._l)); }
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



extern const ExactFloat infty;
inline ExactFloat operator"" _exact(long double lx) { double x=lx; assert(x==lx); return ExactFloat(x); }
inline TwoExp::operator ExactFloat () const { return ExactFloat(this->get_d()); }

inline ExactFloat max(ExactFloat const& x1,  ExactFloat const& x2) { return ExactFloat(max(x1._v,x2._v)); }
inline ExactFloat min(ExactFloat const& x1,  ExactFloat const& x2) { return ExactFloat(min(x1._v,x2._v)); }
inline ExactFloat abs(ExactFloat const& x) { return ExactFloat(abs(x._v)); }

inline ExactFloat nul(ExactFloat const& x) { return ExactFloat(nul(x._v)); }
inline ExactFloat pos(ExactFloat const& x) { return ExactFloat(pos(x._v)); }
inline ExactFloat neg(ExactFloat const& x) { return ExactFloat(neg(x._v)); }
inline ExactFloat half(ExactFloat const& x) { return ExactFloat(half(x._v)); }

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



template<class PR> inline FloatTemplate<Approximate,PR>::FloatTemplate(FloatTemplate<Lower,PR> const& x) : _a(x.raw()) {
}

template<class PR> inline FloatTemplate<Approximate,PR>::FloatTemplate(FloatTemplate<Upper,PR> const& x) : _a(x.raw()) {
}

template<class PR> inline FloatTemplate<Approximate,PR>::FloatTemplate(FloatTemplate<Validated,PR> const& x)
    : _a(half_exact(add_near(x.lower_raw(),x.upper_raw()))) {
}

template<class PR> inline FloatTemplate<Approximate,PR>::FloatTemplate(FloatTemplate<Exact,PR> const& x) : _a(x.raw()) {
}

template<class PR> inline FloatTemplate<Lower,PR>::FloatTemplate(FloatTemplate<Validated,PR> const& x) : _l(x.lower_raw()) {
}

template<class PR> inline FloatTemplate<Lower,PR>::FloatTemplate(FloatTemplate<Exact,PR> const& x) : _l(x.raw()) {
}

template<class PR> inline FloatTemplate<Upper,PR>::FloatTemplate(FloatTemplate<Validated,PR> const& x) : _u(x.upper_raw()) {
}

template<class PR> inline FloatTemplate<Upper,PR>::FloatTemplate(FloatTemplate<Exact,PR> const& x) : _u(x.raw()) {
}

template<class PR> inline FloatTemplate<Validated,PR>::FloatTemplate(FloatTemplate<Exact,PR> const& x) : _l(x.raw()), _u(x.raw()) {
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





template<class PR> inline const FloatTemplate<Exact,PR> FloatTemplate<Validated,PR>::value() const {
    return FloatTemplate<Exact,PR>(med_near(this->_l,this->_u)); }

template<class PR> inline const FloatTemplate<PositiveUpper,PR> FloatTemplate<Validated,PR>::error() const {
    FloatType<PR> _v=med_near(this->_l,this->_u); return FloatTemplate<PositiveUpper,PR>(max(sub_up(this->_u,_v),sub_up(_v,this->_l))); }

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
  template<class A> Void serialize(A& _a, ValidatedFloat& ivl, const Nat version) {
    _a & ivl.lower_raw() & ivl.upper_raw(); }
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

inline ApproximateFloat create_float(Number<Approximate> const& x) { return ApproximateFloat(x); }
inline LowerFloat create_float(LowerFloat const& x) { return LowerFloat(x); }
inline UpperFloat create_float(Number<Upper> const& x) { return UpperFloat(x); }
inline ValidatedFloat create_float(Number<Validated> const& x) { return ValidatedFloat(x); }
inline ValidatedFloat create_float(Number<Effective> const& x) { return ValidatedFloat(x); }
inline ValidatedFloat create_float(Number<Exact> const& x) { return ValidatedFloat(x); }
inline ValidatedFloat create_float(Real const& x) { return ValidatedFloat(x); }

template<class X> struct IsGenericNumber : IsConvertible<X,Real> { };
template<> struct IsGenericNumber<Real> : True { };
template<class P> struct IsGenericNumber<Number<P>> : True { };

template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumber<Y>> =dummy> auto
operator+(X const& x, Y const& y) -> decltype(x+create_float(y)) { return x+create_float(y); }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumber<Y>> =dummy> auto
operator-(X const& x, Y const& y) -> decltype(x-create_float(y)) { return x-create_float(y); }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumber<Y>> =dummy> auto
operator*(X const& x, Y const& y) -> decltype(x*create_float(y)) { return x*create_float(y); }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumber<Y>> =dummy> auto
operator/(X const& x, Y const& y) -> decltype(x/create_float(y)) { return x/create_float(y); }

template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumber<Y>> =dummy> auto
operator+(Y const& y, X const& x) -> decltype(create_float(y,x.raw().precision())+x) { return create_float(y)+x; }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumber<Y>> =dummy> auto
operator-(Y const& y, X const& x) -> decltype(create_float(y)-x) { return create_float(y)-x; }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumber<Y>> =dummy> auto
operator*(Y const& y, X const& x) -> decltype(create_float(y)*x) { return create_float(y)*x; }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumber<Y>> =dummy> auto
operator/(Y const& y, X const& x) -> decltype(create_float(y)/x) { return create_float(y)/x; }


} // namespace Ariadne

#endif
