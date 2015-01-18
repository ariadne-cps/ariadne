/***************************************************************************
 *            numeric/float64.h
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
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file numeric/float64.h
 *  \brief
 */

#ifndef ARIADNE_FLOAT64_H
#define ARIADNE_FLOAT64_H

#include <cmath>
#include <iostream>

#include "utility/metaprogramming.h"
#include "utility/typedefs.h"
#include "numeric/paradigm.h"

#include "numeric/sign.h"
#include "numeric/is_number.h"

#include "numeric/logical.decl.h"
#include "numeric/float.decl.h"
#include "numeric/number.h"
#include "numeric/real.h"
#include "numeric/logical.h"
#include "numeric/flt64.h"

namespace Ariadne {

/************ Selection of floating-point arithmetic for Exact values ********/

struct IsExactTag { };
struct LowerUpperTag { };
struct ValueErrorTag { };

#ifdef USE_VALIDATED_FLOAT_ARITHMETIC
typedef ValueErrorTag DefaultExactArithmeticTag;
#else
typedef LowerUpperTag DefaultExactArithmeticTag;
#endif

/************ Precision class ************************************************/

struct Precision64 {
};

/************ Traits classes *************************************************/

template<class T> struct IsScalar;
template<class P> struct IsScalar<Float64Template<P>> : True { };

/************ Floating-point literals ****************************************/

ExactFloat64 operator"" _x(long double);
ErrorFloat64 operator"" _e(long double);
LowerFloat64 operator"" _l(long double);
UpperFloat64 operator"" _u(long double);
ApproximateFloat64 operator"" _a(long double);

/************ Creation *******************************************************/

class Float64Creator {
  public:
    template<class P> using Type = Float64Type<P>;
    static BoundedFloat64 create(ValidatedNumber, Exact);
    static BoundedFloat64 create(ValidatedNumber, Effective);
    static BoundedFloat64 create(ValidatedNumber, Validated);
    static BoundedFloat64 create(ValidatedNumber, Bounded);
    static MetricFloat64 create(ValidatedNumber, Metric);
    static UpperFloat64 create(UpperNumber, Upper);
    static LowerFloat64 create(LowerNumber, Lower);
    static ApproximateFloat64 create(ApproximateNumber, Approximate);

    template<class P> static auto create(Real r, P p) -> Type<Weaker<Effective,P>> {
        return create(ValidatedNumber(r),p); }
    template<class P> static auto create(unsigned int n, P p) -> Type<P>;
    template<class P> static auto create(int n, P p) -> Type<P>;
    template<class P> static auto create(double d, P p) -> Type<Approximate>;
};

/************ Dispatcher for mixed generic/concrete arithmetic ***************/

struct DispatchFloat64Arithmetic {
    template<class X> using IsNumber = IsConvertible<X,ApproximateNumber>;
    template<class X> using IsConcreteNumber = IsConvertible<X,Real>;
    template<class F> using IsFloatNumber = IsConvertible<F,ApproximateFloat64>;
    template<class X> using IsNonFloatNumber = And<IsNumber<X>,Not<IsFloatNumber<X>>>;
    template<class N> using IsGenericNumber = And<IsNumber<N>,Not<IsConcreteNumber<N>>,Not<IsFloatNumber<N>>>;
    typedef BoundedFloat64 ValidatedFloat64Type;

    // Find the type that a number N needs to be converted to to be added to a float of type F
    static ValidatedFloat64Type conversion_type(ValidatedNumber, ExactFloat64 dummy);
    static MetricFloat64 conversion_type(ValidatedNumber, MetricFloat64 dummy);
    static BoundedFloat64 conversion_type(ValidatedNumber, BoundedFloat64 dummy);
    static UpperFloat64 conversion_type(UpperNumber, UpperFloat64 dummy);
    static LowerFloat64 conversion_type(LowerNumber, LowerFloat64 dummy);
    static ApproximateFloat64 conversion_type(ApproximateNumber, ApproximateFloat64 dummy);
    template<class X, class F> using ConversionType = decltype(DispatchFloat64Arithmetic::conversion_type(declval<X>(),declval<F>()));
    template<class F, class X> static ConversionType<X,F> convert(X x) { return ConversionType<X,F>(x); }

    template<class F> using NegationType = decltype(neg(declval<F>()));

    template<class F, class N> static auto add_number(F f, N n) -> decltype(f+convert<F>(n)) { return f+convert<F>(n); }
    template<class F, class N> static auto add_number(N n, F f) -> decltype(convert<F>(n)+f) { return convert<F>(n)+f; }
    template<class F, class N> static auto sub_number(F f, N n) -> decltype(f-convert<NegationType<F>>(n)) { return f-convert<NegationType<F>>(n); }
    template<class F, class N> static auto sub_number(N n, F f) -> decltype(convert<NegationType<F>>(n)-f) { return convert<NegationType<F>>(n)-f; }
    template<class F, class N> static auto mul_number(F f, N n) -> decltype(f*convert<F>(n)) { return f*convert<F>(n); }
    template<class F, class N> static auto mul_number(N n, F f) -> decltype(convert<F>(n)*f) { return convert<F>(n)*f; }
    template<class F, class N> static auto div_number(F f, N n) -> decltype(f/convert<NegationType<F>>(n)) { return f/convert<NegationType<F>>(n); }
    template<class F, class N> static auto div_number(N n, F f) -> decltype(convert<NegationType<F>>(n)/f) { return convert<NegationType<F>>(n)/f; }

    // Mixed arithmetic with generic number types
    template<class P, class X, EnableIf<IsNonFloatNumber<X>> =dummy> friend auto operator+(Float64Template<P> f, X x) -> decltype(add_number(f,x)) { return add_number(f,x); }
    template<class P, class X, EnableIf<IsNonFloatNumber<X>> =dummy> friend auto operator+(X x, Float64Template<P> f) -> decltype(add_number(x,f)) { return add_number(x,f); }

    template<class P, class X, EnableIf<IsNonFloatNumber<X>> =dummy> friend auto operator-(Float64Template<P> f, X x) -> decltype(sub_number(f,x)) { return sub_number(f,x); }
    template<class P, class X, EnableIf<IsNonFloatNumber<X>> =dummy> friend auto operator-(X x, Float64Template<P> f) -> decltype(sub_number(x,f)) { return sub_number(x,f); }
    template<class P, class X, EnableIf<IsNonFloatNumber<X>> =dummy> friend auto operator*(Float64Template<P> f, X x) -> decltype(mul_number(f,x)) { return mul_number(f,x); }
    template<class P, class X, EnableIf<IsNonFloatNumber<X>> =dummy> friend auto operator*(X x, Float64Template<P> f) -> decltype(mul_number(x,f)) { return mul_number(x,f); }
    template<class P, class X, EnableIf<IsNonFloatNumber<X>> =dummy> friend auto operator/(Float64Template<P> f, X x) -> decltype(div_number(f,x)) { return div_number(f,x); }
    template<class P, class X, EnableIf<IsNonFloatNumber<X>> =dummy> friend auto operator/(X x, Float64Template<P> f) -> decltype(div_number(x,f)) { return div_number(x,f); }

/*
    // Mixed arithmetic with generic number types
    template<class P1, class P2> friend auto operator+(Float64Template<P1> f, Number<P2> n) -> decltype(add_number(f,n)) { return add_number(f,n); }
    template<class P1, class P2> friend auto operator+(Number<P1> n, Float64Template<P2> f) -> decltype(add_number(n,f)) { return add_number(n,f); }

    template<class P1, class P2> friend auto operator-(Float64Template<P1> f, Number<P2> n) -> decltype(sub_number(f,n)) { return sub_number(f,n); }
    template<class P1, class P2> friend auto operator-(Number<P1> n, Float64Template<P2> f) -> decltype(sub_number(n,f)) { return sub_number(n,f); }

    template<class P1, class P2> friend auto operator*(Float64Template<P1> f, Number<P2> n) -> decltype(mul_number(f,n)) { return mul_number(f,n); }
    template<class P1, class P2> friend auto operator*(Number<P1> n, Float64Template<P2> f) -> decltype(mul_number(n,f)) { return mul_number(n,f); }

    template<class P1, class P2> friend auto operator/(Float64Template<P1> f, Number<P2> n) -> decltype(div_number(f,n)) { return div_number(f,n); }
    template<class P1, class P2> friend auto operator/(Number<P1> n, Float64Template<P2> f) -> decltype(div_number(n,f)) { return div_number(n,f); }

    // Mixed arithmetic with concrete real number types
    template<class P> friend auto operator+(Float64Template<P> f, Real n) -> decltype(add_number(f,n)) { return add_number(f,n); }
    template<class P> friend auto operator+(Real n, Float64Template<P> f) -> decltype(add_number(n,f)) { return add_number(n,f); }

    template<class P> friend auto operator-(Float64Template<P> f, Real n) -> decltype(sub_number(f,n)) { return sub_number(f,n); }
    template<class P> friend auto operator-(Real n, Float64Template<P> f) -> decltype(sub_number(n,f)) { return sub_number(n,f); }

    template<class P> friend auto operator*(Float64Template<P> f, Real n) -> decltype(mul_number(f,n)) { return mul_number(f,n); }
    template<class P> friend auto operator*(Real n, Float64Template<P> f) -> decltype(mul_number(n,f)) { return mul_number(n,f); }

    template<class P> friend auto operator/(Float64Template<P> f, Real n) -> decltype(div_number(f,n)) { return div_number(f,n); }
    template<class P> friend auto operator/(Real n, Float64Template<P> f) -> decltype(div_number(n,f)) { return div_number(n,f); }
*/
    // Mixed arithmetic with builtin floating-point type
    using ApF = Float64Template<Approximate>;
    template<class P, class D, EnableIf<IsFloatingPoint<D>> =dummy> friend auto operator+(Float64Template<P> f, D d) -> decltype(add(f,ApF(d))) { return add(f,ApF(d)); }
    template<class P, class D, EnableIf<IsFloatingPoint<D>> =dummy> friend auto operator+(D d, Float64Template<P> f) -> decltype(add(ApF(d),f)) { return add(ApF(d),f); }

    template<class P, class D, EnableIf<IsFloatingPoint<D>> =dummy> friend auto operator-(Float64Template<P> f, D d) -> decltype(sub(f,ApF(d))) { return sub(f,ApF(d)); }
    template<class P, class D, EnableIf<IsFloatingPoint<D>> =dummy> friend auto operator-(D d, Float64Template<P> f) -> decltype(sub(ApF(d),f)) { return sub(ApF(d),f); }

    template<class P, class D, EnableIf<IsFloatingPoint<D>> =dummy> friend auto operator*(Float64Template<P> f, D d) -> decltype(mul(f,ApF(d))) { return mul(f,ApF(d)); }
    template<class P, class D, EnableIf<IsFloatingPoint<D>> =dummy> friend auto operator*(D d, Float64Template<P> f) -> decltype(mul(ApF(d),f)) { return mul(ApF(d),f); }

    template<class P, class D, EnableIf<IsFloatingPoint<D>> =dummy> friend auto operator/(Float64Template<P> f, D d) -> decltype(div(f,ApF(d))) { return div(f,ApF(d)); }
    template<class P, class D, EnableIf<IsFloatingPoint<D>> =dummy> friend auto operator/(D d, Float64Template<P> f) -> decltype(div(ApF(d),f)) { return div(ApF(d),f); }

};

/************ Mixed comparisons ***************************************************/



template<class X, class P, EnableIf<IsFloat<X>> =dummy> inline auto
    operator==(X x, Number<P> y) -> decltype(Number<Paradigm<X>>(x)==y) { return x==x.create(y); }
template<class X, class P, EnableIf<IsFloat<X>> =dummy> inline auto
    operator!=(X x, Number<P> y) -> decltype(Number<Paradigm<X>>(x)!=y) { return x!=x.create(y); }
template<class X, class P, EnableIf<IsFloat<X>> =dummy> inline auto
    operator< (X x, Number<P> y) -> decltype(Number<Paradigm<X>>(x)< y) { return x< x.create(y); }
template<class X, class P, EnableIf<IsFloat<X>> =dummy> inline auto
    operator> (X x, Number<P> y) -> decltype(Number<Paradigm<X>>(x)> y) { return x> x.create(y); }
template<class X, class P, EnableIf<IsFloat<X>> =dummy> inline auto
    operator<=(X x, Number<P> y) -> decltype(Number<Paradigm<X>>(x)<=y) { return x<=x.create(y); }
template<class X, class P, EnableIf<IsFloat<X>> =dummy> inline auto
    operator>=(X x, Number<P> y) -> decltype(Number<Paradigm<X>>(x)>=y) { return x>=x.create(y); }

template<class X, class P, EnableIf<IsFloat<X>> =dummy> inline auto
    operator==(Number<P> y, X x) -> decltype(y==Number<Paradigm<X>>(x)) { return x.create(y)==x; }
template<class X, class P, EnableIf<IsFloat<X>> =dummy> inline auto
    operator!=(Number<P> y, X x) -> decltype(y!=Number<Paradigm<X>>(x)) { return x.create(y)!=x; }
template<class X, class P, EnableIf<IsFloat<X>> =dummy> inline auto
    operator< (Number<P> y, X x) -> decltype(y< Number<Paradigm<X>>(x)) { return x.create(y)< x; }
template<class X, class P, EnableIf<IsFloat<X>> =dummy> inline auto
    operator> (Number<P> y, X x) -> decltype(y> Number<Paradigm<X>>(x)) { return x.create(y)> x; }
template<class X, class P, EnableIf<IsFloat<X>> =dummy> inline auto
    operator<=(Number<P> y, X x) -> decltype(y<=Number<Paradigm<X>>(x)) { return x.create(y)<=x; }
template<class X, class P, EnableIf<IsFloat<X>> =dummy> inline auto
    operator>=(Number<P> y, X x) -> decltype(y>=Number<Paradigm<X>>(x)) { return x.create(y)>=x; }

// Comparisons with concrete real type; don't allow conversions to prevent ambiguous overloads
template<class X, class R, EnableIf<And<IsFloat<X>,IsSame<R,Real>>> =dummy> inline auto
    operator==(X x, R r) -> decltype(Number<Paradigm<X>>()==Number<Paradigm<R>>()) { return x==x.create(r); }
template<class X, class R, EnableIf<And<IsFloat<X>,IsSame<R,Real>>> =dummy> inline auto
    operator!=(X x, R r) -> decltype(Number<Paradigm<X>>()!=Number<Paradigm<R>>()) { return x!=x.create(r); }
template<class X, class R, EnableIf<And<IsFloat<X>,IsSame<R,Real>>> =dummy> inline auto
    operator< (X x, R r) -> decltype(Number<Paradigm<X>>()< Number<Paradigm<R>>()) { return x< x.create(r); }
template<class X, class R, EnableIf<And<IsFloat<X>,IsSame<R,Real>>> =dummy> inline auto
    operator> (X x, R r) -> decltype(Number<Paradigm<X>>()> Number<Paradigm<R>>()) { return x> x.create(r); }
template<class X, class R, EnableIf<And<IsFloat<X>,IsSame<R,Real>>> =dummy> inline auto
    operator<=(X x, R r) -> decltype(Number<Paradigm<X>>()<=Number<Paradigm<R>>()) { return x<=x.create(r); }
template<class X, class R, EnableIf<And<IsFloat<X>,IsSame<R,Real>>> =dummy> inline auto
    operator>=(X x, R r) -> decltype(Number<Paradigm<X>>()>=Number<Paradigm<R>>()) { return x>=x.create(r); }

template<class X, class R, EnableIf<And<IsFloat<X>,IsSame<R,Real>>> =dummy> inline auto
    operator==(R r, X x) -> decltype(Number<Paradigm<R>>()==Number<Paradigm<X>>()) { return x.create(r)==x; }
template<class X, class R, EnableIf<And<IsFloat<X>,IsSame<R,Real>>> =dummy> inline auto
    operator!=(R r, X x) -> decltype(Number<Paradigm<R>>()!=Number<Paradigm<X>>()) { return x.create(r)!=x; }
template<class X, class R, EnableIf<And<IsFloat<X>,IsSame<R,Real>>> =dummy> inline auto
    operator< (R r, X x) -> decltype(Number<Paradigm<R>>()< Number<Paradigm<X>>()) { return x.create(r)< x; }
template<class X, class R, EnableIf<And<IsFloat<X>,IsSame<R,Real>>> =dummy> inline auto
    operator> (R r, X x) -> decltype(Number<Paradigm<R>>()> Number<Paradigm<X>>()) { return x.create(r)> x; }
template<class X, class R, EnableIf<And<IsFloat<X>,IsSame<R,Real>>> =dummy> inline auto
    operator<=(R r, X x) -> decltype(Number<Paradigm<R>>()<=Number<Paradigm<X>>()) { return x.create(r)<=x; }
template<class X, class R, EnableIf<And<IsFloat<X>,IsSame<R,Real>>> =dummy> inline auto
    operator>=(R r, X x) -> decltype(Number<Paradigm<R>>()>=Number<Paradigm<X>>()) { return x.create(r)>=x; }

// Reverse comparisons with concrete rational type

template<class X, EnableIf<IsFloat<X>> =dummy> inline auto
    operator==(Rational const& q, X x) -> decltype(x==q) { return x==q; }
template<class X, EnableIf<IsFloat<X>> =dummy> inline auto
    operator!=(Rational const& q, X x) -> decltype(x!=q) { return x!=q; }
template<class X, EnableIf<IsFloat<X>> =dummy> inline auto
    operator< (Rational const& q, X x) -> decltype(x> q) { return x> q; }
template<class X, EnableIf<IsFloat<X>> =dummy> inline auto
    operator> (Rational const& q, X x) -> decltype(x< q) { return x< q; }
template<class X, EnableIf<IsFloat<X>> =dummy> inline auto
    operator<=(Rational const& q, X x) -> decltype(x>=q) { return x>=q; }
template<class X, EnableIf<IsFloat<X>> =dummy> inline auto
    operator>=(Rational const& q, X x) -> decltype(x<=q) { return x<=q; }


// Comparisons with builtin integral types
template<class X, class N, EnableIf<And<IsFloat<X>,IsIntegral<N>>> =dummy> inline auto
    operator==(X x, N n) -> decltype(x==declval<ExactFloat64>()) { return x==(-x).create(n); }
template<class X, class N, EnableIf<And<IsFloat<X>,IsIntegral<N>>> =dummy> inline auto
    operator!=(X x, N n) -> decltype(x!=declval<ExactFloat64>()) { return x!=(-x).create(n); }
template<class X, class N, EnableIf<And<IsFloat<X>,IsIntegral<N>>> =dummy> inline auto
    operator< (X x, N n) -> decltype(x< declval<ExactFloat64>()) { return x< (-x).create(n); }
template<class X, class N, EnableIf<And<IsFloat<X>,IsIntegral<N>>> =dummy> inline auto
    operator> (X x, N n) -> decltype(x> declval<ExactFloat64>()) { return x> (-x).create(n); }
template<class X, class N, EnableIf<And<IsFloat<X>,IsIntegral<N>>> =dummy> inline auto
    operator<=(X x, N n) -> decltype(x<=declval<ExactFloat64>()) { return x<=(-x).create(n); }
template<class X, class N, EnableIf<And<IsFloat<X>,IsIntegral<N>>> =dummy> inline auto
    operator>=(X x, N n) -> decltype(x>=declval<ExactFloat64>()) { return x>=(-x).create(n); }

// Comparisons of ARBIRTRARY numbers with builtin floating-point types
template<class X, class D, EnableIf<And<IsNumber<X>,IsFloatingPoint<D>>> =dummy> inline auto
    operator==(X x, D d) -> decltype(declval<X>()==declval<ApproximateFloat64>()) { return x==Float64Template<Approximate>(d); }
template<class X, class D, EnableIf<And<IsNumber<X>,IsFloatingPoint<D>>> =dummy> inline auto
    operator!=(X x, D d) -> decltype(declval<X>()!=declval<ApproximateFloat64>()) { return x!=Float64Template<Approximate>(d); }
template<class X, class D, EnableIf<And<IsNumber<X>,IsFloatingPoint<D>>> =dummy> inline auto
    operator< (X x, D d) -> decltype(declval<X>()< declval<ApproximateFloat64>()) { return x< Float64Template<Approximate>(d); }
template<class X, class D, EnableIf<And<IsNumber<X>,IsFloatingPoint<D>>> =dummy> inline auto
    operator> (X x, D d) -> decltype(declval<X>()> declval<ApproximateFloat64>()) { return x> Float64Template<Approximate>(d); }
template<class X, class D, EnableIf<And<IsNumber<X>,IsFloatingPoint<D>>> =dummy> inline auto
    operator<=(X x, D d) -> decltype(declval<X>()<=declval<ApproximateFloat64>()) { return x<=Float64Template<Approximate>(d); }
template<class X, class D, EnableIf<And<IsNumber<X>,IsFloatingPoint<D>>> =dummy> inline auto
    operator>=(X x, D d) -> decltype(declval<X>()>=declval<ApproximateFloat64>()) { return x>=Float64Template<Approximate>(d); }

/************ ExactFloat64 ***************************************************/

extern const ExactFloat64 infty64;

//! \ingroup Float64SubModule
//! \brief Floating-point numbers representing exact values.
template<> class Float64Template<Exact>
    : DispatchFloat64Arithmetic
{
    friend class Float64Template<Bounded>;
    friend class Float64Template<Approximate>;
    volatile double _v;
    static Int output_precision;
  public:
    static Void set_output_precision(Int);
  public:
    typedef Exact Paradigm;

    Float64Template<Exact>();
    template<class N, EnableIf<IsIntegral<N>> =dummy> Float64Template<Exact>(N);
    explicit Float64Template<Exact>(Integer);
    explicit Float64Template<Exact>(Flt64);
    operator Number<Exact>() const;
    operator BoundedFloat64 () const;
    operator MetricFloat64 () const;
    operator LowerFloat64 () const;
    operator UpperFloat64 () const;
    operator ApproximateFloat64 () const;
    Flt64 get_flt() const;
    Flt64 raw() const;
    double get_d() const;

    template<class X> auto create(const X& x) -> decltype(Float64Creator::create(x,Exact())) {
        return Float64Creator::create(x,Exact()); }

    ExactFloat64 create_zero() const;
    MetricFloat64 pm(ErrorFloat64) const;

    friend MetricFloat64 add(ExactFloat64 x1, ExactFloat64 x2, ValueErrorTag);
    friend MetricFloat64 sub(ExactFloat64 x1, ExactFloat64 x2, ValueErrorTag);
    friend MetricFloat64 mul(ExactFloat64 x1, ExactFloat64 x2, ValueErrorTag);
    friend MetricFloat64 div(ExactFloat64 x1, ExactFloat64 x2, ValueErrorTag);

    friend BoundedFloat64 add(ExactFloat64 x1, ExactFloat64 x2, LowerUpperTag);
    friend BoundedFloat64 sub(ExactFloat64 x1, ExactFloat64 x2, LowerUpperTag);
    friend BoundedFloat64 mul(ExactFloat64 x1, ExactFloat64 x2, LowerUpperTag);
    friend BoundedFloat64 div(ExactFloat64 x1, ExactFloat64 x2, LowerUpperTag);

    friend BoundedFloat64 add(BoundedFloat64 x1, BoundedFloat64 x2);
    friend BoundedFloat64 sub(BoundedFloat64 x1, BoundedFloat64 x2);
    friend BoundedFloat64 mul(BoundedFloat64 x1, BoundedFloat64 x2);
    friend BoundedFloat64 div(BoundedFloat64 x1, BoundedFloat64 x2);


    friend ExactFloat64 operator+(ExactFloat64);
    friend ExactFloat64 operator-(ExactFloat64);
    friend ExactFloat64 pos(ExactFloat64);
    friend ExactFloat64 neg(ExactFloat64);
    friend ExactFloat64 half(ExactFloat64);
    friend BoundedFloat64 sqr(ExactFloat64);
    friend BoundedFloat64 rec(ExactFloat64);

    friend BoundedFloat64 pow(ExactFloat64, Nat);
    friend BoundedFloat64 pow(ExactFloat64, Int);
    friend BoundedFloat64 sqrt(ExactFloat64);
    friend BoundedFloat64 exp(ExactFloat64);
    friend BoundedFloat64 log(ExactFloat64);
    friend BoundedFloat64 sin(ExactFloat64);
    friend BoundedFloat64 cos(ExactFloat64);
    friend BoundedFloat64 tan(ExactFloat64);
    friend BoundedFloat64 atan(ExactFloat64);

    friend ExactFloat64 min(ExactFloat64 x1, ExactFloat64 x2);
    friend ExactFloat64 max(ExactFloat64 x1, ExactFloat64 x2);
    friend ExactFloat64 abs(ExactFloat64 x);
    friend ErrorFloat64 mag(ExactFloat64 x);

    friend Boolean operator==(ExactFloat64 x1, ExactFloat64 x2);
    friend Boolean operator!=(ExactFloat64 x1, ExactFloat64 x2);
    friend Boolean operator< (ExactFloat64 x1, ExactFloat64 x2);
    friend Boolean operator> (ExactFloat64 x1, ExactFloat64 x2);
    friend Boolean operator<=(ExactFloat64 x1, ExactFloat64 x2);
    friend Boolean operator>=(ExactFloat64 x1, ExactFloat64 x2);

    friend Bool same(ExactFloat64 x1, ExactFloat64 x2);

    friend OutputStream& operator<<(OutputStream& os, ExactFloat64 const&);

    friend ExactFloat64 operator*(ExactFloat64,TwoExp);
    friend ExactFloat64 operator/(ExactFloat64,TwoExp);
    friend ExactFloat64& operator*=(ExactFloat64&,TwoExp);
    friend ExactFloat64& operator/=(ExactFloat64&,TwoExp);

    // Mixed comparisons with Rational
    friend Boolean operator==(ExactFloat64, Rational const&);
    friend Boolean operator!=(ExactFloat64, Rational const&);
    friend Boolean operator< (ExactFloat64, Rational const&);
    friend Boolean operator> (ExactFloat64, Rational const&);
    friend Boolean operator<=(ExactFloat64, Rational const&);
    friend Boolean operator>=(ExactFloat64, Rational const&);

  private:
    Float64Template<Exact>(Int64);
/*
    friend Boolean operator==(ExactFloat64, Int64);
    friend Boolean operator!=(ExactFloat64, Int64);
    friend Boolean operator< (ExactFloat64, Int64);
    friend Boolean operator> (ExactFloat64, Int64);
    friend Boolean operator<=(ExactFloat64, Int64);
    friend Boolean operator>=(ExactFloat64, Int64);
*/
};

template<class N, EnableIf<IsIntegral<N>>> inline
Float64Template<Exact>::Float64Template(N n) : Float64Template<Exact>(Int64(n)) { assert(this->_v==n); }


/************ ErrorFloat64 ***************************************************/

//! \ingroup Float64SubModule
//! \brief Positive floating-point numbers representing an over-approximation to the error of some quantity.
template<> class Float64Template<Error> {
    volatile double _e;
    static Int output_precision;
  public:
    static Void set_output_precision(Int);
  public:
    typedef Upper Paradigm;
    Float64Template<Error>();
    template<class M, EnableIf<IsSame<M,Nat>> = dummy> Float64Template<Error>(M m);
    template<class X, EnableIf<IsSame<X,double>> = dummy> explicit Float64Template<Error>(X x);
    explicit Float64Template<Error>(double);
    explicit Float64Template<Error>(Flt64);
    explicit Float64Template<Error>(UpperFloat64);
    operator UpperFloat64 () const;
    Flt64 get_flt() const;
    Flt64 raw() const;
    double get_d() const;
    explicit operator double() const;
    friend UpperFloat64 operator+(ErrorFloat64);
    friend LowerFloat64 operator-(ErrorFloat64);
    friend ErrorFloat64 operator+(ErrorFloat64,ErrorFloat64);
    friend ErrorFloat64 operator*(ErrorFloat64,ErrorFloat64);
    friend ErrorFloat64 operator*(ErrorFloat64 x1, UpperFloat64 x2); // TODO: Should we allow this?
    friend ErrorFloat64 operator/(ErrorFloat64 x1, LowerFloat64 x2); // TODO: Should we allow this?
    friend ErrorFloat64& operator+=(ErrorFloat64&,ErrorFloat64);
    friend ErrorFloat64& operator+=(ErrorFloat64&,ExactFloat64);
    friend ErrorFloat64& operator*=(ErrorFloat64&,ErrorFloat64);
    friend UpperFloat64 operator+(UpperFloat64,ErrorFloat64);
    friend LowerFloat64 operator-(LowerFloat64,ErrorFloat64);
    friend ErrorFloat64 pos(ErrorFloat64);
    friend ErrorFloat64 half(ErrorFloat64);
    friend ErrorFloat64 add(ErrorFloat64,ErrorFloat64);
    friend ErrorFloat64 mul(ErrorFloat64,ErrorFloat64);
    friend ErrorFloat64 pow(ErrorFloat64,Nat);
    friend ErrorFloat64 min(ErrorFloat64 x1, ErrorFloat64 x2);
    friend ErrorFloat64 max(ErrorFloat64 x1, ErrorFloat64 x2);
    friend NegSierpinski operator==(UpperFloat64, LowerFloat64);
    friend NegSierpinski operator==(LowerFloat64, UpperFloat64);
    friend Sierpinski operator!=(UpperFloat64, LowerFloat64);
    friend Sierpinski operator!=(LowerFloat64, UpperFloat64);
    friend Sierpinski operator< (UpperFloat64, LowerFloat64);
    friend Sierpinski operator> (LowerFloat64, UpperFloat64);
    friend Sierpinski operator<=(UpperFloat64, LowerFloat64);
    friend Sierpinski operator>=(LowerFloat64, UpperFloat64);
    // Comparison with a raw float since no direct conversion to ApproximateFloat64
    friend Fuzzy operator==(ErrorFloat64, Flt64);
    friend Sierpinski operator<=(ErrorFloat64, Flt64);

    friend Bool same(ErrorFloat64, ErrorFloat64);
    friend Bool refines(ErrorFloat64 x1, ErrorFloat64 x2);

    friend OutputStream& operator<<(OutputStream& os, ErrorFloat64 const&);
};

template<class M, EnableIf<IsSame<M,Nat>>> inline
Float64Template<Error>::Float64Template(M m) : Float64Template<Error>(static_cast<double>(m)) { }

template<class M, EnableIf<IsUnsigned<M>> =dummy> inline
ErrorFloat64 operator/(ErrorFloat64 x1, M m2) {
    return Float64Template<Error>(x1.get_d()/m2); }



/************ MetricFloat64 ***************************************************/

//! \ingroup Float64SubModule
//! \brief Floating-point class representing an approximation to a number with an over-approximation of the error.
template<> class Float64Template<Metric>
    : DispatchFloat64Arithmetic
{
    volatile double _v; volatile double _e;
    static Int output_precision;
  public:
    static Void set_output_precision(Int);
  private:
    Float64Template<Metric>(double d, IsExactTag);
  public:
    typedef Validated Paradigm;
    typedef MetricFloat64 NumericType;

    Float64Template<Metric>();
    explicit Float64Template<Metric>(double);
    explicit Float64Template<Metric>(double,double);
    explicit Float64Template<Metric>(Flt64);
    Float64Template<Metric>(ExactFloat64, ErrorFloat64);
    Float64Template<Metric>(LowerFloat64, UpperFloat64) = delete;

    template<class N, EnableIf<IsIntegral<N>> =dummy> Float64Template<Metric>(N n);
    explicit Float64Template<Metric>(Integer const&);
    explicit Float64Template<Metric>(Rational const&);
    explicit Float64Template<Metric>(Real const&);

    template<class X> auto create(const X& x) -> decltype(Float64Creator::create(x,Validated())) {
        return Float64Creator::create(x,Validated()); };

    MetricFloat64 create_zero() const;
    MetricFloat64 create_constant(Number<Validated>) const;

    operator BoundedFloat64 () const;
    operator LowerFloat64 () const;
    operator UpperFloat64 () const;
    operator ApproximateFloat64 () const;

    operator Number<Validated>() const;

    ExactFloat64 value() const;
    ErrorFloat64 error() const;
    UpperFloat64 upper() const;
    LowerFloat64 lower() const;

    friend MetricFloat64 widen(MetricFloat64 x);

    friend MetricFloat64 operator+(MetricFloat64);
    friend MetricFloat64 operator-(MetricFloat64);
    friend MetricFloat64 operator+(MetricFloat64,MetricFloat64);
    friend MetricFloat64 operator-(MetricFloat64,MetricFloat64);
    friend MetricFloat64 operator*(MetricFloat64,MetricFloat64);
    friend MetricFloat64 operator/(MetricFloat64,MetricFloat64);
    friend MetricFloat64& operator+=(MetricFloat64&,MetricFloat64);
    friend MetricFloat64& operator-=(MetricFloat64&,MetricFloat64);
    friend MetricFloat64& operator*=(MetricFloat64&,MetricFloat64);
    friend MetricFloat64& operator/=(MetricFloat64&,MetricFloat64);

    friend MetricFloat64 add(MetricFloat64 x1, MetricFloat64 x2);
    friend MetricFloat64 sub(MetricFloat64 x1, MetricFloat64 x2);
    friend MetricFloat64 mul(MetricFloat64 x1, MetricFloat64 x2);
    friend MetricFloat64 div(MetricFloat64 x1, MetricFloat64 x2);
    friend MetricFloat64 pos(MetricFloat64);
    friend MetricFloat64 neg(MetricFloat64);
    friend MetricFloat64 sqr(MetricFloat64);
    friend MetricFloat64 rec(MetricFloat64 x);
    friend MetricFloat64 pow(MetricFloat64, Nat);
    friend MetricFloat64 pow(MetricFloat64, Int);
    friend MetricFloat64 min(MetricFloat64 x1, MetricFloat64 x2);
    friend MetricFloat64 max(MetricFloat64 x1, MetricFloat64 x2);
    friend MetricFloat64 abs(MetricFloat64 x);
    friend ErrorFloat64 mag(MetricFloat64);
    friend MetricFloat64 sqrt(MetricFloat64);
    friend MetricFloat64 exp(MetricFloat64);
    friend MetricFloat64 log(MetricFloat64);
    friend MetricFloat64 sin(MetricFloat64);
    friend MetricFloat64 cos(MetricFloat64);
    friend MetricFloat64 tan(MetricFloat64);
    friend MetricFloat64 atan(MetricFloat64);

    friend Tribool operator==(MetricFloat64,MetricFloat64);
    friend Tribool operator!=(MetricFloat64,MetricFloat64);
    friend Tribool operator<=(MetricFloat64,MetricFloat64);
    friend Tribool operator>=(MetricFloat64,MetricFloat64);
    friend Tribool operator< (MetricFloat64,MetricFloat64);
    friend Tribool operator> (MetricFloat64,MetricFloat64);

    // Mixed comparisons with Rational
    friend Tribool operator==(MetricFloat64, Rational const&);
    friend Tribool operator!=(MetricFloat64, Rational const&);
    friend Tribool operator< (MetricFloat64, Rational const&);
    friend Tribool operator> (MetricFloat64, Rational const&);
    friend Tribool operator<=(MetricFloat64, Rational const&);
    friend Tribool operator>=(MetricFloat64, Rational const&);

    friend ExactFloat64 estimate(MetricFloat64 x);
    friend ErrorFloat64 error(MetricFloat64 x);
    friend MetricFloat64 refinement(MetricFloat64 x1, MetricFloat64 x2);
    friend Bool same(MetricFloat64 x1, MetricFloat64 x2);
    friend Bool refines(MetricFloat64 x1, MetricFloat64 x2);
    friend Bool inconsistent(MetricFloat64 x1, MetricFloat64 x2);
    friend Bool represents(MetricFloat64 x, ExactFloat64 y);

    friend OutputStream& operator<<(OutputStream& os, MetricFloat64 const&);

};

template<class N, EnableIf<IsIntegral<N>>> inline
Float64Template<Metric>::Float64Template(N n) : Float64Template<Metric>(double(n),IsExactTag()) { }



/************ BoundedFloat64 ***************************************************/

//! \ingroup Float64SubModule
//! \brief Class providing floating-point lower and upper bounds for a number.
template<> class Float64Template<Bounded>
    : DispatchFloat64Arithmetic
    , public ProvideConvertedOperations<BoundedFloat64,ExactFloat64,Return<BoundedFloat64>>
    , public ProvideConvertedOperations<BoundedFloat64,MetricFloat64,Return<BoundedFloat64>>
    , public ProvideConvertedOperations<BoundedFloat64,Int,Return<BoundedFloat64>>
{
    friend class Float64Template<Approximate>; friend class Float64Template<Metric>;
    volatile double _l; volatile double _u;
    static Int output_precision;
  public:
    static Void set_output_precision(Int);
    static Int get_output_precision();
  public:
    typedef Validated Paradigm;
    typedef BoundedFloat64 NumericType;
    Float64Template<Bounded>();
    explicit Float64Template<Bounded>(double);
    explicit Float64Template<Bounded>(double,double);
    explicit Float64Template<Bounded>(Flt64);
    explicit Float64Template<Bounded>(Flt64,Flt64);
    Float64Template<Bounded>(LowerFloat64, UpperFloat64);
    Float64Template<Bounded>(LowerFloat64, ErrorFloat64) = delete;

    template<class N, EnableIf<IsIntegral<N>> =dummy>
        Float64Template<Bounded>(N n) : Float64Template<Bounded>(double(n),Exact()) { }
    template<class X, EnableIf<IsFloatingPoint<X>> =dummy>
        explicit Float64Template<Bounded>(X x) : Float64Template<Bounded>(double(x),Exact()) { }

    explicit Float64Template<Bounded>(Integer const&);
    explicit Float64Template<Bounded>(Rational const&);
    explicit Float64Template<Bounded>(Real const&);

    operator ApproximateFloat64 () const;
    operator UpperFloat64 () const;
    operator LowerFloat64 () const;
//    explicit operator MetricFloat64 () const;
    operator MetricFloat64 () const;

    operator Number<Validated>() const;

    template<class X> auto create(const X& x) -> decltype(Float64Creator::create(x,Bounded())) {
        return Float64Creator::create(x,Bounded()); };

    BoundedFloat64 create_zero() const;
    BoundedFloat64 create_constant(Number<Validated>) const;

    LowerFloat64 lower() const;
    UpperFloat64 upper() const;
    ExactFloat64 value() const;
    ErrorFloat64 error() const;
    ErrorFloat64 width() const;
    ErrorFloat64 radius() const;

    RawFloat64 upper_raw() const;
    RawFloat64 lower_raw() const;
    double get_d() const;

    friend BoundedFloat64 widen(BoundedFloat64 x);

    friend BoundedFloat64 operator+(BoundedFloat64);
    friend BoundedFloat64 operator-(BoundedFloat64);
    friend BoundedFloat64 operator+(BoundedFloat64,BoundedFloat64);
    friend BoundedFloat64 operator-(BoundedFloat64,BoundedFloat64);
    friend BoundedFloat64 operator*(BoundedFloat64,BoundedFloat64);
    friend BoundedFloat64 operator/(BoundedFloat64,BoundedFloat64);
    friend BoundedFloat64& operator+=(BoundedFloat64&,BoundedFloat64);
    friend BoundedFloat64& operator-=(BoundedFloat64&,BoundedFloat64);
    friend BoundedFloat64& operator*=(BoundedFloat64&,BoundedFloat64);
    friend BoundedFloat64& operator/=(BoundedFloat64&,BoundedFloat64);

    friend BoundedFloat64 add(BoundedFloat64,BoundedFloat64);
    friend BoundedFloat64 sub(BoundedFloat64,BoundedFloat64);
    friend BoundedFloat64 mul(BoundedFloat64,BoundedFloat64);
    friend BoundedFloat64 div(BoundedFloat64,BoundedFloat64);
    friend BoundedFloat64 half(BoundedFloat64);
    friend BoundedFloat64 pos(BoundedFloat64);
    friend BoundedFloat64 neg(BoundedFloat64);
    friend BoundedFloat64 rec(BoundedFloat64);
    friend BoundedFloat64 sqr(BoundedFloat64);
    friend BoundedFloat64 pow(BoundedFloat64, Nat);
    friend BoundedFloat64 pow(BoundedFloat64, Int);
    friend BoundedFloat64 max(BoundedFloat64,BoundedFloat64);
    friend BoundedFloat64 min(BoundedFloat64,BoundedFloat64);
    friend BoundedFloat64 abs(BoundedFloat64);
    friend LowerFloat64 mig(BoundedFloat64);
    friend ErrorFloat64 mag(BoundedFloat64);
    friend BoundedFloat64 sqrt(BoundedFloat64);
    friend BoundedFloat64 exp(BoundedFloat64);
    friend BoundedFloat64 log(BoundedFloat64);
    friend BoundedFloat64 sin(BoundedFloat64);
    friend BoundedFloat64 cos(BoundedFloat64);
    friend BoundedFloat64 tan(BoundedFloat64);
    friend BoundedFloat64 atan(BoundedFloat64);

    friend BoundedFloat64 add(MetricFloat64,BoundedFloat64);
    friend BoundedFloat64 add(BoundedFloat64,MetricFloat64);
    friend BoundedFloat64 sub(MetricFloat64,BoundedFloat64);
    friend BoundedFloat64 sub(BoundedFloat64,MetricFloat64);
    friend BoundedFloat64 mul(MetricFloat64,BoundedFloat64);
    friend BoundedFloat64 mul(BoundedFloat64,MetricFloat64);
    friend BoundedFloat64 div(MetricFloat64,BoundedFloat64);
    friend BoundedFloat64 div(BoundedFloat64,MetricFloat64);

    friend BoundedFloat64 operator+(MetricFloat64,BoundedFloat64);
    friend BoundedFloat64 operator+(BoundedFloat64,MetricFloat64);
    friend BoundedFloat64 operator-(MetricFloat64,BoundedFloat64);
    friend BoundedFloat64 operator-(BoundedFloat64,MetricFloat64);
    friend BoundedFloat64 operator*(MetricFloat64,BoundedFloat64);
    friend BoundedFloat64 operator*(BoundedFloat64,MetricFloat64);
    friend BoundedFloat64 operator/(MetricFloat64,BoundedFloat64);
    friend BoundedFloat64 operator/(BoundedFloat64,MetricFloat64);

    friend Tribool operator==(BoundedFloat64, BoundedFloat64);
    friend Tribool operator!=(BoundedFloat64, BoundedFloat64);
    friend Tribool operator< (BoundedFloat64, BoundedFloat64);
    friend Tribool operator> (BoundedFloat64, BoundedFloat64);
    friend Tribool operator<=(BoundedFloat64, BoundedFloat64);
    friend Tribool operator>=(BoundedFloat64, BoundedFloat64);

    // Mixed comparisons with Rational
    friend Tribool operator==(BoundedFloat64, Rational const&);
    friend Tribool operator!=(BoundedFloat64, Rational const&);
    friend Tribool operator< (BoundedFloat64, Rational const&);
    friend Tribool operator> (BoundedFloat64, Rational const&);
    friend Tribool operator<=(BoundedFloat64, Rational const&);
    friend Tribool operator>=(BoundedFloat64, Rational const&);

    friend ExactFloat64 estimate(BoundedFloat64 x);
    friend ErrorFloat64 error(BoundedFloat64 x);
    friend BoundedFloat64 refinement(BoundedFloat64 x1, BoundedFloat64 x2);
    friend Bool same(BoundedFloat64 x1, BoundedFloat64 x2);
    friend Bool refines(BoundedFloat64 x1, BoundedFloat64 x2);
    friend Bool inconsistent(BoundedFloat64 x1, BoundedFloat64 x2);
    friend Bool represents(BoundedFloat64 x, ExactFloat64 y);

    friend OutputStream& operator<<(OutputStream& os, BoundedFloat64 const&);

  private:
    Float64Template<Bounded>(double,Exact);

};


/************ UpperFloat64 ***************************************************/

//! \ingroup Float64SubModule
//! \brief Floating-point value specifying an upper bound for a number.
template<> class Float64Template<Upper>
    : DispatchFloat64Arithmetic
{
    volatile double _u;
  public:
    typedef Upper Paradigm;

    Float64Template<Upper>();

    template<class N, EnableIf<IsIntegral<N>> =dummy>
        Float64Template<Upper>(N n) : Float64Template<Upper>(double(n),Exact()) { }
    template<class X, EnableIf<IsFloatingPoint<X>> =dummy>
        explicit Float64Template<Upper>(X x) : Float64Template<Upper>(double(x),Exact()) { }

    explicit Float64Template<Upper>(Flt64);
    explicit Float64Template<Upper>(Integer const&);
    explicit Float64Template<Upper>(Rational const&);
    explicit Float64Template<Upper>(Real const&);

    template<class X> auto create(const X& x) -> decltype(Float64Creator::create(x,Upper())) {
        return Float64Creator::create(x,Upper()); }

    operator ApproximateFloat64 () const;
    operator Number<Upper>() const;
    Flt64 raw() const;
    Flt64 get_flt() const;
    double get_d() const;

    friend UpperFloat64 operator+(UpperFloat64);
    friend UpperFloat64 operator-(LowerFloat64);
    friend LowerFloat64 operator-(UpperFloat64);
    friend UpperFloat64 operator+(UpperFloat64,UpperFloat64);
    friend UpperFloat64 operator+(UpperFloat64,ErrorFloat64);
    friend UpperFloat64 operator-(UpperFloat64,LowerFloat64);
    friend LowerFloat64 operator-(LowerFloat64,UpperFloat64);
    friend UpperFloat64 operator*(UpperFloat64,UpperFloat64); // TODO: Should we provide this operator?
    friend UpperFloat64 operator/(UpperFloat64,LowerFloat64); // TODO: Should we provide this operator?
    friend LowerFloat64 operator/(LowerFloat64,UpperFloat64); // TODO: Should we provide this operator?
    friend ApproximateFloat64 operator*(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 operator/(ApproximateFloat64,ApproximateFloat64);
    friend UpperFloat64& operator+=(UpperFloat64&,UpperFloat64);
    friend UpperFloat64& operator-=(UpperFloat64&,LowerFloat64);
    friend UpperFloat64 add(UpperFloat64,UpperFloat64);
    friend UpperFloat64 sub(UpperFloat64,LowerFloat64);
    friend LowerFloat64 sub(LowerFloat64,UpperFloat64);
    friend UpperFloat64 mul(UpperFloat64,UpperFloat64); // TODO: Should we provide this operator?
    friend UpperFloat64 div(UpperFloat64,LowerFloat64); // TODO: Should we provide this operator?
    friend LowerFloat64 div(LowerFloat64,UpperFloat64); // TODO: Should we provide this operator?
    friend UpperFloat64 pos(UpperFloat64);
    friend LowerFloat64 neg(UpperFloat64);
    friend UpperFloat64 neg(LowerFloat64);
    friend UpperFloat64 half(UpperFloat64);
    friend UpperFloat64 sqr(UpperFloat64);
    friend LowerFloat64 rec(UpperFloat64);
    friend UpperFloat64 sqrt(UpperFloat64);
    friend UpperFloat64 exp(UpperFloat64);
    friend UpperFloat64 log(UpperFloat64);
    friend UpperFloat64 atan(UpperFloat64);
    friend UpperFloat64 min(UpperFloat64,UpperFloat64);
    friend UpperFloat64 max(UpperFloat64,UpperFloat64);

    friend ApproximateFloat64 add(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 sub(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 mul(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 div(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 pow(ApproximateFloat64,Int);
    friend ApproximateFloat64 sin(ApproximateFloat64);
    friend ApproximateFloat64 cos(ApproximateFloat64);
    friend ApproximateFloat64 tan(ApproximateFloat64);
    friend ApproximateFloat64 asin(ApproximateFloat64);
    friend ApproximateFloat64 acos(ApproximateFloat64);
    friend ApproximateFloat64 atan(ApproximateFloat64);

    friend OutputStream& operator<<(OutputStream& os, UpperFloat64 const&);

    friend NegSierpinski operator==(UpperFloat64, LowerFloat64);
    friend NegSierpinski operator==(LowerFloat64, UpperFloat64);
    friend Sierpinski operator!=(UpperFloat64, LowerFloat64);
    friend Sierpinski operator!=(LowerFloat64, UpperFloat64);

    friend Sierpinski operator< (UpperFloat64, LowerFloat64);
    friend Sierpinski operator> (LowerFloat64, UpperFloat64);
    friend Sierpinski operator<=(UpperFloat64, LowerFloat64);
    friend Sierpinski operator>=(LowerFloat64, UpperFloat64);
    friend NegSierpinski operator< (LowerFloat64, UpperFloat64);
    friend NegSierpinski operator> (UpperFloat64, LowerFloat64);
    friend NegSierpinski operator<=(LowerFloat64, UpperFloat64);
    friend NegSierpinski operator>=(UpperFloat64, LowerFloat64);

    friend Fuzzy operator==(ApproximateFloat64, ApproximateFloat64);
    friend Fuzzy operator!=(ApproximateFloat64, ApproximateFloat64);
    friend Fuzzy operator< (ApproximateFloat64, ApproximateFloat64);
    friend Fuzzy operator> (ApproximateFloat64, ApproximateFloat64);
    friend Fuzzy operator<=(ApproximateFloat64, ApproximateFloat64);
    friend Fuzzy operator>=(ApproximateFloat64, ApproximateFloat64);

    // Fallback for UpperFloat+LowerFloat, LowerFloat+UpperFloat and UpperFloat-UpperFloat
    friend ApproximateFloat64 operator+(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 operator-(ApproximateFloat64,ApproximateFloat64);


    friend NegSierpinski operator==(UpperFloat64, Rational const&);
    friend Sierpinski operator!=(UpperFloat64, Rational const&);
    friend Sierpinski operator< (UpperFloat64, Rational const&);
    friend NegSierpinski operator> (UpperFloat64, Rational const&);
    friend Sierpinski operator<=(UpperFloat64, Rational const&);
    friend NegSierpinski operator>=(UpperFloat64, Rational const&);

    friend ApproximateFloat64 min(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 max(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 abs(ApproximateFloat64);
  private:
    Float64Template<Upper>(double,Exact);
};

/************ LowerFloat64 ***************************************************/

//! \ingroup Float64SubModule
//! \brief Floating-point value specifying a lower bound for a number.
template<> class Float64Template<Lower>
    : DispatchFloat64Arithmetic
{
    volatile double _l;
  public:
    typedef Lower Paradigm;
    Float64Template<Lower>();

    template<class N, EnableIf<IsIntegral<N>> =dummy>
        Float64Template<Lower>(N n) : Float64Template<Lower>(double(n),Exact()) { }
    template<class X, EnableIf<IsFloatingPoint<X>> =dummy>
        explicit Float64Template<Lower>(X x) : Float64Template<Lower>(double(x),Exact()) { }
    explicit Float64Template<Lower>(Flt64);
    explicit Float64Template<Lower>(Integer const&);
    explicit Float64Template<Lower>(Rational const&);
    explicit Float64Template<Lower>(Real const&);

    template<class X> auto create(const X& x) -> decltype(Float64Creator::create(x,Lower())) {
        return Float64Creator::create(x,Lower()); }

    operator ApproximateFloat64 () const;
    operator Number<Lower>() const;
    Flt64 raw() const;
    Flt64 get_flt() const;
    double get_d() const;

    friend LowerFloat64 operator+(LowerFloat64);
    friend UpperFloat64 operator-(LowerFloat64);
    friend LowerFloat64 operator-(UpperFloat64);
    friend LowerFloat64 operator+(LowerFloat64,LowerFloat64);
    friend LowerFloat64 operator-(LowerFloat64,ErrorFloat64);
    friend LowerFloat64 operator-(LowerFloat64,UpperFloat64);
    friend UpperFloat64 operator-(UpperFloat64,LowerFloat64);
    friend LowerFloat64 operator*(LowerFloat64,LowerFloat64); // TODO: Should we provide this operator?
    friend LowerFloat64 operator/(LowerFloat64,UpperFloat64); // TODO: Should we provide this operator?
    friend UpperFloat64 operator/(UpperFloat64,LowerFloat64); // TODO: Should we provide this operator?
    friend ApproximateFloat64 operator*(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 operator/(ApproximateFloat64,ApproximateFloat64);
    friend LowerFloat64& operator+=(LowerFloat64&,LowerFloat64);
    friend LowerFloat64& operator-=(LowerFloat64&,UpperFloat64);
    friend LowerFloat64 add(LowerFloat64,LowerFloat64);
    friend LowerFloat64 sub(LowerFloat64,UpperFloat64);
    friend UpperFloat64 sub(UpperFloat64,LowerFloat64);
    friend LowerFloat64 mul(LowerFloat64,LowerFloat64); // TODO: Should we provide this operator?
    friend LowerFloat64 div(LowerFloat64,UpperFloat64); // TODO: Should we provide this operator?
    friend UpperFloat64 div(UpperFloat64,LowerFloat64); // TODO: Should we provide this operator?
    friend LowerFloat64 pos(LowerFloat64);
    friend UpperFloat64 neg(LowerFloat64);
    friend LowerFloat64 neg(UpperFloat64);
    friend LowerFloat64 half(LowerFloat64);
    friend LowerFloat64 sqr(LowerFloat64);
    friend UpperFloat64 rec(LowerFloat64);
    friend LowerFloat64 sqrt(LowerFloat64);
    friend LowerFloat64 exp(LowerFloat64);
    friend LowerFloat64 log(LowerFloat64);
    friend LowerFloat64 atan(LowerFloat64);
    friend LowerFloat64 min(LowerFloat64,LowerFloat64);
    friend LowerFloat64 max(LowerFloat64,LowerFloat64);
    friend OutputStream& operator<<(OutputStream& os, LowerFloat64 const&);

    friend ApproximateFloat64 add(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 sub(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 mul(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 div(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 pow(ApproximateFloat64,Int);
    friend ApproximateFloat64 sin(ApproximateFloat64);
    friend ApproximateFloat64 cos(ApproximateFloat64);
    friend ApproximateFloat64 tan(ApproximateFloat64);

    friend NegSierpinski operator==(UpperFloat64, LowerFloat64);
    friend NegSierpinski operator==(LowerFloat64, UpperFloat64);
    friend Sierpinski operator!=(UpperFloat64, LowerFloat64);
    friend Sierpinski operator!=(LowerFloat64, UpperFloat64);

    friend Sierpinski operator< (UpperFloat64, LowerFloat64);
    friend Sierpinski operator> (LowerFloat64, UpperFloat64);
    friend Sierpinski operator<=(UpperFloat64, LowerFloat64);
    friend Sierpinski operator>=(LowerFloat64, UpperFloat64);
    friend NegSierpinski operator< (LowerFloat64, UpperFloat64);
    friend NegSierpinski operator> (UpperFloat64, LowerFloat64);
    friend NegSierpinski operator<=(LowerFloat64, UpperFloat64);
    friend NegSierpinski operator>=(UpperFloat64, LowerFloat64);

    friend Fuzzy operator==(ApproximateFloat64, ApproximateFloat64);
    friend Fuzzy operator!=(ApproximateFloat64, ApproximateFloat64);
    friend Fuzzy operator< (ApproximateFloat64, ApproximateFloat64);
    friend Fuzzy operator> (ApproximateFloat64, ApproximateFloat64);
    friend Fuzzy operator<=(ApproximateFloat64, ApproximateFloat64);
    friend Fuzzy operator>=(ApproximateFloat64, ApproximateFloat64);

    friend NegSierpinski operator==(LowerFloat64, Rational const&);
    friend Sierpinski operator!=(LowerFloat64, Rational const&);
    friend NegSierpinski operator< (LowerFloat64, Rational const&);
    friend Sierpinski operator> (LowerFloat64, Rational const&);
    friend NegSierpinski operator<=(LowerFloat64, Rational const&);
    friend Sierpinski operator>=(LowerFloat64, Rational const&);

    // Fallback for UpperFloat+LowerFloat, LowerFloat+UpperFloat and LowerFloat-LowerFloat
    friend ApproximateFloat64 operator+(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 operator-(ApproximateFloat64,ApproximateFloat64);

    friend ApproximateFloat64 min(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 max(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 abs(ApproximateFloat64);
  private:
    Float64Template<Lower>(double,Exact);
};

ApproximateFloat64 asin(ApproximateFloat64);
ApproximateFloat64 acos(ApproximateFloat64);
ApproximateFloat64 atan(ApproximateFloat64);

/************ ApproximateFloat64 ***************************************************/

//! \ingroup Float64SubModule
//! \brief Floating-point value specifying an approximation for a number, with no guarantees on the quality.
template<> class Float64Template<Approximate>
    : DispatchFloat64Arithmetic
{
    static Int output_precision;
    volatile double _a;
  public:
    static Void set_output_precision(Int);
  public:
    typedef Approximate Paradigm;
    typedef ApproximateFloat64 NumericType;
    Float64Template<Approximate>();
    Float64Template<Approximate>(double);
    Float64Template<Approximate>(Flt64);
    explicit Float64Template<Approximate>(Integer const&);
    explicit Float64Template<Approximate>(Rational const&);
    explicit Float64Template<Approximate>(Real const&);
    RawFloat64 raw() const;

    ApproximateFloat64 create_zero() const;
    ApproximateFloat64 create_constant(Number<Approximate>) const;

    template<class X> auto create(const X& x) -> decltype(Float64Creator::create(x,Approximate())) {
        return Float64Creator::create(x,Approximate()); }

    operator Number<Approximate>() const;

    Flt64 get_flt() const;
    double get_d() const;
    friend ApproximateFloat64 operator+(ApproximateFloat64);
    friend ApproximateFloat64 operator-(ApproximateFloat64);
    friend ApproximateFloat64 operator+(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 operator-(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 operator*(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 operator/(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64& operator+=(ApproximateFloat64&,ApproximateFloat64);
    friend ApproximateFloat64& operator-=(ApproximateFloat64&,ApproximateFloat64);
    friend ApproximateFloat64& operator*=(ApproximateFloat64&,ApproximateFloat64);
    friend ApproximateFloat64& operator/=(ApproximateFloat64&,ApproximateFloat64);
    friend ApproximateFloat64 add(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 sub(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 mul(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 div(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 pos(ApproximateFloat64);
    friend ApproximateFloat64 neg(ApproximateFloat64);
    friend ApproximateFloat64 half(ApproximateFloat64);
    friend ApproximateFloat64 sqr(ApproximateFloat64);
    friend ApproximateFloat64 rec(ApproximateFloat64);
    friend ApproximateFloat64 pow(ApproximateFloat64,Nat);
    friend ApproximateFloat64 pow(ApproximateFloat64,Int);
    friend ApproximateFloat64 min(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 max(ApproximateFloat64,ApproximateFloat64);
    friend ApproximateFloat64 abs(ApproximateFloat64);
    friend ApproximateFloat64 mag(ApproximateFloat64 x);
    friend ApproximateFloat64 sqrt(ApproximateFloat64);
    friend ApproximateFloat64 exp(ApproximateFloat64);
    friend ApproximateFloat64 log(ApproximateFloat64);
    friend ApproximateFloat64 sin(ApproximateFloat64);
    friend ApproximateFloat64 cos(ApproximateFloat64);
    friend ApproximateFloat64 tan(ApproximateFloat64);
    friend ApproximateFloat64 atan(ApproximateFloat64);
    friend Fuzzy operator==(ApproximateFloat64,ApproximateFloat64);
    friend Fuzzy operator!=(ApproximateFloat64,ApproximateFloat64);
    friend Fuzzy operator<=(ApproximateFloat64,ApproximateFloat64);
    friend Fuzzy operator>=(ApproximateFloat64,ApproximateFloat64);
    friend Fuzzy operator< (ApproximateFloat64,ApproximateFloat64);
    friend Fuzzy operator> (ApproximateFloat64,ApproximateFloat64);
    friend Fuzzy operator==(ApproximateFloat64, Rational const&);
    friend Fuzzy operator!=(ApproximateFloat64, Rational const&);
    friend Fuzzy operator< (ApproximateFloat64, Rational const&);
    friend Fuzzy operator> (ApproximateFloat64, Rational const&);
    friend Fuzzy operator<=(ApproximateFloat64, Rational const&);
    friend Fuzzy operator>=(ApproximateFloat64, Rational const&);
    friend Bool same(ApproximateFloat64 x1, ApproximateFloat64 x2);
    friend OutputStream& operator<<(OutputStream& os, ApproximateFloat64 const&);
};



/************ Casting ***************************************************/

MetricFloat64 make_valid(const Flt64&);
BoundedFloat64 make_bound(const Flt64&);
UpperFloat64 make_upper(const Flt64&);
LowerFloat64 make_lower(const Flt64&);

ExactFloat64 make_exact(const ApproximateFloat64&);
ErrorFloat64 make_error(const ApproximateFloat64&);
ApproximateFloat64 make_apprx(const ApproximateFloat64&);

inline Flt64& make_raw(ApproximateFloat64& a) { return reinterpret_cast<Flt64&>(a); }
inline Flt64& make_raw(ErrorFloat64& e) { return reinterpret_cast<Flt64&>(e); }
inline Flt64& make_raw(ExactFloat64& x) { return reinterpret_cast<Flt64&>(x); }
inline Flt64 const& make_raw(ApproximateFloat64 const& a) { return reinterpret_cast<Flt64 const&>(a); }
inline Flt64 const& make_raw(ErrorFloat64 const& e) { return reinterpret_cast<Flt64 const&>(e); }
inline Flt64 const& make_raw(ExactFloat64 const& x) { return reinterpret_cast<Flt64 const&>(x); }

/************ Creation ***************************************************/

template<class P> inline auto Float64Creator::create(unsigned int n, P p) -> Type<P> {
    return Float64Template<Exact>(n); }
template<class P> inline auto Float64Creator::create(int n, P p) -> Type<P> {
    return Float64Template<Exact>(n); }
template<class P> inline auto Float64Creator::create(double d, P p) -> Type<Approximate> {
    return Float64Template<Approximate>(d); }


/************ Validation ***************************************************/

Bool refines(ErrorFloat64 x1, ErrorFloat64 x2);

ErrorFloat64 error(BoundedFloat64 x);
BoundedFloat64 refinement(BoundedFloat64 x1, BoundedFloat64 x2);
Bool refines(BoundedFloat64 x1, BoundedFloat64 x2);
Bool inconsistent(BoundedFloat64 x1, BoundedFloat64 x2);
Bool represents(BoundedFloat64 x1, ExactFloat64 x2);

ErrorFloat64 error(MetricFloat64 x);
MetricFloat64 refinement(MetricFloat64 x1, MetricFloat64 x2);
Bool refines(MetricFloat64 x1, MetricFloat64 x2);
Bool inconsistent(MetricFloat64 x1, MetricFloat64 x2);
Bool represents(MetricFloat64 x1, ExactFloat64 x2);

/************ Exact arithmetic ***************************************************/

inline auto operator+(ExactFloat64 x1, ExactFloat64 x2) -> decltype(add(x1,x2,DefaultExactArithmeticTag())) {
    return add(x1,x2,DefaultExactArithmeticTag()); }
inline auto operator-(ExactFloat64 x1, ExactFloat64 x2) -> decltype(add(x1,x2,DefaultExactArithmeticTag())) {
    return sub(x1,x2,DefaultExactArithmeticTag()); }
inline auto operator*(ExactFloat64 x1, ExactFloat64 x2) -> decltype(add(x1,x2,DefaultExactArithmeticTag())) {
    return mul(x1,x2,DefaultExactArithmeticTag()); }
inline auto operator/(ExactFloat64 x1, ExactFloat64 x2) -> decltype(add(x1,x2,DefaultExactArithmeticTag())) {
    return div(x1,x2,DefaultExactArithmeticTag()); }

/************ Prototype ***************************************************/

template<class T> class Prototype;

template<> class Prototype<ApproximateFloat64> {
    typedef ApproximateFloat64 T;
  public:
    Prototype() { }
    Prototype(T t) { }
    template<class X> T convert(X const& x) { return T(x); }
};

template<> class Prototype<LowerFloat64> {
    typedef LowerFloat64 T;
  public:
    Prototype() { }
    Prototype(T t) { }
    template<class X> T convert(X const& x) { return T(x); }
};

template<> class Prototype<UpperFloat64> {
    typedef UpperFloat64 T;
  public:
    Prototype() { }
    Prototype(T t) { }
    template<class X> T convert(X const& x) { return T(x); }
};

template<> class Prototype<BoundedFloat64> {
    typedef BoundedFloat64 T;
  public:
    Prototype() { }
    Prototype(T t) { }
    template<class X> T convert(X const& x) { return T(x); }
};

template<> class Prototype<MetricFloat64> {
    typedef MetricFloat64 T;
  public:
    Prototype() { }
    Prototype(T t) { }
    template<class X> T convert(X const& x) { return T(x); }
};

template<> class Prototype<ExactFloat64> {
    typedef ExactFloat64 T;
  public:
    Prototype() { }
    Prototype(T t) { }
    template<class X> T convert(X const& x) { return T(x); }
};


} // namespace Ariadne

#endif
