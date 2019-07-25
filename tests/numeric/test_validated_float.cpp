/***************************************************************************
 *            test_validated_float.cpp
 *
 *  Copyright  2006-14  Alberto Casagrande, Pieter Collins
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

#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include "config.hpp"
#include "numeric/builtin.hpp"
#include "numeric/decimal.hpp"
#include "numeric/rational.hpp"
#include "numeric/float.decl.hpp"
#include "numeric/float.hpp"

#include "../test.hpp"

using namespace Ariadne;
using namespace std;

namespace Ariadne {
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator==(FloatDP const& x, Q const& q) { return Rational(x)==q; }
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator!=(FloatDP const& x, Q const& q) { return Rational(x)!=q; }
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator<=(FloatDP const& x, Q const& q) { return Rational(x)<=q; }
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator>=(FloatDP const& x, Q const& q) { return Rational(x)>=q; }
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator< (FloatDP const& x, Q const& q) { return Rational(x)< q; }
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator> (FloatDP const& x, Q const& q) { return Rational(x)> q; }
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator==(FloatMP const& x, Q const& q) { return Rational(x)==q; }
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator!=(FloatMP const& x, Q const& q) { return Rational(x)!=q; }
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator<=(FloatMP const& x, Q const& q) { return Rational(x)<=q; }
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator>=(FloatMP const& x, Q const& q) { return Rational(x)>=q; }
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator< (FloatMP const& x, Q const& q) { return Rational(x)< q; }
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator> (FloatMP const& x, Q const& q) { return Rational(x)> q; }

template<class F> Bool models(LowerBound<F> x, Rational q) { return x.raw() <= q; }
template<class F> Bool models(UpperBound<F> x, Rational q) { return x.raw() >= q; }
template<class F> Bool models(Bounds<F> x, Rational q) { return x.lower_raw() <= q and x.upper_raw() >= q; }
template<class F> Bool models(Ball<F> x, Rational q) { return x.error_raw() >= abs(Rational(x.value_raw())-q); }

template<> String class_name<DoublePrecision>() { return "DoublePrecision"; }
template<> String class_name<MultiplePrecision>() { return "MultiplePrecision"; }
} // namespace Ariadne


template<class PR>
class TestFloats
{
    typedef double Double;
  protected:
    Nat m; Int n; Double ad; ExactDouble ed;
    Integer z; Dyadic w; Decimal d; Rational q; Real r;
  public:
    TestFloats();
  public:
    static Rational to_rational(FloatType<ApproximateTag,PR> x) { return Rational(x.raw()); }
    static Rational to_rational(FloatType<LowerTag,PR> x) { return Rational(x.raw()); }
    static Rational to_rational(FloatType<UpperTag,PR> x) { return Rational(x.raw()); }
    static Rational to_rational(FloatType<ExactTag,PR> x) { return Rational(x.raw()); }
};

template<class PR>
TestFloats<PR>::TestFloats()
    : m(1), n(1), ad(1.0), ed(1.0), z(1), w(1), d(1.0), q(1), r(1)
{
}

template<class PR>
class TestFloatApproximation
    : public TestFloats<PR>
{
    typedef FloatType<ApproximateTag,PR> FloatApproximationType;

    using TestFloats<PR>::m,TestFloats<PR>::n,TestFloats<PR>::ad,TestFloats<PR>::ed;
    using TestFloats<PR>::z,TestFloats<PR>::w,TestFloats<PR>::d,TestFloats<PR>::q;
  private:
    PR precision;
  public:
    TestFloatApproximation(PR prec) : precision(prec) { };
    Void test();
  private:
    Void test_concept();
    Void test_arithmetic();
    Void test_comparison();
};


template<class PR> Void
TestFloatApproximation<PR>::test()
{
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(test_comparison());
}

template<class PR> Void
TestFloatApproximation<PR>::test_concept()
{
    typedef DoublePrecision PRE; PRE pre;
    PR pr=precision;

    FloatApproximation<PR>::set_output_places(17);
    FloatApproximationType ax(pr);
    FloatApproximationType rx(pr);

    // Constructors
    FloatApproximationType(0u,pr); FloatApproximationType(0,pr); FloatApproximationType(0.0,pr);
    FloatApproximationType(Integer(0),pr); FloatApproximationType(Dyadic(0),pr); FloatApproximationType(Rational(0),pr);

    rx=FloatApproximationType(RawFloat<PR>(pr));
    rx=FloatApproximationType(FloatValue<PR>(pr)); FloatApproximationType(FloatBall<PR,PRE>(pr,pre)); FloatApproximationType(FloatBounds<PR>(pr));
    rx=FloatApproximationType(FloatUpperBound<PR>(pr)); rx=FloatApproximationType(FloatLowerBound<PR>(pr));
    rx=FloatApproximationType(FloatApproximation<PR>(pr));

    // Assignment
    rx=m; rx=n; rx=ad; rx=ed; rx=z; rx=w; rx=q;

    // ExactTag operations
    rx=nul(ax); rx=pos(ax); rx=neg(ax); rx=hlf(ax); rx=sqr(ax); rx=rec(ax);

    rx=operator+(ax); rx=operator-(ax); rx=operator+(ax,ax); rx=operator-(ax,ax); rx=operator*(ax,ax); rx=operator/(ax,ax);

    // Arithmetic
    rx=add(ax,ax);
    rx=sub(ax,ax);
    rx=mul(ax,ax);
    rx=div(ax,ax);
    rx=pow(ax,m);
    rx=pow(ax,n);

    // Order
    rx=max(ax,ax); rx=min(ax,ax); rx=abs(ax);

    // Transcendental functions
    rx=sqrt(ax);
    rx=exp(ax);
    rx=log(ax);
    rx=sin(ax);
    rx=cos(ax);
    rx=tan(ax);
    rx=atan(ax);
}

template<class PR> Void
TestFloatApproximation<PR>::test_arithmetic()
{
    RawFloat<PR> third(Rational(1,3),near,precision);
    RawFloat<PR> fifth(Rational(1,5),near,precision);
    RawFloat<PR> seventh(Rational(1,7),near,precision);

    ARIADNE_TEST_SAME(max(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(third));
    ARIADNE_TEST_SAME(min(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(fifth));
    ARIADNE_TEST_SAME(abs(FloatApproximation<PR>(neg(third))),FloatApproximation<PR>(third));

    ARIADNE_TEST_SAME(mag(FloatApproximation<PR>(neg(third))),PositiveFloatApproximation<PR>(third));
    ARIADNE_TEST_SAME(mig(FloatApproximation<PR>(neg(third))),PositiveFloatApproximation<PR>(third));

    ARIADNE_TEST_SAME(operator+(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(add(near,third,fifth)));
    ARIADNE_TEST_SAME(operator-(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(sub(near,third,fifth)));
    ARIADNE_TEST_SAME(operator*(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(mul(near,third,fifth)));
    ARIADNE_TEST_SAME(operator/(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(div(near,third,fifth)));

    ARIADNE_TEST_SAME(nul(FloatApproximation<PR>(third)),FloatApproximation<PR>(nul(third)));
    ARIADNE_TEST_SAME(pos(FloatApproximation<PR>(third)),FloatApproximation<PR>(pos(third)));
    ARIADNE_TEST_SAME(neg(FloatApproximation<PR>(third)),FloatApproximation<PR>(neg(third)));
    ARIADNE_TEST_SAME(hlf(FloatApproximation<PR>(third)),FloatApproximation<PR>(hlf(third)));
    ARIADNE_TEST_SAME(add(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(add(near,third,fifth)));
    ARIADNE_TEST_SAME(sub(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(sub(near,third,fifth)));
    ARIADNE_TEST_SAME(mul(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(mul(near,third,fifth)));
    ARIADNE_TEST_SAME(div(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(div(near,third,fifth)));
    ARIADNE_TEST_SAME(fma(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth),FloatApproximation<PR>(seventh)),
                          FloatApproximation<PR>(fma(near,third,fifth,seventh)));
    ARIADNE_TEST_SAME(pow(FloatApproximation<PR>(fifth),3),FloatApproximation<PR>(pow(near,fifth,3)));
}

inline bool same(LogicalValue l1, LogicalValue l2) { return static_cast<uchar>(l1) == static_cast<uchar>(l2); }

inline bool same(ApproximateKleenean ak1, ApproximateKleenean ak2) {
    return same(reinterpret_cast<LogicalValue const&>(ak1),reinterpret_cast<LogicalValue const&>(ak2)); }

template<class PR> Void
TestFloatApproximation<PR>::test_comparison()
{
    PR pr=precision;
    ARIADNE_TEST_SAME(FloatApproximation<PR>(-2,pr)==FloatApproximation<PR>(-1,pr),ApproximateKleenean(false));
    ARIADNE_TEST_SAME(FloatApproximation<PR>(-2,pr)!=FloatApproximation<PR>(-1,pr),ApproximateKleenean(true));
    ARIADNE_TEST_SAME(FloatApproximation<PR>(-2,pr)< FloatApproximation<PR>(-1,pr),ApproximateKleenean(true));
    ARIADNE_TEST_SAME(FloatApproximation<PR>(-2,pr)> FloatApproximation<PR>(-1,pr),ApproximateKleenean(false));
    ARIADNE_TEST_SAME(FloatApproximation<PR>(-2,pr)<=FloatApproximation<PR>(-1,pr),ApproximateKleenean(true));
    ARIADNE_TEST_SAME(FloatApproximation<PR>(-2,pr)>=FloatApproximation<PR>(-1,pr),ApproximateKleenean(false));

    ARIADNE_TEST_SAME(FloatApproximation<PR>(-1,pr)==FloatApproximation<PR>(-1,pr),ApproximateKleenean(true));
    ARIADNE_TEST_SAME(FloatApproximation<PR>(-1,pr)!=FloatApproximation<PR>(-1,pr),ApproximateKleenean(false));
    ARIADNE_TEST_SAME(FloatApproximation<PR>(-1,pr)< FloatApproximation<PR>(-1,pr),ApproximateKleenean(false));
    ARIADNE_TEST_SAME(FloatApproximation<PR>(-1,pr)> FloatApproximation<PR>(-1,pr),ApproximateKleenean(false));
    ARIADNE_TEST_SAME(FloatApproximation<PR>(-1,pr)<=FloatApproximation<PR>(-1,pr),ApproximateKleenean(true));
    ARIADNE_TEST_SAME(FloatApproximation<PR>(-1,pr)>=FloatApproximation<PR>(-1,pr),ApproximateKleenean(true));
}



template<class PR>
class TestDirectedFloats
    : public TestFloats<PR>
{
    typedef FloatType<ApproximateTag,PR> FloatApproximationType;
    typedef FloatType<LowerTag,PR> FloatLowerBoundType;
    typedef FloatType<UpperTag,PR> FloatUpperBoundType;
    typedef FloatType<BoundedTag,PR> FloatBoundsType;
    typedef FloatType<MetricTag,PR> FloatBallType;
    typedef FloatType<ExactTag,PR> FloatValueType;

    typedef PositiveFloatLowerBound<PR> PositiveFloatLowerBoundType;
    typedef PositiveFloatUpperBound<PR> PositiveFloatUpperBoundType;

  private:
    PR precision;
  public:
    TestDirectedFloats(PR prec) : precision(prec) { };
    Void test();
  private:
    Void test_concept();
    Void test_precision();
    Void test_conversions();
    Void test_validation();
    Void test_rounded_arithmetic();
    Void test_comparison();
};

template<class PR> Void
TestDirectedFloats<PR>::test()
{
    ARIADNE_TEST_CALL(test_precision());
    ARIADNE_TEST_CALL(test_conversions());
    ARIADNE_TEST_CALL(test_validation());
    ARIADNE_TEST_CALL(test_rounded_arithmetic());
    ARIADNE_TEST_CALL(test_comparison());
}

template<class PR> Void
TestDirectedFloats<PR>::test_concept()
{
    Nat m;
    FloatApproximationType ax(1);
    FloatLowerBoundType lx(1);
    FloatUpperBoundType ux(1);
    FloatValueType ex(1);

    PositiveFloatLowerBoundType plx(1u);
    PositiveFloatUpperBoundType pux(1u);

    lx=+lx; lx=-ux; lx=lx+lx; lx=lx-ux; lx=lx*m; lx=lx/m;
    ex=nul(lx); lx=pos(lx); lx=neg(ux); lx=hlf(lx);
    lx=add(lx,lx); lx=sub(lx,ux);
    lx=sqrt(lx); lx=exp(lx); lx=log(lx); lx=atan(lx);
    lx=max(lx,lx); lx=min(lx,lx);

    plx=+plx; plx=plx+plx; plx=plx*plx; plx=plx/pux; plx=plx/m;
    plx=nul(plx); plx=pos(plx); plx=hlf(plx); plx=sqr(plx); plx=rec(pux);
    plx=add(plx,plx); plx=mul(plx,plx); plx=div(plx,pux);
    plx=sqrt(plx); plx=exp(lx); lx=log(plx); plx=atan(plx);
    plx=max(plx,plx); plx=max(plx,lx); plx=max(lx,plx); plx=min(plx,plx);


    ux=+ux; ux=-lx; ux=ux+ux; ux=ux-lx; ux=ux*m; ux=ux/m;
    ex=nul(ux); ux=pos(ux); ux=neg(lx); ux=hlf(ux);
    ux=add(ux,ux); ux=sub(ux,lx);
    ux=sqrt(ux); ux=exp(ux); ux=log(ux);
    ux=max(ux,ux); ux=min(ux,ux);

    pux=+pux; pux=pux+pux; pux=pux*pux; pux=pux/plx; pux=pux/m;
    pux=nul(pux); pux=pos(pux); pux=hlf(pux); pux=sqr(pux); pux=rec(plx);
    pux=add(pux,pux); pux=mul(pux,pux); pux=div(pux,plx);
    pux=sqrt(pux); pux=exp(lx); lx=log(pux); pux=atan(pux);
    pux=max(pux,pux); pux=max(pux,lx); pux=max(lx,pux); pux=min(pux,pux);
}

template<class PR> Void
TestDirectedFloats<PR>::test_precision()
{
    FloatLowerBoundType lx(Rational(1),precision);
    ARIADNE_TEST_EQUALS(max(lx,lx).precision(),precision);
    ARIADNE_TEST_EQUALS(min(lx,lx).precision(),precision);
    ARIADNE_TEST_EQUALS((lx+lx).precision(),precision);
    ARIADNE_TEST_EQUALS(exp(lx).precision(),precision);
    ARIADNE_TEST_EQUALS((lx+2u).precision(),precision);
    ARIADNE_TEST_EQUALS((lx*2u).precision(),precision);
    ARIADNE_TEST_EQUALS((lx/2u).precision(),precision);
    ARIADNE_TEST_EQUALS((lx+2).precision(),precision);
    ARIADNE_TEST_EQUALS((lx-2).precision(),precision);
}

template<class PR> Void
TestDirectedFloats<PR>::test_conversions()
{
    Rational one=1;
    Rational four_thirds=4*one/3;
    Rational five_thirds=5*one/3;
    Rational neg_five_thirds=-5*one/3;

    ValidatedLowerNumber l(four_thirds);
    ValidatedUpperNumber u(five_thirds);

    ARIADNE_TEST_COMPARE(FloatBoundsType(l,u,precision).lower().raw(),<=,four_thirds);
    ARIADNE_TEST_COMPARE(FloatBoundsType(l,u,precision).upper().raw(),>=,five_thirds);

    ARIADNE_TEST_COMPARE(FloatLowerBoundType(five_thirds,precision).raw(),<=,five_thirds);
    ARIADNE_TEST_COMPARE(FloatLowerBoundType(neg_five_thirds,precision).raw(),<=,neg_five_thirds);
    ARIADNE_TEST_COMPARE(FloatUpperBoundType(five_thirds,precision).raw(),>=,five_thirds);
    ARIADNE_TEST_COMPARE(FloatUpperBoundType(neg_five_thirds,precision).raw(),>=,neg_five_thirds);
}

template<class PR> Void
TestDirectedFloats<PR>::test_validation() {
    Rational one=1;
    Rational two_=2;
    ARIADNE_TEST_ASSERT(refines(FloatLowerBoundType(one,precision),FloatLowerBoundType(-one,precision)));
    ARIADNE_TEST_ASSERT(refines(FloatUpperBoundType(-one,precision),FloatUpperBoundType(+one,precision)));
    ARIADNE_TEST_ASSERT(refines(FloatUpperBoundType(-two_,precision),FloatUpperBoundType(-one,precision)));
//    ARIADNE_TEST_ASSERT(refines(rec(FloatUpperBoundType(-two,precision)),rec(FloatUpperBoundType(-one,precision))));
}

template<class PR> Void
TestDirectedFloats<PR>::test_rounded_arithmetic() {
    Rational one=1;
    Rational third=one/3;
    Rational fifth=one/3;
    ARIADNE_TEST_COMPARE((FloatLowerBoundType(third,precision)+FloatLowerBoundType(fifth,precision)).raw(),<=,third+fifth);
    ARIADNE_TEST_COMPARE((FloatLowerBoundType(third,precision)-FloatUpperBoundType(fifth,precision)).raw(),<=,third-fifth);
    ARIADNE_TEST_ASSERT(refines(FloatUpperBoundType(third+fifth,precision),FloatUpperBoundType(third,precision)+FloatUpperBoundType(fifth,precision)));
    ARIADNE_TEST_ASSERT(refines(FloatLowerBoundType(third+fifth,precision),FloatLowerBoundType(third,precision)+FloatLowerBoundType(fifth,precision)));
    ARIADNE_TEST_COMPARE((PositiveFloatLowerBoundType(third,precision)*PositiveFloatLowerBoundType(fifth,precision)).raw(),<=,third*fifth);
    ARIADNE_TEST_COMPARE((PositiveFloatUpperBoundType(third,precision)*PositiveFloatUpperBoundType(fifth,precision)).raw(),>=,third*fifth);
    ARIADNE_TEST_COMPARE((PositiveFloatLowerBoundType(third,precision)/PositiveFloatUpperBoundType(fifth,precision)).raw(),<=,third/fifth);
    ARIADNE_TEST_COMPARE((PositiveFloatUpperBoundType(third,precision)/PositiveFloatLowerBoundType(fifth,precision)).raw(),>=,third/fifth);
}

template<class PR> Void
TestDirectedFloats<PR>::test_comparison() {
    {
        PositiveFloatUpperBoundType one(1u);
        PositiveFloatUpperBoundType third=one/3u;

        ARIADNE_TEST_ASSERT(definitely(third+third < 1u));
        ARIADNE_TEST_ASSERT(definitely(third+third <= 1u));
        ARIADNE_TEST_ASSERT(definitely(5u*third < 2u));
        ARIADNE_TEST_ASSERT(definitely(5u*third <= 2u));
        ARIADNE_TEST_ASSERT(possibly(3u*third > 1u));
        ARIADNE_TEST_ASSERT(possibly(3u*third >= 1u));
    }

    {
        PositiveFloatLowerBoundType one(1u);
        PositiveFloatLowerBoundType two_thirds=one*2u/3u;
        ARIADNE_TEST_ASSERT(possibly(two_thirds+two_thirds < 2u));
        ARIADNE_TEST_ASSERT(possibly(two_thirds+two_thirds <= 2u));
        ARIADNE_TEST_ASSERT(definitely(2u*two_thirds > 1u));
        ARIADNE_TEST_ASSERT(definitely(2u*two_thirds >= 1u));
        ARIADNE_TEST_ASSERT(possibly(3u*two_thirds < 2u));
        ARIADNE_TEST_ASSERT(possibly(3u*two_thirds <= 2u));
    }


}


template<class PR, class PRE=PR>
class TestFloatBall
    : public TestFloats<PR>
{
    typedef RawFloat<PR> RawFloatType;
    typedef FloatType<ApproximateTag,PR> FloatApproximationType;
    typedef FloatType<LowerTag,PR> FloatLowerBoundType;
    typedef FloatType<UpperTag,PR> FloatUpperBoundType;
    typedef FloatType<BoundedTag,PR> FloatBoundsType;
    typedef FloatBall<PR,PRE> FloatBallType;
    typedef FloatType<ExactTag,PR> FloatValueType;
  private:
    PR precision; PRE error_precision;
    static PRE _make_error_precision(PR const& pr) { if constexpr (IsSame<PR,PRE>::value) { return pr; } else { return PRE(); } }
  public:
    TestFloatBall(PR prec) : precision(prec), error_precision(_make_error_precision(prec)) { };
    TestFloatBall(PR prec, PRE err_prec) : precision(prec), error_precision(err_prec) { };
    Void test();
  private:
    using TestFloats<PR>::to_rational;
    Void test_concept();
    Void test_precision();
    Void test_conversions();
    Void test_validation();
    Void test_rounded_arithmetic();
};

template<class PR, class PRE> Void
TestFloatBall<PR,PRE>::test()
{
    ARIADNE_TEST_CALL(test_precision());
    ARIADNE_TEST_CALL(test_conversions());
    ARIADNE_TEST_CALL(test_rounded_arithmetic());
}

template<class PR, class PRE> Void
TestFloatBall<PR,PRE>::test_conversions()
{
    Rational one=1;
    Rational third=one/3;
}

template<class PR, class PRE> Void
TestFloatBall<PR,PRE>::test_precision()
{
    FloatBallType mx(Rational(1),precision);
    ARIADNE_TEST_EQUALS(mx.value().precision(),precision);
    ARIADNE_TEST_EQUALS(mx.error().precision(),error_precision);
    ARIADNE_TEST_EQUALS(max(mx,mx).precision(),precision);
    ARIADNE_TEST_EQUALS(min(mx,mx).precision(),precision);
    ARIADNE_TEST_EQUALS(abs(mx).precision(),precision);
    ARIADNE_TEST_EQUALS((mx+mx).precision(),precision);
    ARIADNE_TEST_EQUALS((mx-mx).precision(),precision);
    ARIADNE_TEST_EQUALS((mx*mx).precision(),precision);
    ARIADNE_TEST_EQUALS((mx/mx).precision(),precision);
    ARIADNE_TEST_EQUALS(sqrt(mx).precision(),precision);
    ARIADNE_TEST_EQUALS(exp(mx).precision(),precision);
    ARIADNE_TEST_EQUALS(log(mx).precision(),precision);
    ARIADNE_TEST_EQUALS(sin(mx).precision(),precision);
    ARIADNE_TEST_EQUALS(cos(mx).precision(),precision);
    ARIADNE_TEST_EQUALS(tan(mx).precision(),precision);
    ARIADNE_TEST_EQUALS(atan(mx).precision(),precision);
    ARIADNE_TEST_EQUALS((mx+2u).precision(),precision);
    ARIADNE_TEST_EQUALS((mx-2u).precision(),precision);
    ARIADNE_TEST_EQUALS((mx*2u).precision(),precision);
    ARIADNE_TEST_EQUALS((mx/2u).precision(),precision);
    ARIADNE_TEST_EQUALS((mx+2).precision(),precision);
    ARIADNE_TEST_EQUALS((mx-2).precision(),precision);
    ARIADNE_TEST_EQUALS((mx*2).precision(),precision);
    ARIADNE_TEST_EQUALS((mx/2).precision(),precision);
}

// TODO: Move to FloatBall if necessary
template<class F, class FE> Bool models(Ball<F,FE> x, Rational y) {
    return abs(Dyadic(x.value())-y) <= Dyadic(x.error_raw()); }
template<class F, class FE> Bool models(Ball<F,FE> x, Value<F> v) {
    return abs(Dyadic(x.value())-Dyadic(v)) <= Dyadic(x.error_raw()); }

template<class PR, class PRE> Void
TestFloatBall<PR,PRE>::test_rounded_arithmetic()
{
    Rational one=1;
    Rational three=3;
    Rational five=5;
    Rational six=6;
    Rational third=one/three;
    Rational fifth=one/five;

    ARIADNE_TEST_BINARY_PREDICATE(models,FloatBallType(third,precision)+FloatBallType(fifth,precision),third+fifth);
    ARIADNE_TEST_BINARY_PREDICATE(models,FloatBallType(third,precision)-FloatBallType(fifth,precision),third-fifth);
    ARIADNE_TEST_BINARY_PREDICATE(models,FloatBallType(third,precision)*FloatBallType(fifth,precision),third*fifth);
    ARIADNE_TEST_BINARY_PREDICATE(models,FloatBallType(third,precision)/FloatBallType(fifth,precision),third/fifth);
//    ARIADNE_TEST_BINARY_PREDICATE(models,1/FloatBallType(three,precision)/FloatBallType(fifth,precision),1/three);

    ARIADNE_TEST_BINARY_PREDICATE(models,pow(FloatBallType(three/five,precision),4u),pow(three/five,4u));
    ARIADNE_TEST_BINARY_PREDICATE(models,pow(FloatBallType(three/five,precision),4),pow(three/five,4));
    ARIADNE_TEST_BINARY_PREDICATE(models,pow(FloatBallType(three/five,precision),-4),pow(three/five,-4));
    ARIADNE_TEST_BINARY_PREDICATE(models,pow(FloatBallType(three/five,precision),-7),pow(three/five,-7));
}



template<class PR>
class TestFloatBounds
{
    typedef RawFloat<PR> RawFloatType;
    typedef FloatType<BoundedTag,PR> FloatBoundsType;
    typedef FloatType<ExactTag,PR> FloatValueType;
  private:
    PR precision;
  public:
    TestFloatBounds(PR prec) : precision(prec) { }
    Void test();
  private:
    Void test_concept();
    Void test_constructors();
    Void test_input();
    Void test_class();
    Void test_comparison();
    Void test_precision();
    Void test_correct_rounded_arithmetic();
    Void test_accurate_rounded_arithmetic();
    Void test_exact_rounded_arithmetic();
    Void test_aliasing();
    Void test_monotone_functions();
    Void test_trigonometric_functions();
    Void regression_tests();
};


template<class PR> Void
TestFloatBounds<PR>::test()
{
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_input());
    ARIADNE_TEST_CALL(test_class());
    ARIADNE_TEST_CALL(test_comparison());
    ARIADNE_TEST_CALL(test_precision());
    ARIADNE_TEST_CALL(test_correct_rounded_arithmetic());
    ARIADNE_TEST_CALL(test_accurate_rounded_arithmetic());
    ARIADNE_TEST_CALL(test_exact_rounded_arithmetic());
    ARIADNE_TEST_CALL(test_monotone_functions());
    ARIADNE_TEST_CALL(test_trigonometric_functions());
    ARIADNE_TEST_CALL(regression_tests());
}

template<class PR> Void
TestFloatBounds<PR>::test_concept()
{
    FloatBoundsType::set_output_places(17);

    Nat m=1;
    Int n=1;
    double d=1;
    FloatValueType x(1);
    RawFloatType a,b;
    FloatBoundsType vx(1);
    FloatBoundsType rx(1);

    // Constructors
    rx=FloatBoundsType(); rx=FloatBoundsType(n); rx=FloatBoundsType(m); rx=FloatBoundsType(d); rx=FloatBoundsType(x); rx=FloatBoundsType(vx);
    rx=FloatBoundsType(n,n); rx=FloatBoundsType(m,m); rx=FloatBoundsType(d,d); rx=FloatBoundsType(a,b);
    rx=FloatBoundsType(n,m); rx=FloatBoundsType(m,d); rx=FloatBoundsType(d,n);
//
    // Assignment
    rx=n; rx=m; rx=x; rx=vx;

    // ExactTag operations
    rx=nul(vx); rx=pos(vx); rx=neg(vx); rx=hlf(vx); rx=sqr(vx); rx=rec(vx);

    rx=operator+(x,x); rx=operator+(x,vx); rx=operator+(vx,x); rx=operator+(vx,vx);
    rx=operator-(x,x); rx=operator-(x,vx); rx=operator-(vx,x); rx=operator-(vx,vx);
    rx=operator*(x,x); rx=operator*(x,vx); rx=operator*(vx,x); rx=operator*(vx,vx);
    rx=operator/(x,x); rx=operator/(x,vx); rx=operator/(vx,x); rx=operator/(vx,vx);

    // Arithmetic
    rx=add(x,x); rx=add(vx,vx);
    rx=sub(x,x); rx=sub(vx,vx);
    rx=mul(x,x); rx=mul(vx,vx);
    rx=div(x,x); rx=div(vx,vx);
    rx=pow(x,m); rx=pow(x,m);
    rx=pow(x,n); rx=pow(x,n);

    // Order
    rx=max(vx,vx); rx=min(vx,vx); rx=abs(vx);

    // Transcendental functions
    rx=sqrt(vx);
    rx=exp(vx);
    rx=log(vx);
    rx=sin(vx);
    rx=cos(vx);
    rx=tan(vx);
    rx=atan(vx);
}

template<class PR> Void
TestFloatBounds<PR>::test_precision()
{

    FloatBoundsType bx(Rational(1),precision);
    ARIADNE_TEST_EQUALS(bx.lower().precision(),precision);
    ARIADNE_TEST_EQUALS(bx.upper().precision(),precision);
    ARIADNE_TEST_EQUALS(max(bx,bx).precision(),precision);
    ARIADNE_TEST_EQUALS(min(bx,bx).precision(),precision);
    ARIADNE_TEST_EQUALS(abs(bx).precision(),precision);
    ARIADNE_TEST_EQUALS((bx+bx).precision(),precision);
    ARIADNE_TEST_EQUALS((bx-bx).precision(),precision);
    ARIADNE_TEST_EQUALS((bx*bx).precision(),precision);
    ARIADNE_TEST_EQUALS((bx/bx).precision(),precision);
    ARIADNE_TEST_EQUALS(sqrt(bx).precision(),precision);
    ARIADNE_TEST_EQUALS(exp(bx).precision(),precision);
    ARIADNE_TEST_EQUALS(log(bx).precision(),precision);
    ARIADNE_TEST_EQUALS(sin(bx).precision(),precision);
    ARIADNE_TEST_EQUALS(cos(bx).precision(),precision);
    ARIADNE_TEST_EQUALS(tan(bx).precision(),precision);
    ARIADNE_TEST_EQUALS(atan(bx).precision(),precision);
    ARIADNE_TEST_EQUALS((bx+2u).precision(),precision);
    ARIADNE_TEST_EQUALS((bx-2u).precision(),precision);
    ARIADNE_TEST_EQUALS((bx*2u).precision(),precision);
    ARIADNE_TEST_EQUALS((bx/2u).precision(),precision);
    ARIADNE_TEST_EQUALS((bx+2).precision(),precision);
    ARIADNE_TEST_EQUALS((bx-2).precision(),precision);
    ARIADNE_TEST_EQUALS((bx*2).precision(),precision);
    ARIADNE_TEST_EQUALS((bx/2).precision(),precision);

    ARIADNE_TEST_EQUALS(abs(FloatBoundsType(-1,2,precision)).lower().precision(),precision);

}



// Test that interval arithmetic is rounded correctly,
// without paying attention to accuracy issues.
template<class PR> Void
TestFloatBounds<PR>::test_correct_rounded_arithmetic()
{
    FloatBoundsType onethird=FloatBoundsType(1)/FloatBoundsType(3);
    ARIADNE_TEST_COMPARE( onethird.lower_raw(), < , onethird.upper_raw() );
    FloatBoundsType one_approx=onethird*FloatBoundsType(3);
    ARIADNE_TEST_COMPARE( one_approx.lower_raw(), < , 1.0 );
    ARIADNE_TEST_COMPARE( one_approx.upper_raw(), > , 1.0 );
}


// Test that interval arithmetic gives the most accurate rounded values
template<class PR> Void
TestFloatBounds<PR>::test_accurate_rounded_arithmetic()
{
    const double min=std::numeric_limits<double>::min();
    const double eps=std::numeric_limits<double>::epsilon();

    ARIADNE_TEST_SAME(FloatBoundsType(1.5)+FloatBoundsType(min),FloatBoundsType(1.5,1.5+eps));
    ARIADNE_TEST_SAME(FloatBoundsType(1.5)-FloatBoundsType(min),FloatBoundsType(1.5-eps,1.5));
    ARIADNE_TEST_SAME(FloatBoundsType(1+eps,1+2*eps)*FloatBoundsType(1+eps,1+3*eps),FloatBoundsType(1+2*eps,1+6*eps));
    ARIADNE_TEST_SAME(FloatBoundsType(1)/FloatBoundsType(3),FloatBoundsType(0.33333333333333331,0.33333333333333337));
    ARIADNE_TEST_SAME(FloatBoundsType(2)/FloatBoundsType(5),FloatBoundsType(0.39999999999999997,0.40000000000000002));

    ARIADNE_TEST_SAME(FloatBoundsType(1.5)+FloatValueType(min),FloatBoundsType(1.5,1.5+eps));
    ARIADNE_TEST_SAME(FloatBoundsType(1.5)-FloatValueType(min),FloatBoundsType(1.5-eps,1.5));
    ARIADNE_TEST_SAME(FloatBoundsType(1+eps,1+2*eps)*FloatValueType(1+eps),FloatBoundsType(1+2*eps,1+4*eps));
    ARIADNE_TEST_SAME(FloatBoundsType(1+3*eps,1+5*eps)/FloatValueType(1+eps),FloatBoundsType(1+eps,1+4*eps));

    ARIADNE_TEST_SAME(FloatValueType(min)-FloatBoundsType(1.5),FloatBoundsType(-1.5,eps-1.5));
    ARIADNE_TEST_SAME(FloatValueType(1+5*eps)/FloatBoundsType(1+2*eps,1+3*eps),FloatBoundsType(1+eps,1+3*eps));

    ARIADNE_TEST_SAME(sqr(FloatBoundsType(1-eps,1+eps)),FloatBoundsType(1-4*eps/2,1+3*eps));

    ARIADNE_TEST_SAME(pow(FloatBoundsType(3,5),-1),FloatBoundsType(0.19999999999999998,0.33333333333333337));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(3,5),-2),FloatBoundsType(0.039999999999999986955,0.11111111111111114658));

    ARIADNE_TEST_SAME(rec(FloatBoundsType(1+2*eps,1+5*eps)),FloatBoundsType(1-10*eps/2,1-3*eps/2));

}


// Test that interval arithmetic gives exact values if possible
template<class PR> Void
TestFloatBounds<PR>::test_exact_rounded_arithmetic()
{
    ARIADNE_TEST_SAME(FloatBoundsType(5,7)+FloatBoundsType(2,4),FloatBoundsType(7,11));
    ARIADNE_TEST_SAME(FloatBoundsType(5,7)-FloatBoundsType(2,6),FloatBoundsType(-1,5));

    ARIADNE_TEST_SAME(FloatBoundsType(5,7)*FloatBoundsType(2,4),FloatBoundsType(10,28));
    ARIADNE_TEST_SAME(FloatBoundsType(5,7)*FloatBoundsType(-2,4),FloatBoundsType(-14,28));
    ARIADNE_TEST_SAME(FloatBoundsType(5,7)*FloatBoundsType(-4,-2),FloatBoundsType(-28,-10));
    ARIADNE_TEST_SAME(FloatBoundsType(-7,5)*FloatBoundsType(2,4),FloatBoundsType(-28,20));
    ARIADNE_TEST_SAME(FloatBoundsType(-7,5)*FloatBoundsType(-2,4),FloatBoundsType(-28,20));
    ARIADNE_TEST_SAME(FloatBoundsType(-7,5)*FloatBoundsType(-4,-2),FloatBoundsType(-20,28));
    ARIADNE_TEST_SAME(FloatBoundsType(-7,-5)*FloatBoundsType(2,4),FloatBoundsType(-28,-10));
    ARIADNE_TEST_SAME(FloatBoundsType(-7,-5)*FloatBoundsType(-2,4),FloatBoundsType(-28,14));
    ARIADNE_TEST_SAME(FloatBoundsType(-7,-5)*FloatBoundsType(-4,-2),FloatBoundsType(10,28));

    ARIADNE_TEST_SAME(FloatBoundsType(5,7)/FloatBoundsType(2,4),FloatBoundsType(1.25,3.50));
    ARIADNE_TEST_SAME(FloatBoundsType(5,7)/FloatBoundsType(-4,-2),FloatBoundsType(-3.50,-1.25));
    ARIADNE_TEST_SAME(FloatBoundsType(-7,5)/FloatBoundsType(2,4),FloatBoundsType(-3.50,2.50));
    ARIADNE_TEST_SAME(FloatBoundsType(-7,5)/FloatBoundsType(-4,-2),FloatBoundsType(-2.50,3.5));
    ARIADNE_TEST_SAME(FloatBoundsType(-7,-5)/FloatBoundsType(2,4),FloatBoundsType(-3.50,-1.25));
    ARIADNE_TEST_SAME(FloatBoundsType(-7,-5)/FloatBoundsType(-4,-2),FloatBoundsType(1.25,3.50));

    ARIADNE_TEST_SAME(pow(FloatBoundsType(5,7),0u),FloatBoundsType(1,1));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-5,7),0u),FloatBoundsType(1,1));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,5),0u),FloatBoundsType(1,1));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,-5),0u),FloatBoundsType(1,1));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(5,7),1u),FloatBoundsType(5,7));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-5,7),1u),FloatBoundsType(-5,7));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,5),1u),FloatBoundsType(-7,5));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,-5),1u),FloatBoundsType(-7,-5));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(5,7),2u),FloatBoundsType(25,49));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-5,7),2u),FloatBoundsType(0,49));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,5),2u),FloatBoundsType(0,49));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,-5),2u),FloatBoundsType(25,49));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(5,7),3u),FloatBoundsType(125,343));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-5,7),3u),FloatBoundsType(-125,343));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,5),3u),FloatBoundsType(-343,125));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,-5),3u),FloatBoundsType(-343,-125));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(5,7),4u),FloatBoundsType(625,2401));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-5,7),4u),FloatBoundsType(0,2401));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,5),4u),FloatBoundsType(0,2401));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,-5),4u),FloatBoundsType(625,2401));

    ARIADNE_TEST_SAME(pow(FloatBoundsType(5,7),0),FloatBoundsType(1,1));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,5),0),FloatBoundsType(1,1));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,-5),0),FloatBoundsType(1,1));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(5,7),1),FloatBoundsType(5,7));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-5,7),1),FloatBoundsType(-5,7));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,5),1),FloatBoundsType(-7,5));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,-5),1),FloatBoundsType(-7,-5));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(5,7),2),FloatBoundsType(25,49));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-5,7),2),FloatBoundsType(0,49));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,5),2),FloatBoundsType(0,49));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,-5),2),FloatBoundsType(25,49));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(5,7),3),FloatBoundsType(125,343));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-5,7),3),FloatBoundsType(-125,343));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,5),3),FloatBoundsType(-343,125));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,-5),3),FloatBoundsType(-343,-125));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(5,7),4),FloatBoundsType(625,2401));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-5,7),4),FloatBoundsType(0,2401));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,5),4),FloatBoundsType(0,2401));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,-5),4),FloatBoundsType(625,2401));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(5,7),5),FloatBoundsType(3125,16807));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-5,7),5),FloatBoundsType(-3125,16807));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,5),5),FloatBoundsType(-16807,3125));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,-5),5),FloatBoundsType(-16807,-3125));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(5,7),7),FloatBoundsType(78125,823543));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-5,7),7),FloatBoundsType(-78125,823543));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,5),7),FloatBoundsType(-823543,78125));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-7,-5),7),FloatBoundsType(-823543,-78125));

    ARIADNE_TEST_SAME(pow(FloatBoundsType(2,4),-1),FloatBoundsType(0.25,0.5));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-4,-2),-1),FloatBoundsType(-0.5,-0.25));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(2,4),-2),FloatBoundsType(0.0625,0.25));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-4,-2),-2),FloatBoundsType(0.0625,0.25));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(2,4),-3),FloatBoundsType(0.015625,0.125));
    ARIADNE_TEST_SAME(pow(FloatBoundsType(-4,-2),-3),FloatBoundsType(-0.125,-0.015625));

    ARIADNE_TEST_SAME(rec(FloatBoundsType(2,4)),FloatBoundsType(0.25,0.50));
    ARIADNE_TEST_SAME(rec(FloatBoundsType(-4,-2)),FloatBoundsType(-0.50,-0.25));
}



template<> Void
TestFloatBounds<DoublePrecision>::test_constructors()
{
    DoublePrecision pr;
    FloatDP zero=0;

    // Construct from pair
    FloatBoundsType xd1(FloatDP(1.125),FloatDP(2.25));
    ARIADNE_TEST_ASSERT(xd1.lower_raw()==1.125); ARIADNE_TEST_ASSERT(xd1.upper_raw()==2.25);

    // Default constructor
    FloatBoundsType xd2;
    if(xd2.lower_raw()>xd2.upper_raw()) {
        ARIADNE_TEST_WARN("FloatBoundsType default constructor returns an empty set.");
    } else {
        ARIADNE_TEST_BINARY_PREDICATE(same,xd2,FloatBoundsType(zero,zero));
    }

    // Constructor from Rational with approximations
    FloatBoundsType xd3(Rational(21,10),pr);
    ARIADNE_TEST_COMPARE(Rational(xd3.lower_raw()),<,Rational(21,10));
    ARIADNE_TEST_COMPARE(Rational(xd3.upper_raw()),>,Rational(21,10));

    // Constructor from rational bounds
    FloatBoundsType xd4(2.1_q,3.2_q,pr);
    ARIADNE_TEST_COMPARE(xd4.lower_raw(),<=,2.1);
    ARIADNE_TEST_COMPARE(xd4.upper_raw(),>=,3.2);

    // Approximate constructor from a single value
    FloatBoundsType xd5(Rational(1,3),pr);
    ARIADNE_TEST_COMPARE(Rational(xd5.lower_raw()),<,Rational(1,3));
    ARIADNE_TEST_COMPARE(Rational(xd5.upper_raw()),>,Rational(1,3));

    // ExactTag constructor from a single value
    FloatBoundsType xd6(FloatDP(1.25));
    ARIADNE_TEST_EQUAL(xd6.lower_raw(),FloatDP(1.25));
    ARIADNE_TEST_EQUAL(xd6.upper_raw(),FloatDP(1.25));
}

template<class PR> Void
TestFloatBounds<PR>::test_constructors()
{
    RawFloatType zero=0;

    // Construct from pair
    FloatBoundsType xd1(RawFloatType(1.125),RawFloatType(2.25));
    ARIADNE_TEST_ASSERT(xd1.lower_raw()==1.125); ARIADNE_TEST_ASSERT(xd1.upper_raw()==2.25);

    // Default constructor
    FloatBoundsType xd2;
    if(xd2.lower_raw()>xd2.upper_raw()) {
        ARIADNE_TEST_WARN("FloatBoundsType default constructor returns an empty set.");
    } else {
        ARIADNE_TEST_BINARY_PREDICATE(same,xd2,FloatBoundsType(zero,zero));
    }

    // Constructor with approximations
    FloatBoundsType xd3(Rational(21,10),precision);
    ARIADNE_TEST_COMPARE(Rational(xd3.lower_raw()),<,Rational(21,10));
    ARIADNE_TEST_COMPARE(Rational(xd3.upper_raw()),>,Rational(21,10));

    // Constructor from approximate values
    FloatBoundsType xd4(2.1_q,3.2_q,precision);
    ARIADNE_TEST_COMPARE(xd4.lower_raw(),<=,Rational(21,10));
    ARIADNE_TEST_COMPARE(xd4.upper_raw(),>=,Rational(32,10));

    // ApproximateTag constructor from a single value
    FloatBoundsType xd5(Rational(1,3),precision);
    ARIADNE_TEST_COMPARE(xd5.lower_raw(),<,Rational(1,3));
    ARIADNE_TEST_COMPARE(xd5.upper_raw(),>,Rational(1,3));

    // ExactTag constructor from a single value
    FloatBoundsType xd6(RawFloatType(1.25));
    ARIADNE_TEST_EQUAL(xd6.lower_raw(),RawFloatType(1.25));
    ARIADNE_TEST_EQUAL(xd6.upper_raw(),RawFloatType(1.25));
}

template<class PR> Void TestFloatBounds<PR>::test_class()
{
    // Test lower, upper, midpoint, radius, width
    RawFloatType one=1.0;
    RawFloatType two_=2.0;

    // Tests for exact operations
    ARIADNE_TEST_EQUAL(FloatBoundsType(-0.25,0.50).lower().raw(),-0.25);
    ARIADNE_TEST_EQUAL(FloatBoundsType(-0.25,0.50).upper().raw(),0.5);
    ARIADNE_TEST_EQUAL(FloatBoundsType(-0.25,0.50).value().raw(),0.125);
    ARIADNE_TEST_EQUAL(FloatBoundsType(-0.25,0.50).error().raw(),0.375)

    // Tests for inexact operations
    ARIADNE_TEST_EQUAL(FloatBoundsType(-1./3,2./3).lower().raw(),-0.33333333333333331483);
    ARIADNE_TEST_EQUAL(FloatBoundsType(-1./3,2./3).upper().raw(),0.66666666666666662966);
    ARIADNE_TEST_EQUAL(FloatBoundsType(-1./3,2./3).value().raw(),0.16666666666666665741);
    ARIADNE_TEST_EQUAL(FloatBoundsType(-1./3,2./3).error().raw(),0.5)

    // Tests for inexact operations
    ARIADNE_TEST_EQUAL(FloatBoundsType(div(down,-one,3),div(up,two_,3)).lower().raw(),-0.33333333333333337034);
    ARIADNE_TEST_EQUAL(FloatBoundsType(div(down,-one,3),div(up,two_,3)).upper().raw(),0.66666666666666674068);
    ARIADNE_TEST_EQUAL(FloatBoundsType(div(down,-one,3),div(up,two_,3)).value().raw(),0.16666666666666668517);
    ARIADNE_TEST_EQUAL(FloatBoundsType(div(down,-one,3),div(up,two_,3)).error().raw(),0.50000000000000011102);
}

template<class PR> Void TestFloatBounds<PR>::test_input()
{
    FloatBoundsType x;

    string input("[1.125,2.75] [0.4,0.6]");
    stringstream iss(input);

    iss >> x;
    ARIADNE_TEST_EQUALS(x.lower_raw(),Rational(9,8))
    ARIADNE_TEST_EQUALS(x.upper_raw(),Rational(11,4))

    iss >> x;
    // TODO: Make stream input produce rounded values
    //ARIADNE_TEST_COMPARE(x.lower_raw(),<,Rational(2,5))
    // ARIADNE_TEST_COMPARE(x.upper_raw(),>,Rational(3,5))
    if(not(x.lower_raw()<=Rational(2,5) and x.upper_raw()>=Rational(3,5))) {
        ARIADNE_TEST_WARN("FloatType<BoundedTag,"<<class_name<PR>()<<"> string constructor returns an approximate interval, not an outwardly rounded interval.");
    }
}

template<class PR> Void TestFloatBounds<PR>::test_comparison() {
    // FIXME: If using Boost style interval tests, uncomment the line below
    // and comment out the line after
    //ARIADNE_TEST_ASSERT(indeterminate(ivld1==ivld2));
    FloatBoundsType ivl1(1.125,2.25);
    FloatBoundsType ivl2=ivl1;

    ARIADNE_TEST_ASSERT(!definitely(ivl1==ivl2));
    ARIADNE_TEST_ASSERT(possibly(ivl1==ivl2));
    FloatBoundsType& ivl1ref=ivl1;
    ivl1ref=FloatBoundsType(5.25,7.375);
    cout << "ivl1ref=" << ivl1ref << endl;
    ARIADNE_TEST_ASSERT(ivl1ref.lower_raw()==RawFloatType(5.25));
}

template<class PR> Void TestFloatBounds<PR>::test_aliasing() {
    FloatValueType ex2(1.5);
    FloatValueType ex3(2.25);

    FloatBoundsType vx1;
    FloatBoundsType vx2(1.5,2.25);
    FloatBoundsType vx3(3.125,4.0625);

    // Check to make sure aliases are handled correctly
    vx1=vx3; vx1=vx2-vx1; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,FloatBoundsType(vx2-vx3));
    vx1=vx3; vx1=vx2*vx1; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,FloatBoundsType(vx2*vx3));
    vx1=vx2; vx1=vx1*vx3; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,FloatBoundsType(vx2*vx3));
    vx1=vx2; vx1=vx1*ex3; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,FloatBoundsType(vx2*ex3));
    vx1=vx3; vx1=ex2*vx1; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,FloatBoundsType(ex2*vx3));
    vx1=vx2; vx1=vx1/vx3; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,FloatBoundsType(vx2/vx3));
    vx1=vx2; vx1=vx1/ex3; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,FloatBoundsType(vx2/ex3));
    vx1=vx3; vx1=ex2/vx1; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,FloatBoundsType(ex2/vx3));
}

template<class PR> Void TestFloatBounds<PR>::test_monotone_functions()
{

    FloatBoundsType two_(2.0);
    FloatBoundsType sqrttwo=sqrt(two_);
    ARIADNE_TEST_PRINT(sqrttwo);
    ARIADNE_TEST_COMPARE(sqrttwo.lower_raw(),<=,1.4142135623730949);
    ARIADNE_TEST_COMPARE(sqrttwo.lower_raw(),> ,1.4142135623730947);
    ARIADNE_TEST_COMPARE(sqrttwo.upper_raw(),>=,1.4142135623730951);
    ARIADNE_TEST_COMPARE(sqrttwo.upper_raw(),< ,1.4142135623730954);

    FloatBoundsType one(1.0);
    FloatBoundsType expone=exp(one);
    ARIADNE_TEST_PRINT(expone);
    ARIADNE_TEST_COMPARE(expone.lower_raw(),<,2.71828182845905);
    ARIADNE_TEST_COMPARE(expone.lower_raw(),>,2.71828182845903);
    ARIADNE_TEST_COMPARE(expone.upper_raw(),>,2.71828182845904);
    ARIADNE_TEST_COMPARE(expone.upper_raw(),<,2.71828182845906);
    ARIADNE_TEST_ASSERT(expone.lower_raw()<expone.upper_raw());

    FloatBoundsType e(2.7182818284590451,2.7182818284590455);
    FloatBoundsType loge=log(e);
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_COMPARE(loge.lower_raw(),<,1.0);
    ARIADNE_TEST_COMPARE(loge.lower_raw(),>,0.9999999999998);
    ARIADNE_TEST_COMPARE(loge.upper_raw(),>,1.0);
    ARIADNE_TEST_COMPARE(loge.upper_raw(),<,1.000000000002);
}

template<class PR> Void TestFloatBounds<PR>::test_trigonometric_functions()
{
    try {
        FloatBoundsType x(6.283185307179586,6.283185307179587);
        FloatBoundsType sinx=sin(x);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_COMPARE(sinx.lower_raw(),<,0.0);
        ARIADNE_TEST_COMPARE(sinx.lower_raw(),>,-1e-14);
        ARIADNE_TEST_COMPARE(sinx.upper_raw(),>,0.0);
        ARIADNE_TEST_COMPARE(sinx.upper_raw(),<,+1e-14);
        ARIADNE_TEST_ASSERT(sinx.lower_raw()<sinx.upper_raw());
    }
    catch(...) { }

    try {
        FloatBoundsType x(7.0685834705770345);
        FloatBoundsType sinx=sin(x);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_COMPARE(sinx.lower_raw(),<,0.7071067811866);
        ARIADNE_TEST_COMPARE(sinx.upper_raw(),>,0.7071067811865);
        ARIADNE_TEST_ASSERT(sinx.lower_raw()<sinx.upper_raw());
    }
    catch(...) { }

}

template<class PR> Void TestFloatBounds<PR>::regression_tests() {
    RawFloatType inf_=RawFloatType::inf(precision);

    // Regression test; fails dramatically on certain types of rounding
    {
        FloatBoundsType x(1.5707963267948966,1.5707963267948968);
        FloatBoundsType cosx=cos(x);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_COMPARE(cosx.lower_raw(),<,0.0);
        ARIADNE_TEST_COMPARE(cosx.lower_raw(),>,-1e-14);
        ARIADNE_TEST_COMPARE(cosx.upper_raw(),>,0.0);
        ARIADNE_TEST_COMPARE(cosx.upper_raw(),<,+1e-14);
        ARIADNE_TEST_ASSERT(cosx.lower_raw()<cosx.upper_raw());
    }

    // Regression test for dividing by interval with lower endpoint -0.0 or upper endpoint +0.0

    ARIADNE_TEST_EQUAL((FloatBoundsType(1.0,2.0)/FloatBoundsType(-0.0,1.0)).upper_raw(),+inf_);
    ARIADNE_TEST_EQUAL((FloatBoundsType(1.0,2.0)/FloatBoundsType(-1.0,+0.0)).lower_raw(),-inf_);

    ARIADNE_TEST_EQUAL(rec(FloatBoundsType(-0.0,+1.0)).upper_raw(),+inf_);
    ARIADNE_TEST_EQUAL(rec(FloatBoundsType(-1.0,+0.0)).lower_raw(),-inf_);
}

template<class PR>
class TestFloatValue
{
    using PRE=DoublePrecision;
    typedef RawFloat<PR> RawFloatType;
    typedef FloatBounds<PR> FloatBoundsType;
    typedef FloatBall<PR,PRE> FloatBallType;
    typedef FloatValue<PR> FloatValueType;
  private:
    PR precision;
  public:
    TestFloatValue(PR prec) : precision(prec) { }
    Void test();
  private:
    Void test_concept();
    Void test_operations();
    Void test_predicates();
};

template<class PR> Void
TestFloatValue<PR>::test()
{
    ARIADNE_TEST_CALL(test_operations());
    ARIADNE_TEST_CALL(test_predicates());
}

template<class PR> Void
TestFloatValue<PR>::test_concept()
{
    FloatValue<DoublePrecision>::set_output_places(17);

    PR pr=precision;
    PRE pre;

    Boolean b;
    Nat m=1u;
    Int n=1;
    Integer z=1;
    Dyadic w=1;
    ExactDouble d(1.0);
    TwoExp t(0);
    RawFloatType f(pr);
    FloatValueType vx(pr);
    FloatValueType rx(pr);
    FloatBoundsType rbx(pr);

    // Constructors
    rx=FloatValueType(m,pr); rx=FloatValueType(n,pr); rx=FloatValueType(z,pr); rx=FloatValueType(w,pr);
    rx=FloatValueType(d,pr); rx=FloatValueType(t,pr);
    rx=FloatValueType(pr); rx=FloatValueType(f);

    // Assignment
    rx=m; rx=n; rx=z; rx=w; rx=d; rx=t;

    // Arithmetic operators
    rx=operator+(vx); rx=operator-(vx); rx=t*vx; rx=vx*t; rx=vx/t;
    rbx=operator+(vx,vx); rbx=operator-(vx,vx); rbx=operator*(vx,vx); rbx=operator/(vx,vx);

    // Exact operations
    rx=nul(vx); rx=pos(vx); rx=neg(vx); rx=hlf(vx);
    rx=mul(vx,t); rx=div(vx,t);

    FloatBall<PR,PRE> rmx=add(vx,vx,pre); rmx=sub(vx,vx,pre); rmx=mul(vx,vx,pre); rmx=div(vx,vx,pre);

    // Arithmetic
    rbx=add(vx,vx); rbx=sub(vx,vx); rbx=mul(vx,vx); rbx=div(vx,vx);
    rbx=sqr(vx); rbx=rec(vx); rbx=pow(vx,m); rbx=pow(vx,n);

    // Order
    rx=max(vx,vx); rx=min(vx,vx); rx=abs(vx);

    // Comparisons
    b=(vx==vx); b=(vx!=vx); b=(vx<=vx); b=(vx>=vx); b=(vx< vx); b=(vx> vx);
}

template<class PR> Void
TestFloatValue<PR>::test_operations()
{
    PR pr=precision;
    DoublePrecision pre;
    FloatValueType vr(pr);

    Nat m(5);
    Int n(-3);
    TwoExp t(-12);
    Dyadic w(-3,1u);
    FloatValueType vx(w,pr);
    Dyadic w1(-3,1u), w2(5,2u);
    FloatValueType vx1(w1,pr), vx2(w2,pr);

    ARIADNE_TEST_EQUALS(vx,w);

    ARIADNE_TEST_EQUALS(FloatValueType(RawFloatType(1.25,pr)),1.25_dy);

    ARIADNE_TEST_EQUALS(FloatValueType(3u,pr),3.0_dy);
    ARIADNE_TEST_EQUALS(FloatValueType(-5,pr),-5.0_dy);
    ARIADNE_TEST_EQUALS(FloatValueType(ExactDouble(1.25),pr),1.25_dy);
    ARIADNE_TEST_EQUALS(FloatValueType(TwoExp(-3),pr),0.125_dy);
    ARIADNE_TEST_EQUALS(FloatValueType(Integer(-23),pr),-23.0_dy);
    ARIADNE_TEST_EQUALS(FloatValueType(Dyadic(-23,3u),pr),-2.875_dy);


    ARIADNE_TEST_EQUALS((vr=3u),3.0_dy);
    ARIADNE_TEST_EQUALS((vr=-5),-5.0_dy);
    ARIADNE_TEST_EQUALS((vr=Integer(-23)),-23.0_dy);
    ARIADNE_TEST_EQUALS((vr=TwoExp(-3)),0.125_dy);
    ARIADNE_TEST_EQUALS((vr=Dyadic(-23,3u)),-2.875_dy);

    ARIADNE_TEST_EQUALS(FloatValueType(w,pr).operator Dyadic(),w);
    ARIADNE_TEST_EQUALS(FloatValueType(w,pr).operator Rational(),w);

    ARIADNE_TEST_EQUALS(nul(vx),nul(w));
    ARIADNE_TEST_EQUALS(pos(vx),pos(w));
    ARIADNE_TEST_EQUALS(neg(vx),neg(w));
    ARIADNE_TEST_EQUALS(hlf(vx),hlf(w));

    ARIADNE_TEST_EQUALS(mul(vx,t),mul(w,t));
    ARIADNE_TEST_EQUALS(div(vx,t),w/t);

    ARIADNE_TEST_BINARY_PREDICATE(models,sqr(vx),sqr(w));
    ARIADNE_TEST_BINARY_PREDICATE(models,rec(vx),rec(w));

//    ARIADNE_TEST_BINARY_PREDICATE(models,sqr(vx,pre),sqr(w));
//    ARIADNE_TEST_BINARY_PREDICATE(models,rec(vx,pre),rec(w));

    ARIADNE_TEST_EQUALS(max(vx1,vx2),max(w1,w2));
    ARIADNE_TEST_EQUALS(min(vx1,vx2),min(w1,w2));
    ARIADNE_TEST_EQUALS(abs(vx),abs(w));
    ARIADNE_TEST_EQUALS(mig(vx).raw(),abs(w));
    ARIADNE_TEST_EQUALS(mag(vx).raw(),abs(w));
    ARIADNE_TEST_SAME(mag(vx),PositiveFloatUpperBound<PR>(abs(w),pr));
    ARIADNE_TEST_SAME(mig(vx),PositiveFloatLowerBound<PR>(abs(w),pr));

    ARIADNE_TEST_EQUALS(t*vx,t*w);
    ARIADNE_TEST_EQUALS(vx*t,w*t);
    ARIADNE_TEST_EQUALS(vx/t,w/t);

    ARIADNE_TEST_BINARY_PREDICATE(models,vx1+vx2,w1+w2);
    ARIADNE_TEST_BINARY_PREDICATE(models,vx1-vx2,w1-w2);
    ARIADNE_TEST_BINARY_PREDICATE(models,vx1*vx2,w1*w2);
    ARIADNE_TEST_BINARY_PREDICATE(models,vx1/vx2,w1/w2);

    ARIADNE_TEST_BINARY_PREDICATE(models,add(vx1,vx2),add(w1,w2));
    ARIADNE_TEST_BINARY_PREDICATE(models,sub(vx1,vx2),sub(w1,w2));
    ARIADNE_TEST_BINARY_PREDICATE(models,mul(vx1,vx2),mul(w1,w2));
    ARIADNE_TEST_BINARY_PREDICATE(models,div(vx1,vx2),div(w1,w2));

    ARIADNE_TEST_BINARY_PREDICATE(models,add(vx1,vx2,pre),add(w1,w2));
    ARIADNE_TEST_BINARY_PREDICATE(models,sub(vx1,vx2,pre),sub(w1,w2));
    ARIADNE_TEST_BINARY_PREDICATE(models,mul(vx1,vx2,pre),mul(w1,w2));
    ARIADNE_TEST_BINARY_PREDICATE(models,div(vx1,vx2,pre),div(w1,w2));


    ARIADNE_TEST_SAME(add(vx1,vx2,pre),FloatBallType(add(w1,w2),pr,pre));

    ARIADNE_TEST_BINARY_PREDICATE(models,pow(vx,m),pow(w,m));
    ARIADNE_TEST_BINARY_PREDICATE(models,pow(vx,n),pow(Rational(w),n));

//    friend Bounds<F> med(Value<F> const& x1, Value<F> const& x2);
//    friend Bounds<F> rad(Value<F> const& x1, Value<F> const& x2);

    ARIADNE_TEST_BINARY_PREDICATE(models,sqr(sqrt(abs(vx))),abs(w));
    ARIADNE_TEST_BINARY_PREDICATE(models,log(exp(vx)),w);
    ARIADNE_TEST_BINARY_PREDICATE(models,exp(log(abs(vx))),abs(w));
    ARIADNE_TEST_BINARY_PREDICATE(models,atan(tan(vx)),w);
    ARIADNE_TEST_BINARY_PREDICATE(models,tan(atan(vx)),w);

    ARIADNE_TEST_SAME(sin(vx),sin(FloatBounds<PR>(vx)));
    ARIADNE_TEST_SAME(cos(vx),cos(FloatBounds<PR>(vx)));
    ARIADNE_TEST_SAME(tan(vx),tan(FloatBounds<PR>(vx)));
    ARIADNE_TEST_SAME(atan(vx),atan(FloatBounds<PR>(vx)));

//    ARIADNE_TEST_SAME(shft(vx,n),shft(w,n));
    ARIADNE_TEST_EQUALS(integer_cast(FloatValue<PR>(Dyadic(2),pr)),Integer(2));
    ARIADNE_TEST_FAIL(integer_cast(FloatValue<PR>(Dyadic(7,2u),pr)));
}

template<class PR> Void
TestFloatValue<PR>::test_predicates()
{
    PR pr=precision;
    ARIADNE_TEST_BINARY_PREDICATE(eq,FloatValue<PR>(-1,pr),FloatValue<PR>(-1,pr));
    ARIADNE_TEST_BINARY_PREDICATE(not eq,FloatValue<PR>(-1,pr),FloatValue<PR>(-2,pr));
    ARIADNE_TEST_BINARY_PREDICATE(lt,FloatValue<PR>(-2,pr),FloatValue<PR>(-1,pr));
    ARIADNE_TEST_BINARY_PREDICATE(not lt,FloatValue<PR>(-2,pr),FloatValue<PR>(-2,pr));
    ARIADNE_TEST_BINARY_PREDICATE(not lt,FloatValue<PR>(-1,pr),FloatValue<PR>(-2,pr));

    Dyadic w(3,1u);
    Integer zl(1), zu(2);
    Dyadic wl(5,2u), wu(7,2u);
    Rational ql(4,3), qu(5,3);

    ARIADNE_TEST_EQUALS(cmp(FloatValue<PR>(w,pr),w),cmp(w,w));
    ARIADNE_TEST_EQUALS(cmp(FloatValue<PR>(w,pr),zl),cmp(w,zl));
    ARIADNE_TEST_EQUALS(cmp(FloatValue<PR>(w,pr),zu),cmp(w,zu));
    ARIADNE_TEST_EQUALS(cmp(FloatValue<PR>(w,pr),wl),cmp(w,wl));
    ARIADNE_TEST_EQUALS(cmp(FloatValue<PR>(w,pr),wu),cmp(w,wu));
    ARIADNE_TEST_EQUALS(cmp(FloatValue<PR>(w,pr),ql),cmp(w,ql));
    ARIADNE_TEST_EQUALS(cmp(FloatValue<PR>(w,pr),qu),cmp(w,qu));
}


Int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    TestFloatApproximation<DoublePrecision>(dp).test();
    TestFloatApproximation<MultiplePrecision>(MultiplePrecision(128)).test();

    TestDirectedFloats<DoublePrecision>(dp).test();
    TestDirectedFloats<MultiplePrecision>(MultiplePrecision(128)).test();

    TestFloatBounds<DoublePrecision>(dp).test();
    TestFloatBounds<MultiplePrecision>(MultiplePrecision(128)).test();

    TestFloatBall<DoublePrecision,DoublePrecision>(dp).test();
    TestFloatBall<MultiplePrecision>(MultiplePrecision(128)).test();
    TestFloatBall<MultiplePrecision,DoublePrecision>(MultiplePrecision(128)).test();

    TestFloatValue<DoublePrecision>(dp).test();
    TestFloatValue<MultiplePrecision>(MultiplePrecision(128)).test();

    return ARIADNE_TEST_FAILURES;
}

