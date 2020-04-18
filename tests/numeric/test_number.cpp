/***************************************************************************
 *            test_number.cpp
 *
 *  Copyright  2016-20  Alberto Casagrande, Pieter Collins
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

#include "numeric/builtin.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/rational.hpp"
#include "numeric/number.hpp"

#include "numeric/floatdp.hpp"
#include "numeric/float_value.hpp"
#include "numeric/float_ball.hpp"
#include "numeric/float_bounds.hpp"
#include "numeric/float_upper_bound.hpp"
#include "numeric/float_lower_bound.hpp"
#include "numeric/float_approximation.hpp"
#include "numeric/float_error.hpp"

#include "../test.hpp"

using namespace std;
using namespace Ariadne;

class TestNumbers
{
  public:
    Void test();
    Void test_float_value_behaviour();
    Void test_operations();
};

void TestNumbers::test()
{
    ARIADNE_TEST_CALL(test_operations());
    ARIADNE_TEST_CALL(test_float_value_behaviour());
}

Void
TestNumbers::test_float_value_behaviour()
{
    ExactNumber y(FloatDPValue(0,dp));
    try {
        y = y+y;
        ARIADNE_TEST_NOTIFY("Binary operations on FloatValue<DP> within ExactNumber yield "<<y.class_name());
    } catch (const DispatchException& e) {
        ARIADNE_TEST_NOTIFY("Binary operations on FloatValue<DP> give error:\n    "<<e.what());
    }
}

Void
TestNumbers::test_operations()
{
    Int n=1; Integer z=1; FloatDPValue v(3); FloatDPBounds b(3);
    ExactNumber yn=n; ExactNumber yz=z; ExactNumber yv=v; ValidatedNumber yb=b;
    ValidatedErrorNumber en=n; ValidatedErrorNumber ev=v;

    Dyadic w2(2);
    Dyadic w3(3);
    Rational q2(2);
    Rational q3(3);
    FloatDPValue x2(2,dp);
    FloatDPValue x3(3,dp);
    ExactNumber y2(x2), y3(x3);
    ARIADNE_TEST_PRINT(w2/w3);
    ARIADNE_TEST_PRINT(x2/x3);
    ARIADNE_TEST_THROWS(y2/y3, DispatchException);
    ARIADNE_TEST_THROWS(q2/y3, DispatchException);
    ARIADNE_TEST_THROWS(y2/q3, DispatchException);

    ARIADNE_TEST_FAIL(add(yv,yv));
    ARIADNE_TEST_PRINT(add(yb,yb));
    ARIADNE_TEST_PRINT(add(yb,yv));
//    ARIADNE_TEST_PRINT(add(b,z));
    ARIADNE_TEST_PRINT(add(b,yb));
    ARIADNE_TEST_PRINT(add(yb,yn));

    ValidatedErrorNumber ym=max(en,ev);
    ARIADNE_TEST_PRINT(max(en,ev));
    ARIADNE_TEST_PRINT(max(en,ev).pointer());
    ARIADNE_TEST_PRINT(max(en,ev).class_name());
    ARIADNE_TEST_PRINT(max(en,ev));

    ARIADNE_TEST_EXECUTE(add(ExactNumber(1),ExactNumber(FloatDPValue(2))));
    ARIADNE_TEST_FAIL(add(ExactNumber(FloatDPValue(1)),ExactNumber(FloatDPValue(2))));

    ARIADNE_TEST_PRINT(max(ExactNumber(1),ExactNumber(2)));
    ARIADNE_TEST_FAIL(max(ExactNumber(FloatDPValue(1)),ExactNumber(FloatDPValue(2))));
    ARIADNE_TEST_PRINT(max(ExactNumber(1),ExactNumber(FloatDPValue(2))));

    ARIADNE_TEST_PRINT(max(ExactNumber(1),ExactNumber(FloatDPValue(2))));

    max(1,FloatDPValue(3));
//    max(1u,FloatDPError(3u));
}

template<class Y> class TestNumber
{
  public:
    Void test();
  private:
    Void test_concept();
    Void test_class();
    Void test_get();
    Void test_operations();
    Void test_comparisons();
};


Int main() {
    TestNumbers().test();

    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    TestNumber<ApproximateNumber>().test();
    TestNumber<ValidatedNumber>().test();
    return ARIADNE_TEST_FAILURES;
    TestNumber<EffectiveNumber>().test();
    TestNumber<ExactNumber>().test();

    std::cerr<<"INCOMPLETE ";
    return ARIADNE_TEST_FAILURES;
//    TestNumber<ValidatedUpperNumber>().test();
//    TestNumber<ValidatedLowerNumber>().test();
}


template<class Y> Void
TestNumber<Y>::test()
{
    ARIADNE_TEST_CALL(test_class());
    ARIADNE_TEST_CALL(test_get());
    ARIADNE_TEST_CALL(test_comparisons());
}

template<class Y> Void
TestNumber<Y>::test_concept()
{
    Int n;
    Y y;
    y=Y(2);
    y=3;
    Y y2(1);

    y=y+y; y=y-y; y=y*y; y=y/y;
    y+=y2; y-=y2; y*=y2; y/=y2;
    y=abs(y); y=max(y,y); y=min(y,y);
    y=add(y,y); y=sub(y,y); y=mul(y,y); y=div(y,y);
    y=sqrt(y); y=exp(y); y=log(y); y=atan(y);
    y=pos(y); y=neg(y); y=sqr(y); y=rec(y);
    y=sin(y); y=cos(y); y=tan(y);

    y==y; y!=y; y<=y; y>=y; y<y; y>y;
    y==n; y!=n; y<=n; y>=n; y<n; y>n;
    cout << y;
}

template<class F, class FE> Bool models(Ball<F,FE> const& x, Rational const& q) {
    return abs(Dyadic(x.value().raw())-q)<=Dyadic(x.error().raw()); }
template<class F> Bool models(Bounds<F> const& x, Rational const& q) {
    return Dyadic(x.lower().raw())<=q && q <= Dyadic(x.upper().raw()); }
template<class F> Bool models(LowerBound<F> const& x, Rational const& q) {
    return Dyadic(x.raw())<=q; }
template<class F> Bool models(UpperBound<F> const& x, Rational const& q) {
    return q<=Dyadic(x.raw()); }


template<class Y> void test_number_get() {
    Rational q=3/5_q;
    Real r(q);
    Y y(r);
    MultiplePrecision mp(192);
    MultiplePrecision mpe(128);

    MetricTag metric;

    ARIADNE_TEST_WITHIN(y.get(ApproximateTag(),dp),q,Dyadic(TwoExp(-52)));
    ARIADNE_TEST_WITHIN(y.get(ApproximateTag(),mp),q,Dyadic(TwoExp(-128)));

    ARIADNE_TEST_BINARY_PREDICATE(models,y.get(LowerTag(),dp),q);
    ARIADNE_TEST_BINARY_PREDICATE(models,y.get(LowerTag(),mp),q);

    ARIADNE_TEST_BINARY_PREDICATE(models,y.get(UpperTag(),dp),q);
    ARIADNE_TEST_BINARY_PREDICATE(models,y.get(UpperTag(),mp),q);

    ARIADNE_TEST_BINARY_PREDICATE(models,y.get(BoundedTag(),dp),q);
    ARIADNE_TEST_BINARY_PREDICATE(models,y.get(BoundedTag(),mp),q);

    ARIADNE_TEST_BINARY_PREDICATE(models,y.get(MetricTag(),dp),q);
    ARIADNE_TEST_BINARY_PREDICATE(models,y.get(MetricTag(),mp),q);

    ARIADNE_TEST_BINARY_PREDICATE(models,y.get(metric,dp,dp),q);
    ARIADNE_TEST_BINARY_PREDICATE(models,y.get(metric,mp,dp),q);
    ARIADNE_TEST_BINARY_PREDICATE(models,y.get(metric,mp,mpe),q);

    ARIADNE_TEST_SAME_TYPE(decltype(y.get(metric,dp,dp)),FloatDPBall);
    ARIADNE_TEST_SAME_TYPE(decltype(y.get(metric,mp,dp)),FloatMPDPBall);
    ARIADNE_TEST_SAME_TYPE(decltype(y.get(metric,mp,mpe)),FloatMPBall);
}

template<> void test_number_get<ExactNumber>() { }

template<> void test_number_get<ApproximateNumber>() {
    Rational q=3/5_q;
    ApproximateNumber y(q);
    MultiplePrecision mp(128);
    ARIADNE_TEST_WITHIN(y.get(ApproximateTag(),dp),q,3e-16);
    ARIADNE_TEST_WITHIN(y.get(ApproximateTag(),mp),q,3e-16);

    // Regression tests for DP/MP conversion.
    Dbl d=-0.6;
    ARIADNE_TEST_ASSIGN(y,d);
    ARIADNE_TEST_EXECUTE(y.get(ApproximateTag(),dp));
    ARIADNE_TEST_EXECUTE(y.get(ApproximateTag(),mp));

    FloatDPApproximation x(-2.5,double_precision);
    ARIADNE_TEST_ASSIGN(y,x);
    ARIADNE_TEST_EXECUTE(y.get(ApproximateTag(),dp));
    ARIADNE_TEST_EXECUTE(y.get(ApproximateTag(),mp));

}

template<class Y> Void
TestNumber<Y>::test_get()
{
    test_number_get<Y>();

}

inline FloatDPValue max(Int y1, FloatDPValue x2) { return max(FloatDPValue(y1,x2.precision()),x2); }

template<class Y> Void
TestNumber<Y>::test_class()
{
    Int n=1; Integer z=1; FloatDPValue v(3); FloatDPBounds b(3);
    ExactNumber yn=n; ExactNumber yz=z; ExactNumber yv=v; ValidatedNumber yb=b;

    ARIADNE_TEST_EQUAL(yz.class_name(),"Integer");
    ARIADNE_TEST_EQUAL(yv.class_name(),"FloatDPValue");
    ARIADNE_TEST_EQUAL(yb.class_name(),"FloatDPBounds");
}


template<class Y> Void
TestNumber<Y>::test_comparisons() {
}

template<> Void
TestNumber<ExactNumber>::test_comparisons() {
    return;
    ARIADNE_TEST_CONSTRUCT(ExactNumber,y1,(Rational(2,3)));
    ARIADNE_TEST_CONSTRUCT(ExactNumber,y2,(Rational(683,1024)));
    ARIADNE_TEST_CONSTRUCT(ExactNumber,pinf,(+ExactDouble::infinity()));
    ARIADNE_TEST_CONSTRUCT(ExactNumber,ninf,(-ExactDouble::infinity()));

/*
    ARIADNE_TEST_BINARY_PREDICATE( operator==,y1,y1);
    ARIADNE_TEST_BINARY_PREDICATE(!operator!=,y1,y1);
    ARIADNE_TEST_BINARY_PREDICATE( operator<=,y1,y1);
    ARIADNE_TEST_BINARY_PREDICATE( operator>=,y1,y1);
    ARIADNE_TEST_BINARY_PREDICATE(!operator< ,y1,y1);
    ARIADNE_TEST_BINARY_PREDICATE(!operator> ,y1,y1);
*/
    ARIADNE_TEST_BINARY_PREDICATE(!operator==,y1,y2);
    ARIADNE_TEST_BINARY_PREDICATE( operator!=,y1,y2);
    ARIADNE_TEST_BINARY_PREDICATE( operator<=,y1,y2);
    ARIADNE_TEST_BINARY_PREDICATE(!operator>=,y1,y2);
    ARIADNE_TEST_BINARY_PREDICATE( operator< ,y1,y2);
    ARIADNE_TEST_BINARY_PREDICATE(!operator> ,y1,y2);

    ARIADNE_TEST_BINARY_PREDICATE(!operator==,y2,y1);
    ARIADNE_TEST_BINARY_PREDICATE( operator!=,y2,y1);
    ARIADNE_TEST_BINARY_PREDICATE(!operator<=,y2,y1);
    ARIADNE_TEST_BINARY_PREDICATE( operator>=,y2,y1);
    ARIADNE_TEST_BINARY_PREDICATE(!operator< ,y2,y1);
    ARIADNE_TEST_BINARY_PREDICATE( operator> ,y2,y1);

    ARIADNE_TEST_BINARY_PREDICATE( operator==,y2,y2);
    ARIADNE_TEST_BINARY_PREDICATE(!operator< ,y2,y2);

    ARIADNE_TEST_BINARY_PREDICATE(operator==,y1,y1);

    ARIADNE_TEST_BINARY_PREDICATE(operator==,pinf,pinf);
    ARIADNE_TEST_BINARY_PREDICATE(operator!=,ninf,pinf);
    ARIADNE_TEST_BINARY_PREDICATE(operator< ,ninf,pinf);
    ARIADNE_TEST_BINARY_PREDICATE(operator!=,y1,pinf);
    ARIADNE_TEST_BINARY_PREDICATE(operator!=,pinf,y1);
    ARIADNE_TEST_BINARY_PREDICATE(operator<,y1,pinf);
    ARIADNE_TEST_BINARY_PREDICATE(operator<,ninf,y1);

    ARIADNE_TEST_ASSERT(definitely(ValidatedUpperNumber(3)<ValidatedLowerNumber(4)));
    ARIADNE_TEST_ASSERT(possibly(ValidatedUpperNumber(4)>ValidatedLowerNumber(3)));
    ARIADNE_TEST_ASSERT(not definitely(ValidatedUpperNumber(4)>ValidatedLowerNumber(3)));

}
