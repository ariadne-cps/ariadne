/***************************************************************************
 *            test_real.cc
 *
 *  Copyright 2013--17  Pieter Collins
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

#include "utility/module.h"
#include "config.h"

#include "numeric/logical.h"
#include "numeric/real.h"

#include "numeric/integer.h"
#include "numeric/dyadic.h"
#include "numeric/rational.h"

#include "numeric/float.h"

#include "test.h"

using namespace std;
using namespace Ariadne;

auto operator==(const Real& x1, double x2) -> decltype(x1==Float64Value(x2)) {
    return Float64Bounds(x1,Precision64())==Float64Bounds(x2);
}

class TestReal
{
  public:
    void test();
  private:
    void test_concept();
    void test_constructors();
    void test_conversions();
    void test_arithmetic();
    void test_transcendental();
    void test_comparison();
    void test_accuracy();
};

void TestReal::test()
{
    Float64Approximation::set_output_places(18);
//    ARIADNE_TEST_CALL(test_concept());
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_conversions());
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(test_transcendental());
    ARIADNE_TEST_CALL(test_comparison());
    ARIADNE_TEST_CALL(test_accuracy());
}

void TestReal::test_concept() {
    Real x,y;
    x=Real(); x=Real(1); x=Real(1.0_exact);
    x=Real(3.14159,3.141593,3.14160);
    x=Real(1);
    y=+x; y=-x; y=x+x; y=x-x; y=x*x; y=x/x;
    y=pow(x,2u); y=pow(x,2);
    y=sqr(x); y=rec(x); y=sqrt(x);
    y=exp(x); y=log(x);
    y=sin(x); y=cos(x); y=tan(x);
}

void TestReal::test_conversions() {
    Real one=1;
    Real pi=4*atan(one);
    auto pi_dp=pi(Precision64());
//    auto pi_mp=pi(PrecisionMP(2));
    ARIADNE_TEST_PRINT(pi_dp);
//    ARIADNE_TEST_PRINT(pi_mp);

}

void TestReal::test_constructors() {
    Precision64 pr;
    ARIADNE_TEST_CONSTRUCT(Real,xv, );
    ARIADNE_TEST_EQUALS(xv.get(pr),0);
    ARIADNE_TEST_EQUALS(xv.lower().get(pr).raw(),0);
    ARIADNE_TEST_EQUALS(xv.upper().get(pr).raw(),0);
    ARIADNE_TEST_CONSTRUCT(Real,xz,(1));
    ARIADNE_TEST_EQUALS(xz.get(pr),1);
    ARIADNE_TEST_CONSTRUCT(Real,xe,(1.5_exact));
    ARIADNE_TEST_EQUALS(xe.get(pr),1.5);
    ARIADNE_TEST_CONSTRUCT(Real,xlau,(3.14159,3.141593,3.14160));
    ARIADNE_TEST_CONSTRUCT(Real,xn,(1.1_q));
    ARIADNE_TEST_COMPARE(Rational(xn.lower().get(pr).raw()),<,Rational(11,10));
    ARIADNE_TEST_COMPARE(Rational(xn.upper().get(pr).raw()),>,Rational(11,10));
    ARIADNE_TEST_CONSTRUCT(Real,xq,(Rational(11,10)));
    ARIADNE_TEST_COMPARE(Rational(xq.lower().get(pr).raw()),<,Rational(11,10));
    ARIADNE_TEST_COMPARE(Rational(xq.upper().get(pr).raw()),>,Rational(11,10));
}

void TestReal::test_arithmetic() {
    Float64Approximation::set_output_places(18);
    Real x(2.5_dyadic);
    Real y(4.0_dyadic);
    ARIADNE_TEST_EQUALS(x, 2.5);
    ARIADNE_TEST_EQUALS(y, 4.0);
    ARIADNE_TEST_EQUALS(+x, 2.5);
    ARIADNE_TEST_EQUALS(-x,-2.5);
    ARIADNE_TEST_EQUALS(x+y, 6.5);
    ARIADNE_TEST_EQUALS(x-y,-1.5);
    ARIADNE_TEST_EQUALS(x*y,10.0);
    ARIADNE_TEST_EQUALS(x/y,0.625);
    ARIADNE_TEST_EQUALS(add(x,y), 6.5);
    ARIADNE_TEST_EQUALS(sub(x,y),-1.5);
    ARIADNE_TEST_EQUALS(mul(x,y),10.0);
    ARIADNE_TEST_EQUALS(div(x,y),0.625);
    ARIADNE_TEST_EQUALS(pow(x,3u),15.625);
    ARIADNE_TEST_EQUALS(pow(x,3),15.625);
    ARIADNE_TEST_EQUALS(pos(x),+2.5);
    ARIADNE_TEST_EQUALS(neg(x),-2.5);
    ARIADNE_TEST_EQUALS(sqr(x),6.25);
    ARIADNE_TEST_EQUALS(rec(y),0.25);
}

void TestReal::test_transcendental() {
    Float64Approximation eps{Float64::eps(Precision64())};
    Real x(2.5_dyadic);
    Float64Approximation ax=x.get(Precision64());
    ARIADNE_TEST_EQUALS(sqrt(Real(4)),2.0);
    ARIADNE_TEST_EQUALS(exp(Real(0)),1.0);
    ARIADNE_TEST_EQUALS(log(Real(1)),0.0);
    ARIADNE_TEST_WITHIN(sqrt(x),sqrt(ax),eps);
    ARIADNE_TEST_WITHIN(exp(x),exp(ax),8*eps);
    ARIADNE_TEST_WITHIN(log(x),log(ax),eps);
    ARIADNE_TEST_WITHIN(sin(x),sin(ax),eps);
    ARIADNE_TEST_WITHIN(cos(x),cos(ax),eps);
    ARIADNE_TEST_WITHIN(tan(x),tan(ax),eps);
    //ARIADNE_TEST_WITHIN(atan(x),atan(ax),eps);
}

void TestReal::test_comparison() {
    Effort effort(0);
    Real one=1;
    Real e=exp(one);
    Real pi=Ariadne::pi;
    Real elogpi=e*log(pi);
//    ARIADNE_TEST_SAME(e,e);
//    ARIADNE_TEST_EQUALS(log(e),one);
    ARIADNE_TEST_ASSERT(possibly(check(log(e)==one,effort)));
    ARIADNE_TEST_ASSERT(possibly((log(e)==one).check(effort)));

    ARIADNE_TEST_PRINT(pi.get(Precision64()));
    ARIADNE_TEST_PRINT(e*log(pi).get(Precision64()));

    ARIADNE_TEST_ASSERT(not check(pi==elogpi,effort));
    ARIADNE_TEST_ASSERT(check(pi!=elogpi,effort));
    ARIADNE_TEST_ASSERT(check(elogpi< pi,effort));
    ARIADNE_TEST_ASSERT(check(elogpi<=pi,effort));
    ARIADNE_TEST_ASSERT(not check(elogpi> pi,effort));
    ARIADNE_TEST_ASSERT(not check(elogpi>=pi,effort));

    Rational zero=0;
    ARIADNE_TEST_ASSERT(check(e*log(pi)-pi<zero,effort));

}


namespace Ariadne {
Bool operator>=(FloatMP const& x1, Float64 x2) { return x1>=x2.get_d(); }
Bool operator<=(FloatMP const& x1, Float64 x2) { return x1<=x2.get_d(); }
Bool operator<=(Float64 x1, FloatMP const& x2) { return x2>=x1.get_d(); }
Bool operator>=(Float64 x1, FloatMP const& x2) { return x2<=x1.get_d(); }
}

void TestReal::test_accuracy() {
    Float64Bounds::set_output_places(18);
    Real one=1;
    Real pi=4*atan(one);

    PrecisionMP mp_high(320);
    RawFloatMP pi_near = FloatMP::pi(mp_high,FloatMP::to_nearest);

    Precision64 dp;
    ARIADNE_TEST_CONSTRUCT(Float64Bounds,pi_dp,(pi.get(dp)));

    PrecisionMP mp(128);
    ARIADNE_TEST_CONSTRUCT(FloatMPBounds,pi_mp,(pi.get(mp)));

    ARIADNE_TEST_ASSERT(pi_dp.lower_raw()<=pi_near);
    ARIADNE_TEST_ASSERT(pi_dp.upper_raw()>=pi_near);

    ARIADNE_TEST_ASSERT(pi_mp.lower_raw()<=pi_near);
    ARIADNE_TEST_ASSERT(pi_mp.upper_raw()>=pi_near);

    Accuracy accuracy{256};
    FloatMPValue error(accuracy.error(),PrecisionMP(320));
    ARIADNE_TEST_CONSTRUCT(FloatMPBounds,pi_met,(pi.evaluate(accuracy)));
    ARIADNE_TEST_PRINT(pi_met.error());
    ARIADNE_TEST_PRINT(abs(sub_up(pi_met.value_raw(),pi_near)));
    ARIADNE_TEST_ASSERT(pi_met.error() <= error);
    ARIADNE_TEST_ASSERT(rad_up(pi_met.value().raw(),pi_near) <= error.raw());
}


int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    TestReal().test();

    return ARIADNE_TEST_FAILURES;
}
