/***************************************************************************
 *            test_real.cpp
 *
 *  Copyright  2013-20  Pieter Collins
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

#include "utility/module.hpp"
#include "config.hpp"

#include "numeric/logical.hpp"
#include "numeric/builtin.hpp"
#include "numeric/real.hpp"

#include "numeric/integer.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/decimal.hpp"
#include "numeric/rational.hpp"

#include "numeric/float_literals.hpp"
#include "numeric/float_bounds.hpp"

namespace Ariadne {

Rational to_rational(String x);
Rational operator"" _q (const char* str, std::size_t);

Rational to_rational(String x) {
    Rational q=0;
    SizeType i=0;
    while(x[i]!=char(0) && x[i]!='.') {
        short d=x[i]-'0';
        q=10*q+d;
        ++i;
    }
    if (x[i]=='.') { ++i; }
    Rational r=1;
    while(x[i]!=char(0)) {
        short d=x[i]-'0';
        r/=10;
        q=q+d*r;
        ++i;
    }
    return q;
}

Rational operator"" _q (const char* str, std::size_t) { return Rational(Decimal(String(str))); }
Decimal operator"" _dec (const char* str, std::size_t) { return Decimal(String(str)); }

//Boolean nondeterministic_greater(Real const& x, Rational const& a, Rational const& b);
//template<class Q> Boolean nondeterministic_greater(Real const& x, Q a, Q b) { return nondeterministic_greater(x,Rational(a),Rational(b)); }

} // namespace Ariadne

#include "numeric/float.hpp"

#include "../test.hpp"

using namespace std;
using namespace Ariadne;

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
    void test_sequence();
};

void TestReal::test()
{
    FloatDPApproximation::set_output_places(18);
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_conversions());
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(test_transcendental());
    ARIADNE_TEST_CALL(test_comparison());
    ARIADNE_TEST_CALL(test_accuracy());
    ARIADNE_TEST_CALL(test_sequence());
}

void TestReal::test_concept() {
    Real x,y;
    x=Real(); x=Real(1); x=Real(1.0_x);
    x=Real(1);
    y=+x; y=-x; y=x+x; y=x-x; y=x*x; y=x/x;
    y=pow(x,2u); y=pow(x,2);
    y=sqr(x); y=rec(x); y=sqrt(x);
    y=exp(x); y=log(x);
    y=sin(x); y=cos(x); y=tan(x);
}

void TestReal::test_conversions() {
    Real one=1;
    Real pi_=4*atan(one);
    auto pi_dp=pi_.get(dp);
//    auto pi_mp=pi(MultiplePrecision(2));
    ARIADNE_TEST_PRINT(pi_dp);
//    ARIADNE_TEST_PRINT(pi_mp);

}

void TestReal::test_constructors() {
    ARIADNE_TEST_CONSTRUCT(Real,xv, );
    ARIADNE_TEST_EQUALS(xv.get(dp),0);
    ARIADNE_TEST_EQUALS(xv.lower().get(dp).raw(),0);
    ARIADNE_TEST_EQUALS(xv.upper().get(dp).raw(),0);
    ARIADNE_TEST_CONSTRUCT(Real,xz,(1));
    ARIADNE_TEST_EQUALS(xz.get(dp),1);
    ARIADNE_TEST_CONSTRUCT(Real,xe,(1.5_x));
    ARIADNE_TEST_EQUALS(xe.get(dp),1.5);
    ARIADNE_TEST_CONSTRUCT(Real,xn,(1.1_q));
    ARIADNE_TEST_COMPARE(Rational(xn.lower().get(dp).raw()),<,Rational(11,10));
    ARIADNE_TEST_COMPARE(Rational(xn.upper().get(dp).raw()),>,Rational(11,10));
    ARIADNE_TEST_CONSTRUCT(Real,xq,(Rational(11,10)));
    ARIADNE_TEST_COMPARE(Rational(xq.lower().get(dp).raw()),<,Rational(11,10));
    ARIADNE_TEST_COMPARE(Rational(xq.upper().get(dp).raw()),>,Rational(11,10));
}

void TestReal::test_arithmetic() {
    FloatDPApproximation::set_output_places(18);
    Real x(2.5_dyadic);
    Real y(4.0_dyadic);
    ARIADNE_TEST_EQUALS(x, 2.5_dy);
    ARIADNE_TEST_EQUALS(y, 4.0_x);
    ARIADNE_TEST_EQUALS(+x, 2.5_x);
    ARIADNE_TEST_EQUALS(-x,-2.5_x);
    ARIADNE_TEST_EQUALS(x+y, 6.5_x);
    ARIADNE_TEST_EQUALS(x-y,-1.5_x);
    ARIADNE_TEST_EQUALS(x*y,10.0_x);
    ARIADNE_TEST_EQUALS(x/y,0.625_dy);
    ARIADNE_TEST_EQUALS(add(x,y), 6.5_dy);
    ARIADNE_TEST_EQUALS(sub(x,y),-1.5_dy);
    ARIADNE_TEST_EQUALS(mul(x,y),10.0_dy);
    ARIADNE_TEST_EQUALS(div(x,y),0.625_dy);
    ARIADNE_TEST_EQUALS(pow(x,3u),15.625_dy);
    ARIADNE_TEST_EQUALS(pow(x,3),15.625_dy);
    ARIADNE_TEST_EQUALS(pos(x),+2.5_dy);
    ARIADNE_TEST_EQUALS(neg(x),-2.5_dy);
    ARIADNE_TEST_EQUALS(hlf(x),1.25_dy);
    ARIADNE_TEST_EQUALS(sqr(x),6.25_dy);
    ARIADNE_TEST_EQUALS(rec(y),0.25_dy);
}

void TestReal::test_transcendental() {
    Dyadic eps{FloatDP::eps(dp)};
    Real x(2.5_dyadic);
    FloatDPApproximation ax=x.get(dp);
    ARIADNE_TEST_EQUALS(sqrt(Real(4)),2.0_dy);
    ARIADNE_TEST_EQUALS(exp(Real(0)),1.0_dy);
    ARIADNE_TEST_EQUALS(log(Real(1)),0.0_dy);
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
    Real elogpi=e*log(pi);
//    ARIADNE_TEST_SAME(e,e);
//    ARIADNE_TEST_EQUALS(log(e),one);
    ARIADNE_TEST_ASSERT(possibly(check(log(e)==one,effort)));
    ARIADNE_TEST_ASSERT(possibly((log(e)==one).check(effort)));

    ARIADNE_TEST_PRINT(pi.get(dp));
    ARIADNE_TEST_PRINT(e*log(pi).get(dp));

    ARIADNE_TEST_ASSERT(not check(pi==elogpi,effort));
    ARIADNE_TEST_ASSERT(check(pi!=elogpi,effort));
    ARIADNE_TEST_ASSERT(check(elogpi< pi,effort));
    ARIADNE_TEST_ASSERT(check(elogpi<=pi,effort));
    ARIADNE_TEST_ASSERT(not check(elogpi> pi,effort));
    ARIADNE_TEST_ASSERT(not check(elogpi>=pi,effort));

    Rational zero=0;
    ARIADNE_TEST_ASSERT(check(e*log(pi)-pi<zero,effort));

}


void TestReal::test_accuracy() {
    FloatDPBounds::set_output_places(18);
    Real one=1;
    Real pi_=4*atan(one);

    MultiplePrecision mp_high(320);
    RawFloatMP pi_near = FloatMP::pi(near,mp_high);

    ARIADNE_TEST_CONSTRUCT(FloatDPBounds,pi_dp,(pi_.get(dp)));

    MultiplePrecision mp(128);
    ARIADNE_TEST_CONSTRUCT(FloatMPBounds,pi_mp,(pi_.get(mp)));

    ARIADNE_TEST_ASSERT(pi_dp.lower_raw()<=pi_near);
    ARIADNE_TEST_ASSERT(pi_dp.upper_raw()>=pi_near);

    ARIADNE_TEST_ASSERT(pi_mp.lower_raw()<=pi_near);
    ARIADNE_TEST_ASSERT(pi_mp.upper_raw()>=pi_near);

    mp=precision(320);
    Effort effort{256};
    ARIADNE_TEST_CONSTRUCT(ValidatedReal,pi_ord,(pi_.compute(effort)));
    Accuracy accuracy{256};
    FloatMPValue error(accuracy.error(),mp);
    ARIADNE_TEST_CONSTRUCT(ValidatedReal,pi_met,(pi_.compute(accuracy)));
    ARIADNE_TEST_CONSTRUCT(FloatMPBounds,pi_met_mp,(pi_met.get(mp)));
    ARIADNE_TEST_PRINT(pi_met_mp.error());
    ARIADNE_TEST_PRINT(abs(sub(up,pi_met_mp.value_raw(),pi_near)));
    ARIADNE_TEST_ASSERT(pi_met_mp.error() <= error);
    ARIADNE_TEST_ASSERT(rad(up,pi_met_mp.value().raw(),pi_near) <= error.raw());

    Dyadic eps(1,1024u);
    ARIADNE_TEST_CONSTRUCT(Real,sin_pi,(sin(pi_)));
    ARIADNE_TEST_PRINT(choose(sin_pi>-eps,sin_pi<eps));
/*
    ARIADNE_TEST_PRINT(nondeterministic_greater(sin_pi,-eps,eps));
    ARIADNE_TEST_PRINT(nondeterministic_greater(sin_pi,-1,eps));
    ARIADNE_TEST_PRINT(nondeterministic_greater(sin_pi,-eps,1));
    ARIADNE_TEST_ASSERT(not nondeterministic_greater(sin_pi,eps,2*eps));
    ARIADNE_TEST_ASSERT(nondeterministic_greater(sin_pi,-3*eps,-2*eps));
*/
}


void TestReal::test_sequence() {
    std::function<Dyadic(Natural)> wfn([&](Natural n){return 1-Dyadic(1,n);});
    FastCauchySequence<Dyadic> wseq(wfn);
    Real wlim=limit(wseq);
    std::cout<<wlim.compute(Accuracy(256))<<"\n";

    std::function<Real(Natural)> rfn([&](Natural n){return exp(Real(-(n+1u)));});
    FastCauchySequence<Real> rseq(rfn);
    Real rlim=limit(rseq);

    ARIADNE_TEST_ASSERT(abs(rlim.compute(Accuracy(256)).get())<Dyadic(1,256u));

    Real pi_=4*atan(Real(1));
    Decimal pi_lower="3.1415926535897932"_dec;
    Decimal pi_upper="3.1415926535897933"_dec;
    ARIADNE_TEST_PRINT(nondeterministic_greater(pi_,"3.1415926535897932"_dec,"3.1415926535897933"_dec));
    ARIADNE_TEST_ASSERT(abs(rlim.compute(Accuracy(256)).get())<Dyadic(1,256u));
    ARIADNE_TEST_ASSERT(nondeterministic_greater(pi_,3.0_dec,3.1_dec));
    ARIADNE_TEST_ASSERT(not nondeterministic_greater(pi_,3.2_dec,3.3_dec));
    ARIADNE_TEST_PRINT(nondeterministic_greater(pi_,"3.1415926535897932"_dec,"3.1415926535897937"_dec));
    ARIADNE_TEST_PRINT(nondeterministic_greater(pi_,"3.1415926535897932"_dec,"3.1415926535897933"_dec));

    ARIADNE_TEST_PRINT(pi_.compute(Accuracy(256)).get(MultiplePrecision(256)));
    ARIADNE_TEST_PRINT(std::boolalpha);
    ARIADNE_TEST_PRINT(choose(pi_>pi_lower,pi_<pi_upper));
    ARIADNE_TEST_PRINT(choose(pi_>pi_lower,pi_<3.2_dec));

    ARIADNE_TEST_PRINT(choose({pi_>pi_lower,Real(4)},{pi_<3.2_dec,Real(3)}));
    ARIADNE_TEST_PRINT(choose({pi_>pi_lower,Real(4)},{pi_<pi_upper,Real(3)}));

    Dyadic pi_upper_bound = pi_.compute(Accuracy(256)).get().upper_raw();
    Real x=pi_-pi_upper_bound;
    ARIADNE_TEST_PRINT((x>=0).check(Effort(5u)));
    ARIADNE_TEST_PRINT((x>=0).check(Effort(6u)));
//    ARIADNE_TEST_PRINT(unify(x>=0,+x,x<=0,-x));
    ARIADNE_TEST_PRINT(when({x>=0,+x},{x<=0,-x}).compute(Effort(5u)));
    ARIADNE_TEST_PRINT(when({x>=0,+x},{x<=0,-x}).compute(Effort(6u)));
    ARIADNE_TEST_PRINT(when({x>=0,+x},{x<=0,-x}).compute(Effort(384u)));
    ARIADNE_TEST_PRINT(when({x>=0,+x},{x<=0,-x}).get(precision(5u)));
    ARIADNE_TEST_PRINT(when({x>=0,+x},{x<=0,-x}).get(precision(6u)));
    ARIADNE_TEST_PRINT(when({x>=0,+x},{x<=0,-x}).get(precision(54)));
    ARIADNE_TEST_PRINT(when({x>=0,+x},{x<=0,-x}).get(precision(64)));
    ARIADNE_TEST_FAIL(when({x>=0,+x},{x<=0,1+x}).compute(Effort(5u)));
    ARIADNE_TEST_PRINT(when({x>=0,+x},{x<=0,1+x}).compute(Effort(6u)));


}


int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    TestReal().test();

    return ARIADNE_TEST_FAILURES;
}
