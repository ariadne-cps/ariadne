/***************************************************************************
 *            test_real.cc
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

#include "utility/module.h"
#include "config.h"

#include "numeric/logical.h"
#include "numeric/real.h"

#include "numeric/integer.h"
#include "numeric/rational.h"

#include "numeric/float.h"
#include "numeric/float-approximate.h"
#include "numeric/float-validated.h"
#include "numeric/float-exact.h"

#include "test.h"

using namespace std;
using namespace Ariadne;

Tribool operator==(const Real& x1, double x2) {
    return BoundFloat(x1)==BoundFloat(x2);
}

class TestReal
{
  public:
    void test();
  private:
    void test_concept();
    void test_constructors();
    void test_arithmetic();
    void test_transcendental();
};

void TestReal::test()
{
    ApproximateFloat::set_output_precision(18);
//    ARIADNE_TEST_CALL(test_concept());
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(test_transcendental());
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

void TestReal::test_constructors() {
    ARIADNE_TEST_CONSTRUCT(Real,xv, );
    ARIADNE_TEST_EQUALS(xv.lower(),0);
    ARIADNE_TEST_EQUALS(xv.approx(),0);
    ARIADNE_TEST_EQUALS(xv.upper(),0);
    ARIADNE_TEST_CONSTRUCT(Real,xz,(1));
    ARIADNE_TEST_EQUALS(xz.lower(),1);
    ARIADNE_TEST_EQUALS(xz.approx(),1);
    ARIADNE_TEST_EQUALS(xz.upper(),1);
    ARIADNE_TEST_CONSTRUCT(Real,xe,(1.5_exact));
    ARIADNE_TEST_EQUALS(xe.lower(),1.5);
    ARIADNE_TEST_EQUALS(xe.approx(),1.5);
    ARIADNE_TEST_EQUALS(xe.upper(),1.5);
    ARIADNE_TEST_CONSTRUCT(Real,xlau,(3.14159,3.141593,3.14160));
    ARIADNE_TEST_CONSTRUCT(Real,xn,(1.1_q));
    ARIADNE_TEST_COMPARE(Rational(make_exact(xn.lower())),<,Rational(11,10));
    ARIADNE_TEST_COMPARE(Rational(xn.upper().raw()),>,Rational(11,10));
    ARIADNE_TEST_CONSTRUCT(Real,xq,(Rational(11,10)));
    ARIADNE_TEST_COMPARE(Rational(xq.lower().raw()),<,Rational(11,10));
    ARIADNE_TEST_COMPARE(Rational(xq.upper().raw()),>,Rational(11,10));
}

void TestReal::test_arithmetic() {
    ApprxFloat::set_output_precision(18);
    Real x(2.5_exact);
    Real y(4.0_exact);
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
    ApproximateFloat eps=std::numeric_limits<double>::epsilon();
    Real x(2.5_exact);
    ApproximateFloat ax=x;
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


int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    TestReal().test();

    return ARIADNE_TEST_FAILURES;
}
