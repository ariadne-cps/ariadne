/***************************************************************************
 *            test_rational.cc
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

#include "config.h"

#include "numeric/rational.h"
#include "numeric/integer.h"
#include "numeric/dyadic.h"
#include "numeric/logical.h"

#include <iomanip>

#include "test.h"

using namespace std;
using namespace Ariadne;


class TestRational
{
  public:
    void test();
  private:
    void test_concept();
    void test_literal();
    void test_conversions();
    void test_arithmetic();
};

void TestRational::test()
{
    ARIADNE_TEST_CALL(test_literal());
    ARIADNE_TEST_CALL(test_conversions());
    ARIADNE_TEST_CALL(test_arithmetic());
}

void TestRational::test_concept() {
    unsigned int m=1; unsigned long int lm=1; int n=-2; long int ln=-2; Integer z=-5;
    Rational q;

    q=Rational(); q=Rational(m); q=Rational(lm); q=Rational(n); q=Rational(ln); q=Rational(z); q=Rational(z);
    q=m; q=lm; q=n; q=ln; q=z; q=q;

    q=+q; q=-q;
    q=q+q; q=q-q; q=q*q; q=q/q;

    q=q+n; q=q-n; q=q*n; q=q/n;
    q=n+q; q=n-q; q=n*q; q=n/q;
    q=q+z; q=q-z; q=q*z; q=q/z;
    q=z+q; q=z-q; q=z*q; q=z/q;

    q=max(q,q); q=min(q,q); q=abs(q);
    q=pos(q); q=neg(q); q=sqr(q); q=rec(q);

    q=1.5_q; q=-1.3_q;

    q==q; q!=q; q<=q; q>=q; q<q; q>q;
    q==n; q!=n; q<=n; q>=n; q<n; q>n;
    n==q; n!=q; n<=q; n>=q; n<q; n>q;
    q==z; q!=z; q<=z; q>=z; q<z; q>z;
    z==q; z!=q; z<=q; z>=q; z<q; z>q;

//    z+1.0;
}

void TestRational::test_literal() {
    ARIADNE_TEST_CONSTRUCT(Rational,q,(3.25_q));
    ARIADNE_TEST_EQUALS(q,Rational(13,4));
    ARIADNE_TEST_EQUALS(3.25_q,Rational(13,4));
    ARIADNE_TEST_EQUALS(-11.375_q,Rational(-91,8));
    ARIADNE_TEST_EQUALS(10.3_q,Rational(103,10));
    ARIADNE_TEST_EQUALS(0.333333333333333333_q,Rational(1,3));
    ARIADNE_TEST_EQUALS(0.2857142857142857_q,Rational(2,7));
    ARIADNE_TEST_FAIL(0.453591850358036834_q);
    ARIADNE_TEST_FAIL(3.1415926535897931_q);
}

void TestRational::test_conversions() {
    ARIADNE_TEST_EQUAL(Rational(Integer(-3)),Rational(-3,1));
    ARIADNE_TEST_EQUAL(Rational(Dyadic(-13)),Rational(-13));
    ARIADNE_TEST_EQUAL(Rational(Dyadic(-13,3u)),Rational(-13,8));
}

void TestRational::test_arithmetic() {
    ARIADNE_TEST_EQUAL(Rational(4,5)+Rational(-2,7),Rational(18,35));
    ARIADNE_TEST_EQUAL(Rational(-4,5)-Rational(-2,7),Rational(-18,35));
    ARIADNE_TEST_EQUAL(Rational(4,5)*Rational(-2,7),Rational(-8,35));
    ARIADNE_TEST_EQUAL(Rational(4,5)/Rational(-2,7),Rational(-14,5));
};



int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    ARIADNE_TEST_CLASS(Rational,TestRational());

    return ARIADNE_TEST_FAILURES;
}
