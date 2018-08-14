/***************************************************************************
 *            test_rational.cpp
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

#include "config.hpp"

#include "numeric/rational.hpp"
#include "numeric/builtin.hpp"
#include "numeric/integer.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/decimal.hpp"
#include "numeric/logical.hpp"

#include <iomanip>

#include "test.hpp"

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
    void test_comparisons();

    void test_decimal();

};

void TestRational::test()
{
    ARIADNE_TEST_CALL(test_concept());
    ARIADNE_TEST_CALL(test_literal());
    ARIADNE_TEST_CALL(test_conversions());
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(test_comparisons());

    ARIADNE_TEST_CALL(test_decimal());
}

void TestRational::test_concept() {
    unsigned int m=1; unsigned long int lm=1; int n=-2; long int ln=-2; Integer z=-5;
    Rational q;

    q=Rational(); q=Rational(m); q=Rational(lm); q=Rational(n); q=Rational(ln); q=Rational(z);
    q=m; q=lm; q=n; q=ln; q=z; q=q;

    q=+q; q=-q;

    q=q/q; q=q+q; q=q-q; q=q*q;

    q=q+n; q=q-n; q=q*n; q=q/n;
    q=n-q; q=n*q; q=n+q; q=n/q;
    q=q+z; q=q-z; q=q*z; q=q/z;
    q=z+q; q=z-q; q=z*q; q=z/q;

    q=max(q,q); q=min(q,q); q=abs(q);
    q=pos(q); q=neg(q); q=sqr(q); q=rec(q);

    q=1.5_q; q=-1.3_q;

    ARIADNE_TEST_ASSERT(q==q);
    ARIADNE_TEST_ASSERT(not(q!=q));
    ARIADNE_TEST_ASSERT(q<=q);
    ARIADNE_TEST_ASSERT(q>=q);
    ARIADNE_TEST_ASSERT(not(q<q));
    ARIADNE_TEST_ASSERT(not(q>q));
    ARIADNE_TEST_ASSERT(not(q==n));
    ARIADNE_TEST_ASSERT(q!=n);
    ARIADNE_TEST_ASSERT(not(q<=n));
    ARIADNE_TEST_ASSERT(q>=n);
    ARIADNE_TEST_ASSERT(not(q<n));
    ARIADNE_TEST_ASSERT(q>n);
    ARIADNE_TEST_ASSERT(not(n==q));
    ARIADNE_TEST_ASSERT(n!=q);
    ARIADNE_TEST_ASSERT(n<=q);
    ARIADNE_TEST_ASSERT(not(n>=q));
    ARIADNE_TEST_ASSERT(n<q);
    ARIADNE_TEST_ASSERT(not(n>q));
    ARIADNE_TEST_ASSERT(not(q==z));
    ARIADNE_TEST_ASSERT(q!=z);
    ARIADNE_TEST_ASSERT(not(q<=z));
    ARIADNE_TEST_ASSERT(q>=z);
    ARIADNE_TEST_ASSERT(not(q<z));
    ARIADNE_TEST_ASSERT(q>z);
    ARIADNE_TEST_ASSERT(not(z==q));
    ARIADNE_TEST_ASSERT(z!=q);
    ARIADNE_TEST_ASSERT(z<=q);
    ARIADNE_TEST_ASSERT(not(z>=q));
    ARIADNE_TEST_ASSERT(z<q);
    ARIADNE_TEST_ASSERT(not(z>q));
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

void TestRational::test_comparisons() {
    ExactDouble inf=ExactDouble::infinity();
    ExactDouble max=ExactDouble(std::numeric_limits<double>::max());
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Rational(-max),Rational(+max));
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Rational(-max),Rational(-4,5));
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Rational(-max),Rational(2,3));
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Rational(-4,5),Rational(-2,7));
    ARIADNE_TEST_BINARY_PREDICATE(operator<,-inf,Rational(18,35));
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Rational(18,35),+inf);
};

void TestRational::test_decimal() {
    ARIADNE_TEST_CONSTRUCT(Decimal,d1,(23,1u));
    ARIADNE_TEST_CONSTRUCT(Decimal,d2,(-42,2u));
    ARIADNE_TEST_EQUALS(Decimal(-3.14),Decimal(-314,2u));
    ARIADNE_TEST_FAIL(Decimal d(0.33333333333));
    ARIADNE_TEST_EQUALS(Decimal(Dyadic(7,3u)),Decimal(875,3u));
    ARIADNE_TEST_EQUALS(Decimal(Dyadic(140,1u)),Decimal(70,0u));
    ARIADNE_TEST_EQUALS(Decimal("-3.14"),Decimal(-314,2u));
    ARIADNE_TEST_EQUALS(Decimal("-3.1400"),Decimal(-314,2u));
    ARIADNE_TEST_EQUALS(Decimal("-0.0031400"),Decimal(-314,5u));
    ARIADNE_TEST_EQUALS(Decimal("3140"),Decimal(3140,0u));
    ARIADNE_TEST_EQUALS(Decimal("3.14159265")*Decimal("3.14159265"),Decimal("9.8696043785340225"));
    ARIADNE_TEST_PRINT(Decimal("3.14159265")*Decimal("3.14159265"));
    ARIADNE_TEST_EQUALS(Rational(d1),Rational(23,10));
    ARIADNE_TEST_EQUALS(Rational(d2),Rational(-21,50));
    ARIADNE_TEST_EQUAL(Rational(+d2),+Rational(d2));
    ARIADNE_TEST_EQUAL(Rational(-d2),-Rational(d2));
    ARIADNE_TEST_EQUAL(Rational(d1+d2),Rational(d1)+Rational(d2));
    ARIADNE_TEST_EQUAL(Rational(d1-d2),Rational(d1)-Rational(d2));
    ARIADNE_TEST_EQUAL(Rational(d1*d2),Rational(d1)*Rational(d2));
    ARIADNE_TEST_EQUAL(d1/d2,Rational(d1)/Rational(d2));
};



int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    ARIADNE_TEST_CLASS(Rational,TestRational());

    return ARIADNE_TEST_FAILURES;
}
