/***************************************************************************
 *            test_dyadic.cpp
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

#include "numeric/twoexp.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/builtin.hpp"
#include "numeric/integer.hpp"
#include "numeric/decimal.hpp"
#include "numeric/logical.hpp"

#include <iomanip>

#include "test.hpp"

using namespace std;
using namespace Ariadne;


class TestDyadic
{
  public:
    void test();
  private:
    void test_concept();
    void test_literal();
    void test_conversions();
    void test_arithmetic();
    void test_comparisons();
    void test_infinity();
};

void TestDyadic::test()
{
    ARIADNE_TEST_CALL(test_literal());
    ARIADNE_TEST_CALL(test_conversions());
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(test_comparisons());
    ARIADNE_TEST_CALL(test_infinity());
}

void TestDyadic::test_concept() {
    unsigned int m=1; unsigned long int lm=1; int n=-2; long int ln=-2; Integer z=-5;
    Dyadic w;

    w=Dyadic(); w=Dyadic(m); w=Dyadic(lm); w=Dyadic(n); w=Dyadic(ln); w=Dyadic(z); w=Dyadic(z);
    w=m; w=lm; w=n; w=ln; w=z; w=w;

    w=+w; w=-w;
    w=w+w; w=w-w; w=w*w;

    w=w+n; w=w-n; w=w*n;
    w=n+w; w=n-w; w=n*w;
    w=w+z; w=w-z; w=w*z;
    w=z+w; w=z-w; w=z*w;

    w=max(w,w); w=min(w,w); w=abs(w);
    w=pos(w); w=neg(w); w=sqr(w); w=hlf(w);

    w=1.5_q2; w=-1.375_q2;

    w==w; w!=w; w<=w; w>=w; w<w; w>w;
    w==n; w!=n; w<=n; w>=n; w<n; w>n;
    n==w; n!=w; n<=w; n>=w; n<w; n>w;
    w==z; w!=z; w<=z; w>=z; w<z; w>z;
    z==w; z!=w; z<=w; z>=w; z<w; z>w;
//    z+1.0;
}

void TestDyadic::test_literal() {
    ARIADNE_TEST_CONSTRUCT(Dyadic,q,(3.25_q2));
    ARIADNE_TEST_EQUALS(q,Dyadic(13,2u));
    ARIADNE_TEST_EQUALS(3.25_q2,Dyadic(13,2u));
    ARIADNE_TEST_EQUALS(-11.375_q2,Dyadic(-91,3u));
    ARIADNE_TEST_EQUALS(0.375_q2,Dyadic(3,3u));

    ARIADNE_TEST_EQUALS(_2^3,Dyadic(8));
    ARIADNE_TEST_EQUALS(_2^-3,Dyadic(1,3u));
    ARIADNE_TEST_EQUALS(+(_2^-4),Dyadic(1,4u));
    ARIADNE_TEST_EQUALS(-(_2^-4),Dyadic(-1,4u));
    ARIADNE_TEST_EQUALS(5*(_2^-3),Dyadic(5,3u));
    ARIADNE_TEST_EQUALS(5/(_2^3),Dyadic(5,3u));
    ARIADNE_TEST_EQUALS(5/(_2^-3),Dyadic(40));
// The following should not compile, since ^ binds less tightly than -,*,/
//    ARIADNE_TEST_EQUALS(-_2^-4,Dyadic(-1,4u));
//    ARIADNE_TEST_EQUALS(5*_2^-4,Dyadic(5,4u));
//    ARIADNE_TEST_EQUALS(5/_2^4,Dyadic(5,4u));
}

void TestDyadic::test_conversions() {
    ARIADNE_TEST_EQUAL(Dyadic(Integer(-3)),Dyadic(-3,0u));
    ARIADNE_TEST_EQUAL(Dyadic(Dyadic(-13)),Dyadic(-13));
    ARIADNE_TEST_EQUAL(Dyadic(Dyadic(-13,3u)),Dyadic(-13,3u));
}

void TestDyadic::test_arithmetic() {
    ARIADNE_TEST_EQUAL(Dyadic(3,2u)+Dyadic(-5,3u),Dyadic(1,3u));
    ARIADNE_TEST_EQUAL(Dyadic(3,2u)-Dyadic(-5,3u),Dyadic(11,3u));
    ARIADNE_TEST_EQUAL(Dyadic(3,2u)*Dyadic(-5,3u),Dyadic(-15,5u));
};

void TestDyadic::test_comparisons() {
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Dyadic(-4,3u),Dyadic(-1,2u));
};

void TestDyadic::test_infinity() {
    double dblmax=std::numeric_limits<double>::max();
    double dblinf=std::numeric_limits<double>::infinity();
    ARIADNE_TEST_CONSTRUCT(Dyadic,inf,(dblinf));
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Dyadic(dblmax),inf);

    ARIADNE_TEST_BINARY_PREDICATE(operator<,Dyadic(4,3u),inf);
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Dyadic(0),inf);
    ARIADNE_TEST_BINARY_PREDICATE(operator<,-inf,Dyadic(0));
    ARIADNE_TEST_BINARY_PREDICATE(operator<,-inf,Dyadic(-4,3u));
    ARIADNE_TEST_BINARY_PREDICATE(operator<,-inf,inf);
    ARIADNE_TEST_BINARY_PREDICATE(operator<=,inf,inf);
    ARIADNE_TEST_BINARY_PREDICATE(operator==,inf,inf);
};


int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    ARIADNE_TEST_CLASS(Dyadic,TestDyadic());

    return ARIADNE_TEST_FAILURES;
}
