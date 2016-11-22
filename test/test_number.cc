/***************************************************************************
 *            test_number.cc
 *
 *  Copyright  2016  Alberto Casagrande, Pieter Collins
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

#include "numeric/builtin.h"
#include "numeric/rational.h"
#include "numeric/number.h"

#include "test.h"

using namespace std;
using namespace Ariadne;

template<class Y> class TestNumber
{
  public:
    Void test();
  private:
    Void test_concept();
    Void test_class();
    Void test_comparisons();
};


Int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    TestNumber<ApproximateNumber>().test();
    TestNumber<ValidatedNumber>().test();
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
    //ARIADNE_TEST_CALL(test_concept());
    ARIADNE_TEST_CALL(test_class());
    ARIADNE_TEST_CALL(test_comparisons());
}

template<class Y> Void
TestNumber<Y>::test_concept()
{
    Int n;
    Y y;
    y=Y(2);
    y=3;

    y=y+y; y=y-y; y=y*y; y=y/y;
    y+=y; y-=y; y*=y; y/=y;
    y=abs(y); y=max(y,y); y=min(y,y);
    y=add(y,y); y=sub(y,y); y=mul(y,y); y=div(y,y);
    y=sqrt(y); y=exp(y); y=log(y); y=atan(y);
    y=pos(y); y=neg(y); y=sqr(y); y=rec(y);
    y=sin(y); y=cos(y); y=tan(y);

    y==y; y!=y; y<=y; y>=y; y<y; y>y;
    y==n; y!=n; y<=n; y>=n; y<n; y>n;
    cout << y;
}


template<class Y> Void
TestNumber<Y>::test_class()
{
}

template<class Y> Void
TestNumber<Y>::test_comparisons() {
}

template<> Void
TestNumber<ExactNumber>::test_comparisons() {
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

}
