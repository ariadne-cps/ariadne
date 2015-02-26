/***************************************************************************
 *            test_interval.cc
 *
 *  Copyright  2006-14  Alberto Casagrande, Pieter Collins
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

#include "numeric/module.h"

#include "numeric/rational.h"
#include "geometry/interval.h"

#include "test.h"

using namespace Ariadne;
using namespace std;


static_assert(IsConstructible<LowerFloat,double>::value,"");
static_assert(IsConstructible<UpperFloat,double>::value,"");
static_assert(!IsConvertible<double,LowerFloat>::value,"");
static_assert(!IsConvertible<double,LowerFloat>::value,"");

static_assert(And<IsConstructible<LowerFloat,double>,IsConstructible<UpperFloat,double>,Not<And<IsConvertible<double,LowerFloat>,IsConvertible<double,UpperFloat>>>>::value,"");



class TestInterval
{
    typedef UpperFloatInterval I;
    typedef Flt R;
  public:
    void test();
  private:
    void test_concept();
    void test_constructors();
    void test_input();
    void test_class();
    void test_comparison();
    void test_geometric_predicates();
};


void
TestInterval::test()
{
    ARIADNE_TEST_CALL(test_concept());
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_input());
    ARIADNE_TEST_CALL(test_class());
    ARIADNE_TEST_CALL(test_comparison());
    ARIADNE_TEST_CALL(test_geometric_predicates());
}

void
TestInterval::test_concept()
{
    int n=1;
    uint m=1;
    ExactFloat x=1;
    ExactFloat a,b;
    UpperFloatInterval i=UpperFloatInterval(1);
    UpperFloatInterval j=UpperFloatInterval(1);


    // Constructors
    j=UpperFloatInterval(); j=UpperFloatInterval(n); j=UpperFloatInterval(m);
    j=UpperFloatInterval(x); j=UpperFloatInterval(i);
    j=UpperFloatInterval(n,n); j=UpperFloatInterval(m,m); j=UpperFloatInterval(x,x); j=UpperFloatInterval(a,b);
    j=UpperFloatInterval(n,m); j=UpperFloatInterval(m,x); j=UpperFloatInterval(x,n);
    // Assignment
    j=n; j=m; j=x; j=i;

}


void
TestInterval::test_constructors()
{

    Flt zero=0;

    // Construct from pair
    UpperFloatInterval ivld1(Flt(1.125),Flt(2.25));
    ARIADNE_TEST_ASSERT(ivld1.lower()==1.125); ARIADNE_TEST_ASSERT(ivld1.upper()==2.25);

    // Default constructor
    UpperFloatInterval ivld2;
    if(ivld2.lower()>ivld2.upper()) {
        ARIADNE_TEST_NOTIFY("UpperFloatInterval default constructor returns an empty set.");
    } else {
        ARIADNE_TEST_ASSERT(decide(ivld2==UpperFloatInterval(zero,zero)));
    }

    // Constructor from approximate values
    UpperFloatInterval ivld3(2.1,3.2);
    ARIADNE_TEST_COMPARE(ivld3.lower(),<=,2.1);
    ARIADNE_TEST_COMPARE(ivld3.upper(),>=,3.2);

/*
    // Constructor with approximations
#ifdef HAVE_GMPXX_H
    UpperFloatInterval ivld4(Rational(21,10),Rational(16,5));
    cout<<ivld3<<std::endl;
    ARIADNE_TEST_COMPARE(Rational(make_exact(ivld4.lower())),<,Rational(21,10));
    ARIADNE_TEST_COMPARE(Rational(make_exact(ivld4.upper())),>,Rational(16,5));
    // Approximate constructor from a single value
    UpperFloatInterval ivld5(Rational(1,3));
    ARIADNE_TEST_COMPARE(Rational(make_exact(ivld5.lower())),<,Rational(1,3));
    ARIADNE_TEST_COMPARE(Rational(make_exact(ivld5.upper())),>,Rational(1,3));
#else
    UpperFloatInterval ivld4(2.1,3.2);
    UpperFloatInterval ivld5(1./3.);
#endif // HAVE_GMPXX_H
*/
    // Exact constructor from a single value
    UpperFloatInterval ivld6(ExactFloat(1.25));
    ARIADNE_TEST_EQUAL(ivld6.lower(),ExactFloat(1.25));
    ARIADNE_TEST_EQUAL(ivld6.upper(),ExactFloat(1.25));

    // Empty interval
    UpperFloatInterval ivld7;
    ARIADNE_TEST_EXECUTE(ivld7=EmptyInterval());
    ARIADNE_TEST_ASSERT(ivld7.lower()==+inf); ARIADNE_TEST_ASSERT(ivld7.upper()==-inf);
}

void TestInterval::test_class()
{
    // Test lower, upper, midpoint, radius, width

    // Tests for exact operations
    ARIADNE_TEST_EQUAL(UpperFloatInterval(-0.25,0.50).lower(),-0.25_x);
    ARIADNE_TEST_EQUAL(UpperFloatInterval(-0.25,0.50).upper(),0.5_x);
    ARIADNE_TEST_EQUAL(UpperFloatInterval(-0.25,0.50).midpoint(),0.125_x);
    ARIADNE_TEST_EQUAL(UpperFloatInterval(-0.25,0.50).radius(),0.375_x)
    ARIADNE_TEST_EQUAL(UpperFloatInterval(-0.25,0.50).width(),0.75_x);

    // Tests for inexact operations
    ARIADNE_TEST_EQUAL(UpperFloatInterval(-1./3,2./3).lower(),-0.33333333333333331483_x);
    ARIADNE_TEST_EQUAL(UpperFloatInterval(-1./3,2./3).upper(),0.66666666666666662966_x);
    ARIADNE_TEST_EQUAL(UpperFloatInterval(-1./3,2./3).midpoint(),0.16666666666666665741_x);
    ARIADNE_TEST_EQUAL(UpperFloatInterval(-1./3,2./3).radius(),0.5)
    ARIADNE_TEST_EQUAL(UpperFloatInterval(-1./3,2./3).width(),1.0);

    // Tests for inexact operations
    ARIADNE_TEST_EQUAL(UpperFloatInterval(div_down(-1,3),div_up(2,3)).lower(),-0.33333333333333337034_x);
    ARIADNE_TEST_EQUAL(UpperFloatInterval(div_down(-1,3),div_up(2,3)).upper(),0.66666666666666674068_x);
    ARIADNE_TEST_EQUAL(UpperFloatInterval(div_down(-1,3),div_up(2,3)).midpoint(),0.16666666666666668517_x);
    ARIADNE_TEST_EQUAL(UpperFloatInterval(div_down(-1,3),div_up(2,3)).radius(),0.50000000000000011102_x)
    ARIADNE_TEST_EQUAL(UpperFloatInterval(div_down(-1,3),div_up(2,3)).width(),1.000000000000000222_x);
}

void TestInterval::test_input()
{
/*
    UpperFloatInterval ivl;
    string input("[-1.125,2.25] [0.4,0.6]");
    stringstream iss(input);

    iss >> ivl;
    ARIADNE_TEST_EQUALS(ivl.lower().get_flt(),-1.125);
    ARIADNE_TEST_EQUALS(ivl.upper().get_flt(),2.25);

    iss >> ivl;
    Flt lower=0.39999999999999997;
    Flt upper=0.60000000000000009;
    if(ivl.lower().get_flt()!=lower || ivl.upper().get_flt()!=upper) {
        ARIADNE_TEST_WARN("UpperFloatInterval string constructor returns an approximate interval, not an outwardly rounded interval.");
    }
*/
}

void TestInterval::test_comparison() {
    // FIXME: If using Boost style interval tests, uncomment the line below
    // and comment out the line after
    //ARIADNE_TEST_ASSERT(indeterminate(ivld1==ivld2));
    ExactFloatInterval ivl1(1.125,2.25);
    ExactFloatInterval ivl2=ivl1;

    ARIADNE_TEST_ASSERT(ivl1==ivl2);
    ExactFloatInterval& ivl1ref=ivl1;
    ivl1ref=ExactFloatInterval(5.25,7.375);
    cout << "ivl1ref=" << ivl1ref << endl;
    ARIADNE_TEST_ASSERT(ivl1ref.lower()==ExactFloat(5.25));
}


void TestInterval::test_geometric_predicates()
{
    ExactFloatInterval empty_interval=EmptyInterval();

    ARIADNE_TEST_PRINT(empty_interval);

    ARIADNE_TEST_BINARY_PREDICATE(intersect,ExactFloatInterval(0.0,1.5),ExactFloatInterval(1.5,3));
    ARIADNE_TEST_BINARY_PREDICATE(!intersect,ExactFloatInterval(0.0,1.5),ExactFloatInterval(1.625,3));

    ARIADNE_TEST_BINARY_PREDICATE(!disjoint,ExactFloatInterval(0.0,1.5),ExactFloatInterval(1.5,3));
    ARIADNE_TEST_BINARY_PREDICATE(disjoint,ExactFloatInterval(0.0,1.5),ExactFloatInterval(1.625,3));

    ARIADNE_TEST_BINARY_PREDICATE(subset,ExactFloatInterval(1.0,1.5),ExactFloatInterval(1.0,2.0));
    ARIADNE_TEST_BINARY_PREDICATE(subset,ExactFloatInterval(1.5,2.0),ExactFloatInterval(1.0,2.0));
    ARIADNE_TEST_BINARY_PREDICATE(subset,ExactFloatInterval(1.0,2.0),ExactFloatInterval(1.0,2.0));

    ARIADNE_TEST_BINARY_PREDICATE(superset,ExactFloatInterval(1.0,2.0),ExactFloatInterval(1.0,1.5));
    ARIADNE_TEST_BINARY_PREDICATE(superset,ExactFloatInterval(1.0,2.0),ExactFloatInterval(1.5,2.0));
    ARIADNE_TEST_BINARY_PREDICATE(superset,ExactFloatInterval(1.0,2.0),ExactFloatInterval(1.0,2.0));

    ARIADNE_TEST_BINARY_PREDICATE(inside,ExactFloatInterval(1.25,1.75),ExactFloatInterval(1.0,2.0));
    ARIADNE_TEST_BINARY_PREDICATE(!inside,ExactFloatInterval(1.00,1.75),ExactFloatInterval(1.0,2.0));
    ARIADNE_TEST_BINARY_PREDICATE(!inside,ExactFloatInterval(1.25,2.00),ExactFloatInterval(1.0,2.0));

    ARIADNE_TEST_BINARY_PREDICATE(covers,ExactFloatInterval(1.0,2.0),ExactFloatInterval(1.25,1.75));
    ARIADNE_TEST_BINARY_PREDICATE(!covers,ExactFloatInterval(1.0,2.0),ExactFloatInterval(1.00,1.75));
    ARIADNE_TEST_BINARY_PREDICATE(!covers,ExactFloatInterval(1.0,2.0),ExactFloatInterval(1.25,2.00));

    ARIADNE_TEST_BINARY_PREDICATE(overlap,ExactFloatInterval(1.0,2.0),ExactFloatInterval(1.75,3.0));
    ARIADNE_TEST_BINARY_PREDICATE(!overlap,ExactFloatInterval(1.0,2.0),ExactFloatInterval(2.00,2.75));
    ARIADNE_TEST_BINARY_PREDICATE(!overlap,ExactFloatInterval(1.0,2.0),ExactFloatInterval(2.25,2.75));

    ARIADNE_TEST_BINARY_PREDICATE(!separated,ExactFloatInterval(1.0,2.0),ExactFloatInterval(1.75,3.0));
    ARIADNE_TEST_BINARY_PREDICATE(!separated,ExactFloatInterval(1.0,2.0),ExactFloatInterval(2.00,2.75));
    ARIADNE_TEST_BINARY_PREDICATE(separated,ExactFloatInterval(1.0,2.0),ExactFloatInterval(2.25,2.75));

    ARIADNE_TEST_BINARY_PREDICATE(disjoint,empty_interval,empty_interval);
    ARIADNE_TEST_BINARY_PREDICATE(!intersect,empty_interval,empty_interval);
    ARIADNE_TEST_BINARY_PREDICATE(superset,empty_interval,empty_interval);
    ARIADNE_TEST_BINARY_PREDICATE(subset,empty_interval,empty_interval);
    ARIADNE_TEST_BINARY_PREDICATE(separated,empty_interval,empty_interval);
    ARIADNE_TEST_BINARY_PREDICATE(!overlap,empty_interval,empty_interval);

    ARIADNE_TEST_NOTIFY("inside(empty_interval,empty_interval) returns "<<inside(empty_interval,empty_interval));
    ARIADNE_TEST_NOTIFY("covers(empty_interval,empty_interval) returns "<<covers(empty_interval,empty_interval));

    ARIADNE_TEST_BINARY_PREDICATE(disjoint,ExactFloatInterval(-1,2),empty_interval);
    ARIADNE_TEST_BINARY_PREDICATE(subset,empty_interval,ExactFloatInterval(0.25,0.75));
    ARIADNE_TEST_BINARY_PREDICATE(covers,ExactFloatInterval(0.25,0.75),empty_interval);
    ARIADNE_TEST_BINARY_PREDICATE(inside,empty_interval,ExactFloatInterval(0.25,0.75));
}




int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    TestInterval().test();

    return ARIADNE_TEST_FAILURES;
}

