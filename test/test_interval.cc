/***************************************************************************
 *            test_interval.cc
 *
 *  Copyright  2006-8  Alberto Casagrande, Pieter Collins
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

#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include "config.h"
#include "numeric/rational.h"
#include "geometry/interval.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

class TestInterval
{
    typedef ExactInterval I;
    typedef Float64 R;
  public:
    Void test();
  private:
    Void test_concept();
    Void test_constructors();
    Void test_input();
    Void test_class();
    Void test_comparison();
    Void test_geometric_predicates();
    Void regression_tests();
};


Void
TestInterval::test()
{
    ARIADNE_TEST_CALL(test_concept());
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_input());
    ARIADNE_TEST_CALL(test_class());
    ARIADNE_TEST_CALL(test_comparison());
    ARIADNE_TEST_CALL(test_geometric_predicates());
    ARIADNE_TEST_CALL(regression_tests());
}

Void
TestInterval::test_concept()
{
    Int n=1;
    Nat m=1;
    double d=1;
    ExactFloat64 x=1;
    Float64 a,b;
    ExactInterval i;
    UpperInterval j;


    // Constructors
    j=I(); j=I(n); j=I(m); j=I(d); j=I(x); j=I(i);
    j=I(n,n); j=I(m,m); j=I(d,d); j=I(a,b);
    j=I(n,m); j=I(m,d); j=I(d,n);
    // Assignment
    j=n; j=m; j=x; j=i;

}


Void
TestInterval::test_constructors()
{

    Float64 zero=0;

    // Construct from pair
    ExactInterval ivld1(Float64(1.125),Float64(2.25));
    ARIADNE_TEST_ASSERT(ivld1.lower().raw()==1.125); ARIADNE_TEST_ASSERT(ivld1.upper().raw()==2.25);

    // Default constructor
    ExactInterval ivld2;
    if(ivld2.lower()>ivld2.upper()) {
        ARIADNE_TEST_WARN("ExactInterval default constructor returns an empty set.");
    } else {
        ARIADNE_TEST_ASSERT((Bool)(ivld2==ExactInterval(zero,zero)));
    }

    // Constructor without approximations
    ExactInterval ivld3(Rational(21,8),Rational(17,4));
    cout<<ivld3<<std::endl;
    ARIADNE_TEST_COMPARE(make_exact(ivld3.lower()),==,Rational(21,8));
    ARIADNE_TEST_COMPARE(make_exact(ivld3.upper()),==,Rational(17,4));

    // Constructor from approximate values
    UpperInterval ivld4(2.1,3.2);
    ARIADNE_TEST_COMPARE(ivld4.lower(),<=,2.1);
    ARIADNE_TEST_COMPARE(ivld4.upper(),>=,3.2);

    // Approximate constructor from a single value
    UpperInterval ivld5(Rational(1,3));
    ARIADNE_TEST_COMPARE(make_exact(ivld5.lower()),<,Rational(1,3));
    ARIADNE_TEST_COMPARE(make_exact(ivld5.upper()),>,Rational(1,3));

    // Exact constructor from a single value
    ExactInterval ivld6(Float64(1.25));
    ARIADNE_TEST_EQUAL(ivld6.lower().raw(),Float64(1.25));
    ARIADNE_TEST_EQUAL(ivld6.upper().raw(),Float64(1.25));

    // Empty interval
    ExactInterval ivld7;
    ARIADNE_TEST_EXECUTE(ivld7.set_empty());
    ARIADNE_TEST_ASSERT(ivld7.lower().raw()==+inf); ARIADNE_TEST_ASSERT(ivld7.upper().raw()==-inf);
}

Void TestInterval::test_class()
{
    // Test lower, upper, midpoint, radius, width

    // Tests for exact operations
    ARIADNE_TEST_EQUAL(ExactInterval(-0.25,0.50).lower().raw(),-0.25);
    ARIADNE_TEST_EQUAL(ExactInterval(-0.25,0.50).upper().raw(),0.5);
    ARIADNE_TEST_EQUAL(ExactInterval(-0.25,0.50).centre().raw(),0.125);
    ARIADNE_TEST_EQUAL(ExactInterval(-0.25,0.50).error().raw(),0.375)
    ARIADNE_TEST_EQUAL(ExactInterval(-0.25,0.50).width().raw(),0.75);

}

Void TestInterval::test_input()
{
    ExactInterval ivl1,ivl2;
    string input("[1.125,2.25] [0.4,0.6]");
    stringstream iss(input);

    iss >> ivl1;
    ivl2=ExactInterval(1.125,2.25);

    cout << "ivl1=" << ivl1 << "  ivl2=" << ivl2 << endl;
    ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,ivl2);
    ARIADNE_TEST_ASSERT(ivl1.lower()==ivl2.lower() && ivl1.upper()==ivl2.upper());
    ARIADNE_TEST_ASSERT(equal(ivl1,ivl2));

    iss >> ivl1;
    ivl2=ExactInterval(0.39999999999999997,0.60000000000000009);
    if(!equal(ivl1,ivl2)) {
        ARIADNE_TEST_WARN("ExactInterval string constructor returns an approximate interval, not an outwardly rounded interval.");
    }
}

Void TestInterval::test_comparison() {
    // FIXME: If using Boost style interval tests, uncomment the line below
    // and comment out the line after
    //ARIADNE_TEST_ASSERT(indeterminate(ivld1==ivld2));
    ExactInterval ivl1(1.125,2.25);
    ExactInterval ivl2=ivl1;

    ARIADNE_TEST_ASSERT(ivl1==ivl2);
    ExactInterval& ivl1ref=ivl1;
    ivl1ref=ExactInterval(5.25,7.375);
    cout << "ivl1ref=" << ivl1ref << endl;
    ARIADNE_TEST_ASSERT(ivl1ref.lower().raw()==Float64(5.25));
}

Void TestInterval::test_geometric_predicates()
{
    ExactInterval empty_interval; empty_interval.set_empty();

    ARIADNE_TEST_PRINT(empty_interval);

    ARIADNE_TEST_BINARY_PREDICATE(intersect,ExactInterval(0.0,1.5),ExactInterval(1.5,3));
    ARIADNE_TEST_BINARY_PREDICATE(!intersect,ExactInterval(0.0,1.5),ExactInterval(1.625,3));

    ARIADNE_TEST_BINARY_PREDICATE(!disjoint,ExactInterval(0.0,1.5),ExactInterval(1.5,3));
    ARIADNE_TEST_BINARY_PREDICATE(disjoint,ExactInterval(0.0,1.5),ExactInterval(1.625,3));

    ARIADNE_TEST_BINARY_PREDICATE(subset,ExactInterval(1.0,1.5),ExactInterval(1.0,2.0));
    ARIADNE_TEST_BINARY_PREDICATE(subset,ExactInterval(1.5,2.0),ExactInterval(1.0,2.0));
    ARIADNE_TEST_BINARY_PREDICATE(subset,ExactInterval(1.0,2.0),ExactInterval(1.0,2.0));

    ARIADNE_TEST_BINARY_PREDICATE(superset,ExactInterval(1.0,2.0),ExactInterval(1.0,1.5));
    ARIADNE_TEST_BINARY_PREDICATE(superset,ExactInterval(1.0,2.0),ExactInterval(1.5,2.0));
    ARIADNE_TEST_BINARY_PREDICATE(superset,ExactInterval(1.0,2.0),ExactInterval(1.0,2.0));

    ARIADNE_TEST_BINARY_PREDICATE(inside,ExactInterval(1.25,1.75),ExactInterval(1.0,2.0));
    ARIADNE_TEST_BINARY_PREDICATE(!inside,ExactInterval(1.00,1.75),ExactInterval(1.0,2.0));
    ARIADNE_TEST_BINARY_PREDICATE(!inside,ExactInterval(1.25,2.00),ExactInterval(1.0,2.0));

    ARIADNE_TEST_BINARY_PREDICATE(covers,ExactInterval(1.0,2.0),ExactInterval(1.25,1.75));
    ARIADNE_TEST_BINARY_PREDICATE(!covers,ExactInterval(1.0,2.0),ExactInterval(1.00,1.75));
    ARIADNE_TEST_BINARY_PREDICATE(!covers,ExactInterval(1.0,2.0),ExactInterval(1.25,2.00));

    ARIADNE_TEST_BINARY_PREDICATE(overlap,ExactInterval(1.0,2.0),ExactInterval(1.75,3.0));
    ARIADNE_TEST_BINARY_PREDICATE(!overlap,ExactInterval(1.0,2.0),ExactInterval(2.00,2.75));
    ARIADNE_TEST_BINARY_PREDICATE(!overlap,ExactInterval(1.0,2.0),ExactInterval(2.25,2.75));

    ARIADNE_TEST_BINARY_PREDICATE(!separated,ExactInterval(1.0,2.0),ExactInterval(1.75,3.0));
    ARIADNE_TEST_BINARY_PREDICATE(!separated,ExactInterval(1.0,2.0),ExactInterval(2.00,2.75));
    ARIADNE_TEST_BINARY_PREDICATE(separated,ExactInterval(1.0,2.0),ExactInterval(2.25,2.75));

    ARIADNE_TEST_BINARY_PREDICATE(disjoint,empty_interval,empty_interval);
    ARIADNE_TEST_BINARY_PREDICATE(!intersect,empty_interval,empty_interval);
    ARIADNE_TEST_BINARY_PREDICATE(superset,empty_interval,empty_interval);
    ARIADNE_TEST_BINARY_PREDICATE(subset,empty_interval,empty_interval);
    ARIADNE_TEST_BINARY_PREDICATE(separated,empty_interval,empty_interval);
    ARIADNE_TEST_BINARY_PREDICATE(!overlap,empty_interval,empty_interval);

    if(!inside(empty_interval,empty_interval)) {
        ARIADNE_TEST_WARN("inside(empty_interval,empty_interval) returns false");
    }

    if(!covers(empty_interval,empty_interval)) {
        ARIADNE_TEST_WARN("covers(empty_interval,empty_interval) returns false");
    }

    ARIADNE_TEST_BINARY_PREDICATE(disjoint,ExactInterval(-1,2),empty_interval);
    ARIADNE_TEST_BINARY_PREDICATE(subset,empty_interval,ExactInterval(0.25,0.75));
    ARIADNE_TEST_BINARY_PREDICATE(covers,ExactInterval(0.25,0.75),empty_interval);
    ARIADNE_TEST_BINARY_PREDICATE(inside,empty_interval,ExactInterval(0.25,0.75));
}

Void TestInterval::regression_tests() {

    // Regression test; fails dramatically on certain types of rounding
    {
        UpperInterval x(1.5707963267948966,1.5707963267948968);
        UpperInterval cosx=cos(x);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_COMPARE(cosx.lower(),<,0.0);
        ARIADNE_TEST_COMPARE(cosx.lower(),>,-1e-14);
        ARIADNE_TEST_COMPARE(cosx.upper(),>,0.0);
        ARIADNE_TEST_COMPARE(cosx.upper(),<,+1e-14);
        ARIADNE_TEST_ASSERT(cosx.lower().raw()<cosx.upper().raw());
    }

}


Int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    TestInterval().test();

    return ARIADNE_TEST_FAILURES;
}

