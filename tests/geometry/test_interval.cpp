/***************************************************************************
 *            test_interval.cpp
 *
 *  Copyright  2006-20  Alberto Casagrande, Pieter Collins
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

#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include "config.hpp"
#include "numeric/rational.hpp"
#include "geometry/interval.hpp"

#include "../test.hpp"

using namespace Ariadne;
using namespace std;

class TestIntervalType
{
    typedef ExactIntervalType I;
    typedef FloatDP R;
  public:
    Void test();
  private:
    Void test_concept();
    Void test_constructors();
    Void test_input();
    Void test_class();
    Void test_comparison();
    Void test_geometric_predicates();
    Void test_arithmetic();
    Void regression_tests();
};


Void
TestIntervalType::test()
{
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_input());
    ARIADNE_TEST_CALL(test_class());
    ARIADNE_TEST_CALL(test_comparison());
    ARIADNE_TEST_CALL(test_geometric_predicates());
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(regression_tests());
}

Void
TestIntervalType::test_concept()
{
    DoublePrecision pr;
    Int n=1;
    Nat m=1;
    double d=1;
    FloatDPValue x={1,pr};
    FloatDP a,b;
    ExactIntervalType xivl;
    UpperIntervalType uivl;

    // Constructors
    uivl=I(); uivl=I(n); uivl=I(m); uivl=I(d); uivl=I(x); uivl=I(xivl);
    uivl=I(n,n); uivl=I(m,m); uivl=I(d,d); uivl=I(a,b);
    uivl=I(n,m); uivl=I(m,d); uivl=I(d,n);
    // Assignment
    uivl=n; uivl=m; uivl=x; uivl=xivl;

}


Void
TestIntervalType::test_constructors()
{

    FloatDP zero=0;

    // Construct from pair
    ExactIntervalType ivld1(FloatDP(1.125),FloatDP(2.25));
    ARIADNE_TEST_ASSERT(ivld1.lower().raw()==1.125); ARIADNE_TEST_ASSERT(ivld1.upper().raw()==2.25);

    // Default constructor
    ExactIntervalType ivld2;
    if(ivld2.lower()>ivld2.upper()) {
        ARIADNE_TEST_WARN("ExactIntervalType default constructor returns an empty set.");
    } else {
        ARIADNE_TEST_ASSERT((Bool)(ivld2==ExactIntervalType(zero,zero)));
    }

    // Constructor without approximations
    RationalInterval ivld3(Rational(21,8),Rational(17,4));
    cout<<ivld3<<std::endl;
    ARIADNE_TEST_COMPARE(ivld3.lower(),==,Rational(21,8));
    ARIADNE_TEST_COMPARE(ivld3.upper(),==,Rational(17,4));

    // Constructor from approximate values
    UpperIntervalType ivld4(2.1,3.2);
    ARIADNE_TEST_COMPARE(ivld4.lower(),<=,2.1);
    ARIADNE_TEST_COMPARE(ivld4.upper(),>=,3.2);

    // ApproximateTag constructor from a single value
    ARIADNE_TEST_WARN("Cannot construct Interval<UpperFloatDP> from Rational.");
//    UpperIntervalType ivld5(Rational(1,3));
//    ARIADNE_TEST_COMPARE(cast_exact(ivld5.lower()),<,Rational(1,3));
//    ARIADNE_TEST_COMPARE(cast_exact(ivld5.upper()),>,Rational(1,3));

    // ExactTag constructor from a single value
    ExactIntervalType ivld6(FloatDP(1.25));
    ARIADNE_TEST_EQUAL(ivld6.lower().raw(),FloatDP(1.25));
    ARIADNE_TEST_EQUAL(ivld6.upper().raw(),FloatDP(1.25));

    // Empty interval
    EmptyInterval empty_interval;
    ExactIntervalType ivld7(empty_interval);
    ARIADNE_TEST_ASSERT(ivld7.lower().raw()==+inf); ARIADNE_TEST_ASSERT(ivld7.upper().raw()==-inf);
}

Void TestIntervalType::test_class()
{
    // Test lower, upper, midpoint, radius, width

    // Tests for exact operations
    ARIADNE_TEST_EQUAL(ExactIntervalType(-0.25,0.50).lower(),-0.25_exact);
    ARIADNE_TEST_EQUAL(ExactIntervalType(-0.25,0.50).upper(),0.5_exact);
    ARIADNE_TEST_EQUAL(ExactIntervalType(-0.25,0.50).midpoint(),0.125_exact);
    ARIADNE_TEST_EQUAL(ExactIntervalType(-0.25,0.50).centre(),0.125_exact);
    ARIADNE_TEST_EQUAL(ExactIntervalType(-0.25,0.50).radius(),0.375_exact)
    ARIADNE_TEST_EQUAL(ExactIntervalType(-0.25,0.50).width(),0.75_exact);

}

Void TestIntervalType::test_input()
{
    ExactIntervalType ivl1,ivl2;
    string input("[1.125,2.25] [0.4,0.6]");
    stringstream iss(input);

    iss >> ivl1;
    ivl2=ExactIntervalType(1.125,2.25);

    cout << "ivl1=" << ivl1 << "  ivl2=" << ivl2 << endl;
    ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,ivl2);
    ARIADNE_TEST_ASSERT(ivl1.lower()==ivl2.lower() && ivl1.upper()==ivl2.upper());
    ARIADNE_TEST_ASSERT(equal(ivl1,ivl2));

    iss >> ivl1;
    ivl2=ExactIntervalType(0.39999999999999997,0.60000000000000009);
    if(!equal(ivl1,ivl2)) {
        ARIADNE_TEST_WARN("ExactIntervalType string constructor returns an approximate interval, not an outwardly rounded interval.");
    }
}

Void TestIntervalType::test_comparison() {
    // FIXME: If using Boost style interval tests, uncomment the line below
    // and comment out the line after
    //ARIADNE_TEST_ASSERT(indeterminate(ivld1==ivld2));
    ExactIntervalType ivl1(1.125,2.25);
    ExactIntervalType ivl2=ivl1;

    ARIADNE_TEST_ASSERT(ivl1==ivl2);
    ExactIntervalType& ivl1ref=ivl1;
    ivl1ref=ExactIntervalType(5.25,7.375);
    cout << "ivl1ref=" << ivl1ref << endl;
    ARIADNE_TEST_ASSERT(ivl1ref.lower().raw()==FloatDP(5.25));
}

Void TestIntervalType::test_geometric_predicates()
{
    ExactIntervalType empty_interval=EmptyInterval();

    ARIADNE_TEST_PRINT(empty_interval);

    ARIADNE_TEST_BINARY_PREDICATE(intersect,ExactIntervalType(0.0,1.5),ExactIntervalType(1.5,3));
    ARIADNE_TEST_BINARY_PREDICATE(!intersect,ExactIntervalType(0.0,1.5),ExactIntervalType(1.625,3));

    ARIADNE_TEST_BINARY_PREDICATE(!disjoint,ExactIntervalType(0.0,1.5),ExactIntervalType(1.5,3));
    ARIADNE_TEST_BINARY_PREDICATE(disjoint,ExactIntervalType(0.0,1.5),ExactIntervalType(1.625,3));

    ARIADNE_TEST_BINARY_PREDICATE(subset,ExactIntervalType(1.0,1.5),ExactIntervalType(1.0,2.0));
    ARIADNE_TEST_BINARY_PREDICATE(subset,ExactIntervalType(1.5,2.0),ExactIntervalType(1.0,2.0));
    ARIADNE_TEST_BINARY_PREDICATE(subset,ExactIntervalType(1.0,2.0),ExactIntervalType(1.0,2.0));

    ARIADNE_TEST_BINARY_PREDICATE(superset,ExactIntervalType(1.0,2.0),ExactIntervalType(1.0,1.5));
    ARIADNE_TEST_BINARY_PREDICATE(superset,ExactIntervalType(1.0,2.0),ExactIntervalType(1.5,2.0));
    ARIADNE_TEST_BINARY_PREDICATE(superset,ExactIntervalType(1.0,2.0),ExactIntervalType(1.0,2.0));

    ARIADNE_TEST_BINARY_PREDICATE(inside,ExactIntervalType(1.25,1.75),ExactIntervalType(1.0,2.0));
    ARIADNE_TEST_BINARY_PREDICATE(!inside,ExactIntervalType(1.00,1.75),ExactIntervalType(1.0,2.0));
    ARIADNE_TEST_BINARY_PREDICATE(!inside,ExactIntervalType(1.25,2.00),ExactIntervalType(1.0,2.0));

    ARIADNE_TEST_BINARY_PREDICATE(covers,ExactIntervalType(1.0,2.0),ExactIntervalType(1.25,1.75));
    ARIADNE_TEST_BINARY_PREDICATE(!covers,ExactIntervalType(1.0,2.0),ExactIntervalType(1.00,1.75));
    ARIADNE_TEST_BINARY_PREDICATE(!covers,ExactIntervalType(1.0,2.0),ExactIntervalType(1.25,2.00));

    ARIADNE_TEST_BINARY_PREDICATE(overlap,ExactIntervalType(1.0,2.0),ExactIntervalType(1.75,3.0));
    ARIADNE_TEST_BINARY_PREDICATE(!overlap,ExactIntervalType(1.0,2.0),ExactIntervalType(2.00,2.75));
    ARIADNE_TEST_BINARY_PREDICATE(!overlap,ExactIntervalType(1.0,2.0),ExactIntervalType(2.25,2.75));

    ARIADNE_TEST_BINARY_PREDICATE(!separated,ExactIntervalType(1.0,2.0),ExactIntervalType(1.75,3.0));
    ARIADNE_TEST_BINARY_PREDICATE(!separated,ExactIntervalType(1.0,2.0),ExactIntervalType(2.00,2.75));
    ARIADNE_TEST_BINARY_PREDICATE(separated,ExactIntervalType(1.0,2.0),ExactIntervalType(2.25,2.75));

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

    ARIADNE_TEST_BINARY_PREDICATE(disjoint,ExactIntervalType(-1,2),empty_interval);
    ARIADNE_TEST_BINARY_PREDICATE(subset,empty_interval,ExactIntervalType(0.25,0.75));
    ARIADNE_TEST_BINARY_PREDICATE(covers,ExactIntervalType(0.25,0.75),empty_interval);
    ARIADNE_TEST_BINARY_PREDICATE(inside,empty_interval,ExactIntervalType(0.25,0.75));
}

Void TestIntervalType::test_arithmetic() {
    UpperIntervalType e=EmptyInterval();
    ARIADNE_TEST_SAME(sqrt(UpperIntervalType(-4,-1)),e);
    ARIADNE_TEST_SAME(sqrt(UpperIntervalType(-4,0)),UpperIntervalType(0,0));
    ARIADNE_TEST_BINARY_PREDICATE(refines,UpperIntervalType(-3,3),sqrt(UpperIntervalType(-4,9)));
    ARIADNE_TEST_SAME(log(UpperIntervalType(-4,-1)),e);
    ARIADNE_TEST_SAME(log(UpperIntervalType(-4,0)),e);
    ARIADNE_TEST_SAME(log(UpperIntervalType(-4,1)),UpperIntervalType(-inf,0));
    ARIADNE_TEST_SAME(log(UpperIntervalType( 0,1)),UpperIntervalType(-inf,0));
}

Void TestIntervalType::regression_tests() {

    // Regression test; fails dramatically on certain types of rounding
    {
        UpperIntervalType x(1.5707963267948966,1.5707963267948968);
        UpperIntervalType cosx=cos(x);
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

    TestIntervalType().test();

    return ARIADNE_TEST_FAILURES;
}

