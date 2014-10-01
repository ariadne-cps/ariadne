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
#include "rational.h"
#include "interval.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

namespace Ariadne {
inline ExactFloatType const& make_exact(ExactFloatType const& x) { return x; }
}

class TestInterval
{
    typedef Interval I;
    typedef Float R;
  public:
    void test();
  private:
    void test_concept();
    void test_constructors();
    void test_input();
    void test_class();
    void test_comparison();
    void test_correct_rounded_arithmetic();
    void test_accurate_rounded_arithmetic();
    void test_exact_rounded_arithmetic();
    void test_aliasing();
    void test_monotone_functions();
    void test_trigonometric_functions();
    void test_geometric_predicates();
    void regression_tests();
};


void
TestInterval::test()
{
    ARIADNE_TEST_CALL(test_concept());
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_input());
    ARIADNE_TEST_CALL(test_class());
    ARIADNE_TEST_CALL(test_comparison());
    ARIADNE_TEST_CALL(test_correct_rounded_arithmetic());
    ARIADNE_TEST_CALL(test_accurate_rounded_arithmetic());
    ARIADNE_TEST_CALL(test_exact_rounded_arithmetic());
    ARIADNE_TEST_CALL(test_monotone_functions());
    ARIADNE_TEST_CALL(test_trigonometric_functions());
    ARIADNE_TEST_CALL(test_geometric_predicates());
    ARIADNE_TEST_CALL(regression_tests());
}

void
TestInterval::test_concept()
{
    int n=1;
    uint m=1;
    double d=1;
    ExactFloat x=1;
    Float a,b;
    Interval i=1;
    Interval j=1;


    // Constructors
    j=I(); j=I(n); j=I(m); j=I(d); j=I(x); j=I(i);
    j=I(n,n); j=I(m,m); j=I(d,d); j=I(a,b);
    j=I(n,m); j=I(m,d); j=I(d,n);
    // Assignment
    j=n; j=m; j=d; j=x; j=i;

    // Exact operations
    j=abs(i); j=neg(i); j=rec(i);

    // Arithmetic
    //j=add(x,x); j=add(x,i); j=add(i,x); j=add(i,i);
    //j=sub(x,x); j=sub(x,i); j=sub(i,x); j=sub(i,i);
    //j=mul(x,x); j=mul(x,i); j=mul(i,x); j=mul(i,i);
    //j=div(x,x); j=div(x,i); j=div(i,x); j=div(i,i);
    //j=pow(x,n); j=pow(x,n);
    //j=pow(x,m); j=pow(x,m);

    // Transcendental functions
    j=sqrt(i);
    j=exp(i);
    j=log(i);
    /*
      j=sin(i);
      j=cos(i);
      j=tan(i);
      j=asin(i);
      j=acos(i);
      j=atan(i);

      j=sinh(i);
      j=cosh(i);
      j=tanh(i);
      j=asinh(i);
      j=acosh(i);
      j=atanh(i);
    */
}



// Test that interval arithmetic is rounded correctly,
// without paying attention to accuracy issues.
void
TestInterval::test_correct_rounded_arithmetic()
{
    Interval onethird=Interval(1)/Interval(3);
    ARIADNE_TEST_COMPARE( onethird.lower(), < , onethird.upper() );
    Interval one_approx=onethird*Interval(3);
    ARIADNE_TEST_COMPARE( one_approx.lower(), < , 1.0 );
    ARIADNE_TEST_COMPARE( one_approx.upper(), > , 1.0 );
}


// Test that interval arithmetic gives the most accurate rounded values
void
TestInterval::test_accurate_rounded_arithmetic()
{
    const double min=std::numeric_limits<double>::min();
    const double eps=std::numeric_limits<double>::epsilon();

    ARIADNE_TEST_EQUAL(Interval(1.5)+Interval(min),Interval(1.5,1.5+eps));
    ARIADNE_TEST_EQUAL(Interval(1.5)-Interval(min),Interval(1.5-eps,1.5));
    ARIADNE_TEST_EQUAL(Interval(1+eps,1+2*eps)*Interval(1+eps,1+3*eps),Interval(1+2*eps,1+6*eps));
    ARIADNE_TEST_EQUAL(Interval(1)/Interval(3),Interval(0.33333333333333331,0.33333333333333337));
    ARIADNE_TEST_EQUAL(Interval(2)/Interval(5),Interval(0.39999999999999997,0.40000000000000002));

    ARIADNE_TEST_EQUAL(Interval(1.5)+Float(min),Interval(1.5,1.5+eps));
    ARIADNE_TEST_EQUAL(Interval(1.5)-Float(min),Interval(1.5-eps,1.5));
    ARIADNE_TEST_EQUAL(Interval(1+eps,1+2*eps)*Float(1+eps),Interval(1+2*eps,1+4*eps));
    ARIADNE_TEST_EQUAL(Interval(1+3*eps,1+5*eps)/Float(1+eps),Interval(1+eps,1+4*eps));

    ARIADNE_TEST_EQUAL(Float(min)-Interval(1.5),Interval(-1.5,eps-1.5));
    ARIADNE_TEST_EQUAL(Float(1+5*eps)/Interval(1+2*eps,1+3*eps),Interval(1+eps,1+3*eps));

    ARIADNE_TEST_EQUAL(sqr(Interval(1-eps,1+eps)),Interval(1-4*eps/2,1+3*eps));

    ARIADNE_TEST_EQUAL(pow(Interval(3,5),-1),Interval(0.19999999999999998,0.33333333333333337));
    ARIADNE_TEST_EQUAL(pow(Interval(3,5),-2),Interval(0.039999999999999986955,0.11111111111111114658));

    ARIADNE_TEST_EQUAL(rec(Interval(1+2*eps,1+5*eps)),Interval(1-10*eps/2,1-3*eps/2));

}


// Test that interval arithmetic gives exact values if possible
void
TestInterval::test_exact_rounded_arithmetic()
{
    ARIADNE_TEST_EQUAL(Interval(5,7)+Interval(2,4),Interval(7,11));
    ARIADNE_TEST_EQUAL(Interval(5,7)-Interval(2,6),Interval(-1,5));

    ARIADNE_TEST_EQUAL(Interval(5,7)*Interval(2,4),Interval(10,28));
    ARIADNE_TEST_EQUAL(Interval(5,7)*Interval(-2,4),Interval(-14,28));
    ARIADNE_TEST_EQUAL(Interval(5,7)*Interval(-4,-2),Interval(-28,-10));
    ARIADNE_TEST_EQUAL(Interval(-7,5)*Interval(2,4),Interval(-28,20));
    ARIADNE_TEST_EQUAL(Interval(-7,5)*Interval(-2,4),Interval(-28,20));
    ARIADNE_TEST_EQUAL(Interval(-7,5)*Interval(-4,-2),Interval(-20,28));
    ARIADNE_TEST_EQUAL(Interval(-7,-5)*Interval(2,4),Interval(-28,-10));
    ARIADNE_TEST_EQUAL(Interval(-7,-5)*Interval(-2,4),Interval(-28,14));
    ARIADNE_TEST_EQUAL(Interval(-7,-5)*Interval(-4,-2),Interval(10,28));

    ARIADNE_TEST_EQUAL(Interval(5,7)/Interval(2,4),Interval(1.25,3.50));
    ARIADNE_TEST_EQUAL(Interval(5,7)/Interval(-4,-2),Interval(-3.50,-1.25));
    ARIADNE_TEST_EQUAL(Interval(-7,5)/Interval(2,4),Interval(-3.50,2.50));
    ARIADNE_TEST_EQUAL(Interval(-7,5)/Interval(-4,-2),Interval(-2.50,3.5));
    ARIADNE_TEST_EQUAL(Interval(-7,-5)/Interval(2,4),Interval(-3.50,-1.25));
    ARIADNE_TEST_EQUAL(Interval(-7,-5)/Interval(-4,-2),Interval(1.25,3.50));

    ARIADNE_TEST_EQUAL(pow(Interval(5,7),0u),Interval(1,1));
    ARIADNE_TEST_EQUAL(pow(Interval(-5,7),0u),Interval(1,1));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,5),0u),Interval(1,1));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,-5),0u),Interval(1,1));
    ARIADNE_TEST_EQUAL(pow(Interval(5,7),1u),Interval(5,7));
    ARIADNE_TEST_EQUAL(pow(Interval(-5,7),1u),Interval(-5,7));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,5),1u),Interval(-7,5));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,-5),1u),Interval(-7,-5));
    ARIADNE_TEST_EQUAL(pow(Interval(5,7),2u),Interval(25,49));
    ARIADNE_TEST_EQUAL(pow(Interval(-5,7),2u),Interval(0,49));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,5),2u),Interval(0,49));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,-5),2u),Interval(25,49));
    ARIADNE_TEST_EQUAL(pow(Interval(5,7),3u),Interval(125,343));
    ARIADNE_TEST_EQUAL(pow(Interval(-5,7),3u),Interval(-125,343));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,5),3u),Interval(-343,125));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,-5),3u),Interval(-343,-125));
    ARIADNE_TEST_EQUAL(pow(Interval(5,7),4u),Interval(625,2401));
    ARIADNE_TEST_EQUAL(pow(Interval(-5,7),4u),Interval(0,2401));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,5),4u),Interval(0,2401));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,-5),4u),Interval(625,2401));

    ARIADNE_TEST_EQUAL(pow(Interval(5,7),0),Interval(1,1));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,5),0),Interval(1,1));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,-5),0),Interval(1,1));
    ARIADNE_TEST_EQUAL(pow(Interval(5,7),1),Interval(5,7));
    ARIADNE_TEST_EQUAL(pow(Interval(-5,7),1),Interval(-5,7));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,5),1),Interval(-7,5));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,-5),1),Interval(-7,-5));
    ARIADNE_TEST_EQUAL(pow(Interval(5,7),2),Interval(25,49));
    ARIADNE_TEST_EQUAL(pow(Interval(-5,7),2),Interval(0,49));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,5),2),Interval(0,49));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,-5),2),Interval(25,49));
    ARIADNE_TEST_EQUAL(pow(Interval(5,7),3),Interval(125,343));
    ARIADNE_TEST_EQUAL(pow(Interval(-5,7),3),Interval(-125,343));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,5),3),Interval(-343,125));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,-5),3),Interval(-343,-125));
    ARIADNE_TEST_EQUAL(pow(Interval(5,7),4),Interval(625,2401));
    ARIADNE_TEST_EQUAL(pow(Interval(-5,7),4),Interval(0,2401));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,5),4),Interval(0,2401));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,-5),4),Interval(625,2401));
    ARIADNE_TEST_EQUAL(pow(Interval(5,7),5),Interval(3125,16807));
    ARIADNE_TEST_EQUAL(pow(Interval(-5,7),5),Interval(-3125,16807));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,5),5),Interval(-16807,3125));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,-5),5),Interval(-16807,-3125));
    ARIADNE_TEST_EQUAL(pow(Interval(5,7),7),Interval(78125,823543));
    ARIADNE_TEST_EQUAL(pow(Interval(-5,7),7),Interval(-78125,823543));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,5),7),Interval(-823543,78125));
    ARIADNE_TEST_EQUAL(pow(Interval(-7,-5),7),Interval(-823543,-78125));

    ARIADNE_TEST_EQUAL(pow(Interval(2,4),-1),Interval(0.25,0.5));
    ARIADNE_TEST_EQUAL(pow(Interval(-4,-2),-1),Interval(-0.5,-0.25));
    ARIADNE_TEST_EQUAL(pow(Interval(2,4),-2),Interval(0.0625,0.25));
    ARIADNE_TEST_EQUAL(pow(Interval(-4,-2),-2),Interval(0.0625,0.25));
    ARIADNE_TEST_EQUAL(pow(Interval(2,4),-3),Interval(0.015625,0.125));
    ARIADNE_TEST_EQUAL(pow(Interval(-4,-2),-3),Interval(-0.125,-0.015625));

    ARIADNE_TEST_EQUAL(rec(Interval(2,4)),Interval(0.25,0.50));
    ARIADNE_TEST_EQUAL(rec(Interval(-4,-2)),Interval(-0.50,-0.25));
}



void
TestInterval::test_constructors()
{
    typedef Interval I;


    Float zero=0;

    // Construct from pair
    Interval ivld1(Float(1.125),Float(2.25));
    ARIADNE_TEST_ASSERT(ivld1.lower()==1.125); ARIADNE_TEST_ASSERT(ivld1.upper()==2.25);

    // Default constructor
    Interval ivld2;
    if(ivld2.lower()>ivld2.upper()) {
        ARIADNE_TEST_WARN("Interval default constructor returns an empty set.");
    } else {
        ARIADNE_TEST_ASSERT((bool)(ivld2==Interval(zero,zero)));
    }

    // Constructor with approximations
#ifdef HAVE_GMPXX_H
    Interval ivld3(Rational(21,10),Rational(16,5));
    cout<<ivld3<<std::endl;
    ARIADNE_TEST_COMPARE(make_exact(ivld3.lower()),<,Rational(21,10));
    ARIADNE_TEST_COMPARE(make_exact(ivld3.upper()),>,Rational(16,5));
#else
    Interval ivld3(2.1,3.2);
#endif // HAVE_GMPXX_H

    // Constructor from approximate values
    Interval ivld4(2.1,3.2);
    ARIADNE_TEST_COMPARE(ivld4.lower(),<=,2.1);
    ARIADNE_TEST_COMPARE(ivld4.upper(),>=,3.2);

#ifdef HAVE_GMPXX_H
    // Approximate constructor from a single value
    Interval ivld5(Rational(1,3));
    ARIADNE_TEST_COMPARE(make_exact(ivld5.lower()),<,Rational(1,3));
    ARIADNE_TEST_COMPARE(make_exact(ivld5.upper()),>,Rational(1,3));
#else
    Interval ivld5(1./3.);
#endif // HAVE_GMPXX_H

    // Exact constructor from a single value
    Interval ivld6(Float(1.25));
    ARIADNE_TEST_EQUAL(ivld6.lower().raw(),Float(1.25));
    ARIADNE_TEST_EQUAL(ivld6.upper().raw(),Float(1.25));

    // Empty interval
    Interval ivld7;
    ARIADNE_TEST_EXECUTE(ivld7.set_empty());
    ARIADNE_TEST_ASSERT(ivld7.lower()==+inf); ARIADNE_TEST_ASSERT(ivld7.upper()==-inf);
}

void TestInterval::test_class()
{
    // Test lower, upper, midpoint, radius, width

    // Tests for exact operations
    ARIADNE_TEST_EQUAL(Interval(-0.25,0.50).lower().raw(),-0.25);
    ARIADNE_TEST_EQUAL(Interval(-0.25,0.50).upper().raw(),0.5);
    ARIADNE_TEST_EQUAL(Interval(-0.25,0.50).centre().raw(),0.125);
    ARIADNE_TEST_EQUAL(Interval(-0.25,0.50).error().raw(),0.375)
    ARIADNE_TEST_EQUAL(Interval(-0.25,0.50).width().raw(),0.75);

    // Tests for inexact operations
    ARIADNE_TEST_EQUAL(Interval(-1./3,2./3).lower().raw(),-0.33333333333333331483);
    ARIADNE_TEST_EQUAL(Interval(-1./3,2./3).upper().raw(),0.66666666666666662966);
    ARIADNE_TEST_EQUAL(Interval(-1./3,2./3).centre().raw(),0.16666666666666665741);
    ARIADNE_TEST_EQUAL(Interval(-1./3,2./3).error().raw(),0.5)
    ARIADNE_TEST_EQUAL(Interval(-1./3,2./3).width().raw(),1.0);

    // Tests for inexact operations
    ARIADNE_TEST_EQUAL(Interval(div_down(-1,3),div_up(2,3)).lower().raw(),-0.33333333333333337034);
    ARIADNE_TEST_EQUAL(Interval(div_down(-1,3),div_up(2,3)).upper().raw(),0.66666666666666674068);
    ARIADNE_TEST_EQUAL(Interval(div_down(-1,3),div_up(2,3)).centre().raw(),0.16666666666666668517);
    ARIADNE_TEST_EQUAL(Interval(div_down(-1,3),div_up(2,3)).error().raw(),0.50000000000000011102)
    ARIADNE_TEST_EQUAL(Interval(div_down(-1,3),div_up(2,3)).width().raw(),1.000000000000000222);
}

void TestInterval::test_input()
{
    Interval ivl1,ivl2;
    string input("[1.125,2.25] [0.4,0.6]");
    stringstream iss(input);

    iss >> ivl1;
    ivl2=Interval(1.125,2.25);

    cout << "ivl1=" << ivl1 << "  ivl2=" << ivl2 << endl;
    ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,ivl2);
    ARIADNE_TEST_ASSERT(ivl1.lower()==ivl2.lower() && ivl1.upper()==ivl2.upper());
    ARIADNE_TEST_ASSERT(equal(ivl1,ivl2));

    iss >> ivl1;
    ivl2=Interval(0.39999999999999997,0.60000000000000009);
    if(!equal(ivl1,ivl2)) {
        ARIADNE_TEST_WARN("Interval string constructor returns an approximate interval, not an outwardly rounded interval.");
    }
}

void TestInterval::test_comparison() {
    // FIXME: If using Boost style interval tests, uncomment the line below
    // and comment out the line after
    //ARIADNE_TEST_ASSERT(indeterminate(ivld1==ivld2));
    Interval ivl1(1.125,2.25);
    Interval ivl2=ivl1;

    ARIADNE_TEST_ASSERT(ivl1==ivl2);
    Interval& ivl1ref=ivl1;
    ivl1ref=Interval(5.25,7.375);
    cout << "ivl1ref=" << ivl1ref << endl;
    ARIADNE_TEST_ASSERT(ivl1ref.lower().raw()==Float(5.25));
}

void TestInterval::test_aliasing() {

    Float x2=1.5;
    Float x3=2.25;

    Interval ivl1;
    Interval ivl2(1.5,2.25);
    Interval ivl3(3.125,4.0625);

    // Check to make sure aliases are handled correctly
    ivl1=ivl3; ivl1=ivl2-ivl1; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,Interval(ivl2-ivl3));
    ivl1=ivl3; ivl1=ivl2*ivl1; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,Interval(ivl2*ivl3));
    ivl1=ivl2; ivl1=ivl1*ivl3; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,Interval(ivl2*ivl3));
    ivl1=ivl2; ivl1=ivl1*x3; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,Interval(ivl2*x3));
    ivl1=ivl3; ivl1=x2*ivl1; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,Interval(x2*ivl3));
    ivl1=ivl2; ivl1=ivl1/ivl3; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,Interval(ivl2/ivl3));
    ivl1=ivl2; ivl1=ivl1/x3; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,Interval(ivl2/x3));
    ivl1=ivl3; ivl1=x2/ivl1; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,Interval(x2/ivl3));
}

void TestInterval::test_monotone_functions()
{

    Interval two(2.0);
    Interval sqrttwo=sqrt(two);
    ARIADNE_TEST_PRINT(sqrttwo);
    ARIADNE_TEST_COMPARE(sqrttwo.lower(),<=,1.4142135623730949);
    ARIADNE_TEST_COMPARE(sqrttwo.lower(),> ,1.4142135623730947);
    ARIADNE_TEST_COMPARE(sqrttwo.upper(),>=,1.4142135623730951);
    ARIADNE_TEST_COMPARE(sqrttwo.upper(),< ,1.4142135623730954);

    Interval one(1.0);
    Interval expone=exp(one);
    ARIADNE_TEST_PRINT(expone);
    ARIADNE_TEST_COMPARE(expone.lower(),<,2.71828182845905);
    ARIADNE_TEST_COMPARE(expone.lower(),>,2.71828182845903);
    ARIADNE_TEST_COMPARE(expone.upper(),>,2.71828182845904);
    ARIADNE_TEST_COMPARE(expone.upper(),<,2.71828182845906);
    ARIADNE_TEST_ASSERT(expone.lower()<expone.upper());

    Interval e(2.7182818284590451,2.7182818284590455);
    Interval loge=log(e);
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_COMPARE(loge.lower(),<,1);
    ARIADNE_TEST_COMPARE(loge.lower(),>,0.9999999999998);
    ARIADNE_TEST_COMPARE(loge.upper(),>,1);
    ARIADNE_TEST_COMPARE(loge.upper(),<,1.000000000002);
}

void TestInterval::test_trigonometric_functions()
{
    try {
        Interval x(6.283185307179586,6.283185307179587);
        Interval sinx=sin(x);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_COMPARE(sinx.lower(),<,0.0);
        ARIADNE_TEST_COMPARE(sinx.lower(),>,-1e-14);
        ARIADNE_TEST_COMPARE(sinx.upper(),>,0.0);
        ARIADNE_TEST_COMPARE(sinx.upper(),<,+1e-14);
        ARIADNE_TEST_ASSERT(sinx.lower()<sinx.upper());
    }
    catch(...) { }

    try {
        Interval x(7.0685834705770345);
        Interval sinx=sin(x);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_COMPARE(sinx.lower(),<,0.7071067811866);
        ARIADNE_TEST_COMPARE(sinx.upper(),>,0.7071067811865);
        ARIADNE_TEST_ASSERT(sinx.lower()<sinx.upper());
    }
    catch(...) { }

}

void TestInterval::test_geometric_predicates()
{
    Interval empty_interval; empty_interval.set_empty();

    ARIADNE_TEST_PRINT(empty_interval);

    ARIADNE_TEST_BINARY_PREDICATE(intersect,Interval(0.0,1.5),Interval(1.5,3));
    ARIADNE_TEST_BINARY_PREDICATE(!intersect,Interval(0.0,1.5),Interval(1.625,3));

    ARIADNE_TEST_BINARY_PREDICATE(!disjoint,Interval(0.0,1.5),Interval(1.5,3));
    ARIADNE_TEST_BINARY_PREDICATE(disjoint,Interval(0.0,1.5),Interval(1.625,3));

    ARIADNE_TEST_BINARY_PREDICATE(subset,Interval(1.0,1.5),Interval(1.0,2.0));
    ARIADNE_TEST_BINARY_PREDICATE(subset,Interval(1.5,2.0),Interval(1.0,2.0));
    ARIADNE_TEST_BINARY_PREDICATE(subset,Interval(1.0,2.0),Interval(1.0,2.0));

    ARIADNE_TEST_BINARY_PREDICATE(superset,Interval(1.0,2.0),Interval(1.0,1.5));
    ARIADNE_TEST_BINARY_PREDICATE(superset,Interval(1.0,2.0),Interval(1.5,2.0));
    ARIADNE_TEST_BINARY_PREDICATE(superset,Interval(1.0,2.0),Interval(1.0,2.0));

    ARIADNE_TEST_BINARY_PREDICATE(inside,Interval(1.25,1.75),Interval(1.0,2.0));
    ARIADNE_TEST_BINARY_PREDICATE(!inside,Interval(1.00,1.75),Interval(1.0,2.0));
    ARIADNE_TEST_BINARY_PREDICATE(!inside,Interval(1.25,2.00),Interval(1.0,2.0));

    ARIADNE_TEST_BINARY_PREDICATE(covers,Interval(1.0,2.0),Interval(1.25,1.75));
    ARIADNE_TEST_BINARY_PREDICATE(!covers,Interval(1.0,2.0),Interval(1.00,1.75));
    ARIADNE_TEST_BINARY_PREDICATE(!covers,Interval(1.0,2.0),Interval(1.25,2.00));

    ARIADNE_TEST_BINARY_PREDICATE(overlap,Interval(1.0,2.0),Interval(1.75,3.0));
    ARIADNE_TEST_BINARY_PREDICATE(!overlap,Interval(1.0,2.0),Interval(2.00,2.75));
    ARIADNE_TEST_BINARY_PREDICATE(!overlap,Interval(1.0,2.0),Interval(2.25,2.75));

    ARIADNE_TEST_BINARY_PREDICATE(!separated,Interval(1.0,2.0),Interval(1.75,3.0));
    ARIADNE_TEST_BINARY_PREDICATE(!separated,Interval(1.0,2.0),Interval(2.00,2.75));
    ARIADNE_TEST_BINARY_PREDICATE(separated,Interval(1.0,2.0),Interval(2.25,2.75));

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

    ARIADNE_TEST_BINARY_PREDICATE(disjoint,Interval(-1,2),empty_interval);
    ARIADNE_TEST_BINARY_PREDICATE(subset,empty_interval,Interval(0.25,0.75));
    ARIADNE_TEST_BINARY_PREDICATE(covers,Interval(0.25,0.75),empty_interval);
    ARIADNE_TEST_BINARY_PREDICATE(inside,empty_interval,Interval(0.25,0.75));
}

void TestInterval::regression_tests() {

    // Regression test; fails dramatically on certain types of rounding
    {
        Interval x(1.5707963267948966,1.5707963267948968);
        Interval cosx=cos(x);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_COMPARE(cosx.lower(),<,0.0);
        ARIADNE_TEST_COMPARE(cosx.lower(),>,-1e-14);
        ARIADNE_TEST_COMPARE(cosx.upper(),>,0.0);
        ARIADNE_TEST_COMPARE(cosx.upper(),<,+1e-14);
        ARIADNE_TEST_ASSERT(cosx.lower()<cosx.upper());
    }

}


int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    TestInterval().test();

    return ARIADNE_TEST_FAILURES;
}

