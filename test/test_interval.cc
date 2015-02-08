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
    Void test_correct_rounded_arithmetic();
    Void test_accurate_rounded_arithmetic();
    Void test_exact_rounded_arithmetic();
    Void test_aliasing();
    Void test_monotone_functions();
    Void test_trigonometric_functions();
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
    ARIADNE_TEST_CALL(test_correct_rounded_arithmetic());
    ARIADNE_TEST_CALL(test_accurate_rounded_arithmetic());
    ARIADNE_TEST_CALL(test_exact_rounded_arithmetic());
    ARIADNE_TEST_CALL(test_monotone_functions());
    ARIADNE_TEST_CALL(test_trigonometric_functions());
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

    // Exact operations
    //j=abs(i); j=neg(i); j=rec(i);

    // Arithmetic
    //j=add(x,x); j=add(x,i); j=add(i,x); j=add(i,i);
    //j=sub(x,x); j=sub(x,i); j=sub(i,x); j=sub(i,i);
    //j=mul(x,x); j=mul(x,i); j=mul(i,x); j=mul(i,i);
    //j=div(x,x); j=div(x,i); j=div(i,x); j=div(i,i);
    //j=pow(x,n); j=pow(x,n);
    //j=pow(x,m); j=pow(x,m);

    // Transcendental functions
    /*
      j=sqrt(i);
      j=exp(i);
      j=log(i);
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
Void
TestInterval::test_correct_rounded_arithmetic()
{
    UpperInterval onethird=UpperInterval(1)/UpperInterval(3);
    ARIADNE_TEST_COMPARE( onethird.lower(), < , onethird.upper() );
    UpperInterval one_approx=onethird*UpperInterval(3);
    ARIADNE_TEST_COMPARE( one_approx.lower(), < , 1.0 );
    ARIADNE_TEST_COMPARE( one_approx.upper(), > , 1.0 );
}


// Test that interval arithmetic gives the most accurate rounded values
Void
TestInterval::test_accurate_rounded_arithmetic()
{
    const Float64 min=Float64::min();
    const Float64 eps=Float64::eps();
    const Float64 x=1.5;

    ARIADNE_TEST_EQUAL(UpperInterval(x)+UpperInterval(min),UpperInterval(x,x+eps));
    ARIADNE_TEST_EQUAL(UpperInterval(x)-UpperInterval(min),UpperInterval(x-eps,x));
    ARIADNE_TEST_EQUAL(UpperInterval(1+eps,1+2*eps)*UpperInterval(1+eps,1+3*eps),UpperInterval(1+2*eps,1+6*eps));
    ARIADNE_TEST_EQUAL(UpperInterval(1)/UpperInterval(3),UpperInterval(0.33333333333333331,0.33333333333333337));
    ARIADNE_TEST_EQUAL(UpperInterval(2)/UpperInterval(5),UpperInterval(0.39999999999999997,0.40000000000000002));

    ARIADNE_TEST_EQUAL(UpperInterval(x)+Float64(min),UpperInterval(x,x+eps));
    ARIADNE_TEST_EQUAL(UpperInterval(x)-Float64(min),UpperInterval(x-eps,x));
    ARIADNE_TEST_EQUAL(UpperInterval(1+eps,1+2*eps)*Float64(1+eps),UpperInterval(1+2*eps,1+4*eps));
    ARIADNE_TEST_EQUAL(UpperInterval(1+3*eps,1+5*eps)/Float64(1+eps),UpperInterval(1+eps,1+4*eps));

    ARIADNE_TEST_EQUAL(Float64(min)-UpperInterval(x),UpperInterval(-x,eps-x));
    ARIADNE_TEST_EQUAL(Float64(1+5*eps)/UpperInterval(1+2*eps,1+3*eps),UpperInterval(1+eps,1+3*eps));

    ARIADNE_TEST_EQUAL(sqr(UpperInterval(1-eps,1+eps)),UpperInterval(1-4*eps/2,1+3*eps));

    ARIADNE_TEST_EQUAL(pow(UpperInterval(3,5),-1),UpperInterval(0.19999999999999998,0.33333333333333337));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(3,5),-2),UpperInterval(0.039999999999999986955,0.11111111111111114658));

    ARIADNE_TEST_EQUAL(rec(UpperInterval(1+2*eps,1+5*eps)),UpperInterval(1-10*eps/2,1-3*eps/2));

}


// Test that interval arithmetic gives exact values if possible
Void
TestInterval::test_exact_rounded_arithmetic()
{
    ARIADNE_TEST_EQUAL(UpperInterval(5,7)+UpperInterval(2,4),UpperInterval(7,11));
    ARIADNE_TEST_EQUAL(UpperInterval(5,7)-UpperInterval(2,6),UpperInterval(-1,5));

    ARIADNE_TEST_EQUAL(UpperInterval(5,7)*UpperInterval(2,4),UpperInterval(10,28));
    ARIADNE_TEST_EQUAL(UpperInterval(5,7)*UpperInterval(-2,4),UpperInterval(-14,28));
    ARIADNE_TEST_EQUAL(UpperInterval(5,7)*UpperInterval(-4,-2),UpperInterval(-28,-10));
    ARIADNE_TEST_EQUAL(UpperInterval(-7,5)*UpperInterval(2,4),UpperInterval(-28,20));
    ARIADNE_TEST_EQUAL(UpperInterval(-7,5)*UpperInterval(-2,4),UpperInterval(-28,20));
    ARIADNE_TEST_EQUAL(UpperInterval(-7,5)*UpperInterval(-4,-2),UpperInterval(-20,28));
    ARIADNE_TEST_EQUAL(UpperInterval(-7,-5)*UpperInterval(2,4),UpperInterval(-28,-10));
    ARIADNE_TEST_EQUAL(UpperInterval(-7,-5)*UpperInterval(-2,4),UpperInterval(-28,14));
    ARIADNE_TEST_EQUAL(UpperInterval(-7,-5)*UpperInterval(-4,-2),UpperInterval(10,28));

    ARIADNE_TEST_EQUAL(UpperInterval(5,7)/UpperInterval(2,4),UpperInterval(1.25,3.50));
    ARIADNE_TEST_EQUAL(UpperInterval(5,7)/UpperInterval(-4,-2),UpperInterval(-3.50,-1.25));
    ARIADNE_TEST_EQUAL(UpperInterval(-7,5)/UpperInterval(2,4),UpperInterval(-3.50,2.50));
    ARIADNE_TEST_EQUAL(UpperInterval(-7,5)/UpperInterval(-4,-2),UpperInterval(-2.50,3.5));
    ARIADNE_TEST_EQUAL(UpperInterval(-7,-5)/UpperInterval(2,4),UpperInterval(-3.50,-1.25));
    ARIADNE_TEST_EQUAL(UpperInterval(-7,-5)/UpperInterval(-4,-2),UpperInterval(1.25,3.50));

    ARIADNE_TEST_EQUAL(pow(UpperInterval(5,7),0u),UpperInterval(1,1));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-5,7),0u),UpperInterval(1,1));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,5),0u),UpperInterval(1,1));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,-5),0u),UpperInterval(1,1));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(5,7),1u),UpperInterval(5,7));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-5,7),1u),UpperInterval(-5,7));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,5),1u),UpperInterval(-7,5));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,-5),1u),UpperInterval(-7,-5));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(5,7),2u),UpperInterval(25,49));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-5,7),2u),UpperInterval(0,49));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,5),2u),UpperInterval(0,49));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,-5),2u),UpperInterval(25,49));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(5,7),3u),UpperInterval(125,343));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-5,7),3u),UpperInterval(-125,343));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,5),3u),UpperInterval(-343,125));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,-5),3u),UpperInterval(-343,-125));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(5,7),4u),UpperInterval(625,2401));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-5,7),4u),UpperInterval(0,2401));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,5),4u),UpperInterval(0,2401));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,-5),4u),UpperInterval(625,2401));

    ARIADNE_TEST_EQUAL(pow(UpperInterval(5,7),0),UpperInterval(1,1));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,5),0),UpperInterval(1,1));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,-5),0),UpperInterval(1,1));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(5,7),1),UpperInterval(5,7));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-5,7),1),UpperInterval(-5,7));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,5),1),UpperInterval(-7,5));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,-5),1),UpperInterval(-7,-5));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(5,7),2),UpperInterval(25,49));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-5,7),2),UpperInterval(0,49));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,5),2),UpperInterval(0,49));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,-5),2),UpperInterval(25,49));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(5,7),3),UpperInterval(125,343));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-5,7),3),UpperInterval(-125,343));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,5),3),UpperInterval(-343,125));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,-5),3),UpperInterval(-343,-125));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(5,7),4),UpperInterval(625,2401));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-5,7),4),UpperInterval(0,2401));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,5),4),UpperInterval(0,2401));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,-5),4),UpperInterval(625,2401));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(5,7),5),UpperInterval(3125,16807));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-5,7),5),UpperInterval(-3125,16807));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,5),5),UpperInterval(-16807,3125));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,-5),5),UpperInterval(-16807,-3125));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(5,7),7),UpperInterval(78125,823543));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-5,7),7),UpperInterval(-78125,823543));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,5),7),UpperInterval(-823543,78125));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-7,-5),7),UpperInterval(-823543,-78125));

    ARIADNE_TEST_EQUAL(pow(UpperInterval(2,4),-1),UpperInterval(0.25,0.5));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-4,-2),-1),UpperInterval(-0.5,-0.25));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(2,4),-2),UpperInterval(0.0625,0.25));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-4,-2),-2),UpperInterval(0.0625,0.25));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(2,4),-3),UpperInterval(0.015625,0.125));
    ARIADNE_TEST_EQUAL(pow(UpperInterval(-4,-2),-3),UpperInterval(-0.125,-0.015625));

    ARIADNE_TEST_EQUAL(rec(UpperInterval(2,4)),UpperInterval(0.25,0.50));
    ARIADNE_TEST_EQUAL(rec(UpperInterval(-4,-2)),UpperInterval(-0.50,-0.25));
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

    // Tests for inexact operations
    ARIADNE_TEST_EQUAL(ExactInterval(-1./3,2./3).lower().raw(),-0.33333333333333331483);
    ARIADNE_TEST_EQUAL(ExactInterval(-1./3,2./3).upper().raw(),0.66666666666666662966);
    ARIADNE_TEST_EQUAL(ExactInterval(-1./3,2./3).centre().raw(),0.16666666666666665741);
    ARIADNE_TEST_EQUAL(ExactInterval(-1./3,2./3).error().raw(),0.5)
    ARIADNE_TEST_EQUAL(ExactInterval(-1./3,2./3).width().raw(),1.0);

    // Tests for inexact operations
    ARIADNE_TEST_EQUAL(ExactInterval(div_down(-1,3),div_up(2,3)).lower().raw(),-0.33333333333333337034);
    ARIADNE_TEST_EQUAL(ExactInterval(div_down(-1,3),div_up(2,3)).upper().raw(),0.66666666666666674068);
    ARIADNE_TEST_EQUAL(ExactInterval(div_down(-1,3),div_up(2,3)).centre().raw(),0.16666666666666668517);
    ARIADNE_TEST_EQUAL(ExactInterval(div_down(-1,3),div_up(2,3)).error().raw(),0.50000000000000011102)
    ARIADNE_TEST_EQUAL(ExactInterval(div_down(-1,3),div_up(2,3)).width().raw(),1.000000000000000222);
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

Void TestInterval::test_aliasing() {

    Float64 x2=1.5;
    Float64 x3=2.25;

    UpperInterval ivl1;
    UpperInterval ivl2(1.5,2.25);
    UpperInterval ivl3(3.125,4.0625);

    // Check to make sure aliases are handled correctly
    ivl1=ivl3; ivl1=ivl2-ivl1; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,UpperInterval(ivl2-ivl3));
    ivl1=ivl3; ivl1=ivl2*ivl1; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,UpperInterval(ivl2*ivl3));
    ivl1=ivl2; ivl1=ivl1*ivl3; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,UpperInterval(ivl2*ivl3));
    ivl1=ivl2; ivl1=ivl1*x3; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,UpperInterval(ivl2*x3));
    ivl1=ivl3; ivl1=x2*ivl1; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,UpperInterval(x2*ivl3));
    ivl1=ivl2; ivl1=ivl1/ivl3; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,UpperInterval(ivl2/ivl3));
    ivl1=ivl2; ivl1=ivl1/x3; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,UpperInterval(ivl2/x3));
    ivl1=ivl3; ivl1=x2/ivl1; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,UpperInterval(x2/ivl3));
}

Void TestInterval::test_monotone_functions()
{

    UpperInterval two(2.0);
    UpperInterval sqrttwo=sqrt(two);
    ARIADNE_TEST_PRINT(sqrttwo);
    ARIADNE_TEST_COMPARE(sqrttwo.lower(),<=,1.4142135623730949);
    ARIADNE_TEST_COMPARE(sqrttwo.lower(),> ,1.4142135623730947);
    ARIADNE_TEST_COMPARE(sqrttwo.upper(),>=,1.4142135623730951);
    ARIADNE_TEST_COMPARE(sqrttwo.upper(),< ,1.4142135623730954);

    UpperInterval one(1.0);
    UpperInterval expone=exp(one);
    ARIADNE_TEST_PRINT(expone);
    ARIADNE_TEST_COMPARE(expone.lower(),<,2.71828182845905);
    ARIADNE_TEST_COMPARE(expone.lower(),>,2.71828182845903);
    ARIADNE_TEST_COMPARE(expone.upper(),>,2.71828182845904);
    ARIADNE_TEST_COMPARE(expone.upper(),<,2.71828182845906);
    ARIADNE_TEST_ASSERT(expone.lower()<expone.upper());

    UpperInterval e(2.7182818284590451,2.7182818284590455);
    UpperInterval loge=log(e);
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_COMPARE(loge.lower(),<,1);
    ARIADNE_TEST_COMPARE(loge.lower(),>,0.9999999999998);
    ARIADNE_TEST_COMPARE(loge.upper(),>,1);
    ARIADNE_TEST_COMPARE(loge.upper(),<,1.000000000002);
}

Void TestInterval::test_trigonometric_functions()
{
    try {
        UpperInterval x(6.283185307179586,6.283185307179587);
        UpperInterval sinx=sin(x);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_COMPARE(sinx.lower(),<,0.0);
        ARIADNE_TEST_COMPARE(sinx.lower(),>,-1e-14);
        ARIADNE_TEST_COMPARE(sinx.upper(),>,0.0);
        ARIADNE_TEST_COMPARE(sinx.upper(),<,+1e-14);
        ARIADNE_TEST_ASSERT(sinx.lower()<sinx.upper());
    }
    catch(...) { }

    try {
        UpperInterval x(7.0685834705770345);
        UpperInterval sinx=sin(x);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_COMPARE(sinx.lower(),<,0.7071067811866);
        ARIADNE_TEST_COMPARE(sinx.upper(),>,0.7071067811865);
        ARIADNE_TEST_ASSERT(sinx.lower()<sinx.upper());
    }
    catch(...) { }

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
        ARIADNE_TEST_ASSERT(cosx.lower()<cosx.upper());
    }

}


Int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    TestInterval().test();

    return ARIADNE_TEST_FAILURES;
}

