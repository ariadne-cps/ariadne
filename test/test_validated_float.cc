/***************************************************************************
 *            test_validated_float.cc
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

#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include "config.h"
#include "numeric/rational.h"
#include "numeric/float-user.h"

#include "test.h"

using namespace Ariadne;
using namespace std;


class TestValidatedFloat
{
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
    Void regression_tests();
};


Void
TestValidatedFloat::test()
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
    ARIADNE_TEST_CALL(regression_tests());
}

Void
TestValidatedFloat::test_concept()
{
    ValidatedFloat64::set_output_precision(17);

    Int n=1;
    Nat m=1;
    double d=1;
    ExactFloat64 x(1);
    Float64 a,b;
    ValidatedFloat64 i(1);
    ValidatedFloat64 j(1);


    // Constructors
    j=ValidatedFloat64(); j=ValidatedFloat64(n); j=ValidatedFloat64(m); j=ValidatedFloat64(d); j=ValidatedFloat64(x); j=ValidatedFloat64(i);
    j=ValidatedFloat64(n,n); j=ValidatedFloat64(m,m); j=ValidatedFloat64(d,d); j=ValidatedFloat64(a,b);
    j=ValidatedFloat64(n,m); j=ValidatedFloat64(m,d); j=ValidatedFloat64(d,n);
    // Assignment
    j=n; j=m; j=x; j=i;

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

    j=sin(i);
    j=cos(i);
    j=tan(i);
    //j=asin(i);
    //j=acos(i);
    //j=atan(i);

    //j=sinh(i);
    //j=cosh(i);
    //j=tanh(i);
    //j=asinh(i);
    //j=acosh(i);
    //j=atanh(i);
}



// Test that interval arithmetic is rounded correctly,
// without paying attention to accuracy issues.
Void
TestValidatedFloat::test_correct_rounded_arithmetic()
{
    ValidatedFloat64 onethird=ValidatedFloat64(1)/ValidatedFloat64(3);
    ARIADNE_TEST_COMPARE( onethird.lower_raw(), < , onethird.upper_raw() );
    ValidatedFloat64 one_approx=onethird*ValidatedFloat64(3);
    ARIADNE_TEST_COMPARE( one_approx.lower_raw(), < , 1.0 );
    ARIADNE_TEST_COMPARE( one_approx.upper_raw(), > , 1.0 );
}


// Test that interval arithmetic gives the most accurate rounded values
Void
TestValidatedFloat::test_accurate_rounded_arithmetic()
{
    const double min=std::numeric_limits<double>::min();
    const double eps=std::numeric_limits<double>::epsilon();

    ARIADNE_TEST_SAME(ValidatedFloat64(1.5)+ValidatedFloat64(min),ValidatedFloat64(1.5,1.5+eps));
    ARIADNE_TEST_SAME(ValidatedFloat64(1.5)-ValidatedFloat64(min),ValidatedFloat64(1.5-eps,1.5));
    ARIADNE_TEST_SAME(ValidatedFloat64(1+eps,1+2*eps)*ValidatedFloat64(1+eps,1+3*eps),ValidatedFloat64(1+2*eps,1+6*eps));
    ARIADNE_TEST_SAME(ValidatedFloat64(1)/ValidatedFloat64(3),ValidatedFloat64(0.33333333333333331,0.33333333333333337));
    ARIADNE_TEST_SAME(ValidatedFloat64(2)/ValidatedFloat64(5),ValidatedFloat64(0.39999999999999997,0.40000000000000002));

    ARIADNE_TEST_SAME(ValidatedFloat64(1.5)+ExactFloat64(min),ValidatedFloat64(1.5,1.5+eps));
    ARIADNE_TEST_SAME(ValidatedFloat64(1.5)-ExactFloat64(min),ValidatedFloat64(1.5-eps,1.5));
    ARIADNE_TEST_SAME(ValidatedFloat64(1+eps,1+2*eps)*ExactFloat64(1+eps),ValidatedFloat64(1+2*eps,1+4*eps));
    ARIADNE_TEST_SAME(ValidatedFloat64(1+3*eps,1+5*eps)/ExactFloat64(1+eps),ValidatedFloat64(1+eps,1+4*eps));

    ARIADNE_TEST_SAME(ExactFloat64(min)-ValidatedFloat64(1.5),ValidatedFloat64(-1.5,eps-1.5));
    ARIADNE_TEST_SAME(ExactFloat64(1+5*eps)/ValidatedFloat64(1+2*eps,1+3*eps),ValidatedFloat64(1+eps,1+3*eps));

    ARIADNE_TEST_SAME(sqr(ValidatedFloat64(1-eps,1+eps)),ValidatedFloat64(1-4*eps/2,1+3*eps));

    ARIADNE_TEST_SAME(pow(ValidatedFloat64(3,5),-1),ValidatedFloat64(0.19999999999999998,0.33333333333333337));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(3,5),-2),ValidatedFloat64(0.039999999999999986955,0.11111111111111114658));

    ARIADNE_TEST_SAME(rec(ValidatedFloat64(1+2*eps,1+5*eps)),ValidatedFloat64(1-10*eps/2,1-3*eps/2));

}


// Test that interval arithmetic gives exact values if possible
Void
TestValidatedFloat::test_exact_rounded_arithmetic()
{
    ARIADNE_TEST_SAME(ValidatedFloat64(5,7)+ValidatedFloat64(2,4),ValidatedFloat64(7,11));
    ARIADNE_TEST_SAME(ValidatedFloat64(5,7)-ValidatedFloat64(2,6),ValidatedFloat64(-1,5));

    ARIADNE_TEST_SAME(ValidatedFloat64(5,7)*ValidatedFloat64(2,4),ValidatedFloat64(10,28));
    ARIADNE_TEST_SAME(ValidatedFloat64(5,7)*ValidatedFloat64(-2,4),ValidatedFloat64(-14,28));
    ARIADNE_TEST_SAME(ValidatedFloat64(5,7)*ValidatedFloat64(-4,-2),ValidatedFloat64(-28,-10));
    ARIADNE_TEST_SAME(ValidatedFloat64(-7,5)*ValidatedFloat64(2,4),ValidatedFloat64(-28,20));
    ARIADNE_TEST_SAME(ValidatedFloat64(-7,5)*ValidatedFloat64(-2,4),ValidatedFloat64(-28,20));
    ARIADNE_TEST_SAME(ValidatedFloat64(-7,5)*ValidatedFloat64(-4,-2),ValidatedFloat64(-20,28));
    ARIADNE_TEST_SAME(ValidatedFloat64(-7,-5)*ValidatedFloat64(2,4),ValidatedFloat64(-28,-10));
    ARIADNE_TEST_SAME(ValidatedFloat64(-7,-5)*ValidatedFloat64(-2,4),ValidatedFloat64(-28,14));
    ARIADNE_TEST_SAME(ValidatedFloat64(-7,-5)*ValidatedFloat64(-4,-2),ValidatedFloat64(10,28));

    ARIADNE_TEST_SAME(ValidatedFloat64(5,7)/ValidatedFloat64(2,4),ValidatedFloat64(1.25,3.50));
    ARIADNE_TEST_SAME(ValidatedFloat64(5,7)/ValidatedFloat64(-4,-2),ValidatedFloat64(-3.50,-1.25));
    ARIADNE_TEST_SAME(ValidatedFloat64(-7,5)/ValidatedFloat64(2,4),ValidatedFloat64(-3.50,2.50));
    ARIADNE_TEST_SAME(ValidatedFloat64(-7,5)/ValidatedFloat64(-4,-2),ValidatedFloat64(-2.50,3.5));
    ARIADNE_TEST_SAME(ValidatedFloat64(-7,-5)/ValidatedFloat64(2,4),ValidatedFloat64(-3.50,-1.25));
    ARIADNE_TEST_SAME(ValidatedFloat64(-7,-5)/ValidatedFloat64(-4,-2),ValidatedFloat64(1.25,3.50));

    ARIADNE_TEST_SAME(pow(ValidatedFloat64(5,7),0u),ValidatedFloat64(1,1));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-5,7),0u),ValidatedFloat64(1,1));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,5),0u),ValidatedFloat64(1,1));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,-5),0u),ValidatedFloat64(1,1));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(5,7),1u),ValidatedFloat64(5,7));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-5,7),1u),ValidatedFloat64(-5,7));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,5),1u),ValidatedFloat64(-7,5));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,-5),1u),ValidatedFloat64(-7,-5));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(5,7),2u),ValidatedFloat64(25,49));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-5,7),2u),ValidatedFloat64(0,49));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,5),2u),ValidatedFloat64(0,49));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,-5),2u),ValidatedFloat64(25,49));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(5,7),3u),ValidatedFloat64(125,343));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-5,7),3u),ValidatedFloat64(-125,343));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,5),3u),ValidatedFloat64(-343,125));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,-5),3u),ValidatedFloat64(-343,-125));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(5,7),4u),ValidatedFloat64(625,2401));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-5,7),4u),ValidatedFloat64(0,2401));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,5),4u),ValidatedFloat64(0,2401));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,-5),4u),ValidatedFloat64(625,2401));

    ARIADNE_TEST_SAME(pow(ValidatedFloat64(5,7),0),ValidatedFloat64(1,1));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,5),0),ValidatedFloat64(1,1));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,-5),0),ValidatedFloat64(1,1));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(5,7),1),ValidatedFloat64(5,7));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-5,7),1),ValidatedFloat64(-5,7));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,5),1),ValidatedFloat64(-7,5));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,-5),1),ValidatedFloat64(-7,-5));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(5,7),2),ValidatedFloat64(25,49));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-5,7),2),ValidatedFloat64(0,49));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,5),2),ValidatedFloat64(0,49));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,-5),2),ValidatedFloat64(25,49));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(5,7),3),ValidatedFloat64(125,343));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-5,7),3),ValidatedFloat64(-125,343));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,5),3),ValidatedFloat64(-343,125));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,-5),3),ValidatedFloat64(-343,-125));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(5,7),4),ValidatedFloat64(625,2401));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-5,7),4),ValidatedFloat64(0,2401));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,5),4),ValidatedFloat64(0,2401));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,-5),4),ValidatedFloat64(625,2401));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(5,7),5),ValidatedFloat64(3125,16807));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-5,7),5),ValidatedFloat64(-3125,16807));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,5),5),ValidatedFloat64(-16807,3125));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,-5),5),ValidatedFloat64(-16807,-3125));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(5,7),7),ValidatedFloat64(78125,823543));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-5,7),7),ValidatedFloat64(-78125,823543));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,5),7),ValidatedFloat64(-823543,78125));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-7,-5),7),ValidatedFloat64(-823543,-78125));

    ARIADNE_TEST_SAME(pow(ValidatedFloat64(2,4),-1),ValidatedFloat64(0.25,0.5));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-4,-2),-1),ValidatedFloat64(-0.5,-0.25));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(2,4),-2),ValidatedFloat64(0.0625,0.25));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-4,-2),-2),ValidatedFloat64(0.0625,0.25));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(2,4),-3),ValidatedFloat64(0.015625,0.125));
    ARIADNE_TEST_SAME(pow(ValidatedFloat64(-4,-2),-3),ValidatedFloat64(-0.125,-0.015625));

    ARIADNE_TEST_SAME(rec(ValidatedFloat64(2,4)),ValidatedFloat64(0.25,0.50));
    ARIADNE_TEST_SAME(rec(ValidatedFloat64(-4,-2)),ValidatedFloat64(-0.50,-0.25));
}



Void
TestValidatedFloat::test_constructors()
{
    Float64 zero=0;

    // Construct from pair
    ValidatedFloat64 xd1(Float64(1.125),Float64(2.25));
    ARIADNE_TEST_ASSERT(xd1.lower_raw()==1.125); ARIADNE_TEST_ASSERT(xd1.upper_raw()==2.25);

    // Default constructor
    ValidatedFloat64 xd2;
    if(xd2.lower_raw()>xd2.upper_raw()) {
        ARIADNE_TEST_WARN("ValidatedFloat64 default constructor returns an empty set.");
    } else {
        ARIADNE_TEST_BINARY_PREDICATE(same,xd2,ValidatedFloat64(zero,zero));
    }

    // Constructor with approximations
    ValidatedFloat64 xd3(Rational(21,10));
    cout<<xd3<<std::endl;
    ARIADNE_TEST_COMPARE(Rational(xd3.lower_raw()),<,Rational(21,10));
    ARIADNE_TEST_COMPARE(Rational(xd3.upper_raw()),>,Rational(21,10));

    // Constructor from approximate values
    ValidatedFloat64 xd4(2.1,3.2);
    ARIADNE_TEST_COMPARE(xd4.lower_raw(),<=,2.1);
    ARIADNE_TEST_COMPARE(xd4.upper_raw(),>=,3.2);

    // Approximate constructor from a single value
    ValidatedFloat64 xd5(Rational(1,3));
    ARIADNE_TEST_COMPARE(Rational(xd5.lower_raw()),<,Rational(1,3));
    ARIADNE_TEST_COMPARE(Rational(xd5.upper_raw()),>,Rational(1,3));

    // Exact constructor from a single value
    ValidatedFloat64 xd6(Float64(1.25));
    ARIADNE_TEST_EQUAL(xd6.lower_raw(),Float64(1.25));
    ARIADNE_TEST_EQUAL(xd6.upper_raw(),Float64(1.25));
}

Void TestValidatedFloat::test_class()
{
    // Test lower, upper, midpoint, radius, width
    Float64 one=1.0;
    Float64 two=2.0;

    // Tests for exact operations
    ARIADNE_TEST_EQUAL(ValidatedFloat64(-0.25,0.50).lower().raw(),-0.25);
    ARIADNE_TEST_EQUAL(ValidatedFloat64(-0.25,0.50).upper().raw(),0.5);
    ARIADNE_TEST_EQUAL(ValidatedFloat64(-0.25,0.50).value().raw(),0.125);
    ARIADNE_TEST_EQUAL(ValidatedFloat64(-0.25,0.50).error().raw(),0.375)

    // Tests for inexact operations
    ARIADNE_TEST_EQUAL(ValidatedFloat64(-1./3,2./3).lower().raw(),-0.33333333333333331483);
    ARIADNE_TEST_EQUAL(ValidatedFloat64(-1./3,2./3).upper().raw(),0.66666666666666662966);
    ARIADNE_TEST_EQUAL(ValidatedFloat64(-1./3,2./3).value().raw(),0.16666666666666665741);
    ARIADNE_TEST_EQUAL(ValidatedFloat64(-1./3,2./3).error().raw(),0.5)

    // Tests for inexact operations
    ARIADNE_TEST_EQUAL(ValidatedFloat64(div_down(-one,3),div_up(two,3)).lower().raw(),-0.33333333333333337034);
    ARIADNE_TEST_EQUAL(ValidatedFloat64(div_down(-one,3),div_up(two,3)).upper().raw(),0.66666666666666674068);
    ARIADNE_TEST_EQUAL(ValidatedFloat64(div_down(-one,3),div_up(two,3)).value().raw(),0.16666666666666668517);
    ARIADNE_TEST_EQUAL(ValidatedFloat64(div_down(-one,3),div_up(two,3)).error().raw(),0.50000000000000011102);
}

Void TestValidatedFloat::test_input()
{
    ValidatedFloat64 x1,x2;
    string input("[1.125,2.25] [0.4,0.6]");
    stringstream iss(input);

    iss >> x1;
    x2=ValidatedFloat64(1.125,2.25);

    cout << "x1=" << x1 << "  x2=" << x2 << endl;
    ARIADNE_TEST_BINARY_PREDICATE(same,x1,x2);
    ARIADNE_TEST_ASSERT(x1.lower_raw()==x2.lower_raw() && x1.upper_raw()==x2.upper_raw());
    ARIADNE_TEST_ASSERT(same(x1,x2));

    iss >> x1;
    x2=ValidatedFloat64(0.39999999999999997,0.60000000000000009);
    if(!same(x1,x2)) {
        ARIADNE_TEST_WARN("ValidatedFloat64 string constructor returns an approximate interval, not an outwardly rounded interval.");
    }
}

Void TestValidatedFloat::test_comparison() {
    // FIXME: If using Boost style interval tests, uncomment the line below
    // and comment out the line after
    //ARIADNE_TEST_ASSERT(indeterminate(ivld1==ivld2));
    ValidatedFloat64 ivl1(1.125,2.25);
    ValidatedFloat64 ivl2=ivl1;

    ARIADNE_TEST_ASSERT(!definitely(ivl1==ivl2));
    ARIADNE_TEST_ASSERT(possibly(ivl1==ivl2));
    ValidatedFloat64& ivl1ref=ivl1;
    ivl1ref=ValidatedFloat64(5.25,7.375);
    cout << "ivl1ref=" << ivl1ref << endl;
    ARIADNE_TEST_ASSERT(ivl1ref.lower_raw()==Float64(5.25));
}

Void TestValidatedFloat::test_aliasing() {
    ExactFloat64 ex2(1.5);
    ExactFloat64 ex3(2.25);

    ValidatedFloat64 vx1;
    ValidatedFloat64 vx2(1.5,2.25);
    ValidatedFloat64 vx3(3.125,4.0625);

    // Check to make sure aliases are handled correctly
    vx1=vx3; vx1=vx2-vx1; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,ValidatedFloat64(vx2-vx3));
    vx1=vx3; vx1=vx2*vx1; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,ValidatedFloat64(vx2*vx3));
    vx1=vx2; vx1=vx1*vx3; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,ValidatedFloat64(vx2*vx3));
    vx1=vx2; vx1=vx1*ex3; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,ValidatedFloat64(vx2*ex3));
    vx1=vx3; vx1=ex2*vx1; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,ValidatedFloat64(ex2*vx3));
    vx1=vx2; vx1=vx1/vx3; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,ValidatedFloat64(vx2/vx3));
    vx1=vx2; vx1=vx1/ex3; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,ValidatedFloat64(vx2/ex3));
    vx1=vx3; vx1=ex2/vx1; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,ValidatedFloat64(ex2/vx3));
}

Void TestValidatedFloat::test_monotone_functions()
{

    ValidatedFloat64 two(2.0);
    ValidatedFloat64 sqrttwo=sqrt(two);
    ARIADNE_TEST_PRINT(sqrttwo);
    ARIADNE_TEST_COMPARE(sqrttwo.lower_raw(),<=,1.4142135623730949);
    ARIADNE_TEST_COMPARE(sqrttwo.lower_raw(),> ,1.4142135623730947);
    ARIADNE_TEST_COMPARE(sqrttwo.upper_raw(),>=,1.4142135623730951);
    ARIADNE_TEST_COMPARE(sqrttwo.upper_raw(),< ,1.4142135623730954);

    ValidatedFloat64 one(1.0);
    ValidatedFloat64 expone=exp(one);
    ARIADNE_TEST_PRINT(expone);
    ARIADNE_TEST_COMPARE(expone.lower_raw(),<,2.71828182845905);
    ARIADNE_TEST_COMPARE(expone.lower_raw(),>,2.71828182845903);
    ARIADNE_TEST_COMPARE(expone.upper_raw(),>,2.71828182845904);
    ARIADNE_TEST_COMPARE(expone.upper_raw(),<,2.71828182845906);
    ARIADNE_TEST_ASSERT(expone.lower_raw()<expone.upper_raw());

    ValidatedFloat64 e(2.7182818284590451,2.7182818284590455);
    ValidatedFloat64 loge=log(e);
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_COMPARE(loge.lower_raw(),<,1);
    ARIADNE_TEST_COMPARE(loge.lower_raw(),>,0.9999999999998);
    ARIADNE_TEST_COMPARE(loge.upper_raw(),>,1);
    ARIADNE_TEST_COMPARE(loge.upper_raw(),<,1.000000000002);
}

Void TestValidatedFloat::test_trigonometric_functions()
{
    try {
        ValidatedFloat64 x(6.283185307179586,6.283185307179587);
        ValidatedFloat64 sinx=sin(x);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_COMPARE(sinx.lower_raw(),<,0.0);
        ARIADNE_TEST_COMPARE(sinx.lower_raw(),>,-1e-14);
        ARIADNE_TEST_COMPARE(sinx.upper_raw(),>,0.0);
        ARIADNE_TEST_COMPARE(sinx.upper_raw(),<,+1e-14);
        ARIADNE_TEST_ASSERT(sinx.lower_raw()<sinx.upper_raw());
    }
    catch(...) { }

    try {
        ValidatedFloat64 x(7.0685834705770345);
        ValidatedFloat64 sinx=sin(x);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_COMPARE(sinx.lower_raw(),<,0.7071067811866);
        ARIADNE_TEST_COMPARE(sinx.upper_raw(),>,0.7071067811865);
        ARIADNE_TEST_ASSERT(sinx.lower_raw()<sinx.upper_raw());
    }
    catch(...) { }

}

Void TestValidatedFloat::regression_tests() {

    // Regression test; fails dramatically on certain types of rounding
    {
        ValidatedFloat64 x(1.5707963267948966,1.5707963267948968);
        ValidatedFloat64 cosx=cos(x);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_COMPARE(cosx.lower_raw(),<,0.0);
        ARIADNE_TEST_COMPARE(cosx.lower_raw(),>,-1e-14);
        ARIADNE_TEST_COMPARE(cosx.upper_raw(),>,0.0);
        ARIADNE_TEST_COMPARE(cosx.upper_raw(),<,+1e-14);
        ARIADNE_TEST_ASSERT(cosx.lower_raw()<cosx.upper_raw());
    }

    // Regression test for dividing by interval with lower endpoint -0.0 or upper endpoint +0.0

    ARIADNE_TEST_EQUAL((ValidatedFloat64(1.0,2.0)/ValidatedFloat64(-0.0,1.0)).upper_raw(),+inf);
    ARIADNE_TEST_EQUAL((ValidatedFloat64(1.0,2.0)/ValidatedFloat64(-1.0,+0.0)).lower_raw(),-inf);

    ARIADNE_TEST_EQUAL(rec(ValidatedFloat64(-0.0,+1.0)).upper_raw(),+inf);
    ARIADNE_TEST_EQUAL(rec(ValidatedFloat64(-1.0,+0.0)).lower_raw(),-inf);
}

Int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    TestValidatedFloat().test();

    return ARIADNE_TEST_FAILURES;
}

