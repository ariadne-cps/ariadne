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
    ValidatedFloat::set_output_precision(17);

    Int n=1;
    Nat m=1;
    double d=1;
    ExactFloat x(1);
    Float a,b;
    ValidatedFloat i(1);
    ValidatedFloat j(1);


    // Constructors
    j=ValidatedFloat(); j=ValidatedFloat(n); j=ValidatedFloat(m); j=ValidatedFloat(d); j=ValidatedFloat(x); j=ValidatedFloat(i);
    j=ValidatedFloat(n,n); j=ValidatedFloat(m,m); j=ValidatedFloat(d,d); j=ValidatedFloat(a,b);
    j=ValidatedFloat(n,m); j=ValidatedFloat(m,d); j=ValidatedFloat(d,n);
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
    ValidatedFloat onethird=ValidatedFloat(1)/ValidatedFloat(3);
    ARIADNE_TEST_COMPARE( onethird.lower_raw(), < , onethird.upper_raw() );
    ValidatedFloat one_approx=onethird*ValidatedFloat(3);
    ARIADNE_TEST_COMPARE( one_approx.lower_raw(), < , 1.0 );
    ARIADNE_TEST_COMPARE( one_approx.upper_raw(), > , 1.0 );
}


// Test that interval arithmetic gives the most accurate rounded values
Void
TestValidatedFloat::test_accurate_rounded_arithmetic()
{
    const double min=std::numeric_limits<double>::min();
    const double eps=std::numeric_limits<double>::epsilon();

    ARIADNE_TEST_EQUAL(ValidatedFloat(1.5)+ValidatedFloat(min),ValidatedFloat(1.5,1.5+eps));
    ARIADNE_TEST_EQUAL(ValidatedFloat(1.5)-ValidatedFloat(min),ValidatedFloat(1.5-eps,1.5));
    ARIADNE_TEST_EQUAL(ValidatedFloat(1+eps,1+2*eps)*ValidatedFloat(1+eps,1+3*eps),ValidatedFloat(1+2*eps,1+6*eps));
    ARIADNE_TEST_EQUAL(ValidatedFloat(1)/ValidatedFloat(3),ValidatedFloat(0.33333333333333331,0.33333333333333337));
    ARIADNE_TEST_EQUAL(ValidatedFloat(2)/ValidatedFloat(5),ValidatedFloat(0.39999999999999997,0.40000000000000002));

    ARIADNE_TEST_EQUAL(ValidatedFloat(1.5)+ExactFloat(min),ValidatedFloat(1.5,1.5+eps));
    ARIADNE_TEST_EQUAL(ValidatedFloat(1.5)-ExactFloat(min),ValidatedFloat(1.5-eps,1.5));
    ARIADNE_TEST_EQUAL(ValidatedFloat(1+eps,1+2*eps)*ExactFloat(1+eps),ValidatedFloat(1+2*eps,1+4*eps));
    ARIADNE_TEST_EQUAL(ValidatedFloat(1+3*eps,1+5*eps)/ExactFloat(1+eps),ValidatedFloat(1+eps,1+4*eps));

    ARIADNE_TEST_EQUAL(ExactFloat(min)-ValidatedFloat(1.5),ValidatedFloat(-1.5,eps-1.5));
    ARIADNE_TEST_EQUAL(ExactFloat(1+5*eps)/ValidatedFloat(1+2*eps,1+3*eps),ValidatedFloat(1+eps,1+3*eps));

    ARIADNE_TEST_EQUAL(sqr(ValidatedFloat(1-eps,1+eps)),ValidatedFloat(1-4*eps/2,1+3*eps));

    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(3,5),-1),ValidatedFloat(0.19999999999999998,0.33333333333333337));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(3,5),-2),ValidatedFloat(0.039999999999999986955,0.11111111111111114658));

    ARIADNE_TEST_EQUAL(rec(ValidatedFloat(1+2*eps,1+5*eps)),ValidatedFloat(1-10*eps/2,1-3*eps/2));

}


// Test that interval arithmetic gives exact values if possible
Void
TestValidatedFloat::test_exact_rounded_arithmetic()
{
    ARIADNE_TEST_EQUAL(ValidatedFloat(5,7)+ValidatedFloat(2,4),ValidatedFloat(7,11));
    ARIADNE_TEST_EQUAL(ValidatedFloat(5,7)-ValidatedFloat(2,6),ValidatedFloat(-1,5));

    ARIADNE_TEST_EQUAL(ValidatedFloat(5,7)*ValidatedFloat(2,4),ValidatedFloat(10,28));
    ARIADNE_TEST_EQUAL(ValidatedFloat(5,7)*ValidatedFloat(-2,4),ValidatedFloat(-14,28));
    ARIADNE_TEST_EQUAL(ValidatedFloat(5,7)*ValidatedFloat(-4,-2),ValidatedFloat(-28,-10));
    ARIADNE_TEST_EQUAL(ValidatedFloat(-7,5)*ValidatedFloat(2,4),ValidatedFloat(-28,20));
    ARIADNE_TEST_EQUAL(ValidatedFloat(-7,5)*ValidatedFloat(-2,4),ValidatedFloat(-28,20));
    ARIADNE_TEST_EQUAL(ValidatedFloat(-7,5)*ValidatedFloat(-4,-2),ValidatedFloat(-20,28));
    ARIADNE_TEST_EQUAL(ValidatedFloat(-7,-5)*ValidatedFloat(2,4),ValidatedFloat(-28,-10));
    ARIADNE_TEST_EQUAL(ValidatedFloat(-7,-5)*ValidatedFloat(-2,4),ValidatedFloat(-28,14));
    ARIADNE_TEST_EQUAL(ValidatedFloat(-7,-5)*ValidatedFloat(-4,-2),ValidatedFloat(10,28));

    ARIADNE_TEST_EQUAL(ValidatedFloat(5,7)/ValidatedFloat(2,4),ValidatedFloat(1.25,3.50));
    ARIADNE_TEST_EQUAL(ValidatedFloat(5,7)/ValidatedFloat(-4,-2),ValidatedFloat(-3.50,-1.25));
    ARIADNE_TEST_EQUAL(ValidatedFloat(-7,5)/ValidatedFloat(2,4),ValidatedFloat(-3.50,2.50));
    ARIADNE_TEST_EQUAL(ValidatedFloat(-7,5)/ValidatedFloat(-4,-2),ValidatedFloat(-2.50,3.5));
    ARIADNE_TEST_EQUAL(ValidatedFloat(-7,-5)/ValidatedFloat(2,4),ValidatedFloat(-3.50,-1.25));
    ARIADNE_TEST_EQUAL(ValidatedFloat(-7,-5)/ValidatedFloat(-4,-2),ValidatedFloat(1.25,3.50));

    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(5,7),0u),ValidatedFloat(1,1));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-5,7),0u),ValidatedFloat(1,1));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,5),0u),ValidatedFloat(1,1));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,-5),0u),ValidatedFloat(1,1));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(5,7),1u),ValidatedFloat(5,7));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-5,7),1u),ValidatedFloat(-5,7));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,5),1u),ValidatedFloat(-7,5));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,-5),1u),ValidatedFloat(-7,-5));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(5,7),2u),ValidatedFloat(25,49));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-5,7),2u),ValidatedFloat(0,49));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,5),2u),ValidatedFloat(0,49));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,-5),2u),ValidatedFloat(25,49));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(5,7),3u),ValidatedFloat(125,343));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-5,7),3u),ValidatedFloat(-125,343));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,5),3u),ValidatedFloat(-343,125));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,-5),3u),ValidatedFloat(-343,-125));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(5,7),4u),ValidatedFloat(625,2401));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-5,7),4u),ValidatedFloat(0,2401));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,5),4u),ValidatedFloat(0,2401));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,-5),4u),ValidatedFloat(625,2401));

    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(5,7),0),ValidatedFloat(1,1));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,5),0),ValidatedFloat(1,1));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,-5),0),ValidatedFloat(1,1));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(5,7),1),ValidatedFloat(5,7));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-5,7),1),ValidatedFloat(-5,7));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,5),1),ValidatedFloat(-7,5));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,-5),1),ValidatedFloat(-7,-5));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(5,7),2),ValidatedFloat(25,49));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-5,7),2),ValidatedFloat(0,49));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,5),2),ValidatedFloat(0,49));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,-5),2),ValidatedFloat(25,49));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(5,7),3),ValidatedFloat(125,343));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-5,7),3),ValidatedFloat(-125,343));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,5),3),ValidatedFloat(-343,125));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,-5),3),ValidatedFloat(-343,-125));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(5,7),4),ValidatedFloat(625,2401));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-5,7),4),ValidatedFloat(0,2401));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,5),4),ValidatedFloat(0,2401));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,-5),4),ValidatedFloat(625,2401));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(5,7),5),ValidatedFloat(3125,16807));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-5,7),5),ValidatedFloat(-3125,16807));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,5),5),ValidatedFloat(-16807,3125));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,-5),5),ValidatedFloat(-16807,-3125));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(5,7),7),ValidatedFloat(78125,823543));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-5,7),7),ValidatedFloat(-78125,823543));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,5),7),ValidatedFloat(-823543,78125));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-7,-5),7),ValidatedFloat(-823543,-78125));

    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(2,4),-1),ValidatedFloat(0.25,0.5));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-4,-2),-1),ValidatedFloat(-0.5,-0.25));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(2,4),-2),ValidatedFloat(0.0625,0.25));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-4,-2),-2),ValidatedFloat(0.0625,0.25));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(2,4),-3),ValidatedFloat(0.015625,0.125));
    ARIADNE_TEST_EQUAL(pow(ValidatedFloat(-4,-2),-3),ValidatedFloat(-0.125,-0.015625));

    ARIADNE_TEST_EQUAL(rec(ValidatedFloat(2,4)),ValidatedFloat(0.25,0.50));
    ARIADNE_TEST_EQUAL(rec(ValidatedFloat(-4,-2)),ValidatedFloat(-0.50,-0.25));
}



Void
TestValidatedFloat::test_constructors()
{
    Float zero=0;

    // Construct from pair
    ValidatedFloat xd1(Float(1.125),Float(2.25));
    ARIADNE_TEST_ASSERT(xd1.lower_raw()==1.125); ARIADNE_TEST_ASSERT(xd1.upper_raw()==2.25);

    // Default constructor
    ValidatedFloat xd2;
    if(xd2.lower_raw()>xd2.upper_raw()) {
        ARIADNE_TEST_WARN("ValidatedFloat default constructor returns an empty set.");
    } else {
        ARIADNE_TEST_ASSERT((Bool)(xd2==ValidatedFloat(zero,zero)));
    }

    // Constructor with approximations
    ValidatedFloat xd3(Rational(21,10),Rational(16,5));
    cout<<xd3<<std::endl;
    ARIADNE_TEST_COMPARE(Rational(xd3.lower_raw()),<,Rational(21,10));
    ARIADNE_TEST_COMPARE(Rational(xd3.upper_raw()),>,Rational(16,5));

    // Constructor from approximate values
    ValidatedFloat xd4(2.1,3.2);
    ARIADNE_TEST_COMPARE(xd4.lower_raw(),<=,2.1);
    ARIADNE_TEST_COMPARE(xd4.upper_raw(),>=,3.2);

    // Approximate constructor from a single value
    ValidatedFloat xd5(Rational(1,3));
    ARIADNE_TEST_COMPARE(Rational(xd5.lower_raw()),<,Rational(1,3));
    ARIADNE_TEST_COMPARE(Rational(xd5.upper_raw()),>,Rational(1,3));

    // Exact constructor from a single value
    ValidatedFloat xd6(Float(1.25));
    ARIADNE_TEST_EQUAL(xd6.lower_raw(),Float(1.25));
    ARIADNE_TEST_EQUAL(xd6.upper_raw(),Float(1.25));
}

Void TestValidatedFloat::test_class()
{
    // Test lower, upper, midpoint, radius, width

    // Tests for exact operations
    ARIADNE_TEST_EQUAL(ValidatedFloat(-0.25,0.50).lower().raw(),-0.25);
    ARIADNE_TEST_EQUAL(ValidatedFloat(-0.25,0.50).upper().raw(),0.5);
    ARIADNE_TEST_EQUAL(ValidatedFloat(-0.25,0.50).value().raw(),0.125);
    ARIADNE_TEST_EQUAL(ValidatedFloat(-0.25,0.50).error().raw(),0.375)

    // Tests for inexact operations
    ARIADNE_TEST_EQUAL(ValidatedFloat(-1./3,2./3).lower().raw(),-0.33333333333333331483);
    ARIADNE_TEST_EQUAL(ValidatedFloat(-1./3,2./3).upper().raw(),0.66666666666666662966);
    ARIADNE_TEST_EQUAL(ValidatedFloat(-1./3,2./3).value().raw(),0.16666666666666665741);
    ARIADNE_TEST_EQUAL(ValidatedFloat(-1./3,2./3).error().raw(),0.5)

    // Tests for inexact operations
    ARIADNE_TEST_EQUAL(ValidatedFloat(div_down(-1,3),div_up(2,3)).lower().raw(),-0.33333333333333337034);
    ARIADNE_TEST_EQUAL(ValidatedFloat(div_down(-1,3),div_up(2,3)).upper().raw(),0.66666666666666674068);
    ARIADNE_TEST_EQUAL(ValidatedFloat(div_down(-1,3),div_up(2,3)).value().raw(),0.16666666666666668517);
    ARIADNE_TEST_EQUAL(ValidatedFloat(div_down(-1,3),div_up(2,3)).error().raw(),0.50000000000000011102)
}

Void TestValidatedFloat::test_input()
{
    ValidatedFloat x1,x2;
    string input("[1.125,2.25] [0.4,0.6]");
    stringstream iss(input);

    iss >> x1;
    x2=ValidatedFloat(1.125,2.25);

    cout << "x1=" << x1 << "  x2=" << x2 << endl;
    ARIADNE_TEST_BINARY_PREDICATE(same,x1,x2);
    ARIADNE_TEST_ASSERT(x1.lower_raw()==x2.lower_raw() && x1.upper_raw()==x2.upper_raw());
    ARIADNE_TEST_ASSERT(same(x1,x2));

    iss >> x1;
    x2=ValidatedFloat(0.39999999999999997,0.60000000000000009);
    if(!same(x1,x2)) {
        ARIADNE_TEST_WARN("ValidatedFloat string constructor returns an approximate interval, not an outwardly rounded interval.");
    }
}

Void TestValidatedFloat::test_comparison() {
    // FIXME: If using Boost style interval tests, uncomment the line below
    // and comment out the line after
    //ARIADNE_TEST_ASSERT(indeterminate(ivld1==ivld2));
    ValidatedFloat ivl1(1.125,2.25);
    ValidatedFloat ivl2=ivl1;

    ARIADNE_TEST_ASSERT(ivl1==ivl2);
    ValidatedFloat& ivl1ref=ivl1;
    ivl1ref=ValidatedFloat(5.25,7.375);
    cout << "ivl1ref=" << ivl1ref << endl;
    ARIADNE_TEST_ASSERT(ivl1ref.lower_raw()==Float(5.25));
}

Void TestValidatedFloat::test_aliasing() {
    ExactFloat ex2(1.5);
    ExactFloat ex3(2.25);

    ValidatedFloat vx1;
    ValidatedFloat vx2(1.5,2.25);
    ValidatedFloat vx3(3.125,4.0625);

    // Check to make sure aliases are handled correctly
    vx1=vx3; vx1=vx2-vx1; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,ValidatedFloat(vx2-vx3));
    vx1=vx3; vx1=vx2*vx1; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,ValidatedFloat(vx2*vx3));
    vx1=vx2; vx1=vx1*vx3; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,ValidatedFloat(vx2*vx3));
    vx1=vx2; vx1=vx1*ex3; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,ValidatedFloat(vx2*ex3));
    vx1=vx3; vx1=ex2*vx1; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,ValidatedFloat(ex2*vx3));
    vx1=vx2; vx1=vx1/vx3; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,ValidatedFloat(vx2/vx3));
    vx1=vx2; vx1=vx1/ex3; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,ValidatedFloat(vx2/ex3));
    vx1=vx3; vx1=ex2/vx1; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,ValidatedFloat(ex2/vx3));
}

Void TestValidatedFloat::test_monotone_functions()
{

    ValidatedFloat two(2.0);
    ValidatedFloat sqrttwo=sqrt(two);
    ARIADNE_TEST_PRINT(sqrttwo);
    ARIADNE_TEST_COMPARE(sqrttwo.lower_raw(),<=,1.4142135623730949);
    ARIADNE_TEST_COMPARE(sqrttwo.lower_raw(),> ,1.4142135623730947);
    ARIADNE_TEST_COMPARE(sqrttwo.upper_raw(),>=,1.4142135623730951);
    ARIADNE_TEST_COMPARE(sqrttwo.upper_raw(),< ,1.4142135623730954);

    ValidatedFloat one(1.0);
    ValidatedFloat expone=exp(one);
    ARIADNE_TEST_PRINT(expone);
    ARIADNE_TEST_COMPARE(expone.lower_raw(),<,2.71828182845905);
    ARIADNE_TEST_COMPARE(expone.lower_raw(),>,2.71828182845903);
    ARIADNE_TEST_COMPARE(expone.upper_raw(),>,2.71828182845904);
    ARIADNE_TEST_COMPARE(expone.upper_raw(),<,2.71828182845906);
    ARIADNE_TEST_ASSERT(expone.lower_raw()<expone.upper_raw());

    ValidatedFloat e(2.7182818284590451,2.7182818284590455);
    ValidatedFloat loge=log(e);
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_COMPARE(loge.lower_raw(),<,1);
    ARIADNE_TEST_COMPARE(loge.lower_raw(),>,0.9999999999998);
    ARIADNE_TEST_COMPARE(loge.upper_raw(),>,1);
    ARIADNE_TEST_COMPARE(loge.upper_raw(),<,1.000000000002);
}

Void TestValidatedFloat::test_trigonometric_functions()
{
    try {
        ValidatedFloat x(6.283185307179586,6.283185307179587);
        ValidatedFloat sinx=sin(x);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_COMPARE(sinx.lower_raw(),<,0.0);
        ARIADNE_TEST_COMPARE(sinx.lower_raw(),>,-1e-14);
        ARIADNE_TEST_COMPARE(sinx.upper_raw(),>,0.0);
        ARIADNE_TEST_COMPARE(sinx.upper_raw(),<,+1e-14);
        ARIADNE_TEST_ASSERT(sinx.lower_raw()<sinx.upper_raw());
    }
    catch(...) { }

    try {
        ValidatedFloat x(7.0685834705770345);
        ValidatedFloat sinx=sin(x);
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
        ValidatedFloat x(1.5707963267948966,1.5707963267948968);
        ValidatedFloat cosx=cos(x);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_COMPARE(cosx.lower_raw(),<,0.0);
        ARIADNE_TEST_COMPARE(cosx.lower_raw(),>,-1e-14);
        ARIADNE_TEST_COMPARE(cosx.upper_raw(),>,0.0);
        ARIADNE_TEST_COMPARE(cosx.upper_raw(),<,+1e-14);
        ARIADNE_TEST_ASSERT(cosx.lower_raw()<cosx.upper_raw());
    }

}

Int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    TestValidatedFloat().test();

    return ARIADNE_TEST_FAILURES;
}

