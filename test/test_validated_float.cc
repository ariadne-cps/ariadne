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
#include "rational.h"
#include "float-validated.h"

#include "test.h"

using namespace Ariadne;
using namespace std;


class TestValidatedFloat
{
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
    void regression_tests();
};


void
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

void
TestValidatedFloat::test_concept()
{
    ValidatedFloat::set_output_precision(17);

    int n=1;
    uint m=1;
    double d=1;
    ExactFloat x=1;
    Float a,b;
    ValidatedFloat i=1;
    ValidatedFloat j=1;


    // Constructors
    j=ValidatedFloat(); j=ValidatedFloat(n); j=ValidatedFloat(m); j=ValidatedFloat(d); j=ValidatedFloat(x); j=ValidatedFloat(i);
    j=ValidatedFloat(n,n); j=ValidatedFloat(m,m); j=ValidatedFloat(d,d); j=ValidatedFloat(a,b);
    j=ValidatedFloat(n,m); j=ValidatedFloat(m,d); j=ValidatedFloat(d,n);
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
void
TestValidatedFloat::test_correct_rounded_arithmetic()
{
    ValidatedFloat onethird=ValidatedFloat(1)/ValidatedFloat(3);
    ARIADNE_TEST_COMPARE( onethird.lower_value(), < , onethird.upper_value() );
    ValidatedFloat one_approx=onethird*ValidatedFloat(3);
    ARIADNE_TEST_COMPARE( one_approx.lower_value(), < , 1.0 );
    ARIADNE_TEST_COMPARE( one_approx.upper_value(), > , 1.0 );
}


// Test that interval arithmetic gives the most accurate rounded values
void
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
void
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



void
TestValidatedFloat::test_constructors()
{
    Float zero=0;

    // Construct from pair
    ValidatedFloat ivld1(Float(1.125),Float(2.25));
    ARIADNE_TEST_ASSERT(ivld1.lower_value()==1.125); ARIADNE_TEST_ASSERT(ivld1.upper_value()==2.25);

    // Default constructor
    ValidatedFloat ivld2;
    if(ivld2.lower_value()>ivld2.upper_value()) {
        ARIADNE_TEST_WARN("ValidatedFloat default constructor returns an empty set.");
    } else {
        ARIADNE_TEST_ASSERT((bool)(ivld2==ValidatedFloat(zero,zero)));
    }

    // Constructor with approximations
#ifdef HAVE_GMPXX_H
    ValidatedFloat ivld3(Rational(21,10),Rational(16,5));
    cout<<ivld3<<std::endl;
    ARIADNE_TEST_COMPARE(make_exact(ivld3.lower_value()),<,Rational(21,10));
    ARIADNE_TEST_COMPARE(make_exact(ivld3.upper_value()),>,Rational(16,5));
#else
    ValidatedFloat ivld3(2.1,3.2);
#endif // HAVE_GMPXX_H

    // Constructor from approximate values
    ValidatedFloat ivld4(2.1,3.2);
    ARIADNE_TEST_COMPARE(ivld4.lower_value(),<=,2.1);
    ARIADNE_TEST_COMPARE(ivld4.upper_value(),>=,3.2);

#ifdef HAVE_GMPXX_H
    // Approximate constructor from a single value
    ValidatedFloat ivld5(Rational(1,3));
    ARIADNE_TEST_COMPARE(make_exact(ivld5.lower_value()),<,Rational(1,3));
    ARIADNE_TEST_COMPARE(make_exact(ivld5.upper_value()),>,Rational(1,3));
#else
    ValidatedFloat ivld5(1./3.);
#endif // HAVE_GMPXX_H

    // Exact constructor from a single value
    ValidatedFloat ivld6(Float(1.25));
    ARIADNE_TEST_EQUAL(ivld6.lower_value(),Float(1.25));
    ARIADNE_TEST_EQUAL(ivld6.upper_value(),Float(1.25));

    // Empty interval
    ValidatedFloat ivld7;
    ARIADNE_TEST_EXECUTE(ivld7.set_empty());
    ARIADNE_TEST_ASSERT(ivld7.lower_value()==+inf); ARIADNE_TEST_ASSERT(ivld7.upper_value()==-inf);
}

void TestValidatedFloat::test_class()
{
    // Test lower, upper, midpoint, radius, width

    // Tests for exact operations
    ARIADNE_TEST_EQUAL(ValidatedFloat(-0.25,0.50).lower().raw(),-0.25);
    ARIADNE_TEST_EQUAL(ValidatedFloat(-0.25,0.50).upper().raw(),0.5);
    ARIADNE_TEST_EQUAL(ValidatedFloat(-0.25,0.50).midpoint().raw(),0.125);
    ARIADNE_TEST_EQUAL(ValidatedFloat(-0.25,0.50).radius().raw(),0.375)
    ARIADNE_TEST_EQUAL(ValidatedFloat(-0.25,0.50).width().raw(),0.75);

    // Tests for inexact operations
    ARIADNE_TEST_EQUAL(ValidatedFloat(-1./3,2./3).lower().raw(),-0.33333333333333331483);
    ARIADNE_TEST_EQUAL(ValidatedFloat(-1./3,2./3).upper().raw(),0.66666666666666662966);
    ARIADNE_TEST_EQUAL(ValidatedFloat(-1./3,2./3).midpoint().raw(),0.16666666666666665741);
    ARIADNE_TEST_EQUAL(ValidatedFloat(-1./3,2./3).radius().raw(),0.5)
    ARIADNE_TEST_EQUAL(ValidatedFloat(-1./3,2./3).width().raw(),1.0);

    // Tests for inexact operations
    ARIADNE_TEST_EQUAL(ValidatedFloat(div_down(-1,3),div_up(2,3)).lower().raw(),-0.33333333333333337034);
    ARIADNE_TEST_EQUAL(ValidatedFloat(div_down(-1,3),div_up(2,3)).upper().raw(),0.66666666666666674068);
    ARIADNE_TEST_EQUAL(ValidatedFloat(div_down(-1,3),div_up(2,3)).midpoint().raw(),0.16666666666666668517);
    ARIADNE_TEST_EQUAL(ValidatedFloat(div_down(-1,3),div_up(2,3)).radius().raw(),0.50000000000000011102)
    ARIADNE_TEST_EQUAL(ValidatedFloat(div_down(-1,3),div_up(2,3)).width().raw(),1.000000000000000222);
}

void TestValidatedFloat::test_input()
{
    ValidatedFloat ivl1,ivl2;
    string input("[1.125,2.25] [0.4,0.6]");
    stringstream iss(input);

    iss >> ivl1;
    ivl2=ValidatedFloat(1.125,2.25);

    cout << "ivl1=" << ivl1 << "  ivl2=" << ivl2 << endl;
    ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,ivl2);
    ARIADNE_TEST_ASSERT(ivl1.lower_value()==ivl2.lower_value() && ivl1.upper_value()==ivl2.upper_value());
    ARIADNE_TEST_ASSERT(equal(ivl1,ivl2));

    iss >> ivl1;
    ivl2=ValidatedFloat(0.39999999999999997,0.60000000000000009);
    if(!equal(ivl1,ivl2)) {
        ARIADNE_TEST_WARN("ValidatedFloat string constructor returns an approximate interval, not an outwardly rounded interval.");
    }
}

void TestValidatedFloat::test_comparison() {
    // FIXME: If using Boost style interval tests, uncomment the line below
    // and comment out the line after
    //ARIADNE_TEST_ASSERT(indeterminate(ivld1==ivld2));
    ValidatedFloat ivl1(1.125,2.25);
    ValidatedFloat ivl2=ivl1;

    ARIADNE_TEST_ASSERT(ivl1==ivl2);
    ValidatedFloat& ivl1ref=ivl1;
    ivl1ref=ValidatedFloat(5.25,7.375);
    cout << "ivl1ref=" << ivl1ref << endl;
    ARIADNE_TEST_ASSERT(ivl1ref.lower_value()==Float(5.25));
}

void TestValidatedFloat::test_aliasing() {

    ExactFloat x2(1.5);
    ExactFloat x3(2.25);

    ValidatedFloat ivl1;
    ValidatedFloat ivl2(1.5,2.25);
    ValidatedFloat ivl3(3.125,4.0625);

    // Check to make sure aliases are handled correctly
    ivl1=ivl3; ivl1=ivl2-ivl1; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,ValidatedFloat(ivl2-ivl3));
    ivl1=ivl3; ivl1=ivl2*ivl1; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,ValidatedFloat(ivl2*ivl3));
    ivl1=ivl2; ivl1=ivl1*ivl3; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,ValidatedFloat(ivl2*ivl3));
    ivl1=ivl2; ivl1=ivl1*x3; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,ValidatedFloat(ivl2*x3));
    ivl1=ivl3; ivl1=x2*ivl1; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,ValidatedFloat(x2*ivl3));
    ivl1=ivl2; ivl1=ivl1/ivl3; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,ValidatedFloat(ivl2/ivl3));
    ivl1=ivl2; ivl1=ivl1/x3; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,ValidatedFloat(ivl2/x3));
    ivl1=ivl3; ivl1=x2/ivl1; ARIADNE_TEST_BINARY_PREDICATE(equal,ivl1,ValidatedFloat(x2/ivl3));
}

void TestValidatedFloat::test_monotone_functions()
{

    ValidatedFloat two(2.0);
    ValidatedFloat sqrttwo=sqrt(two);
    ARIADNE_TEST_PRINT(sqrttwo);
    ARIADNE_TEST_COMPARE(sqrttwo.lower_value(),<=,1.4142135623730949);
    ARIADNE_TEST_COMPARE(sqrttwo.lower_value(),> ,1.4142135623730947);
    ARIADNE_TEST_COMPARE(sqrttwo.upper_value(),>=,1.4142135623730951);
    ARIADNE_TEST_COMPARE(sqrttwo.upper_value(),< ,1.4142135623730954);

    ValidatedFloat one(1.0);
    ValidatedFloat expone=exp(one);
    ARIADNE_TEST_PRINT(expone);
    ARIADNE_TEST_COMPARE(expone.lower_value(),<,2.71828182845905);
    ARIADNE_TEST_COMPARE(expone.lower_value(),>,2.71828182845903);
    ARIADNE_TEST_COMPARE(expone.upper_value(),>,2.71828182845904);
    ARIADNE_TEST_COMPARE(expone.upper_value(),<,2.71828182845906);
    ARIADNE_TEST_ASSERT(expone.lower_value()<expone.upper_value());

    ValidatedFloat e(2.7182818284590451,2.7182818284590455);
    ValidatedFloat loge=log(e);
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_COMPARE(loge.lower_value(),<,1);
    ARIADNE_TEST_COMPARE(loge.lower_value(),>,0.9999999999998);
    ARIADNE_TEST_COMPARE(loge.upper_value(),>,1);
    ARIADNE_TEST_COMPARE(loge.upper_value(),<,1.000000000002);
}

void TestValidatedFloat::test_trigonometric_functions()
{
    try {
        ValidatedFloat x(6.283185307179586,6.283185307179587);
        ValidatedFloat sinx=sin(x);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_COMPARE(sinx.lower_value(),<,0.0);
        ARIADNE_TEST_COMPARE(sinx.lower_value(),>,-1e-14);
        ARIADNE_TEST_COMPARE(sinx.upper_value(),>,0.0);
        ARIADNE_TEST_COMPARE(sinx.upper_value(),<,+1e-14);
        ARIADNE_TEST_ASSERT(sinx.lower_value()<sinx.upper_value());
    }
    catch(...) { }

    try {
        ValidatedFloat x(7.0685834705770345);
        ValidatedFloat sinx=sin(x);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_COMPARE(sinx.lower_value(),<,0.7071067811866);
        ARIADNE_TEST_COMPARE(sinx.upper_value(),>,0.7071067811865);
        ARIADNE_TEST_ASSERT(sinx.lower_value()<sinx.upper_value());
    }
    catch(...) { }

}

void TestValidatedFloat::regression_tests() {

    // Regression test; fails dramatically on certain types of rounding
    {
        ValidatedFloat x(1.5707963267948966,1.5707963267948968);
        ValidatedFloat cosx=cos(x);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_COMPARE(cosx.lower_value(),<,0.0);
        ARIADNE_TEST_COMPARE(cosx.lower_value(),>,-1e-14);
        ARIADNE_TEST_COMPARE(cosx.upper_value(),>,0.0);
        ARIADNE_TEST_COMPARE(cosx.upper_value(),<,+1e-14);
        ARIADNE_TEST_ASSERT(cosx.lower_value()<cosx.upper_value());
    }

}

int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    TestValidatedFloat().test();

    return ARIADNE_TEST_FAILURES;
}

