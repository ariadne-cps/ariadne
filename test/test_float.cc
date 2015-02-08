/***************************************************************************
 *            test_float.cc
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
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdexcept>
#include <fenv.h>

#include "config.h"
#include "numeric/rounding.h"
#include "numeric/float.h"
#include "numeric/numeric.h"

#include "test.h"

using namespace std;
using namespace Ariadne;


class TestFloat
{
  public:
    Void test();
  private:
    Void test_concept();
    Void test_class();
    Void test_conversion();
    Void test_stream();
    Void test_comparison();
    Void test_rounding();
    Void test_arithmetic();
    Void test_cosine();
    Void test_function();
};


Int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    TestFloat().test();

    return ARIADNE_TEST_FAILURES;
}


Void
TestFloat::test()
{
    //ARIADNE_TEST_CALL(test_concept());
    ARIADNE_TEST_CALL(test_class());
    ARIADNE_TEST_CALL(test_conversion());
    ARIADNE_TEST_CALL(test_stream());
    ARIADNE_TEST_CALL(test_comparison());
    ARIADNE_TEST_CALL(test_rounding());
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(test_cosine());
    ARIADNE_TEST_CALL(test_function());
}




// Test that the type implements all operations of
// the Float64 concept without testing correctness
Void
TestFloat::test_concept()
{
    Bool b=true;
    Int n=1;
    Nat m=1;
    double d=1;
    Float64 x=1;

    // Constructors
    x=Float64(); x=Float64(n); x=Float64(m); x=Float64(d); x=Float64(x);

    // Assignment
    x=n; x=m; x=d; x=x;

    // Conversion
    d=x.get_d();

    // Maximum and minimum and absolute value
    x=max(x,x); x=min(x,x); x=abs(x);


    // Exact operations
    x=nul(x); x=pos(x); x=neg(x); x=half(x);

    // Rounded arithmetic
    x=add_near(x,x); x=add_approx(x,x); x=add_down(x,x); x=add_up(x,x); // x=add_chop(x,x);
    x=sub_near(x,x); x=add_approx(x,x); x=sub_down(x,x); x=sub_up(x,x); // x=sub_chop(x,x);
    x=mul_near(x,x); x=add_approx(x,x); x=mul_down(x,x); x=mul_up(x,x); // x=mul_chop(x,x);
    x=div_near(x,x); x=add_approx(x,x); x=div_down(x,x); x=div_up(x,x); // x=div_chop(x,x);
    x=fma_approx(x,x,x); x=fma_down(x,x,x); x=fma_up(x,x,x);
    x=pow_approx(x,n); x=pow_down(x,n); x=pow_up(x,n); // x=pow_chop(x,n);
    x=pow_approx(x,m); x=pow_down(x,m); x=pow_up(x,m); // x=pow_chop(x,m);

    x=med_near(x,x); x=rad_up(x,x);

    // Mixed Float64/Int arithmetic
    x=mul_approx(n,x); x=mul_down(n,x); x=mul_up(n,x); // x=mul_chop(n,x);
    x=mul_approx(m,x); x=mul_down(m,x); x=mul_up(m,x); // x=mul_chop(m,x);
    x=mul_approx(x,n); x=mul_down(x,n); x=mul_up(x,n); // x=mul_chop(x,n);
    x=mul_approx(x,m); x=mul_down(x,m); x=mul_up(x,m); // x=mul_chop(x,m);
    x=div_approx(x,n); x=div_down(x,n); x=div_up(x,n); // x=div_chop(x,n);
    x=div_approx(x,m); x=div_down(x,m); x=div_up(x,m); // x=div_chop(x,m);

    // Mixed Float64/double arithmetic
    x=mul_approx(d,x); x=mul_approx(x,d); x=div_approx(x,d);

    // Reset x to 1
    x=1; x=1.0;

    // Operators in rounding mode
    x=+x; x=-x;
    x=x+x; x=x-x; x=x*x; x=x/x;
    x+x; x-=x; x*=x; x/=x;

    // Comparisons
    b=(x==n); b=(x!=n); b=(x<=n); b=(x>=n); b=(x<n); b=(x>n);
    b=(n==x); b=(n!=x); b=(n<=x); b=(n>=x); b=(n<x); b=(n>x);
    b=(x==m); b=(x!=m); b=(x<=m); b=(x>=m); b=(x<m); b=(x>m);
    b=(m==x); b=(m!=x); b=(m<=x); b=(m>=x); b=(m<x); b=(m>x);
    b=(x==d); b=(x!=d); b=(x<=d); b=(x>=d); b=(x<d); b=(x>d);
    b=(d==x); b=(d!=x); b=(d<=x); b=(d>=x); b=(d<x); b=(d>x);
    b=(x==x); b=(x!=x); b=(x<=x); b=(x>=x); b=(x<x); b=(x>x);

    // Rounding mode
    Float64::set_rounding_to_nearest();
    Float64::set_rounding_downward();
    Float64::set_rounding_upward();
    Float64::set_rounding_toward_zero();

    Float64::set_rounding_mode(Float64::to_nearest);
    Float64::set_rounding_mode(Float64::downward);
    Float64::set_rounding_mode(Float64::upward);
    Float64::set_rounding_mode(Float64::toward_zero);

    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64::set_rounding_mode(rnd);

    // Precision64
    Float64::PrecisionType pr=Float64::get_default_precision();
    pr=x.precision();
    x.set_precision(pr);

}


Void
TestFloat::test_class()
{
    cout << __PRETTY_FUNCTION__ << endl;
    // Construct from an Int
    Float64 f1(2);
    ARIADNE_TEST_ASSERT(f1==2);
    // Construct from a double
    Float64 f2(1.25);
    ARIADNE_TEST_ASSERT(f2==1.25);
    // Copy constructor
    Float64 f3(f2);
    ARIADNE_TEST_ASSERT(f3==f2);

    // Assign from an Int
    f1=3;
    ARIADNE_TEST_ASSERT(f1==3);
    // Assign from a double
    f2=2.25;
    ARIADNE_TEST_ASSERT(f2==2.25);
    // Copy assignment
    f3=f2;
    ARIADNE_TEST_ASSERT(f3==f2);

}


Void
TestFloat::test_conversion()
{
    cout << __PRETTY_FUNCTION__ << endl;

    // Convert from integers
    Int n;
    n=std::numeric_limits<Int>::min();
    ARIADNE_TEST_ASSERT(Float64(n)==n);
    n=std::numeric_limits<Int>::max();
    ARIADNE_TEST_ASSERT(Float64(n)==n);
    n=std::numeric_limits<Nat>::max();
    ARIADNE_TEST_ASSERT(Float64(n)==n);

    // Convert to a rational
    ARIADNE_TEST_EQUAL(Rational(Float64(1.0)),Rational(1));
    ARIADNE_TEST_EQUAL(Rational(Float64(2.25)),Rational(9,4));

    // Convert from a rational
    Int num=1; Int den=3;
    Rational q(1,3);
    ARIADNE_TEST_COMPARE(Rational(Float64(q,Float64::downward)),<,q);
    ARIADNE_TEST_COMPARE(Rational(Float64(Rational(num,den),Float64::upward)),>,q);
    ARIADNE_TEST_COMPARE(Float64(Rational(num,den),Float64::downward),<,Float64(Rational(num,den),Float64::upward));

    // Convert from a negative rational
    num=-2; den=5;
    ARIADNE_TEST_COMPARE(Rational(Float64(q,Float64::downward)),<,q);
    ARIADNE_TEST_COMPARE(Rational(Float64(q,Float64::upward)),>,q);
    ARIADNE_TEST_COMPARE(Float64(q,Float64::downward),<,Float64(q,Float64::upward));

};


Void
TestFloat::test_stream()
{
    cout << __PRETTY_FUNCTION__ << endl;

    stringstream ss("1.25 -2.25 42:2 375e1 2.35e1");
    Float64 f1,f2,f3,f4,f5;
    ss >> f1;
    cout << f1 << endl;
    ARIADNE_TEST_ASSERT(f1==1.25);
    ss >> f2;
    cout << f2 << endl;
    ARIADNE_TEST_ASSERT(f2==-2.25);
    ss >> f3;
    cout << f3 << endl;
    ARIADNE_TEST_ASSERT(f3==42.0);
    try {
        ss >> f4;
        cout << f4 << endl;
        if(f4!=23.75) {
            ARIADNE_TEST_WARN("Cannot create float from string literal in exponential form 2.375e1");
        }
    }
    catch(std::exception& e) {
        cerr << e.what() << endl;;
    }

    try {
        ss >> f5;
        cout << f5 << endl;
        if(f4!=23.5) {
            ARIADNE_TEST_WARN("Cannot create float from string literal in exponential form 2.35e1");
        }
    }
    catch(std::exception& e) {
        cerr << e.what() << endl;
    }

}


Void
TestFloat::test_comparison()
{
    cout << __PRETTY_FUNCTION__ << endl;

    Float64 f1(1.25); Float64 f2(-1.25); Float64 f3(-2.25); Float64 f4(1.25);

    // Test comparison of two equal numbers
    ARIADNE_TEST_ASSERT(f1==f4); ARIADNE_TEST_ASSERT(!(f1!=f4));
    ARIADNE_TEST_ASSERT(f1<=f4); ARIADNE_TEST_ASSERT(!(f1> f4));
    ARIADNE_TEST_ASSERT(f1>=f4); ARIADNE_TEST_ASSERT(!(f1< f4));


    // Test comparison of two different numbers
    ARIADNE_TEST_ASSERT(!(f1==f2)); ARIADNE_TEST_ASSERT(f1!=f2);
    ARIADNE_TEST_ASSERT(!(f1<=f2)); ARIADNE_TEST_ASSERT(f1> f2);
    ARIADNE_TEST_ASSERT(f1>=f2); ARIADNE_TEST_ASSERT(!(f1< f2));

    // Test comparison of two negative numbers
    ARIADNE_TEST_ASSERT(!(f2==f3)); ARIADNE_TEST_ASSERT(f2!=f3);
    ARIADNE_TEST_ASSERT(!(f2<=f3)); ARIADNE_TEST_ASSERT(f2> f3);
    ARIADNE_TEST_ASSERT(f2>=f3); ARIADNE_TEST_ASSERT(!(f2< f3));

    // Test comparison with in integer
    Int i2=1;
    ARIADNE_TEST_ASSERT(!(f1==i2)); ARIADNE_TEST_ASSERT(f1!=i2);
    ARIADNE_TEST_ASSERT(!(f1<=i2)); ARIADNE_TEST_ASSERT(f1> i2);
    ARIADNE_TEST_ASSERT(f1>=i2); ARIADNE_TEST_ASSERT(!(f1< i2));

    Int i1=1;
    ARIADNE_TEST_ASSERT(!(i1==f2)); ARIADNE_TEST_ASSERT(i1!=f2);
    ARIADNE_TEST_ASSERT(!(i1<=f2)); ARIADNE_TEST_ASSERT(i1> f2);
    ARIADNE_TEST_ASSERT(i1>=f2); ARIADNE_TEST_ASSERT(!(i1< f2));

    // Test comparison with a double
    double x2=1.0;
    ARIADNE_TEST_ASSERT(!(f1==x2)); ARIADNE_TEST_ASSERT(f1!=x2);
    ARIADNE_TEST_ASSERT(!(f1<=x2)); ARIADNE_TEST_ASSERT(f1> x2);
    ARIADNE_TEST_ASSERT(f1>=x2); ARIADNE_TEST_ASSERT(!(f1< x2));

    double x1=1.0;
    ARIADNE_TEST_ASSERT(!(x1==f2)); ARIADNE_TEST_ASSERT(x1!=f2);
    ARIADNE_TEST_ASSERT(!(x1<=f2)); ARIADNE_TEST_ASSERT(x1> f2);
    ARIADNE_TEST_ASSERT(x1>=f2); ARIADNE_TEST_ASSERT(!(x1< f2));

    // Test comparison with a rational
    //Rational q2=1;
    //ARIADNE_TEST_ASSERT(!(f1==q2)); ARIADNE_TEST_ASSERT(f1!=q2);
    //ARIADNE_TEST_ASSERT(!(f1<=q2)); ARIADNE_TEST_ASSERT(f1> q2);
    //ARIADNE_TEST_ASSERT(f1>=q2); ARIADNE_TEST_ASSERT(!(f1< q2));

    //Rational q1=Rational(-5,4);
    //ARIADNE_TEST_ASSERT(q1==f2)); ARIADNE_TEST_ASSERT(!(q1!=f2));
    //ARIADNE_TEST_ASSERT(q1<=f2)); ARIADNE_TEST_ASSERT(!(q1> f2));
    //ARIADNE_TEST_ASSERT(!(q1>=f2)); ARIADNE_TEST_ASSERT(q1< f2);

}

Void
TestFloat::test_rounding()
{
    volatile double one   = 1;
    volatile double two   = 2;
    volatile double three = 3;
    volatile double five  = 5;
    const double onethirddown    = 0.33333333333333331483;
    const double onethirdup      = 0.33333333333333337034;
    const double onethirdchop    = 0.33333333333333331483;
    const double onethirdnearest = 0.33333333333333331483;
    const double twofifthsdown   = 0.39999999999999996669;
    const double twofifthsup     = 0.40000000000000002220;
    const double twofifthschop   = 0.39999999999999996669;
    const double twofifthsnearest= 0.40000000000000002220;

    Ariadne::set_rounding_mode(Ariadne::downward);
    double onethirdrounddown=one/three;
    ARIADNE_TEST_EQUAL(onethirdrounddown, onethirddown);
    Ariadne::set_rounding_mode(Ariadne::upward);
    double onethirdroundup=one/three;
    ARIADNE_TEST_EQUAL(onethirdroundup, onethirdup);
    Ariadne::set_rounding_mode(Ariadne::toward_zero);
    double onethirdroundchop=one/three;
    ARIADNE_TEST_EQUAL(onethirdroundchop, onethirdchop);
    Ariadne::set_rounding_mode(Ariadne::to_nearest);
    double onethirdroundnearest=one/three;
    ARIADNE_TEST_EQUAL(onethirdroundnearest, onethirdnearest);

    Ariadne::set_rounding_downward();
    double twofifthsrounddown=two/five;
    ARIADNE_TEST_EQUAL(twofifthsrounddown, twofifthsdown);
    Ariadne::set_rounding_upward();
    double twofifthsroundup=two/five;
    ARIADNE_TEST_EQUAL(twofifthsroundup, twofifthsup);
    Ariadne::set_rounding_toward_zero();
    double twofifthsroundchop=two/five;
    ARIADNE_TEST_EQUAL(twofifthsroundchop, twofifthschop);
    Ariadne::set_rounding_to_nearest();
    double twofifthsroundnearest=two/five;
    ARIADNE_TEST_EQUAL(twofifthsroundnearest, twofifthsnearest);
}

Void
TestFloat::test_arithmetic()
{
    cout << __PRETTY_FUNCTION__ << endl;

    static const double eps = std::numeric_limits<double>::epsilon();

    // Set next_up some variables
    Float64 f1(1.25); Float64 f2(2.25); Float64 f3(-3.25); Float64 f4; Float64 f5;

    // Minimum (this should always remain exact)
    f4=min(f1,f2);
    cout << "min(" << f1 << "," << f2 << ") = " << f4 << endl;
    ARIADNE_TEST_ASSERT(f4==f1);
    f4=min(f1,f3);
    cout << "min(" << f1 << "," << f3 << ") = " << f4 << endl;
    ARIADNE_TEST_ASSERT(f4==f3);

    // Maximum (this should always remain exact)
    f4=max(f1,f2);
    cout << "max(" << f1 << "," << f2 << ") = " << f4 << endl;
    ARIADNE_TEST_ASSERT(f4==f2);
    f4=max(f1,f3);
    cout << "max(" << f1 << "," << f3 << ") = " << f4 << endl;
    ARIADNE_TEST_ASSERT(f4==f1);

    // Absolute value (this should always remain exact)
    f4=abs(f1);
    cout << "abs(" << f1 << ") = " << f4 << endl;
    ARIADNE_TEST_ASSERT(f4==1.25);
    f5=abs(f3);
    cout << "abs(" << f3 << ") = " << f5 << endl;
    ARIADNE_TEST_ASSERT(f5==3.25);

    // Median (this should remain exact here)
    f3=med_near(f1,f2);
    cout << f1 << " <= med(" << f1 << "," << f2 << ")=" << f3 << " <= " << f2 << endl;
    ARIADNE_TEST_ASSERT(f1<=f3); ARIADNE_TEST_ASSERT(f3<=f2);
    ARIADNE_TEST_ASSERT(f3==1.75);

    // Negation (this should always remain exact)
    f3=neg(f1);
    cout << "neg(" << f1 << ") = " << f3 << endl;
    f3=-f1;
    cout << "- " << f1 << " = " << f3 << endl;
    ARIADNE_TEST_ASSERT(f3==-1.25);

    // Rounding
    f3=next_down(f1);
    f4=next_up(f1);
    cout << f3 << " < " << f1 << " < " << f4 << endl;
    ARIADNE_TEST_COMPARE(f3,<,f1); ARIADNE_TEST_COMPARE(f4,>,f1);

    // Addition (this should remain exact here)
    f3=add_down(f1,f2);
    f4=add_up(f1,f2);
    cout << f3 << " <= " << f1 << " + " << f2 << " <= " << f4 << endl;
    ARIADNE_TEST_ASSERT(f3<=3.5); ARIADNE_TEST_ASSERT(f4>=3.5);
    // Addition
    ARIADNE_TEST_COMPARE(add_up(Float64(1.0),Float64(eps/2)),>,1.0);
    ARIADNE_TEST_COMPARE(add_down(Float64(1.0),Float64(eps)),>,1.0);
    if(add_down(Float64(1.0),Float64(eps/2)) != Float64(1.0)) {
        ARIADNE_TEST_WARN("Results of floating-point operations stored to higher-precision in registers than memory.");
        ARIADNE_TEST_COMPARE(add_down(Float64(1.0),Float64(eps/(1<<11))),>,Float64(1.0));
    }
    ARIADNE_TEST_COMPARE(add_down(Float64(1.0),Float64(eps/(1<<12))),==,Float64(1.0));
    Float64 one_add_down_half_epsilon=add_down(Float64(1.0),Float64(eps/2));
    ARIADNE_TEST_COMPARE(one_add_down_half_epsilon,==,Float64(1.0));
    ARIADNE_TEST_EQUAL(Float64(1.0),1.0);
    ARIADNE_TEST_COMPARE(Float64(1.0),==,1.0);


    // Subtraction (this should remain exact here)
    f3=sub_down(f1,f2);
    f4=sub_up(f1,f2);
    cout << f3 << " <= " << f1 << " - " << f2 << " <= " << f4 << endl;
    ARIADNE_TEST_ASSERT(f3<=-1); ARIADNE_TEST_ASSERT(f4>=-1);

    // Multiplication (this should remain exact here)
    f3=mul_down(f1,f2);
    f4=mul_up(f1,f2);
    cout << f3 << " <= " << f1 << " * " << f2 << " <= " << f4 << endl;
    ARIADNE_TEST_ASSERT(f3<=2.8125); ARIADNE_TEST_ASSERT(f4>=2.8125);

    // Division (not exact; should catch errors here)
    f3=div_down(f1,f2);
    f4=div_up(f1,f2);
    f5=div_approx(f1,f2);
    cout << f3 << " <= " << f1 << " / " << f2 << " <= " << f4 << endl;
    Float64 five = 5;
    Float64 nine = 9;

    Float64 expected_five_ninths_up=0.55555555555555558023;
    Float64::set_rounding_downward();
    ARIADNE_TEST_COMPARE(five/nine,<,expected_five_ninths_up);
    ARIADNE_TEST_COMPARE(div(five,nine),<,expected_five_ninths_up);
    Float64 five_divby_nine_down=five;
    five_divby_nine_down/=nine;
    ARIADNE_TEST_COMPARE(five_divby_nine_down,<,expected_five_ninths_up);

    Float64::set_rounding_to_nearest();
    ARIADNE_TEST_COMPARE(div_down(Float64(5.0),Float64(9.0)),<,div_up(Float64(5.0),Float64(9.0)));
    ARIADNE_TEST_COMPARE(div_down(Float64(5),Float64(9)),<,div_up(Float64(5),Float64(9)));
    ARIADNE_TEST_COMPARE(div_down(five,nine),<,div_up(five,nine));
    Float64 five_ninths_down = div_down(five,nine);
    Float64 five_ninths_up = div_up(five,nine);
    ARIADNE_TEST_COMPARE(five_ninths_down, < , five_ninths_up);
    ARIADNE_TEST_COMPARE(mul_down(five_ninths_down,nine), < , five);
    ARIADNE_TEST_COMPARE(mul_up(five_ninths_up,nine), > , five);


    // Power (not exact; should catch errors here)
    f3=pow_down(f1,3);
    f4=pow_up(f1,3);
    cout << f3 << " <= pow(" << f1 << ",3) <= " << f4 << endl;
    ARIADNE_TEST_ASSERT(f3<=1.953125); ARIADNE_TEST_ASSERT(f4>=1.953125);

    f3=pow_down(f1,-2);
    f4=pow_up(f1,-2);
    cout << f3 << " <= pow(" << f1 << ",-2) <= " << f4 << endl;
    //ARIADNE_TEST_ASSERT(Rational(f3)<Rational(16,25));
    //ARIADNE_TEST_ASSERT(Rational(f4)>Rational(16,25));

    // Floor and ceiling
    f2=Float64(-3.25); f3=Float64(-2);

    ARIADNE_TEST_ASSERT(floor(f1)==1); ARIADNE_TEST_ASSERT(ceil(f1)==2);
    ARIADNE_TEST_ASSERT(floor(f2)==-4); ARIADNE_TEST_ASSERT(ceil(f2)==-3);
    ARIADNE_TEST_ASSERT(floor(f3)==-2); ARIADNE_TEST_ASSERT(ceil(f3)==-2);

    // Conversion to integer types
    Int i3,i4;
    i3=numeric_cast<Int>(floor(f1));
    i4=numeric_cast<Int>(ceil(f1));
    cout << i3 << " < " << f1 << " < " << i4 << endl;
    ARIADNE_TEST_ASSERT(i3==1); ARIADNE_TEST_ASSERT(i4==2);
    i3=numeric_cast<Int>(floor(f2));
    i4=numeric_cast<Int>(ceil(f2));
    cout << i3 << " < " << f2 << " < " << i4 << endl;
    ARIADNE_TEST_ASSERT(i3==-4); ARIADNE_TEST_ASSERT(i4==-3);

    // Check interval conversions
    Float64 z(0); Float64 o(1); Float64 t(3);

    Float64 odtd=div_down(o,t);
    Float64 odtu=div_up(o,t);
    cout << odtd << " <= 1/3 <= " << odtu << endl;
    ARIADNE_TEST_COMPARE(odtd,<,odtu);
    // Regression test to catch errors when Float64 result is not assigned to a variable
    cout << div_down(o,t) << " <= 1/3 <= " << div_up(o,t) << endl;
    ARIADNE_TEST_COMPARE(div_down(o,t),<,div_up(o,t));

    ARIADNE_TEST_COMPARE(med_near(Float64(2),Float64(3)),==,2.5);
    ARIADNE_TEST_COMPARE(rad_up(Float64(2),Float64(3)),>=,0.5);
    ARIADNE_TEST_COMPARE(rad_up(Float64(2),Float64(3)),<=,0.5000000000000002);

    // The following line should not compile
    // f5=f1+f2;

}


Void
TestFloat::test_function()
{
    cout << __PRETTY_FUNCTION__ << endl;

    cout << setprecision(20);

    // Set up some variables
    Float64 x; Float64 ra; Float64 rl; Float64 ru;

    x=1;
    ra=exp(x);
    rl=next_down(ra);
    ru=next_up(ra);
    ARIADNE_TEST_PRINT(rl);
    ARIADNE_TEST_PRINT(ru);
    ARIADNE_TEST_ASSERT(rl<ru);
    ARIADNE_TEST_ASSERT(2.71<rl);
    ARIADNE_TEST_ASSERT(ru<2.72);


    //The following don't work as rounded operators not exported.
    //test_inverse_pair("sin",&sin_down,&sin_up,&asin_down,&asin_up);
    //The following don't work as acos is decreasing
    //test_inverse_pair("cos",&cos_down,&cos_up,&acos_down,&acos_up);
    //test_inverse_pair("cos",&cos_up,&cos_down,&acos_down,&acos_up);
    //test_inverse_pair("tan",&tan_down,&tan_up,&atan_down,&atan_up);
    //test_inverse_pair("sinh",&sinh_down,&sinh_up,&asinh_down,&asinh_up);
    //test_inverse_pair("cosh",&cosh_down,&cosh_up,&acosh_down,&acosh_up);
    //test_inverse_pair("tanh",&tanh_down,&tanh_up,&atanh_down,&atanh_up);
}

Void
TestFloat::test_cosine()

{
    //3.14159265358979323846264338327950288419716939937510
    static const Float64 pi_down=3.1415926535897931;
    //static const double pi_approx=3.1415926535897931;
    static const Float64 pi_up=3.1415926535897936;

    static const Float64 third_pi_down=1.0471975511965976;
    static const Float64 third_pi_up=1.0471975511965979;

    Float64::set_rounding_mode(Float64::upward);
    ARIADNE_TEST_EQUAL(cos(Float64(0.0)),1.0);
    ARIADNE_TEST_COMPARE(cos(third_pi_down),>,0.5);
    ARIADNE_TEST_COMPARE(sqr(cos(pi_down/4)),>,0.5);
    ARIADNE_TEST_COMPARE(cos(pi_down/2),>,0.0);
    ARIADNE_TEST_COMPARE(cos(pi_up/2),<=,0.0);
    ARIADNE_TEST_COMPARE(cos(pi_down),>,-1.0);
    ARIADNE_TEST_COMPARE(cos(pi_up),>,-1.0);
    ARIADNE_TEST_COMPARE(cos(2*pi_down),==,1.0);
    ARIADNE_TEST_COMPARE(cos(2*pi_up),==,1.0);
    ARIADNE_TEST_COMPARE(cos(3*pi_down),>,-1.0);
    ARIADNE_TEST_COMPARE(cos(3*pi_up),>,-1.0);

    Float64::set_rounding_mode(Float64::downward);
    ARIADNE_TEST_EQUAL(cos(Float64(0.0)),1.0);
    ARIADNE_TEST_COMPARE(cos(third_pi_up),<,0.5);
    ARIADNE_TEST_COMPARE(sqr(cos(pi_up/4)),<,0.5);
    ARIADNE_TEST_COMPARE(cos(pi_down/2),>=,0.0);
    ARIADNE_TEST_COMPARE(cos(pi_up/2),<,0.0);
    ARIADNE_TEST_COMPARE(cos(pi_down),==,-1.0);
    ARIADNE_TEST_COMPARE(cos(pi_up),==,-1.0);
    ARIADNE_TEST_COMPARE(cos(2*pi_down),<,1.0);
    ARIADNE_TEST_COMPARE(cos(2*pi_up),<,1.0);
    ARIADNE_TEST_COMPARE(cos(3*pi_down),==,-1.0);
    ARIADNE_TEST_COMPARE(cos(3*pi_up),==,-1.0);

    static const Float64 one=1.0;
    Float64::set_rounding_mode(Float64::upward);
    Float64 sin_rnd_up_one=sin(one);
    Float64::set_rounding_mode(Float64::downward);
    Float64 sin_rnd_down_one=sin(one);
    Float64::set_rounding_mode(Float64::to_nearest);
    Float64 sin_rnd_approx_one=sin(one);
    ARIADNE_TEST_COMPARE(sin_rnd_down_one,<=,sin_rnd_approx_one);
    ARIADNE_TEST_COMPARE(sin_rnd_approx_one,<=,sin_rnd_up_one);
    ARIADNE_TEST_COMPARE(sin_rnd_down_one,< ,sin_rnd_up_one);
}
