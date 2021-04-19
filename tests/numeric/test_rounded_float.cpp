/***************************************************************************
 *            test_rounded_float.cpp
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
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdexcept>
#include <fenv.h>

#include "config.hpp"
#include "numeric/rounding.hpp"
#include "numeric/float.hpp"
#include "numeric/floatdp.hpp"
#include "numeric/floatmp.hpp"
#include "numeric/rounded_float.hpp"
#include "numeric/numeric.hpp"
#include "numeric/dyadic.hpp"

#include "../test.hpp"

using namespace std;
using namespace Ariadne;

template<class FLT>
class TestRounded
{
    using FloatType = FLT;
    using RoundedFloatType = Rounded<FLT>;
    using PR=typename FLT::PrecisionType;
    PR precision;

  public:
    TestRounded(PR prec);
    Void test();
  private:
    Void test_concept();
    Void test_class();
    Void test_conversion_from_to();
    Void test_comparison();
    Void test_arithmetic();
    Void test_cosine();
    Void test_arctan();
    Void test_function();
};


Int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    return ARIADNE_TEST_FAILURES;
}


template<class FLT>
TestRounded<FLT>::TestRounded(PR pr)
    : precision(pr)
{
}

template<class FLT> Void
TestRounded<FLT>::test()
{
    ARIADNE_TEST_CALL(test_class());
    ARIADNE_TEST_CALL(test_conversion_from_to());
    ARIADNE_TEST_CALL(test_comparison());
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(test_function());
    ARIADNE_TEST_CALL(test_cosine());
    ARIADNE_TEST_CALL(test_arctan());
}


// Test that the type implements all operations of
// the RoundedFloatType concept without testing correctness
template<class FLT> Void
TestRounded<FLT>::test_concept()
{
    Bool b=true;
    Int n=1;
    Nat m=1;
    double d=1;
    RoundedFloatType x=1, x2=1;

    // Constructors
    x=RoundedFloatType(); x=RoundedFloatType(n); x=RoundedFloatType(m); x=RoundedFloatType(d); x=RoundedFloatType(x);

    // Assignment
    x=n; x=m; x=d; x=x2;

    // Conversion
    d=x.get_d();

    // Maximum and minimum and absolute value
    x=max(x,x); x=min(x,x); x=abs(x);


    // ExactTag operations
    x=nul(x); x=pos(x); x=neg(x); x=hlf(x);

    // Rounded arithmetic operations
    x=add(x,x);
    x=sub(x,x);
    x=mul(x,x);
    x=div(x,x);
    x=fma(x,x,x);
    x=pow(x,n);
    x=pow(x,m);

    // Non-exact operations
    x=sqr(x);
    x=rec(x);
    x=sqrt(x);
    x=exp(x);
    x=log(x);
    x=sin(x);
    x=cos(x);
    x=tan(x);
    x=asin(x);
    x=acos(x);
    x=atan(x);

    x=med(x,x); x=rad(x,x);

    // Mixed RoundedFloatType/Int arithmetic
    x=mul(n,x);
    x=mul(m,x);
    x=mul(x,n);
    x=mul(x,m);
    x=div(x,n);
    x=div(x,m);

    // Mixed RoundedFloatType/double arithmetic
    x=mul(d,x);

    // Reset x to zero
    x=0; x=0.0;

    // Reset x to 1
    x=1; x=1.0;

    // Operators in rounding mode
    x=+x; x=-x;
    x=x+x; x=x-x; x=x*x; x=x/x;
    x+x; x-=x2; x*=x2; x/=x2;

    // Comparisons
    b=(x==n); b=(x!=n); b=(x<=n); b=(x>=n); b=(x<n); b=(x>n);
    b=(n==x); b=(n!=x); b=(n<=x); b=(n>=x); b=(n<x); b=(n>x);
    b=(x==m); b=(x!=m); b=(x<=m); b=(x>=m); b=(x<m); b=(x>m);
    b=(m==x); b=(m!=x); b=(m<=x); b=(m>=x); b=(m<x); b=(m>x);
    b=(x==d); b=(x!=d); b=(x<=d); b=(x>=d); b=(x<d); b=(x>d);
    b=(d==x); b=(d!=x); b=(d<=x); b=(d>=x); b=(d<x); b=(d>x);
    b=(x==x); b=(x!=x); b=(x<=x); b=(x>=x); b=(x<x); b=(x>x);

    // Rounded mode
    RoundedFloatType::set_rounding_to_nearest();
    RoundedFloatType::set_rounding_downward();
    RoundedFloatType::set_rounding_upward();
    RoundedFloatType::set_rounding_toward_zero();

    RoundedFloatType::set_rounding_mode(to_nearest);
    RoundedFloatType::set_rounding_mode(downward);
    RoundedFloatType::set_rounding_mode(upward);
    RoundedFloatType::set_rounding_mode(toward_zero);

    RoundedFloatType::set_rounding_mode(near);
    RoundedFloatType::set_rounding_mode(down);
    RoundedFloatType::set_rounding_mode(up);

    typename RoundedFloatType::RoundedModeType rnd=RoundedFloatType::get_rounding_mode();
    RoundedFloatType::set_rounding_mode(rnd);

    // DoublePrecision
    typename RoundedFloatType::PrecisionType pr=RoundedFloatType::get_default_precision();
    pr=x.precision();
    x.set_precision(pr);

}


template<class FLT> Void
TestRounded<FLT>::test_class()
{
    cout << __PRETTY_FUNCTION__ << endl;
    // Construct from an Int
    ARIADNE_TEST_CONSTRUCT(RoundedFloatType,f1,(2));
    ARIADNE_TEST_EQUALS(f1,2);
    // Construct from a double
    ARIADNE_TEST_CONSTRUCT(RoundedFloatType,f2,(1.25));
    ARIADNE_TEST_EQUALS(f2,1.25_x);
    // Copy constructor
    ARIADNE_TEST_CONSTRUCT(RoundedFloatType,f3,(f2));
    ARIADNE_TEST_EQUAL(f3,f2);


    // Assign from an Int
    ARIADNE_TEST_EXECUTE(f1=3);
    ARIADNE_TEST_EQUALS(f1,3);
    // Assign from a double
    ARIADNE_TEST_EXECUTE(f2=FloatType(2.25,precision));
    ARIADNE_TEST_EQUALS(f2,2.25_x);
    // Copy assignment
    ARIADNE_TEST_EXECUTE(f3=f2);
    ARIADNE_TEST_EQUAL(f3,f2);
    // Self-assignment
    ARIADNE_TEST_EXECUTE(f3=f3);
    ARIADNE_TEST_EQUAL(f3,f2);

    if(not (precision==FloatType::get_default_precision())) {
        PR low_precision=min(precision,FloatType::get_default_precision());
        PR high_precision=max(precision,FloatType::get_default_precision());
        ARIADNE_TEST_CONSTRUCT(RoundedFloatType,f4,(2.25,high_precision));
        ARIADNE_TEST_CONSTRUCT(RoundedFloatType,f5,(3.75,low_precision));
        ARIADNE_TEST_EXECUTE(f5=f4);
        ARIADNE_TEST_EQUAL(f5.precision(),high_precision);
        ARIADNE_TEST_EQUAL(f5,f4);
        ARIADNE_TEST_CONSTRUCT(RoundedFloatType,f6,(low_precision));
        ARIADNE_TEST_EXECUTE(f6=std::move(f4));
        ARIADNE_TEST_EQUAL(f6.precision(),high_precision);
        ARIADNE_TEST_EQUAL(f6,f5);
    }


}



template<class FLT> Void
TestRounded<FLT>::test_conversion_from_to()
{
    // Convert from dyadic
    RoundedFloatType f1(Dyadic(5,2u),precision);
    ARIADNE_TEST_EQUAL(f1,1.25_x);

    // Convert from self
    RoundedFloatType f2(f1,precision);
    ARIADNE_TEST_EQUAL(f1,f2);
    RoundedFloatType f3(f1,precision);
    ARIADNE_TEST_EQUAL(f1,f3);

    // Convert from integers
    Int n;
    n=std::numeric_limits<Int>::min();
    ARIADNE_TEST_ASSERT(RoundedFloatType(n)==n);
    n=std::numeric_limits<Int>::max();
    ARIADNE_TEST_ASSERT(RoundedFloatType(n)==n);
    Nat nn = std::numeric_limits<Nat>::max();
    ARIADNE_TEST_ASSERT(RoundedFloatType(nn)==nn);

    // Convert to a rational
    ARIADNE_TEST_EQUAL(Rational(RoundedFloatType(1.0).raw()),Rational(1));
    ARIADNE_TEST_EQUAL(Rational(RoundedFloatType(2.25).raw()),Rational(9,4));

    // Convert from a rational
    Int num=1; Int den=3;
    Rational q(1,3);
    RoundedFloatType::set_rounding_downward();
    ARIADNE_TEST_COMPARE(Rational(RoundedFloatType(q,precision).raw()),<,q);
    RoundedFloatType::set_rounding_upward();
    ARIADNE_TEST_COMPARE(Rational(RoundedFloatType(Rational(num,den),precision).raw()),>,q);

    // Convert from a negative rational
    num=-2; den=5;
    RoundedFloatType::set_rounding_downward();
    ARIADNE_TEST_COMPARE(Rational(RoundedFloatType(q,precision).raw()),<,q);
    RoundedFloatType::set_rounding_upward();
    ARIADNE_TEST_COMPARE(Rational(RoundedFloatType(q,precision).raw()),>,q);
}

template<class FLT> Void
TestRounded<FLT>::test_comparison()
{
    cout << __PRETTY_FUNCTION__ << endl;

    ARIADNE_TEST_CONSTRUCT(RoundedFloatType,f0,(3.00));
    ARIADNE_TEST_CONSTRUCT(RoundedFloatType,f1,(1.25));
    ARIADNE_TEST_CONSTRUCT(RoundedFloatType,f2,(-1.25));
    ARIADNE_TEST_CONSTRUCT(RoundedFloatType,f3,(-2.25));
    ARIADNE_TEST_CONSTRUCT(RoundedFloatType,f4,(1.25));

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

    // Test comparison with in int
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

    // Test comparison with an integer
    ARIADNE_TEST_CONSTRUCT(Integer,z0,(-3));
    ARIADNE_TEST_ASSERT(!(f1==z0)); ARIADNE_TEST_ASSERT(f1!=z0);
    ARIADNE_TEST_ASSERT(!(f1<=z0)); ARIADNE_TEST_ASSERT(f1> z0);
    ARIADNE_TEST_ASSERT(f1>=z0); ARIADNE_TEST_ASSERT(!(f1< z0));
    ARIADNE_TEST_ASSERT(!(z0==f1)); ARIADNE_TEST_ASSERT(z0<f1);

    // Test comparison with a dyadic
    ARIADNE_TEST_CONSTRUCT(Dyadic,w2,(-5,2u));
    ARIADNE_TEST_ASSERT(!(f1==w2)); ARIADNE_TEST_ASSERT(f1!=w2);
    ARIADNE_TEST_ASSERT(!(w2==f1)); ARIADNE_TEST_ASSERT(w2<f1);
    ARIADNE_TEST_ASSERT(!(f1<=w2)); ARIADNE_TEST_ASSERT(f1> w2);
    ARIADNE_TEST_ASSERT(f1>=w2); ARIADNE_TEST_ASSERT(!(f1< w2));
    ARIADNE_TEST_ASSERT(f2==w2); ARIADNE_TEST_ASSERT(!(f2!=w2));

    // Test comparison with a rational
    ARIADNE_TEST_CONSTRUCT(Rational,q2,(-5,4));
    ARIADNE_TEST_ASSERT(!(f1==q2)); ARIADNE_TEST_ASSERT(f1!=q2);
    ARIADNE_TEST_ASSERT(!(f1<=q2)); ARIADNE_TEST_ASSERT(f1> q2);
    ARIADNE_TEST_ASSERT(f1>=q2); ARIADNE_TEST_ASSERT(!(f1< q2));
    ARIADNE_TEST_ASSERT(f2==q2); ARIADNE_TEST_ASSERT(!(f2!=q2));
}

template<class FLT> Void
TestRounded<FLT>::test_arithmetic()
{
    cout << __PRETTY_FUNCTION__ << endl;

    static bool full_precision_warning = false;
    RoundedFloatType::set_rounding_to_nearest();
    const RoundedFloatType pi_near=RoundedFloatType::pi(precision);
    const RoundedFloatType four   =RoundedFloatType(4,precision);
    RoundedFloatType::set_rounding_upward();
    if( pi_near * 4.0_x != mul(pi_near,four) ) {
        if(not full_precision_warning) {
            full_precision_warning=true;
            ARIADNE_TEST_WARN("Mixed RoundedFloatType/BuiltinTag operations may not use full precision.");
        }
    }

    static const double eps = std::numeric_limits<double>::epsilon();

    // Set next_up some variables
    RoundedFloatType f1(1.25); RoundedFloatType f2(2.25); RoundedFloatType f3(-3.25); RoundedFloatType f4; RoundedFloatType f5;

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
    ARIADNE_TEST_ASSERT(f4==1.25_x);
    f5=abs(f3);
    cout << "abs(" << f3 << ") = " << f5 << endl;
    ARIADNE_TEST_ASSERT(f5==3.25_x);

    // Median (this should remain exact here)
    f3=med(f1,f2);
    cout << f1 << " <= med(" << f1 << "," << f2 << ")=" << f3 << " <= " << f2 << endl;
    ARIADNE_TEST_ASSERT(f1<=f3); ARIADNE_TEST_ASSERT(f3<=f2);
    ARIADNE_TEST_ASSERT(f3==1.75_x);

    // Negation (this should always remain exact)
    f3=neg(f1);
    cout << "neg(" << f1 << ") = " << f3 << endl;
    f3=-f1;
    cout << "- " << f1 << " = " << f3 << endl;
    ARIADNE_TEST_ASSERT(f3==-1.25_x);

    // Addition (this should remain exact here)
    f3=add(f1,f2);
    f4=add(f1,f2);
    cout << f3 << " <= " << f1 << " + " << f2 << " <= " << f4 << endl;
    ARIADNE_TEST_ASSERT(f3<=3.5_x); ARIADNE_TEST_ASSERT(f4>=3.5_x);
    // Addition
    ARIADNE_TEST_EXECUTE(RoundedFloatType::set_rounding_upward());
    ARIADNE_TEST_COMPARE(add(RoundedFloatType(1.0),RoundedFloatType(eps/2)),>,1.0);
    ARIADNE_TEST_EXECUTE(RoundedFloatType::set_rounding_downward());
    ARIADNE_TEST_COMPARE(add(RoundedFloatType(1.0),RoundedFloatType(eps)),>,1.0);
    if(add(RoundedFloatType(1.0),RoundedFloatType(eps/2)) != RoundedFloatType(1.0)) {
        ARIADNE_TEST_WARN("Results of floating-point operations stored to higher-precision in registers than memory.");
        ARIADNE_TEST_COMPARE(add(RoundedFloatType(1.0),RoundedFloatType(eps/(1<<11))),>,RoundedFloatType(1.0));
    }
    ARIADNE_TEST_COMPARE(add(RoundedFloatType(1.0),RoundedFloatType(eps/(1<<12))),==,RoundedFloatType(1.0));
    ARIADNE_TEST_CONSTRUCT(RoundedFloatType,one_add_down_half_epsilon,(add(RoundedFloatType(1.0),RoundedFloatType(eps/2))));
    ARIADNE_TEST_COMPARE(one_add_down_half_epsilon,==,RoundedFloatType(1.0));
    ARIADNE_TEST_EQUAL(RoundedFloatType(1.0),1.0);
    ARIADNE_TEST_COMPARE(RoundedFloatType(1.0),==,1.0);


    // Subtraction (this should remain exact here)
    RoundedFloatType::set_rounding_downward();
    f3=sub(f1,f2);
    RoundedFloatType::set_rounding_upward();
    f4=sub(f1,f2);
    cout << f3 << " <= " << f1 << " - " << f2 << " <= " << f4 << endl;
    ARIADNE_TEST_ASSERT(f3<=-1); ARIADNE_TEST_ASSERT(f4>=-1);

    // Multiplication (this should remain exact here)
    RoundedFloatType::set_rounding_downward();
    f3=mul(f1,f2);
    RoundedFloatType::set_rounding_upward();
    f4=mul(f1,f2);
    cout << f3 << " <= " << f1 << " * " << f2 << " <= " << f4 << endl;
    ARIADNE_TEST_ASSERT(f3<=2.8125_x); ARIADNE_TEST_ASSERT(f4>=2.8125_x);

    // Division (not exact; should catch errors here)
    RoundedFloatType::set_rounding_downward();
    f3=div(f1,f2);
    RoundedFloatType::set_rounding_upward();
    f4=div(f1,f2);
    RoundedFloatType::set_rounding_to_nearest();
    f5=div(f1,f2);
    cout << f3 << " <= " << f1 << " / " << f2 << " <= " << f4 << endl;
    RoundedFloatType five(5,precision);
    RoundedFloatType nine(9,precision);

    if (Same<FloatType,FloatDP>) {
        RoundedFloatType expected_five_ninths_up(0.55555555555555558023);
        RoundedFloatType::set_rounding_downward();
        ARIADNE_TEST_COMPARE(five/nine,<,expected_five_ninths_up);
        ARIADNE_TEST_COMPARE(div(five,nine),<,expected_five_ninths_up);
        RoundedFloatType five_divby_nine_down=five;
        five_divby_nine_down/=nine;
        ARIADNE_TEST_COMPARE(five_divby_nine_down,<,expected_five_ninths_up);
    }

    RoundedFloatType::set_rounding_downward();
    RoundedFloatType five_ninths_down = div(five,nine);
    RoundedFloatType five_down = mul(five_ninths_down,nine);
    RoundedFloatType::set_rounding_upward();
    RoundedFloatType five_ninths_up = div(five,nine);
    RoundedFloatType five_up = mul(five_ninths_up,nine);
    ARIADNE_TEST_COMPARE(five_ninths_down, < , five_ninths_up);
    ARIADNE_TEST_COMPARE(five_down, < , five);
    ARIADNE_TEST_COMPARE(five_up, > , five);


    // Power (not exact; should catch errors here)
    RoundedFloatType::set_rounding_downward();
    f3=pow(f1,3);
    RoundedFloatType::set_rounding_upward();
    f4=pow(f1,3);
    cout << f3 << " <= pow(" << f1 << ",3) <= " << f4 << endl;
    ARIADNE_TEST_ASSERT(f3<=1.953125_x); ARIADNE_TEST_ASSERT(f4>=1.953125_x);

    RoundedFloatType::set_rounding_downward();
    f3=pow(f1,-2);
    RoundedFloatType::set_rounding_upward();
    f4=pow(f1,-2);
    cout << f3 << " <= pow(" << f1 << ",-2) <= " << f4 << endl;
    //ARIADNE_TEST_ASSERT(Rational(f3)<Rational(16,25));
    //ARIADNE_TEST_ASSERT(Rational(f4)>Rational(16,25));

    // Floor and ceiling
    f2=RoundedFloatType(-3.25,precision); f3=RoundedFloatType(-2,precision);

    ARIADNE_TEST_ASSERT(floor(f1)==1); ARIADNE_TEST_ASSERT(ceil(f1)==2);
    ARIADNE_TEST_ASSERT(floor(f2)==-4); ARIADNE_TEST_ASSERT(ceil(f2)==-3);
    ARIADNE_TEST_ASSERT(floor(f3)==-2); ARIADNE_TEST_ASSERT(ceil(f3)==-2);

    // Check interval conversions
    RoundedFloatType z(0); RoundedFloatType o(1); RoundedFloatType t(3);

    RoundedFloatType::set_rounding_downward();
    RoundedFloatType odtd=div(o,t);
    RoundedFloatType::set_rounding_upward();
    RoundedFloatType odtu=div(o,t);
    cout << odtd << " <= 1/3 <= " << odtu << endl;
    ARIADNE_TEST_COMPARE(odtd,<,odtu);
    // Regression test to catch errors when RoundedFloatType result is not assigned to a variable

    ARIADNE_TEST_COMPARE(med(RoundedFloatType(2),RoundedFloatType(3)),==,2.5_x);
    ARIADNE_TEST_COMPARE(rad(RoundedFloatType(2),RoundedFloatType(3)),==,0.5_x);
}


template<class FLT> Void
TestRounded<FLT>::test_function()
{
    cout << __PRETTY_FUNCTION__ << endl;

    cout << setprecision(20);

    // Set up some variables
    RoundedFloatType x; RoundedFloatType ra; RoundedFloatType rl; RoundedFloatType ru;

    x=RoundedFloatType(1.0_x,precision);
    ra=exp(x);

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

template<class FLT> Void
TestRounded<FLT>::test_cosine()
{
    //3.14159265358979323846264338327950288419716939937510
    const RoundedFloatType pi_down=RoundedFloatType(FloatType::pi(down,precision));
    const RoundedFloatType pi_up  =RoundedFloatType(FloatType::pi(up,precision));
    const RoundedFloatType three  =RoundedFloatType(3,precision);

    RoundedFloatType::set_rounding_mode(downward);
    static const RoundedFloatType third_pi_down=div(pi_down,three);
    RoundedFloatType::set_rounding_mode(upward);
    static const RoundedFloatType third_pi_up  =div(pi_up  ,three);

    RoundedFloatType::set_rounding_mode(upward);
    ARIADNE_TEST_COMPARE(cos(hlf(pi_down)),>,0.0);
    ARIADNE_TEST_EQUAL(cos(RoundedFloatType(0.0)),1.0);
    ARIADNE_TEST_COMPARE(cos(third_pi_down),>,0.5_x);
    ARIADNE_TEST_COMPARE(sqr(cos(pi_down/4)),>,0.5_x);
    ARIADNE_TEST_COMPARE(cos(pi_down/2),>,0.0);
    ARIADNE_TEST_COMPARE(cos(pi_up/2),<=,0.0);
    ARIADNE_TEST_COMPARE(cos(pi_down),>,-1.0);
    ARIADNE_TEST_COMPARE(cos(pi_up),>,-1.0);
    ARIADNE_TEST_COMPARE(cos(2*pi_down),==,1.0);
    ARIADNE_TEST_COMPARE(cos(2*pi_up),==,1.0);
    ARIADNE_TEST_COMPARE(cos(3*pi_down),>,-1.0);
    ARIADNE_TEST_COMPARE(cos(3*pi_up),>,-1.0);

    RoundedFloatType::set_rounding_mode(downward);
    ARIADNE_TEST_EQUAL(cos(RoundedFloatType(0.0)),1.0);
    ARIADNE_TEST_COMPARE(cos(third_pi_up),<,0.5_x);
    ARIADNE_TEST_COMPARE(sqr(cos(pi_up/4)),<,0.5_x);
    ARIADNE_TEST_COMPARE(cos(pi_down/2),>=,0.0);
    ARIADNE_TEST_COMPARE(cos(pi_up/2),<,0.0);
    ARIADNE_TEST_COMPARE(cos(pi_down),==,-1.0);
    ARIADNE_TEST_COMPARE(cos(pi_up),==,-1.0);
    ARIADNE_TEST_COMPARE(cos(2*pi_down),<,1.0);
    ARIADNE_TEST_COMPARE(cos(2*pi_up),<,1.0);
    ARIADNE_TEST_COMPARE(cos(3*pi_down),==,-1.0);
    ARIADNE_TEST_COMPARE(cos(3*pi_up),==,-1.0);

    static const RoundedFloatType one(1.0,precision);
    RoundedFloatType::set_rounding_mode(upward);
    RoundedFloatType sin_rnd_up_one=sin(one);
    RoundedFloatType::set_rounding_mode(downward);
    RoundedFloatType sin_rnd_down_one=sin(one);
    RoundedFloatType::set_rounding_mode(to_nearest);
    RoundedFloatType sin_rnd_approx_one=sin(one);
    ARIADNE_TEST_COMPARE(sin_rnd_down_one,<=,sin_rnd_approx_one);
    ARIADNE_TEST_COMPARE(sin_rnd_approx_one,<=,sin_rnd_up_one);
    ARIADNE_TEST_COMPARE(sin_rnd_down_one,< ,sin_rnd_up_one);
}

template<> Void
TestRounded<FloatDP>::test_arctan()
{
    PR pr;

    static const RoundedFloatType pi_down(3.1415926535897931_x,pr);
    //static const RoundedFloatType pi_near(3.1415926535897931_x,pr);
    static const RoundedFloatType pi_up  (3.1415926535897936_x,pr);

    static const RoundedFloatType atan_quarter_down(0.244978663126864143_x,pr);
    //static const RoundedFloatType atan_quarter_near(0.244978663126864154_x,pr);
    static const RoundedFloatType atan_quarter_up  (0.244978663126864171_x,pr);

    static const RoundedFloatType sqrt_three_down(1.732050807568877193_x,pr);
    //static const RoundedFloatType sqrt_three_near(1.732050807568877294_x,pr);
    static const RoundedFloatType sqrt_three_up  (1.732050807568877415_x,pr);

    static const RoundedFloatType zero(0.0_x,pr);
    static const RoundedFloatType one(1.0_x,pr);
    static const RoundedFloatType eps(FloatType::eps(pr));

    RoundedFloatType::set_rounding_mode(upward);
    ARIADNE_TEST_EQUAL(atan(zero),zero);
    ARIADNE_TEST_COMPARE(atan(one)*4,>=,pi_up);
    ARIADNE_TEST_COMPARE(atan(-one)*4,>=,-pi_down);
    ARIADNE_TEST_COMPARE(atan(sqrt_three_up)*3,>=,pi_up);
    ARIADNE_TEST_COMPARE(atan(-sqrt_three_down)*3,>=,-pi_down);
    ARIADNE_TEST_COMPARE(atan(1/sqrt_three_down)*6,>=,pi_up);
    ARIADNE_TEST_COMPARE(atan(-1/sqrt_three_up)*6,>=,-pi_down);
    ARIADNE_TEST_COMPARE(atan(one/4),>=,atan_quarter_up);
    ARIADNE_TEST_COMPARE(atan(-one/4),>=,-atan_quarter_down);

    RoundedFloatType::set_rounding_mode(downward);
    ARIADNE_TEST_EQUAL(atan(zero),zero);
    ARIADNE_TEST_COMPARE(atan(+one)*4,<=,pi_down);
    ARIADNE_TEST_COMPARE(atan(-one)*4,<=,-pi_up);
    ARIADNE_TEST_COMPARE(atan(+sqrt_three_down)*3,<=,pi_down);
    ARIADNE_TEST_COMPARE(atan(-sqrt_three_up)*3,<=,-pi_up);
    ARIADNE_TEST_COMPARE(atan(1/sqrt_three_up)*6,<=,pi_down);
    ARIADNE_TEST_COMPARE(atan(-1/sqrt_three_down)*6,<=,-pi_up);
    ARIADNE_TEST_COMPARE(atan(+one/4),<=,atan_quarter_down);
    ARIADNE_TEST_COMPARE(atan(-one/4),<=,-atan_quarter_up);

    ARIADNE_TEST_COMPARE(atan(one)*4-pi_up,<=,4*eps);
    ARIADNE_TEST_COMPARE(atan(sqrt_three_up)*3-pi_up,<=,4*eps);
    ARIADNE_TEST_COMPARE(atan(1/sqrt_three_down)*6-pi_up,<=,4*eps);
    ARIADNE_TEST_COMPARE(atan(one/4)-atan_quarter_up,<=,eps/2);

}


template<class FLT> Void
TestRounded<FLT>::test_arctan()
{
    //pi~=3.14159265358979323846264338327950288419716939937510
    const RoundedFloatType pi_down(FloatType::pi(down,precision));
    const RoundedFloatType pi_near(FloatType::pi(near,precision));
    const RoundedFloatType pi_up  (FloatType::pi(up,precision));

    if constexpr(std::is_same<FloatType,FloatDP>::value) {
        ARIADNE_TEST_EQUAL(pi_down.raw(),3.1415926535897931);
        ARIADNE_TEST_EQUAL(pi_near.raw(),3.1415926535897931);
        ARIADNE_TEST_EQUAL(pi_up.raw(),  3.1415926535897936);
    }

    RoundedFloatType three=RoundedFloatType(3,precision);
    RoundedFloatType::set_rounding_mode(downward);
    const RoundedFloatType sqrt_three_down=sqrt(three);
    RoundedFloatType::set_rounding_mode(upward);
    const RoundedFloatType sqrt_three_up  =sqrt(three);

    const RoundedFloatType zero(0,precision);
    const RoundedFloatType one(1,precision);
    const RoundedFloatType four(4,precision);
    const RoundedFloatType six(6,precision);
    const RoundedFloatType eps(FloatType::eps(precision));

    RoundedFloatType::set_rounding_mode(upward);
    ARIADNE_TEST_EQUAL(atan(zero),zero);
    ARIADNE_TEST_COMPARE(mul(atan(+one),four),>=,pi_up);
    ARIADNE_TEST_COMPARE(mul(atan(-one),four),>=,-pi_down);
    ARIADNE_TEST_COMPARE(mul(atan(+sqrt_three_up),three),>=,pi_up);
    ARIADNE_TEST_COMPARE(mul(atan(-sqrt_three_down),three),>=,-pi_down);
    ARIADNE_TEST_COMPARE(mul(atan(+1/sqrt_three_down),six),>=,pi_up);
    ARIADNE_TEST_COMPARE(mul(atan(-1/sqrt_three_up),six),>=,-pi_down);

    RoundedFloatType::set_rounding_mode(downward);
    ARIADNE_TEST_EQUAL(atan(zero),zero);
    ARIADNE_TEST_COMPARE(mul(atan(+one),four),<=,pi_down);
    ARIADNE_TEST_COMPARE(mul(atan(-one),four),<=,-pi_up);
    ARIADNE_TEST_COMPARE(mul(atan(+sqrt_three_down),three),<=,pi_down);
    ARIADNE_TEST_COMPARE(mul(atan(-sqrt_three_up),three),<=,-pi_up);
    ARIADNE_TEST_COMPARE(mul(atan(+1/sqrt_three_up),six),<=,pi_down);
    ARIADNE_TEST_COMPARE(mul(atan(-1/sqrt_three_down),six),<=,-pi_up);

    ARIADNE_TEST_COMPARE(sub(mul(atan(one),four),pi_up),<=,4*eps);
    ARIADNE_TEST_COMPARE(sub(mul(atan(sqrt_three_up),three),pi_up),<=,4*eps);
    ARIADNE_TEST_COMPARE(sub(mul(atan(1/sqrt_three_down),six),pi_up),<=,4*eps);

    ARIADNE_TEST_COMPARE(atan(one)*4-pi_up,<=,4*eps);
    ARIADNE_TEST_COMPARE(atan(sqrt_three_up)*3-pi_up,<=,4*eps);
    ARIADNE_TEST_COMPARE(atan(1/sqrt_three_down)*6-pi_up,<=,4*eps);

    // Regression test
    ARIADNE_TEST_COMPARE(pi_down/4,<,pi_up/4);

}
