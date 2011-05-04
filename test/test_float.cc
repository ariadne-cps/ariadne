/***************************************************************************
 *            test_float.cc
 *
 *  Copyright  2006-11  Alberto Casagrande, Pieter Collins
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

#include "rounding.h"
#include "numeric.h"

#include "test.h"

using namespace std;
using namespace Ariadne;


class TestFloat 
{
  public:
    void test();
  private:
    void test_concept();
    void test_class();
    void test_conversion();
    void test_stream();
    void test_comparison();
    void test_rounding();
    void test_arithmetic();
    void test_cosine();
    void test_function();
};


int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    TestFloat().test();

    return ARIADNE_TEST_FAILURES;
}


void
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
// the Float concept without testing correctness
void
TestFloat::test_concept()
{
    bool b=true;
    int n=1; 
    uint m=1;
    double d=1;
    Float x=1;

    // Constructors
    x=Float(); x=Float(n); x=Float(m); x=Float(d); x=Float(x);
  
    // Assignment
    x=n; x=m; x=d; x=x; 
  
    // Conversion
    d=x;

    // Maximum and minimum
    x=max(x,x); x=min(x,x);


    // Exact operations
    x=abs(x); x=neg(x); x=+x; x=-x;

    // Rounded arithmetic
    x=add_approx(x,x); x=add_down(x,x); x=add_up(x,x); // x=add_chop(x,x);
    x=sub_approx(x,x); x=sub_down(x,x); x=sub_up(x,x); // x=sub_chop(x,x);
    x=mul_approx(x,x); x=mul_down(x,x); x=mul_up(x,x); // x=mul_chop(x,x);
    x=div_approx(x,x); x=div_down(x,x); x=div_up(x,x); // x=div_chop(x,x);
    pow_approx(x,n); x=pow_down(x,n); x=pow_up(x,n); // x=pow_chop(x,n);
    pow_approx(x,m); x=pow_down(x,m); x=pow_up(x,m); // x=pow_chop(x,m);

    x=med_approx(x,x); x=rad_up(x,x);

    // Mixed Float/int arithmetic
    x=mul_approx(n,x); x=mul_down(n,x); x=mul_up(n,x); // x=mul_chop(n,x);
    x=mul_approx(m,x); x=mul_down(m,x); x=mul_up(m,x); // x=mul_chop(m,x);
    x=mul_approx(x,n); x=mul_down(x,n); x=mul_up(x,n); // x=mul_chop(x,n);
    x=mul_approx(x,m); x=mul_down(x,m); x=mul_up(x,m); // x=mul_chop(x,m);
    x=div_approx(x,n); x=div_down(x,n); x=div_up(x,n); // x=div_chop(x,n);
    x=div_approx(x,m); x=div_down(x,m); x=div_up(x,m); // x=div_chop(x,m);

    // Mixed Float/double arithmetic
    x=mul_approx(d,x); x=mul_approx(x,d); x=div_approx(x,d); 

    // Reset x to 1
    x=1;

    // Comparisons
    b=(x==n); b=(x!=n); b=(x<=n); b=(x>=n); b=(x<n); b=(x>n);
    b=(n==x); b=(n!=x); b=(n<=x); b=(n>=x); b=(n<x); b=(n>x);
    b=(x==m); b=(x!=m); b=(x<=m); b=(x>=m); b=(x<m); b=(x>m);
    b=(m==x); b=(m!=x); b=(m<=x); b=(m>=x); b=(m<x); b=(m>x);
    b=(x==d); b=(x!=d); b=(x<=d); b=(x>=d); b=(x<d); b=(x>d);
    b=(d==x); b=(d!=x); b=(d<=x); b=(d>=x); b=(d<x); b=(d>x);
    b=(x==x); b=(x!=x); b=(x<=x); b=(x>=x); b=(x<x); b=(x>x);

    // Rounding mode
    set_rounding_to_nearest();
    set_rounding_downward();
    set_rounding_upward();
    set_rounding_toward_zero();

    set_rounding_mode(to_nearest);
    set_rounding_mode(downward);
    set_rounding_mode(upward);
    set_rounding_mode(toward_zero);

    rounding_mode_t rnd=get_rounding_mode();
    set_rounding_mode(rnd);
}


void
TestFloat::test_class()
{
    cout << __PRETTY_FUNCTION__ << endl;
    // Construct from an int
    Float f1(2);
    ARIADNE_TEST_ASSERT(f1==2);
    // Construct from a double
    Float f2(1.25);
    ARIADNE_TEST_ASSERT(f2==1.25);
    // Copy constructor
    Float f3(f2);
    ARIADNE_TEST_ASSERT(f3==f2);
  
    // Assign from an int
    f1=3;
    ARIADNE_TEST_ASSERT(f1==3);
    // Assign from a double
    f2=2.25;
    ARIADNE_TEST_ASSERT(f2==2.25);
    // Copy assignment
    f3=f2;
    ARIADNE_TEST_ASSERT(f3==f2);

}


void
TestFloat::test_conversion()
{
    cout << __PRETTY_FUNCTION__ << endl;
 
    // Convert from integers
    int n;
    n=std::numeric_limits<int>::min();
    ARIADNE_TEST_ASSERT(Float(n)==n);
    n=std::numeric_limits<int>::max();
    ARIADNE_TEST_ASSERT(Float(n)==n);
    n=std::numeric_limits<unsigned int>::max();
    ARIADNE_TEST_ASSERT(Float(n)==n);

    // Conversion from long integers cannot be performed exactly, so is banned

    /*
    // Convert to a rational
    Rational q;
    q=Rational(R(3));
    ARIADNE_TEST_ASSERT(q==3);

    // Assign to a rational
    q=Rational(R(2.25));
    ARIADNE_TEST_ASSERT(q==2.25);
  
    // Convert from a rational
    q=Rational(1,3); 
    R f1=R(q,round_down); 
    R f2=R(q,round_up);
    cout << f1 << " <= " << q << " <= " << f2 << endl;
    ARIADNE_TEST_ASSERT(f1<=q); ARIADNE_TEST_ASSERT(f2>=q); ARIADNE_TEST_ASSERT(f1<f2);
  
    // Convert from a negative rational
    q=Rational(-2,5); cout << q << endl;
    f1.set(q,round_down); cout << f1 << endl;
    f2.set(q,round_up);
    cout << f1 << " <= " << q << " <= " << f2 << endl;
    ARIADNE_TEST_ASSERT(f1<=q); ARIADNE_TEST_ASSERT(f2>=q); ARIADNE_TEST_ASSERT(f1<f2);
    */
};


void
TestFloat::test_stream()
{
    cout << __PRETTY_FUNCTION__ << endl;

    stringstream ss("1.25 -2.25 42 375e1 2.35e1");
    Float f1,f2,f3,f4,f5;

    ss >> f1;
    cout << f1 << endl;
    ARIADNE_TEST_ASSERT(f1==1.25);
    ss >> f2;
    cout << f2 << endl;
    ARIADNE_TEST_ASSERT(f2==-2.25);
    ss >> f3;
    cout << f3 << endl;
    ARIADNE_TEST_ASSERT(f3==42);
    ss >> f4;
    cout << f4 << endl;
    ARIADNE_TEST_ASSERT(f4==3750);
    ss >> f5;
    cout << f5 << endl;
    ARIADNE_TEST_ASSERT(f5==23.5);
}


void
TestFloat::test_comparison()
{
    cout << __PRETTY_FUNCTION__ << endl;
  
    Float f1(1.25); Float f2(-1.25); Float f3(-2.25); Float f4(1.25);

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
    int i2=1;
    ARIADNE_TEST_ASSERT(!(f1==i2)); ARIADNE_TEST_ASSERT(f1!=i2); 
    ARIADNE_TEST_ASSERT(!(f1<=i2)); ARIADNE_TEST_ASSERT(f1> i2);
    ARIADNE_TEST_ASSERT(f1>=i2); ARIADNE_TEST_ASSERT(!(f1< i2));
  
    int i1=1;
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
  
void
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
    
    set_rounding_mode(downward);
    double onethirdrounddown=one/three;
    ARIADNE_TEST_EQUAL(onethirdrounddown, onethirddown);
    set_rounding_mode(upward);
    double onethirdroundup=one/three;
    ARIADNE_TEST_EQUAL(onethirdroundup, onethirdup);
    set_rounding_mode(toward_zero);
    double onethirdroundchop=one/three;
    ARIADNE_TEST_EQUAL(onethirdroundchop, onethirdchop);
    set_rounding_mode(to_nearest);
    double onethirdroundnearest=one/three;
    ARIADNE_TEST_EQUAL(onethirdroundnearest, onethirdnearest);

    set_rounding_downward();
    double twofifthsrounddown=two/five;
    ARIADNE_TEST_EQUAL(twofifthsrounddown, twofifthsdown);
    set_rounding_upward();
    double twofifthsroundup=two/five;
    ARIADNE_TEST_EQUAL(twofifthsroundup, twofifthsup);
    set_rounding_toward_zero();
    double twofifthsroundchop=two/five;
    ARIADNE_TEST_EQUAL(twofifthsroundchop, twofifthschop);
    set_rounding_to_nearest();
    double twofifthsroundnearest=two/five;
    ARIADNE_TEST_EQUAL(twofifthsroundnearest, twofifthsnearest);
}

void
TestFloat::test_arithmetic()
{
    cout << __PRETTY_FUNCTION__ << endl;
  
    // Set up some variables
    Float f1(1.25); Float f2(2.25); Float f3(-3.25); Float f4; Float f5;
 
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
    f3=med_approx(f1,f2);
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
    f3=down(f1);
    f4=up(f1);
    cout << f3 << " < " << f1 << " < " << f4 << endl;
    ARIADNE_TEST_ASSERT(f3<f1); ARIADNE_TEST_ASSERT(f4>f1);

    // Addition (this should remain exact here)
    f3=add_down(f1,f2);
    f4=add_up(f1,f2);
    cout << f3 << " <= " << f1 << " + " << f2 << " <= " << f4 << endl;
    ARIADNE_TEST_ASSERT(f3<=3.5); ARIADNE_TEST_ASSERT(f4>=3.5);
  
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
    ARIADNE_TEST_ASSERT(f3<f4); ARIADNE_TEST_ASSERT(f3<=f5); ARIADNE_TEST_ASSERT(f4>=f5);
    //ARIADNE_TEST_ASSERT(Rational(f3)<=Rational(5,9));
    //ARIADNE_TEST_ASSERT(Rational(f4)>=Rational(5,9));
  
    // Check division my multipyling back
    cout << mul_down(f3,f2) << " <= (" << f1 << "/" << f2 << ")*" << f2 << " <= " << mul_up(f4,f2) << endl;
    ARIADNE_TEST_ASSERT(mul_down(f3,f2)<f1);
    ARIADNE_TEST_ASSERT(mul_up(f4,f2)>f1);
  
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
    f2=Float(-3.25); f3=Float(-2);
  
    ARIADNE_TEST_ASSERT(floor(f1)==1); ARIADNE_TEST_ASSERT(ceil(f1)==2);
    ARIADNE_TEST_ASSERT(floor(f2)==-4); ARIADNE_TEST_ASSERT(ceil(f2)==-3);
    ARIADNE_TEST_ASSERT(floor(f3)==-2); ARIADNE_TEST_ASSERT(ceil(f3)==-2);
  
    // Conversion to integer types
    int i3,i4;
    i3=int(floor(f1));
    i4=int(ceil(f1));
    cout << i3 << " < " << f1 << " < " << i4 << endl;
    ARIADNE_TEST_ASSERT(i3==1); ARIADNE_TEST_ASSERT(i4==2);
    i3=int(floor(f2));
    i4=int(ceil(f2));
    cout << i3 << " < " << f2 << " < " << i4 << endl;
    ARIADNE_TEST_ASSERT(i3==-4); ARIADNE_TEST_ASSERT(i4==-3);
  
    // Check interval conversions
    Float z(0); Float o(1); Float t(3);
    Interval io(1); Interval it(3);
  
    cout << "o/t=" << Interval(o/t) << endl;

    Interval iao=(o/t)*t;
    cout << "(o/t)*t=" << iao << endl;
    cout << iao.lower() << " " << iao.upper() << endl;

    cout << div_down(o,t) << " <= 1/3 <= " << div_up(o,t) << endl;
    ARIADNE_TEST_ASSERT(div_down(o,t)<div_up(o,t));
    cout << Interval(o/t) << endl; 
    ARIADNE_TEST_ASSERT(contains(iao,o));
    Interval iaz=iao-io;
    Interval iz=Float(1)-Float(1);
    cout << iaz << endl;
    ARIADNE_TEST_ASSERT(contains(iaz,z)); 
    ARIADNE_TEST_ASSERT(!bool(!subset(iz,iaz)));
    cout << endl;

    ARIADNE_TEST_ASSERT(med_approx(Float(2),Float(3))==2.5);
    ARIADNE_TEST_ASSERT(rad_up(Float(2),Float(3))>=0.5);
    ARIADNE_TEST_ASSERT(rad_up(Float(2),Float(3))<=0.5000000000000002);

    // The following line should not compile
    // f5=f1+f2;
  
}


void
TestFloat::test_function()
{
    cout << __PRETTY_FUNCTION__ << endl;
  
    cout << setprecision(20);
  
    // Set up some variables
    Float x; Float ra; Float rl; Float ru;
  
    x=1;
    ra=Ariadne::exp(x);
    rl=down(ra);
    ru=up(ra);
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

void
TestFloat::test_cosine()

{   
    //3.14159265358979323846264338327950288419716939937510
    static const double pi_down=3.1415926535897931;
    //static const double pi_approx=3.1415926535897931;
    static const double pi_up=3.1415926535897936;
     
    static const double third_pi_down=1.0471975511965976;
    static const double third_pi_up=1.0471975511965979;
    
    static const double sqrt_half_down=0.70710678118654746;
    //static const double sqrt_half_approx=0.70710678118654757;
    static const double sqrt_half_up=0.70710678118654757;

    cout << __PRETTY_FUNCTION__ << endl;

    set_rounding_mode(upward);
    ARIADNE_TEST_EQUAL(cos_rnd(0.0),1.0);
    ARIADNE_TEST_EQUAL(cos_rnd(0.0),1.0);
    ARIADNE_TEST_COMPARE(cos_rnd(third_pi_down),>,0.5);
    ARIADNE_TEST_COMPARE(cos_rnd(pi_down/4),>,sqrt_half_up);
    ARIADNE_TEST_COMPARE(cos_rnd(pi_down/2),>,0.0);
    ARIADNE_TEST_COMPARE(cos_rnd(pi_up/2),<=,0.0);
    ARIADNE_TEST_COMPARE(cos_rnd(pi_down),>,-1.0);
    ARIADNE_TEST_COMPARE(cos_rnd(pi_up),>,-1.0);
    ARIADNE_TEST_COMPARE(cos_rnd(2*pi_down),==,1.0);
    ARIADNE_TEST_COMPARE(cos_rnd(2*pi_up),==,1.0);
    ARIADNE_TEST_COMPARE(cos_rnd(3*pi_down),>,-1.0);
    ARIADNE_TEST_COMPARE(cos_rnd(3*pi_up),>,-1.0);

    set_rounding_mode(downward);
    ARIADNE_TEST_EQUAL(cos_rnd(0.0),1.0);
    ARIADNE_TEST_EQUAL(cos_rnd(0.0),1.0);
    ARIADNE_TEST_COMPARE(cos_rnd(third_pi_up),<,0.5);
    ARIADNE_TEST_COMPARE(cos_rnd(pi_up/4),<,sqrt_half_down);
    ARIADNE_TEST_COMPARE(cos_rnd(pi_down/2),>=,0.0);
    ARIADNE_TEST_COMPARE(cos_rnd(pi_up/2),<,0.0);
    ARIADNE_TEST_COMPARE(cos_rnd(pi_down),==,-1.0);
    ARIADNE_TEST_COMPARE(cos_rnd(pi_up),==,-1.0);
    ARIADNE_TEST_COMPARE(cos_rnd(2*pi_down),<,1.0);
    ARIADNE_TEST_COMPARE(cos_rnd(2*pi_up),<,1.0);
    ARIADNE_TEST_COMPARE(cos_rnd(3*pi_down),==,-1.0);
    ARIADNE_TEST_COMPARE(cos_rnd(3*pi_up),==,-1.0);

    set_rounding_mode(to_nearest);

    double x=1.23; double y=1.23;
    ARIADNE_TEST_COMPARE(x,<=,y);

}
