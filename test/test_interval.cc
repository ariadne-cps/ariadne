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
#include "numeric.h"

#include "test.h"

using namespace Ariadne;
using namespace std;




class TestInterval
{
    typedef Interval I;
    typedef Float R;
  public:
    void test();
  private:
    void test_concept();
    void test_misc();
    void test_class();
    void test_conversion();
};


void
TestInterval::test()
{
    test_concept();
    test_misc();
}

void
TestInterval::test_concept()
{
    int n=1; 
    uint m=1;
    double d=1;
    Float x=1;
    Interval i=1;
    Interval j=1;
 
  
    // Constructors
    j=I(); j=I(n); j=I(m); j=I(d); j=I(x); j=I(i);
    j=I(n,n); j=I(m,m); j=I(d,d); j=I(x,x);
    j=I(n,m); j=I(m,d); j=I(d,n); 
    // Assignment
    j=n; j=m; j=d; j=x; j=i;

    // Exact operations
    j=abs(i); j=neg(i); j=rec(i);

    // Arithmetic
    j=add(x,x); j=add(x,i); j=add(i,x); j=add(i,i);
    j=sub(x,x); j=sub(x,i); j=sub(i,x); j=sub(i,i);
    j=mul(x,x); j=mul(x,i); j=mul(i,x); j=mul(i,i);
    j=div(x,x); j=div(x,i); j=div(i,x); j=div(i,i);
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


void
TestInterval::test_misc()
{
    cout << __PRETTY_FUNCTION__ << endl;
    typedef Interval I;


    Float zero=0;
  
    // Construct from pair
    Interval ivld1(Float(1.125),Float(2.25));
    ARIADNE_TEST_ASSERT(ivld1.lower()==1.125); ARIADNE_TEST_ASSERT(ivld1.upper()==2.25);
    // Default constructor
    Interval ivld2;
    if(ivld2.lower()>ivld2.upper()) { 
        cerr << "Warning: Interval default constructor returns an empty set\n"; 
    } else {
        ARIADNE_TEST_ASSERT((bool)(ivld2==Interval(zero,zero)));
    }
    // Constructor with approximations
#ifdef HAVE_GMPXX_H
    Interval ivld3(Rational(21,10),Rational(16,5));
    cout<<ivld3<<std::endl;
    ARIADNE_TEST_ASSERT(ivld3.lower()<Rational(21,10));
    ARIADNE_TEST_ASSERT(ivld3.upper()>Rational(16,5));
#else
    Interval ivld3(2.1,1.6);
#endif // HAVE_GMPXX_H

    // Constructor from approximate values
    Interval ivld4(2.1,3.2);
    ARIADNE_TEST_ASSERT(ivld4.lower()<=2.1);
    ARIADNE_TEST_ASSERT(ivld4.upper()>=3.2);

#ifdef HAVE_GMPXX_H
    // Approximate constructor from a single value
    Interval ivld5(Rational(1,3));
    ARIADNE_TEST_ASSERT(ivld5.lower()<Rational(1,3));
    ARIADNE_TEST_ASSERT(ivld5.upper()>Rational(1,3));
#else
    Interval ivld5(1./3.);
#endif // HAVE_GMPXX_H

    // Exact constructor from a single value
    Interval ivld6(Float(1.25));
    ARIADNE_TEST_ASSERT(ivld6.lower()==Float(1.25));
    ARIADNE_TEST_ASSERT(ivld6.upper()==Float(1.25));
  
    Interval ivlq1(1.1,2.2);
    Interval ivlq2(1.125,1.125);
    Interval ivlq3(2.125,3.25);
  
    Float f0=-2.25;
    Float f1=1.5;
    Float f2=2.25;
    Float f3=3.125;
    Float f4=4.0625;
    Interval ivlf1;
    Interval ivlf2(f1,f2);
    Interval ivlf3(f3,f4);
  
    Float r1=1.5;
    Float r2=2.25;
    Float r3=3.125;
    Float r4=4.0625;
    Interval ivlr1;
    Interval ivlr2(r1,r2);
    Interval ivlr3(r3,r4);
    Interval ivlr4(r1,r4);

    cout << "ivlr1=" << ivlr1 << ", ivlr2=" << ivlr2 << ", ivlr3=" << ivlr3 
         << ", ivlr4=" << ivlr4 << endl;
  
    Interval ivlf0(f2,f1);
    cout << "ivlf0=" << ivlf0 << endl;
    cout << "ivlf0.empty()=" << ivlf0.empty() << ", ivlf1.empty()=" << ivlf1.empty() << endl;
  
    ivlf1=Interval(f1,f3);
    ivlf2=Interval(f0,f1);
    ivlf3=Interval(f2,f4);
    cout << "ivlf1=" << ivlf1 << ", ivlf2=" << ivlf2 << ", ivlf3=" << ivlf3 << endl;
    cout << "max(ivlf1,ivlf3)=" << max(ivlf1,ivlf3) << endl;
    cout << "min(ivlf1,ivlf3)=" << min(ivlf1,ivlf3) << endl;
    cout << "lower(" << ivlf1 << ")=" << ivlf1.lower() << endl;
    cout << "upper(" << ivlf1 << ")=" << ivlf1.upper() << endl;
    cout << "midpoint(" << ivlf1 << ")=" << ivlf1.midpoint() << endl;
    cout << "radius(" << ivlf2 << ")=" << radius(ivlf2) << endl;
    cout << "width(" << ivlf2 << ")=" << width(ivlf2) << endl;
    ARIADNE_TEST_ASSERT(ivlf1.midpoint()==Float(2.3125));
  
    try {
        string input("[1.125,2.25] ");
        stringstream iss(input);
    
        iss >> ivld2;
        cout << "ivld1=" << ivld1 << "  ivld2=" << ivld2 << endl;
        ARIADNE_TEST_ASSERT(equal(ivld1,ivld2));
        ARIADNE_TEST_ASSERT(ivld1.lower()==ivld2.lower() && ivld1.upper()==ivld2.upper());
    
        ARIADNE_TEST_ASSERT(equal(ivld1,ivld2));
    
        // FIXME: If using Boost style interval tests, uncomment the line below
        // and comment out the line after
        //ARIADNE_TEST_ASSERT(indeterminate(ivld1==ivld2));
        ARIADNE_TEST_ASSERT(ivld1==ivld2);
    
        Interval& ivlf1ref=ivlf1;
        ivlf1ref=Interval(5.25,7.375);
        cout << "ivlf1ref=" << ivlf1ref << endl;
        ARIADNE_TEST_ASSERT(ivlf1ref.lower()==Float(5.25));
    
        ivld1 = ivld2+ivld3;
        cout << ivld2 << " + " << ivld3 << " = " << ivld1 << endl;
        ivld1 = ivld2-ivld3;
        cout << ivld2 << " - " << ivld3 << " = " << ivld1 << endl;
        ivld1 = ivld2*ivld3;
        cout << ivld2 << " * " << ivld3 << " = " << ivld1 << endl;
        ivld1 = ivld2/ivld3;
        cout << ivld2 << " / " << ivld3 << " = " << ivld1 << endl;
        cout << endl;
        //ivld1 = cos(ivld2);
    
        ivlf1 = ivlf2+ivlf3;
        cout << ivlf2 << " + " << ivlf3 << " = " << ivlf1 << endl;
        ivlf1 = ivlf2-ivlf3;
        cout << ivlf2 << " - " << ivlf3 << " = " << ivlf1 << endl;
        ivlf1 = ivlf2*ivlf3;
        cout << ivlf2 << " * " << ivlf3 << " = " << ivlf1 << endl;
        ivlf1 = ivlf2/ivlf3;
        cout << ivlf2 << " / " << ivlf3 << " = " << ivlf1 << endl;

        ivlf1 = ivlf2+f3;
        cout << ivlf2 << " + " << f3 << " = " << ivlf1 << endl;
        ivlf1 = f2+ivlf3;
        cout << f2 << " + " << ivlf3 << " = " << ivlf1 << endl;
        ivlf1 = ivlf2-f3;
        cout << ivlf2 << " - " << f3 << " = " << ivlf1 << endl;
        ivlf1 = f2-ivlf3;
        cout << f2 << " - " << ivlf3 << " = " << ivlf1 << endl;
        ivlf1 = ivlf2*f3;
        cout << ivlf2 << " * " << f3 << " = " << ivlf1 << endl;
        ivlf1 = f2*ivlf3;
        cout << f2 << " * " << ivlf3 << " = " << ivlf1 << endl;
        ivlf1 = ivlf2/f3;
        cout << ivlf2 << " / " << f3 << " = " << ivlf1 << endl;
        ivlf1 = f2/ivlf3;
        cout << f2 << " / " << ivlf3 << " = " << ivlf1 << endl;
        cout << endl;
    
        // Check to make sure aliases are handled correctly
        ivlf1=ivlf3; ivlf1=ivlf2-ivlf1; ARIADNE_TEST_ASSERT(equal(ivlf1,I(ivlf2-ivlf3)));
        ivlf1=ivlf3; ivlf1=ivlf2*ivlf1; ARIADNE_TEST_ASSERT(equal(ivlf1,I(ivlf2*ivlf3)));
        ivlf1=ivlf2; ivlf1=ivlf1*ivlf3; ARIADNE_TEST_ASSERT(equal(ivlf1,I(ivlf2*ivlf3)));
        ivlf1=ivlf2; ivlf1=ivlf1*f3; ARIADNE_TEST_ASSERT(equal(ivlf1,I(ivlf2*f3)));
        ivlf1=ivlf3; ivlf1=f2*ivlf1; ARIADNE_TEST_ASSERT(equal(ivlf1,I(f2*ivlf3)));
        ivlf1=ivlf2; ivlf1=ivlf1/ivlf3; ARIADNE_TEST_ASSERT(equal(ivlf1,I(ivlf2/ivlf3)));
        ivlf1=ivlf2; ivlf1=ivlf1/f3; ARIADNE_TEST_ASSERT(equal(ivlf1,I(ivlf2/f3)));
        ivlf1=ivlf3; ivlf1=f2/ivlf1; ARIADNE_TEST_ASSERT(equal(ivlf1,I(f2/ivlf3)));
    

        ivlq1 = ivlq2+ivlq3;
        cout << ivlq2 << " + " << ivlq3 << " = " << ivlq1 << endl;
        ivlq1 = ivlq2-ivlq3;
        cout << ivlq2 << " - " << ivlq3 << " = " << ivlq1 << endl;
        ivlq1 = ivlq2*ivlq3;
        cout << ivlq2 << " * " << ivlq3 << " = " << ivlq1 << endl;
        ivlq1 = ivlq2/ivlq3; 
        cout << ivlq2 << " / " << ivlq3 << " =" << ivlq1 << endl;
        cout << endl;
        //ivlr1 = sin(ivlr2);

        // ensure proper rounding using ARIADNE_TEST_ASSERTions
        Interval ivlo(1.0);
        Interval ivlt(3.0);
        Interval ivlodt=ivlo/ivlt;
        Interval ivloa=ivlodt*ivlt;
        cout << ivlo << " / " << ivlt << " = " << ivlodt << endl;
        cout << ivlo << " in " << Interval(ivlodt * ivlt) << endl;
        ARIADNE_TEST_ASSERT(ivlodt.lower() < ivlodt.upper());
        ARIADNE_TEST_ASSERT(ivloa.lower() < ivlo.lower() && ivloa.upper() > ivlo.upper());
    
    

        ivlf1=Interval(-13,-7);
        ivlf2=Interval(-3,2);
        ivlf3=Interval(5,11);
    
        // Fixme: The explicit cast should not be needed in the output
        cout << ivlf1 << " * " << ivlf1 << " = " << Interval(ivlf1*ivlf1) << endl;
        cout << ivlf1 << " * " << ivlf2 << " = " << Interval(ivlf1*ivlf2) << endl;
        cout << ivlf1 << " * " << ivlf3 << " = " << Interval(ivlf1*ivlf3) << endl;
        cout << ivlf2 << " * " << ivlf1 << " = " << Interval(ivlf2*ivlf1) << endl;
        cout << ivlf2 << " * " << ivlf2 << " = " << Interval(ivlf2*ivlf2) << endl;
        cout << ivlf2 << " * " << ivlf3 << " = " << Interval(ivlf2*ivlf3) << endl;
        cout << ivlf3 << " * " << ivlf1 << " = " << Interval(ivlf3*ivlf1) << endl;
        cout << ivlf3 << " * " << ivlf2 << " = " << Interval(ivlf3*ivlf2) << endl;
        cout << ivlf3 << " * " << ivlf3 << " = " << Interval(ivlf3*ivlf3) << endl;

        cout << ivlf1 << " / " << ivlf1 << " = " << Interval(ivlf1/ivlf1) << endl;
        cout << ivlf1 << " / " << ivlf3 << " = " << Interval(ivlf1/ivlf3) << endl;
        cout << ivlf2 << " / " << ivlf1 << " = " << Interval(ivlf2/ivlf1) << endl;
        cout << ivlf2 << " / " << ivlf3 << " = " << Interval(ivlf2/ivlf3) << endl;
        cout << ivlf3 << " / " << ivlf1 << " = " << Interval(ivlf3/ivlf1) << endl;
        cout << ivlf3 << " / " << ivlf3 << " = " << Interval(ivlf3/ivlf3) << endl;


    }
    catch(exception& e) {
        cout << "EXCEPTION " << e.what() << "\n";
        throw e;
    }
  
    // e=2.718281828459045235360287471352662497[78]
    // 2*pi = 6.283185307179586476925286766559005768[3,5]

    // Regression test; fails dramatically on certain types of rounding
    {
        Interval one(1.0);
        Interval expone=exp(one);
        ARIADNE_TEST_PRINT(one);
        ARIADNE_TEST_COMPARE(expone.lower(),<,2.71828182845905);
        ARIADNE_TEST_COMPARE(expone.lower(),>,2.71828182845903);
        ARIADNE_TEST_COMPARE(expone.upper(),>,2.71828182845904);
        ARIADNE_TEST_COMPARE(expone.upper(),<,2.71828182845906);
        ARIADNE_TEST_ASSERT(expone.lower()<expone.upper());
    }

    {
        Interval e(2.7182818284590451,2.7182818284590455);
        Interval loge=log(e);
        ARIADNE_TEST_PRINT(e);
        ARIADNE_TEST_COMPARE(loge.lower(),<,1);
        ARIADNE_TEST_COMPARE(loge.lower(),>,0.9999999999998);
        ARIADNE_TEST_COMPARE(loge.upper(),>,1);
        ARIADNE_TEST_COMPARE(loge.upper(),<,1.000000000002);
    }


  
    // Regression test; fails dramatically on certain types of rounding
    {
        Interval x(6.283185307179586,6.283185307179587);
        Interval sinx=sin(x);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_COMPARE(sinx.lower(),<,0.0);
        ARIADNE_TEST_COMPARE(sinx.lower(),>,-1e-14);
        ARIADNE_TEST_COMPARE(sinx.upper(),>,0.0);
        ARIADNE_TEST_COMPARE(sinx.upper(),<,+1e-14);
        ARIADNE_TEST_ASSERT(sinx.lower()<sinx.upper());
    }

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

    {
        Interval x(7.0685834705770345);
        Interval sinx=sin(x);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_COMPARE(sinx.lower(),<,0.7071067811866);
        ARIADNE_TEST_COMPARE(sinx.upper(),>,0.7071067811865);
        ARIADNE_TEST_ASSERT(sinx.lower()<sinx.upper());
    }

  

}


int main() {
    TestInterval().test();
   
    cerr << "INCOMPLETE ";
    return ARIADNE_TEST_FAILURES;
}

