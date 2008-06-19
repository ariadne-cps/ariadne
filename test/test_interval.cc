/***************************************************************************
 *            test_interval.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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

#include "numeric/rational.h"
#include "numeric/float.h"
#include "numeric/interval.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne;
using namespace std;

template class Interval<Rational>;
#ifdef ENABLE_FLOAT64
template class Interval<Float64>;
#endif
#ifdef ENABLE_FLOATMP
template class Interval<FloatMP>;
#endif

template<class R> void test_interval();
template<> void test_interval<Rational>();



template<class R> class TestInterval;

template<class T> 
class TestInterval< Float<T> >
{
 public:
  void test();
 private:
  void test_concept();
};

template<> 
class TestInterval< Rational >
{
 public:
  void test();
 private:
  void test_concept();
};


template<class T> void
TestInterval< Float<T> >::test()
{
  test_concept();
}

template<class T> void
TestInterval< Float<T> >::test_concept()
{
  typedef Float<T> R;
  typedef Interval< Float<T> > I;

  int n=0; 
  uint m=0;
  double d=0;
  Integer z=0;
  Rational q=0;
  Float<T> x=0;
  Interval< Float<T> > i=0;
 
  // Constructors
  i=I(); i=I(n); i=I(m); i=I(d); i=I(z); i=I(x); i=I(q); i=I(i);
  i=I(n,n); i=I(m,m); i=I(d,d); i=I(z,z); i=I(x,x); i=I(q,q);
  i=I(n,m); i=I(m,d); i=I(d,n); 
  i=I(z,x); i=I(x,q); i=I(q,z); 
  // Assignment
  i=n; i=m; i=d; i=z; i=x; i=q; i=i;

  // Exact operations
  i=abs(i); i=pos(i); i=neg(i);

  // Arithmetic
  i=add(x,x); i=add(x,i); i=add(i,x); i=add(i,i);
  i=sub(x,x); i=sub(x,i); i=sub(i,x); i=sub(i,i);
  i=mul(x,x); i=mul(x,i); i=mul(i,x); i=mul(i,i);
  i=div(x,x); i=div(x,i); i=div(i,x); i=div(i,i);
  //i=pow(x,n); i=pow(x,n);
  //i=pow(x,m); i=pow(x,m);

  // Transcendental functions
  i=sqrt(x); i=sqrt(i);
  i=exp(x); i=exp(i);
  i=log(x); i=log(i);
  i=sin(x); i=sin(i);
  i=cos(x); i=cos(i);
  i=tan(x); i=tan(i);
  i=asin(x); i=asin(i);
  i=acos(x); i=acos(i);
  i=atan(x); i=atan(i);
  
  /*
  i=sinh(x); i=sinh(i);
  i=cosh(x); i=cosh(i);
  i=tanh(x); i=tanh(i);
  i=asinh(x); i=asinh(i);
  i=acosh(x); i=acosh(i);
  i=atanh(x); i=atanh(i);
  */
}


template<class R>
void
test_interval()
{
  cout << __PRETTY_FUNCTION__ << endl;
  typedef Interval<R> I;


  R zero=0;
  
  // Construct from pair
  Interval<R> ivld1(R(1.125),R(2.25));
  ARIADNE_TEST_ASSERT(ivld1.lower()==1.125); ARIADNE_TEST_ASSERT(ivld1.upper()==2.25);
  // Default constructor
  Interval<R> ivld2;
  if(ivld2.empty()) { 
    cerr << "Warning: Interval<R> default constructor returns an empty set\n"; 
  } else {
    ARIADNE_TEST_ASSERT((bool)(ivld2==Interval<R>(zero,zero)));
  }
  // Constructor with approximations
  Interval<R> ivld3(Rational(21,10),Rational(16,5));
  cout<<ivld3<<std::endl;
  ARIADNE_TEST_ASSERT(ivld3.lower()<Rational(21,10));
  ARIADNE_TEST_ASSERT(ivld3.upper()>Rational(16,5));
  // Constructor from approximate values
  Interval<R> ivld4(2.1,3.2);
  ARIADNE_TEST_ASSERT(ivld4.lower()<=2.1);
  ARIADNE_TEST_ASSERT(ivld4.upper()>=3.2);
  // Approximate constructor from a single value
  Interval<R> ivld5(Rational(1,3));
  ARIADNE_TEST_ASSERT(ivld5.lower()<Rational(1,3));
  ARIADNE_TEST_ASSERT(ivld5.upper()>Rational(1,3));
  // Exact constructor from a single value
  Interval<R> ivld6(R(1.25));
  ARIADNE_TEST_ASSERT(ivld6.lower()==R(1.25));
  ARIADNE_TEST_ASSERT(ivld6.upper()==R(1.25));
  
  Interval<R> ivlq1(1.1,2.2);
  Interval<R> ivlq2(1.125,1.125);
  Interval<R> ivlq3(2.125,3.25);
  
  R f0=-2.25;
  R f1=1.5;
  R f2=2.25;
  R f3=3.125;
  R f4=4.0625;
  Interval<R> ivlf1;
  Interval<R> ivlf2(f1,f2);
  Interval<R> ivlf3(f3,f4);
  
  R r1=1.5;
  R r2=2.25;
  R r3=3.125;
  R r4=4.0625;
  Interval<R> ivlr1;
  Interval<R> ivlr2(r1,r2);
  Interval<R> ivlr3(r3,r4);
  Interval<R> ivlr4(r1,r4);

  cout << "ivlr1=" << ivlr1 << ", ivlr2=" << ivlr2 << ", ivlr3=" << ivlr3 
       << ", ivlr4=" << ivlr4 << endl;
  
  Interval<R> ivlf0(f2,f1);
  cout << "ivlf0=" << ivlf0 << endl;
  cout << "ivlf0.empty()=" << ivlf0.empty() << ", ivlf1.empty()=" << ivlf1.empty() << endl;
  
  ivlf1=Interval<R>(f1,f3);
  ivlf2=Interval<R>(f0,f1);
  ivlf3=Interval<R>(f2,f4);
  cout << "ivlf1=" << ivlf1 << ", ivlf2=" << ivlf2 << ", ivlf3=" << ivlf3 << endl;
  cout << "max(ivlf1,ivlf3)=" << max(ivlf1,ivlf3) << endl;
  cout << "min(ivlf1,ivlf3)=" << min(ivlf1,ivlf3) << endl;
  cout << "lower(" << ivlf1 << ")=" << ivlf1.lower() << endl;
  cout << "upper(" << ivlf1 << ")=" << ivlf1.upper() << endl;
  cout << "midpoint(" << ivlf1 << ")=" << ivlf1.midpoint() << endl;
  cout << "radius(" << ivlf2 << ")=" << radius(ivlf2) << endl;
  cout << "width(" << ivlf2 << ")=" << width(ivlf2) << endl;
  ARIADNE_TEST_ASSERT(ivlf1.midpoint()==R(2.3125));
  
  try {
    string input("[1.125,2.25] ");
    stringstream iss(input);
    
    iss >> ivld2;
    cout << "ivld1=" << ivld1 << "  ivld2=" << ivld2 << endl;
    ARIADNE_TEST_ASSERT(equal(ivld1,ivld2));
    ARIADNE_TEST_ASSERT(ivld1.lower()==ivld2.lower() && ivld1.upper()==ivld2.upper());
    
    ARIADNE_TEST_ASSERT(equal(ivld1,ivld2));
    
    ARIADNE_TEST_ASSERT(indeterminate(ivld1==ivld2));
    
    Interval<R>& ivlf1ref=ivlf1;
    ivlf1ref=Interval<R>(5.25,7.375);
    cout << "ivlf1ref=" << ivlf1ref << endl;
    ARIADNE_TEST_ASSERT(ivlf1ref.lower()==R(5.25));
    
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
    Interval<R> ivlo(1.0);
    Interval<R> ivlt(3.0);
    Interval<R> ivlodt=ivlo/ivlt;
    Interval<R> ivloa=ivlodt*ivlt;
    cout << ivlo << " / " << ivlt << " = " << ivlodt << endl;
    cout << ivlo << " in " << Interval<R>(ivlodt * ivlt) << endl;
    ARIADNE_TEST_ASSERT(ivlodt.lower() < ivlodt.upper());
    ARIADNE_TEST_ASSERT(ivloa.lower() < ivlo.lower() && ivloa.upper() > ivlo.upper());
    
    

    ivlf1=Interval<R>(-13,-7);
    ivlf2=Interval<R>(-3,2);
    ivlf3=Interval<R>(5,11);
    
    // Fixme: The explicit cast should not be needed in the output
    cout << ivlf1 << " * " << ivlf1 << " = " << Interval<R>(ivlf1*ivlf1) << endl;
    cout << ivlf1 << " * " << ivlf2 << " = " << Interval<R>(ivlf1*ivlf2) << endl;
    cout << ivlf1 << " * " << ivlf3 << " = " << Interval<R>(ivlf1*ivlf3) << endl;
    cout << ivlf2 << " * " << ivlf1 << " = " << Interval<R>(ivlf2*ivlf1) << endl;
    cout << ivlf2 << " * " << ivlf2 << " = " << Interval<R>(ivlf2*ivlf2) << endl;
    cout << ivlf2 << " * " << ivlf3 << " = " << Interval<R>(ivlf2*ivlf3) << endl;
    cout << ivlf3 << " * " << ivlf1 << " = " << Interval<R>(ivlf3*ivlf1) << endl;
    cout << ivlf3 << " * " << ivlf2 << " = " << Interval<R>(ivlf3*ivlf2) << endl;
    cout << ivlf3 << " * " << ivlf3 << " = " << Interval<R>(ivlf3*ivlf3) << endl;

    cout << ivlf1 << " / " << ivlf1 << " = " << Interval<R>(ivlf1/ivlf1) << endl;
    cout << ivlf1 << " / " << ivlf3 << " = " << Interval<R>(ivlf1/ivlf3) << endl;
    cout << ivlf2 << " / " << ivlf1 << " = " << Interval<R>(ivlf2/ivlf1) << endl;
    cout << ivlf2 << " / " << ivlf3 << " = " << Interval<R>(ivlf2/ivlf3) << endl;
    cout << ivlf3 << " / " << ivlf1 << " = " << Interval<R>(ivlf3/ivlf1) << endl;
    cout << ivlf3 << " / " << ivlf3 << " = " << Interval<R>(ivlf3/ivlf3) << endl;


  }
  catch(exception& e) {
    cout << "EXCEPTION " << e.what() << "\n";
    throw e;
  }
  
  // e=2.718281828459045235360287471352662497[78]
  // 2*pi = 6.283185307179586476925286766559005768[3,5]

  // Regression test; fails dramatically on certain types of rounding
  {
    Interval<R> one(1.0);
    Interval<R> expone=exp(one);
    ARIADNE_TEST_PRINT(one);
    ARIADNE_TEST_COMPARE(expone.lower(),<,2.71828182845905);
    ARIADNE_TEST_COMPARE(expone.lower(),>,2.71828182845903);
    ARIADNE_TEST_COMPARE(expone.upper(),>,2.71828182845904);
    ARIADNE_TEST_COMPARE(expone.upper(),<,2.71828182845906);
    ARIADNE_TEST_ASSERT(expone.lower()<expone.upper());
  }

  {
    Interval<R> e(Rational(2.71828182845904),Rational(2.71828182845905));
    Interval<R> loge=log(e);
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_COMPARE(loge.lower(),<,1);
    ARIADNE_TEST_COMPARE(loge.lower(),>,0.9999999999998);
    ARIADNE_TEST_COMPARE(loge.upper(),>,1);
    ARIADNE_TEST_COMPARE(loge.upper(),<,1.000000000002);
 }


  // Regression test; fails dramatically on certain types of rounding
  {
    Interval<R> x(6.283185307179586,6.283185307179587);
    Interval<R> sinx=sin(x);
    ARIADNE_TEST_PRINT(x);
    ARIADNE_TEST_COMPARE(sinx.lower(),<,0.0);
    ARIADNE_TEST_COMPARE(sinx.lower(),>,-2e-15);
    ARIADNE_TEST_COMPARE(sinx.upper(),>,0.0);
    ARIADNE_TEST_COMPARE(sinx.upper(),<,+2e-15);
    ARIADNE_TEST_ASSERT(sinx.lower()<sinx.upper());
  }

  // Regression test; fails dramatically on certain types of rounding
  {
    Interval<R> x(1.5707963267948966,1.5707963267948968);
    Interval<R> cosx=cos(x);
    ARIADNE_TEST_PRINT(x);
    ARIADNE_TEST_COMPARE(cosx.lower(),<,0.0);
    ARIADNE_TEST_COMPARE(cosx.lower(),>,-2e-15);
    ARIADNE_TEST_COMPARE(cosx.upper(),>,0.0);
    ARIADNE_TEST_COMPARE(cosx.upper(),<,+2e-15);
    ARIADNE_TEST_ASSERT(cosx.lower()<cosx.upper());
  }

  {
    Interval<R> x(7.0685834705770345);
    Interval<R> sinx=sin(x);
    ARIADNE_TEST_PRINT(x);
    ARIADNE_TEST_COMPARE(sinx.lower(),<,0.7071067811866);
    ARIADNE_TEST_COMPARE(sinx.upper(),>,0.7071067811865);
    ARIADNE_TEST_ASSERT(sinx.lower()<sinx.upper());
  }
}

template<>
void
test_interval<Rational>()
{
  typedef Rational R;
  
  cout << "test_interval<" << name<R>() << ">" << endl;
  
  Interval<R> ivld1(R(1.125),R(2.25));
  cout << "ivld1=" << ivld1 << endl;
  Interval<R> ivld2;
  cout << "ivld2=" << ivld2 << endl;
  Interval<R> ivld3(2.1,3.2);
  cout << "ivld3=" << ivld3 << endl;

  Interval<R> ivlq1(1.1,2.2);
  Interval<R> ivlq2(1.125,1.125);
  Interval<R> ivlq3(2.125,3.25);
  
  R f0=-2.25;
  R f1=1.5;
  R f2=2.25;
  R f3=3.125;
  R f4=4.0625;
  Interval<R> ivlf1;
  Interval<R> ivlf2(f1,f2);
  Interval<R> ivlf3(f3,f4);
  
  R r1=1.5;
  R r2=2.25;
  R r3=3.125;
  R r4=4.0625;
  Interval<R> ivlr1;
  Interval<R> ivlr2(r1,r2);
  Interval<R> ivlr3(r3,r4);
  Interval<R> ivlr4(r1,r4);

  cout << "ivlr1=" << ivlr1 << ", ivlr2=" << ivlr2 << ", ivlr3=" << ivlr3 
       << ", ivlr4=" << ivlr4 << endl;
  
  Interval<R> ivlf0(f2,f1);
  cout << "ivlf0=" << ivlf0 << endl;
  cout << "ivlf0.empty()=" << ivlf0.empty() << ", ivlf1.empty()=" << ivlf1.empty() << endl;
  
  ivlf1=Interval<R>(f1,f3);
  ivlf2=Interval<R>(f0,f1);
  ivlf3=Interval<R>(f2,f4);
  cout << "min(ivlf1,ivlf2)=" << min(ivlf1,ivlf3) << endl;
  cout << "max(ivlf1,ivlf2)=" << max(ivlf1,ivlf3) << endl;
  cout << "abs(" << ivlf2 << ")=" << abs(ivlf2) << endl;
   
  try {
    string input("[5/4,9/4] ");
    stringstream iss(input);
    
    ivld1=Interval<R>(1.25,2.25);
    iss >> ivld2;
    cout << "ivld1=" << ivld1 << "  ivld2=" << ivld2 << endl;
    ARIADNE_TEST_ASSERT(equal(ivld1,ivld2));
    ARIADNE_TEST_ASSERT(indeterminate(ivld1==ivld2));
    
    Interval<R>& ivlf1ref=ivlf1;
    ivlf1ref=Interval<R>(5.25,7.375);
    cout << "ivlf1ref=" << ivlf1ref << endl;
    ARIADNE_TEST_ASSERT(ivlf1ref.lower()==R(5.25));
    
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

  }
  catch(exception& e) {
    cout << "EXCEPTION " << e.what() << "\n";
    throw e;
  }
  
}


int main() {
  initialise();

  cout << setprecision(20);
  mpf_set_default_prec (8);

  test_interval<Rational>();

#ifdef ENABLE_FLOAT64
  test_interval<Float64>();
  TestInterval<Float64>().test();
#endif
   
#ifdef ENABLE_FLOATMP
  test_interval<FloatMP>();
  TestInterval<FloatMP>().test();
#endif
   
  cerr << "INCOMPLETE ";
  return ARIADNE_TEST_FAILURES;
}

