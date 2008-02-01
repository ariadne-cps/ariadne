/***************************************************************************
 *            test_float.cc
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
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

#include "numeric/integer.h"
#include "numeric/rational.h"
#include "numeric/float.h"
#include "numeric/interval.h"

#include "test.h"

using namespace std;
using namespace Ariadne;
using namespace Ariadne::Numeric;


template<class R> int test_float();
template<class R> int test_class();
template<class R> int test_conversion();
template<class R> int test_stream();
template<class R> int test_comparison();
template<class R> int test_arithmetic();
template<class R> int test_function();




int main() {

  cout << setprecision(20);
  mpf_set_default_prec (8);

#ifdef ENABLE_FLOAT64
  test_float<Float64>();
  cout << endl;
#endif

#ifdef ENABLE_FLOATMP
  test_float<FloatMP>();
  cout << endl;
#endif

  return 0;
}


template<class R>
int
test_float() {
  test_class<R>();
  test_conversion<R>();
  test_stream<R>();
  test_comparison<R>();
  test_arithmetic<R>();
  test_function<R>();
  return 0;
}

 




template<class R>
int
test_inverse_pair(
  std::string name, 
  R(*fnl)(const R&),
  R(*fnu)(const R&),
  R(*ifnl)(const R&),
  R(*ifnu)(const R&) )
{
  cout << name << endl;
  R o=1;
  R iml=fnl(o);
  R imu=fnu(o);
  R ol=ifnl(iml);
  R ou=ifnu(imu);
  cout << iml << " <= " << name << "(1) <= " << imu << endl;
  cout << ol << " <=    1   <= " << ou << endl;
  ARIADNE_TEST_ASSERT(iml<imu);
  ARIADNE_TEST_ASSERT(ol<=o);
  ARIADNE_TEST_ASSERT(o<=ou);
  return 0;
}
  
 

template<class R>
int
test_class()
{
  cout << __PRETTY_FUNCTION__ << endl;
  // Construct from an int
  R f1(2);
  ARIADNE_TEST_ASSERT(f1==2);
  // Construct from a double
  R f2(1.25);
  ARIADNE_TEST_ASSERT(f2==1.25);
  // Copy constructor
  R f3(f2);
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

  return 0;
}


template<class R>
int
test_conversion()
{
  cout << __PRETTY_FUNCTION__ << endl;
 
  // Convert from integers
  int n;
  n=std::numeric_limits<int>::min();
  ARIADNE_TEST_ASSERT(Rational(R(n))==n);
  n=std::numeric_limits<int>::max();
  ARIADNE_TEST_ASSERT(Rational(R(n))==n);
  n=std::numeric_limits<unsigned int>::max();
  ARIADNE_TEST_ASSERT(Rational(R(n))==n);

  // Conversion from long integers cannot be performed exactly, so is banned

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
  set_(f1,q,round_down); cout << f1 << endl;
  set_(f2,q,round_up);
  cout << f1 << " <= " << q << " <= " << f2 << endl;
  ARIADNE_TEST_ASSERT(f1<=q); ARIADNE_TEST_ASSERT(f2>=q); ARIADNE_TEST_ASSERT(f1<f2);
  
  return 0;
};


template<class R>
int
test_stream()
{
  cout << __PRETTY_FUNCTION__ << endl;

  stringstream ss("1.25 -2.25 42:2 375e1 2.35e1");
  R f1,f2,f3,f4,f5;
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
      throw std::runtime_error("WARNING: cannot create float from string literal in exponential form 2.375e1");
    }
  } 
  catch(std::exception& e) {
    cerr << e.what() << endl;;
  }
  
  try {
    ss >> f5;
    cout << f5 << endl;
    if(f4!=23.5) {
      throw std::runtime_error("WARNING: cannot create float from string literal in exponential form 2.35e1");
    }
  } 
  catch(std::exception& e) {
    cerr << e.what() << endl;
  }
  
  return 0;
}


template<class R>
int
test_comparison()
{
  cout << __PRETTY_FUNCTION__ << endl;
  
  // Test comparison of two equal numbers
  R f1(1.25); R f2(-1.25); R f3(-2.25); R f4(1.25);
  ARIADNE_TEST_ASSERT(!(f1==f2)); ARIADNE_TEST_ASSERT(f1!=f2); 
  ARIADNE_TEST_ASSERT(!(f1<=f2)); ARIADNE_TEST_ASSERT(f1> f2);
  ARIADNE_TEST_ASSERT(f1>=f2); ARIADNE_TEST_ASSERT(!(f1< f2));
  
  // Test comparison of two different numbers
  ARIADNE_TEST_ASSERT(f1==f4); ARIADNE_TEST_ASSERT(!(f1!=f4));
  ARIADNE_TEST_ASSERT(f1<=f4); ARIADNE_TEST_ASSERT(!(f1> f4));
  ARIADNE_TEST_ASSERT(f1>=f4); ARIADNE_TEST_ASSERT(!(f1< f4));
  
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
  Rational q2=1;
  ARIADNE_TEST_ASSERT(!(f1==q2)); ARIADNE_TEST_ASSERT(f1!=q2); 
  ARIADNE_TEST_ASSERT(!(f1<=q2)); ARIADNE_TEST_ASSERT(f1> q2);
  ARIADNE_TEST_ASSERT(f1>=q2); ARIADNE_TEST_ASSERT(!(f1< q2));
  
  //Rational q1=Rational(-5,4);
  //ARIADNE_TEST_ASSERT(q1==f2)); ARIADNE_TEST_ASSERT(!(q1!=f2)); 
  //ARIADNE_TEST_ASSERT(q1<=f2)); ARIADNE_TEST_ASSERT(!(q1> f2));
  //ARIADNE_TEST_ASSERT(!(q1>=f2)); ARIADNE_TEST_ASSERT(q1< f2);
  
  return 0;
}
  
template<class R>
int
test_arithmetic()
{
  cout << __PRETTY_FUNCTION__ << endl;
  
  // Set up some variables
  R f1(1.25); R f2(2.25); R f3(-3.25); R f4; R f5;
 
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
  f3=next_down(f1);
  f4=next_up(f1);
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
  ARIADNE_TEST_ASSERT(Rational(f3)<=Rational(5,9));
  ARIADNE_TEST_ASSERT(Rational(f4)>=Rational(5,9));
  
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
  ARIADNE_TEST_ASSERT(Rational(f3)<Rational(16,25)); 
  ARIADNE_TEST_ASSERT(Rational(f4)>Rational(16,25));
  
  // Floor and ceiling
  f2=R(-3.25); f3=R(-2);
  
  ARIADNE_TEST_ASSERT(floor(f1)==1); ARIADNE_TEST_ASSERT(ceil(f1)==2);
  ARIADNE_TEST_ASSERT(floor(f2)==-4); ARIADNE_TEST_ASSERT(ceil(f2)==-3);
  ARIADNE_TEST_ASSERT(floor(f3)==-2); ARIADNE_TEST_ASSERT(ceil(f3)==-2);
  
  // Conversion to integer types
  int i3,i4;
  i3=floor(f1);
  i4=ceil(f1);
  cout << i3 << " < " << f1 << " < " << i4 << endl;
  ARIADNE_TEST_ASSERT(i3==1); ARIADNE_TEST_ASSERT(i4==2);
  i3=floor(f2);
  i4=ceil(f2);
  cout << i3 << " < " << f2 << " < " << i4 << endl;
  ARIADNE_TEST_ASSERT(i3==-4); ARIADNE_TEST_ASSERT(i4==-3);
  
  // Check interval conversions
  R z(0); R o(1); R t(3);
  Interval<R> io(1); Interval<R> it(3);
  
  cout << "o/t=" << Interval<R>(o/t) << endl;

  Interval<R> iao=(o/t)*t;
  cout << "(o/t)*t=" << iao << endl;
  cout << iao.lower() << " " << iao.upper() << endl;

  cout << div_down(o,t) << " <= 1/3 <= " << div_up(o,t) << endl;
  ARIADNE_TEST_ASSERT(div_down(o,t)<div_up(o,t));
  cout << Interval<R>(o/t) << endl; 
  ARIADNE_TEST_ASSERT(encloses(iao,o));
  Interval<R> iaz=iao-io;
  Interval<R> iz=R(1)-R(1);
  cout << iaz << endl;
  ARIADNE_TEST_ASSERT(encloses(iaz,z)); 
  ARIADNE_TEST_ASSERT(!bool(!subset(iz,iaz)));
  cout << endl;

  ARIADNE_TEST_ASSERT(med_approx(R(2),R(3))==2.5);
  ARIADNE_TEST_ASSERT(rad_up(R(2),R(3))==0.5);

  // The following line should not compile
  // f5=f1+f2;
  
  return 0;
}


template<class R>
int
test_function()
{
  cout << __PRETTY_FUNCTION__ << endl;
  
  cout << setprecision(20);
  mpf_set_default_prec (128);
  
  //test_inverse_pair("exp",&exp_down<R>,&exp_up<R>,&log_down<R>,&log_up<R>);

  //The following don't work as rounded operators not exported.
  //test_inverse_pair("sin",&sin_down<R>,&sin_up<R>,&asin_down<R>,&asin_up<R>);
  //The following don't work as acos is decreasing
  //test_inverse_pair("cos",&cos_down<R>,&cos_up<R>,&acos_down<R>,&acos_up<R>);
  //test_inverse_pair("cos",&cos_up<R>,&cos_down<R>,&acos_down<R>,&acos_up<R>);
  //test_inverse_pair("tan",&tan_down<R>,&tan_up<R>,&atan_down<R>,&atan_up<R>);
  //test_inverse_pair("sinh",&sinh_down<R>,&sinh_up<R>,&asinh_down<R>,&asinh_up<R>);
  //test_inverse_pair("cosh",&cosh_down<R>,&cosh_up<R>,&acosh_down<R>,&acosh_up<R>);
  //test_inverse_pair("tanh",&tanh_down<R>,&tanh_up<R>,&atanh_down<R>,&atanh_up<R>);
  return 0;
}

