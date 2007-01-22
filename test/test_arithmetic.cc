/***************************************************************************
 *            test_arithmetic.cc
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

#include <fstream>
#include <sstream>
#include <iomanip>

#include <gmpxx.h>
#include <mpfr.h>
#include <boost/numeric/interval.hpp>

#include "numeric/interval.h"
#include "numeric/float64.h"
#include "numeric/mpfloat.h"
#include "numeric/rational.h"
#include "numeric/arithmetic.h"

#include "test.h"

using namespace std;
using namespace Ariadne::Numeric;
using Ariadne::name;
using Ariadne::Rational;
using Ariadne::MPFloat;
using Ariadne::Float64;

int test_boost_rounding(); 

template<class R> int test_class();
template<class R> int test_conversion();
template<class R> int test_stream();
template<> int test_stream<Rational>();
template<class R> int test_comparison();
template<class R> int test_arithmetic();
template<> int test_arithmetic<Rational>();


int main() {

  cout << setprecision(20);
  mpf_set_default_prec (8);

  test_class<Rational>();
  test_class<MPFloat>();
  test_class<Float64>();
  cout << endl;
  
  test_conversion<MPFloat>();
  test_conversion<Float64>();
  cout << endl;

  test_stream<Rational>();
  test_stream<MPFloat>();
  test_stream<Float64>();
  cout << endl;
  
  test_comparison<Rational>();
  test_comparison<MPFloat>();
  test_comparison<Float64>();
  cout << endl;
  
  test_arithmetic<Rational>();
  test_arithmetic<MPFloat>();
  test_arithmetic<Float64>();
 
  test_boost_rounding();
  
  return 0;
}

int
test_boost_rounding() 
{
  cout << __PRETTY_FUNCTION__ << endl;

  double x=1;
  double y=3;
  double zl,zu;
  
  { 
    boost::numeric::interval_lib::rounded_arith_std<double> rnd;
    zl=rnd.div_down(x,y);
    zu=rnd.div_up(x,y);
  }
  cout << zl << " <= " << x << "/" << y << " <= " << zu << endl;
  if(!(zl<zu)) {
    cerr << "Warning: boost::numeric::interval_lib::rounded_arith_std<double> does not round correctly\n";
  }
  
  return 0;
}


template<class R>
int
test_class()
{
  cout << __PRETTY_FUNCTION__ << endl;
  // Construct from an int
  R f1(2);
  assert(f1==2);
  // Construct from a double
  R f2(1.25);
  assert(f2==1.25);
  // Copy constructor
  R f3(f2);
  assert(f3==f2);
  
  // Assign from an int
  f1=3;
  assert(f1==3);
  // Assign from a double
  f2=2.25;
  assert(f2==2.25);
  // Copy assignment
  f3=f2;
  assert(f3==f2);

  return 0;
}

template<class R>
int
test_conversion()
{
  cout << __PRETTY_FUNCTION__ << endl;

  R f1=3;
  R f2=2.25;
  
  // Convert to a rational
  Rational q(f1);
  assert(q==3);
  // Assign to a rational
  q=f2;
  assert(q==2.25);
  
  // Convert from a rational
  q=Rational(1,3);
  f1=conv_down<R>(q);
  f2=conv_up<R>(q);
  cout << f1 << " <= " << q << " <= " << f2 << endl;
  assert(f1<=q); assert(f2>=q); assert(f1<f2);
  
  // Convert from a negative rational
  q=Rational(-2,5);
  f1=conv_down<R>(q);
  f2=conv_up<R>(q);
  cout << f1 << " <= " << q << " <= " << f2 << endl;
  assert(f1<=q); assert(f2>=q); assert(f1<f2);
  
  return 0;
};

template<class R>
int
test_stream()
{
  cout << __PRETTY_FUNCTION__ << endl;

  stringstream ss("1.25 -2.25 42 2.35e1");
  R f1,f2,f3,f4;
  ss >> f1;
  assert(f1==1.25);
  ss >> f2;
  assert(f2==-2.25);
  ss >> f3;
  assert(f3==42.0);
  ss >> f4;
  assert(f4==23.5);
  
  return 0;
}

template<>
int
test_stream<Rational>()
{
  cout << __PRETTY_FUNCTION__ << endl;

  typedef Rational R;
  
  stringstream ss("42 -23/5 8/2 1.25");
  R f1,f2,f3,f4;
  ss >> f1;
  cout << "f1 = " << f3 << endl;
  assert(f1==42);
  ss >> f2;
  cout << "f2 = " << f3 << endl;
  assert(f2==Rational(-23,5));
  ss >> f3;
  cout << "f3 = " << f3 << endl;
  assert(f3==Rational(8,2));
  assert(Rational(12,3)==Rational(8,2));
  assert(f3==Rational(4));
  ss >> f4;
  cout << "f4 = " << f4 << endl;
  if(f4!=1.25) {
    std::cerr << "Warning: Rational class cannot handle decimal input\n";
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
  assert(!(f1==f2)); assert(f1!=f2); 
  assert(!(f1<=f2)); assert(f1> f2);
  assert(f1>=f2); assert(!(f1< f2));
  
  // Test comparison of two different numbers
  assert(f1==f4); assert(!(f1!=f4));
  assert(f1<=f4); assert(!(f1> f4));
  assert(f1>=f4); assert(!(f1< f4));
  
  // Test comparison with in integer
  int i2=1;
  assert(!(f1==i2)); assert(f1!=i2); 
  assert(!(f1<=i2)); assert(f1> i2);
  assert(f1>=i2); assert(!(f1< i2));
  
  int i1=1;
  assert(!(i1==f2)); assert(i1!=f2); 
  assert(!(i1<=f2)); assert(i1> f2);
  assert(i1>=f2); assert(!(i1< f2));
  
  // Test comparison with a double
  double x2=1.0;
  assert(!(f1==x2)); assert(f1!=x2); 
  assert(!(f1<=x2)); assert(f1> x2);
  assert(f1>=x2); assert(!(f1< x2));
  
  double x1=1.0;
  assert(!(x1==f2)); assert(x1!=f2); 
  assert(!(x1<=f2)); assert(x1> f2);
  assert(x1>=f2); assert(!(x1< f2));
  
  // Test comparison with a rational
  Rational q2=1;
  assert(!(f1==q2)); assert(f1!=q2); 
  assert(!(f1<=q2)); assert(f1> q2);
  assert(f1>=q2); assert(!(f1< q2));
  
  //Rational q1=Rational(-5,4);
  //assert(q1==f2)); assert(!(q1!=f2)); 
  //assert(q1<=f2)); assert(!(q1> f2));
  //assert(!(q1>=f2)); assert(q1< f2);
  
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
  f4=min_exact(f1,f2);
  cout << "min(" << f1 << "," << f2 << ") = " << f4 << endl; 
  assert(f4==f1);
  f4=min_exact(f1,f3);
  cout << "min(" << f1 << "," << f3 << ") = " << f4 << endl; 
  assert(f4==f3);
  
  // Maximum (this should always remain exact)
  f4=max_exact(f1,f2);
  cout << "max(" << f1 << "," << f2 << ") = " << f4 << endl; 
  assert(f4==f2);
  f4=max_exact(f1,f3);
  cout << "max(" << f1 << "," << f3 << ") = " << f4 << endl; 
  assert(f4==f1);
  
  // Absolute value (this should always remain exact)
  f4=abs_exact(f1);
  cout << "abs(" << f1 << ") = " << f4 << endl; 
  assert(f4==1.25);
  f5=abs_exact(f3);
  cout << "abs(" << f3 << ") = " << f5 << endl; 
  assert(f5==3.25);
   
  // Median (this should remain exact here)
  f3=med_approx(f1,f2);
  cout << f1 << " <= ( " << f1 << " + " << f2 << " ) / 2 <= " << f2 << endl;
  assert(f1<=f3); assert(f3<=f2); 
  assert(f3==1.75);
  
  // Negation (this should always remain exact)
  f3=neg_exact(f1);
  cout << "-" << f1 << " = " << f3 << endl; 
  assert(f3==-1.25);
  
  // Addition (this should remain exact here)
  f3=add_down(f1,f2);
  f4=add_up(f1,f2);
  cout << f3 << " <= " << f1 << " + " << f2 << " <= " << f4 << endl;
  assert(f3<=3.5); assert(f4>=3.5);
  
  // Subtraction (this should remain exact here)
  f3=sub_down(f1,f2);
  f4=sub_up(f1,f2);
  cout << f3 << " <= " << f1 << " - " << f2 << " <= " << f4 << endl;
  assert(f3<=-1); assert(f4>=-1);
  
  // Multiplication (this should remain exact here)
  f3=mul_down(f1,f2);
  f4=mul_up(f1,f2);
  cout << f3 << " <= " << f1 << " * " << f2 << " <= " << f4 << endl;
  assert(f3<=2.8125); assert(f4>=2.8125);
  
  // Division (not exact; should catch errors here)
  f3=div_down(f1,f2);
  f4=div_up(f1,f2);
  f5=div_approx(f1,f2);
  cout << f3 << " <= " << f1 << " / " << f2 << " <= " << f4 << endl;
  assert(f3<f4); assert(f3<=f5); assert(f4>=f5);
  assert(Rational(f3)<=Rational(5,9));
  assert(Rational(f4)>=Rational(5,9));
  
  // Check division my multipyling back
  cout << mul_down(f3,f2) << " <= (" << f1 << "/" << f2 << ")*" << f2 << " <= " << mul_up(f4,f2) << endl;
  assert(mul_down(f3,f2)<f1);
  assert(mul_up(f4,f2)>f1);
  
  // Power (not exact; should catch errors here)
  f3=pow_down(f1,3);
  f4=pow_up(f1,3);
  cout << f3 << " <= pow(" << f1 << ",3) <= " << f4 << endl;
  assert(f3<=1.953125); assert(f4>=1.953125);
  
  f3=pow_down(f1,-2);
  f4=pow_up(f1,-2);
  cout << f3 << " <= pow(" << f1 << ",-2) <= " << f4 << endl;
  assert(Rational(f3)<Rational(16,25)); assert(Rational(f4)>Rational(16,25));
  
  // Floor and ceiling
  f2=R(-3.25); f3=R(-2);
  
  assert(Ariadne::Numeric::floor(f1)==1); assert(Ariadne::Numeric::ceil(f1)==2);
  assert(floor(f2)==-4); assert(ceil(f2)==-3);
  assert(floor(f3)==-2); assert(ceil(f3)==-2);
  
  // Conversion to integer types
  int i3,i4;
  i3=int_down<int>(f1);
  i4=int_up<int>(f1);
  cout << i3 << " < " << f1 << " < " << i4 << endl;
  assert(i3==1); assert(i4==2);
  i3=int_down<int>(f2);
  i4=int_up<int>(f2);
  cout << i3 << " < " << f2 << " < " << i4 << endl;
  assert(i3==-4); assert(i4==-3);
  
  // Check interval conversions
  R z(0); R o(1); R t(3);
  Interval<R> io(1); Interval<R> it(3);
  
  Interval<R> iao=(o/t)*t;
  cout << iao << endl;
  cout << div_down(o,t) << " <= 1/3 <= " << div_up(o,t) << endl;
  assert(div_down(o,t)<div_up(o,t));
  cout << o/t << endl;
  assert(contains_value(iao,o));
  Interval<R> iaz=iao-io;
  Interval<R> iz=R(1)-R(1);
  cout << iaz << endl;
  assert(contains_value(iaz,z)); 
  assert(!bool(!subset(iz,iaz)));
  cout << endl;

  // The following line should not compile
  // f5=f1+f2;
  
  return 0;
}


template<>
int
test_arithmetic<Rational>()
{
  cout << __PRETTY_FUNCTION__ << endl;
  
  typedef Rational R;

  R f1(1.25);
  R f2(2.25);
  R f3;

  f3=add(f1,f2);
  cout << f1 << " + " << f2 << " = " << f3 << endl;
  assert(f3==R(7,2));
  f3=sub(f1,f2);
  cout << f1 << " - " << f2 << " = " << f3 << endl;
  assert(f3==R(-1,1));
  f3=mul(f1,f2);
  cout << f1 << " * " << f2 << " = " << f3 << endl;
  assert(f3==R(45,16));
  f3=div(f1,f2);
  cout << f1 << " / " << f2 << " = " << f3 << endl;
  assert(f3==R(5,9));

  cout << endl;
  
  return 0;
}
