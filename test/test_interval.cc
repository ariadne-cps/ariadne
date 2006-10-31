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

#include "numeric/interval.h"

#include "numeric/float64.h"
#include "numeric/mpfloat.h"
#include "numeric/rational.h"

#include "test.h"

using namespace Ariadne;
using std::exception;
using std::ostream; using std::ofstream; using std::stringstream;
using std::cout; using std::cerr;
using std::endl; using std::flush; using std::setprecision;
using std::string;

template class Interval<Float64>;
template class Interval<MPFloat>;
template class Interval<Rational>;

template<class R> int test_interval();
template<> int test_interval<Rational>();

int main() {
  cout << setprecision(20);
  mpf_set_default_prec (8);

  test_interval<Float64>();
  test_interval<MPFloat>();
  test_interval<Rational>();
  
  cerr << "INCOMPLETE ";
  return 0;
}



template<class R>
int
test_interval()
{
  cout << __PRETTY_FUNCTION__ << endl;
  
  // Construct from pair
  Interval<R> ivld1(R(1.125),R(2.25));
  assert(ivld1.lower()==1.125); assert(ivld1.upper()==2.25);
  // Default constructor
  Interval<R> ivld2;
  assert(ivld1.empty());
  // Constructor with approximations
  Interval<R> ivld3(Rational(21,10),Rational(16,5));
  assert(ivld3.lower()<Rational(21,10));
  assert(ivld3.upper()>Rational(16,5));
  // Constructor from approximate values
  Interval<R> ivld4(2.1,3.2);
  assert(ivld4.lower()<=2.1);
  assert(ivld4.upper()>=3.2);
  // Approximate constructor from a single value
  Interval<R> ivld5(Rational(1,3));
  assert(ivld5.lower()<Rational(1,3));
  assert(ivld5.upper()>Rational(1,3));
  // Exact constructor from a single value
  Interval<R> ivld6(R(1.25));
  assert(ivld6.lower()==R(1.25));
  assert(ivld6.upper()==R(1.25));
  
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
  cout << "abs(" << ivlf2 << ")=" << abs(ivlf2) << endl;
  cout << "centre(" << ivlf1 << ")=" << ivlf1.centre() << endl;
  assert(ivlf1.centre()==R(2.3125));
  
  try {
    string input("[1.125,2.25] ");
    stringstream iss(input);
    
    iss >> ivld2;
    cout << "ivld1=" << ivld1 << "  ivld2=" << ivld2 << endl;
    assert(equal(ivld1,ivld2));
    assert(ivld1.lower()==ivld2.lower() && ivld1.upper()==ivld2.upper());
    
    assert(equal(ivld1,ivld2));
    
    assert(indeterminate(ivld1==ivld2));
    
    Interval<R>& ivlf1ref=ivlf1;
    ivlf1ref=Interval<R>(5.25,7.375);
    cout << "ivlf1ref=" << ivlf1ref << endl;
    assert(ivlf1ref.lower()==R(5.25));
    
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

    // ensure proper rounding using assertions
    Interval<R> ivlo(1.0);
    Interval<R> ivlt(3.0);
    Interval<R> ivlodt=ivlo/ivlt;
    Interval<R> ivloa=ivlodt*ivlt;
    cout << ivlo << " / " << ivlt << " = " << ivlodt << endl;
    cout << ivlo << " in " << ivlodt * ivlt << endl;
    assert(ivlodt.lower()< ivlodt.upper());
    assert(ivloa.lower() < ivlo.lower() && ivloa.upper() > ivlo.upper());
    
    

    ivlf1=Interval<R>(-13,-7);
    ivlf2=Interval<R>(-3,2);
    ivlf3=Interval<R>(5,11);
    
    cout << ivlf1 << " * " << ivlf1 << " = " << ivlf1*ivlf1 << endl;
    cout << ivlf1 << " * " << ivlf2 << " = " << ivlf1*ivlf2 << endl;
    cout << ivlf1 << " * " << ivlf3 << " = " << ivlf1*ivlf3 << endl;
    cout << ivlf2 << " * " << ivlf1 << " = " << ivlf2*ivlf1 << endl;
    cout << ivlf2 << " * " << ivlf2 << " = " << ivlf2*ivlf2 << endl;
    cout << ivlf2 << " * " << ivlf3 << " = " << ivlf2*ivlf3 << endl;
    cout << ivlf3 << " * " << ivlf1 << " = " << ivlf3*ivlf1 << endl;
    cout << ivlf3 << " * " << ivlf2 << " = " << ivlf3*ivlf2 << endl;
    cout << ivlf3 << " * " << ivlf3 << " = " << ivlf3*ivlf3 << endl;

    cout << ivlf1 << " / " << ivlf1 << " = " << ivlf1/ivlf1 << endl;
    cout << ivlf1 << " / " << ivlf3 << " = " << ivlf1/ivlf3 << endl;
    cout << ivlf2 << " / " << ivlf1 << " = " << ivlf2/ivlf1 << endl;
    cout << ivlf2 << " / " << ivlf3 << " = " << ivlf2/ivlf3 << endl;
    cout << ivlf3 << " / " << ivlf1 << " = " << ivlf3/ivlf1 << endl;
    cout << ivlf3 << " / " << ivlf3 << " = " << ivlf3/ivlf3 << endl;


  }
  catch(exception& e) {
    cout << "EXCEPTION " << e.what() << "\n";
    throw e;
  }
  
  ivlf1=Interval<R>(R(1),R(1));
  cout << "one=" << ivlf1 << endl;
  ivlf2=exp(ivlf1);
  cout << "exp(" << ivlf1 << ")=" << ivlf2 << endl;
  ivlf3=log(ivlf1);
  cout << "log(" << ivlf2 << ")=" << ivlf3 << endl;
 
  return 0;
}

template<>
int
test_interval<Rational>()
{
  typedef Rational R;
  
  cout << "test_interval<" << name<R>() << ">" << endl;
  
  Interval<R> ivld1(R(1.125),R(2.25));
  Interval<R> ivld2;
  Interval<R> ivld3(2.1,3.2);
  
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
    assert(equal(ivld1,ivld2));
    assert(indeterminate(ivld1==ivld2));
    
    Interval<R>& ivlf1ref=ivlf1;
    ivlf1ref=Interval<R>(5.25,7.375);
    cout << "ivlf1ref=" << ivlf1ref << endl;
    assert(ivlf1ref.lower()==R(5.25));
    
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
  
  return 0;
}
