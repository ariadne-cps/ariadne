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

#include <iostream>
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
using std::endl; using std::flush;
using std::string;

template class Interval<Float64>;
template class Interval<MPFloat>;
template class Interval<Rational>;

template<typename R> int test_interval(std::ostream&);
template<> int test_interval<Rational>(std::ostream&);

int main() {
  cout << "test_interval: " << flush;
  ofstream clog("test_interval.log");
  test_interval<Float64>(clog);
  test_interval<MPFloat>(clog);
  test_interval<Rational>(clog);
//  int dbi=test_interval<Float64>(clog);
//  int fli=test_interval<MPFloat>(clog);
//  int qi=test_interval<Rational>(clog);
  
  clog.close();
  cout << "INCOMPLETE\n";
  
  return 0;
//  return dbi || fli || qi;
}


template<typename Real>
int
test_interval(ostream& clog)
{
  clog << "Test interval<" << name<Real>() << ">" << endl;
  
  Interval<Real> ivld1(Real(1.125),Real(2.25));
  Interval<Real> ivld2;
  Interval<Real> ivld3(2.1,3.2);
  
  Interval<Real> ivlq1(1.1,2.2);
  Interval<Real> ivlq2(1.125,1.125);
  Interval<Real> ivlq3(2.125,3.25);
  
  Real f0=-2.25;
  Real f1=1.5;
  Real f2=2.25;
  Real f3=3.125;
  Real f4=4.0625;
  Interval<Real> ivlf1;
  Interval<Real> ivlf2(f1,f2);
  Interval<Real> ivlf3(f3,f4);
  
  Real r1=1.5;
  Real r2=2.25;
  Real r3=3.125;
  Real r4=4.0625;
  Interval<Real> ivlr1;
  Interval<Real> ivlr2(r1,r2);
  Interval<Real> ivlr3(r3,r4);
  Interval<Real> ivlr4(r1,r4);

  clog << "ivlr1=" << ivlr1 << ", ivlr2=" << ivlr2 << ", ivlr3=" << ivlr3 
       << ", ivlr4=" << ivlr4 << endl;
  
  Interval<Real> ivlf0(f2,f1);
  clog << "ivlf0=" << ivlf0 << endl;
  clog << "ivlf0.empty()=" << ivlf0.empty() << ", ivlf1.empty()=" << ivlf1.empty() << endl;
  
  ivlf1=Interval<Real>(f1,f3);
  ivlf2=Interval<Real>(f0,f1);
  ivlf3=Interval<Real>(f2,f4);
  clog << "ivlf1=" << ivlf1 << ", ivlf2=" << ivlf2 << ", ivlf3=" << ivlf3 << endl;
  clog << "max(ivlf1,ivlf3)=" << max(ivlf1,ivlf3) << endl;
  clog << "min(ivlf1,ivlf3)=" << min(ivlf1,ivlf3) << endl;
  clog << "abs(" << ivlf2 << ")=" << abs(ivlf2) << endl;
  clog << "centre(" << ivlf1 << ")=" << ivlf1.centre() << endl;
  assert(ivlf1.centre()==Real(2.3125));
  
  ivlf1=Interval<Real>(Real(1),Real(1));
  clog << "one=" << ivlf1 << endl;
  ivlf2=exp(ivlf1);
  clog << "exp(" << ivlf1 << ")=" << ivlf2 << endl;
  ivlf3=log(ivlf1);
  clog << "log(" << ivlf2 << ")=" << ivlf3 << endl;
 
  try {
    string input("[1.125,2.25] ");
    stringstream iss(input);
    
    iss >> ivld2;
    clog << "ivld1=" << ivld1 << "  ivld2=" << ivld2 << endl;
    test_assert(ivld1.lower()==ivld2.lower() && ivld1.upper()==ivld2.upper(),
                "construction from stream");
    
    try {
      test_assert(equal(ivld1,ivld2), "equal");
    }
    catch(...) { }
    
    try {
      test_assert(ivld1==ivld2, "operator==");
    }
    catch(...) { }
    
    Interval<Real>& ivlf1ref=ivlf1;
    ivlf1ref=Interval<Real>(5.25,7.375);
    clog << "ivlf1ref=" << ivlf1ref << endl;
    test_assert(ivlf1ref.lower()==Real(5.25),"copy assignment");
    
    ivld1 = ivld2+ivld3;
    clog << ivld2 << " + " << ivld3 << " = " << ivld1 << endl;
    ivld1 = ivld2-ivld3;
    clog << ivld2 << " - " << ivld3 << " = " << ivld1 << endl;
    ivld1 = ivld2*ivld3;
    clog << ivld2 << " * " << ivld3 << " = " << ivld1 << endl;
    ivld1 = ivld2/ivld3;
    clog << ivld2 << " / " << ivld3 << " = " << ivld1 << endl;
    clog << endl;
    //ivld1 = cos(ivld2);
    
    ivlf1 = ivlf2+ivlf3;
    clog << ivlf2 << " + " << ivlf3 << " = " << ivlf1 << endl;
    ivlf1 = ivlf2-ivlf3;
    clog << ivlf2 << " - " << ivlf3 << " = " << ivlf1 << endl;
    ivlf1 = ivlf2*ivlf3;
    clog << ivlf2 << " * " << ivlf3 << " = " << ivlf1 << endl;
    ivlf1 = ivlf2/ivlf3;
    clog << ivlf2 << " / " << ivlf3 << " = " << ivlf1 << endl;

    ivlf1 = ivlf2+f3;
    clog << ivlf2 << " + " << f3 << " = " << ivlf1 << endl;
    ivlf1 = f2+ivlf3;
    clog << f2 << " + " << ivlf3 << " = " << ivlf1 << endl;
    ivlf1 = ivlf2-f3;
    clog << ivlf2 << " - " << f3 << " = " << ivlf1 << endl;
    ivlf1 = f2-ivlf3;
    clog << f2 << " - " << ivlf3 << " = " << ivlf1 << endl;
    ivlf1 = ivlf2*f3;
    clog << ivlf2 << " * " << f3 << " = " << ivlf1 << endl;
    ivlf1 = f2*ivlf3;
    clog << f2 << " * " << ivlf3 << " = " << ivlf1 << endl;
    ivlf1 = ivlf2/f3;
    clog << ivlf2 << " / " << f3 << " = " << ivlf1 << endl;
    ivlf1 = f2/ivlf3;
    clog << f2 << " / " << ivlf3 << " = " << ivlf1 << endl;
    clog << endl;
    

    ivlq1 = ivlq2+ivlq3;
    clog << ivlq2 << " + " << ivlq3 << " = " << ivlq1 << endl;
    ivlq1 = ivlq2-ivlq3;
    clog << ivlq2 << " - " << ivlq3 << " = " << ivlq1 << endl;
    ivlq1 = ivlq2*ivlq3;
    clog << ivlq2 << " * " << ivlq3 << " = " << ivlq1 << endl;
    ivlq1 = ivlq2/ivlq3; 
    clog << ivlq2 << " / " << ivlq3 << " =" << ivlq1 << endl;
    clog << endl;
    //ivlr1 = sin(ivlr2);

  }
  catch(exception& e) {
    cout << "EXCEPTION " << e.what() << "\n";
    return 1;
  }
  
  return 0;
}

template<>
int
test_interval<Rational>(ostream& clog)
{
  typedef Rational Real;
  
  clog << "test_interval<" << name<Real>() << ">" << endl;
  
  Interval<Real> ivld1(Real(1.125),Real(2.25));
  Interval<Real> ivld2;
  Interval<Real> ivld3(2.1,3.2);
  
  Interval<Real> ivlq1(1.1,2.2);
  Interval<Real> ivlq2(1.125,1.125);
  Interval<Real> ivlq3(2.125,3.25);
  
  Real f0=-2.25;
  Real f1=1.5;
  Real f2=2.25;
  Real f3=3.125;
  Real f4=4.0625;
  Interval<Real> ivlf1;
  Interval<Real> ivlf2(f1,f2);
  Interval<Real> ivlf3(f3,f4);
  
  Real r1=1.5;
  Real r2=2.25;
  Real r3=3.125;
  Real r4=4.0625;
  Interval<Real> ivlr1;
  Interval<Real> ivlr2(r1,r2);
  Interval<Real> ivlr3(r3,r4);
  Interval<Real> ivlr4(r1,r4);

  clog << "ivlr1=" << ivlr1 << ", ivlr2=" << ivlr2 << ", ivlr3=" << ivlr3 
       << ", ivlr4=" << ivlr4 << endl;
  
  Interval<Real> ivlf0(f2,f1);
  clog << "ivlf0=" << ivlf0 << endl;
  clog << "ivlf0.empty()=" << ivlf0.empty() << ", ivlf1.empty()=" << ivlf1.empty() << endl;
  
  ivlf1=Interval<Real>(f1,f3);
  ivlf2=Interval<Real>(f0,f1);
  ivlf3=Interval<Real>(f2,f4);
  clog << "min(ivlf1,ivlf2)=" << min(ivlf1,ivlf3) << endl;
  clog << "max(ivlf1,ivlf2)=" << max(ivlf1,ivlf3) << endl;
  clog << "abs(" << ivlf2 << ")=" << abs(ivlf2) << endl;
   
  try {
    string input("[1.125,2.25] ");
    stringstream iss(input);
    
    iss >> ivld2;
    clog << "ivld1=" << ivld1 << "  ivld2=" << ivld2 << endl;
    test_assert(ivld1.lower()==ivld2.lower() && ivld1.upper()==ivld2.upper(),
                "construction from stream");
    
    try {
      test_assert(equal(ivld1,ivld2), "equal");
    }
    catch(...) { }
    
    try {
      test_assert(ivld1==ivld2, "operator==");
    }
    catch(...) { }
    
    Interval<Real>& ivlf1ref=ivlf1;
    ivlf1ref=Interval<Real>(5.25,7.375);
    clog << "ivlf1ref=" << ivlf1ref << endl;
    test_assert(ivlf1ref.lower()==Real(5.25),"copy assignment");
    
    ivld1 = ivld2+ivld3;
    clog << ivld2 << " + " << ivld3 << " = " << ivld1 << endl;
    ivld1 = ivld2-ivld3;
    clog << ivld2 << " - " << ivld3 << " = " << ivld1 << endl;
    ivld1 = ivld2*ivld3;
    clog << ivld2 << " * " << ivld3 << " = " << ivld1 << endl;
    ivld1 = ivld2/ivld3;
    clog << ivld2 << " / " << ivld3 << " = " << ivld1 << endl;
    clog << endl;
    //ivld1 = cos(ivld2);
    
    ivlf1 = ivlf2+ivlf3;
    clog << ivlf2 << " + " << ivlf3 << " = " << ivlf1 << endl;
    ivlf1 = ivlf2-ivlf3;
    clog << ivlf2 << " - " << ivlf3 << " = " << ivlf1 << endl;
    ivlf1 = ivlf2*ivlf3;
    clog << ivlf2 << " * " << ivlf3 << " = " << ivlf1 << endl;
    ivlf1 = ivlf2/ivlf3;
    clog << ivlf2 << " / " << ivlf3 << " = " << ivlf1 << endl;

    ivlf1 = ivlf2+f3;
    clog << ivlf2 << " + " << f3 << " = " << ivlf1 << endl;
    ivlf1 = f2+ivlf3;
    clog << f2 << " + " << ivlf3 << " = " << ivlf1 << endl;
    ivlf1 = ivlf2-f3;
    clog << ivlf2 << " - " << f3 << " = " << ivlf1 << endl;
    ivlf1 = f2-ivlf3;
    clog << f2 << " - " << ivlf3 << " = " << ivlf1 << endl;
    ivlf1 = ivlf2*f3;
    clog << ivlf2 << " * " << f3 << " = " << ivlf1 << endl;
    ivlf1 = f2*ivlf3;
    clog << f2 << " * " << ivlf3 << " = " << ivlf1 << endl;
    ivlf1 = ivlf2/f3;
    clog << ivlf2 << " / " << f3 << " = " << ivlf1 << endl;
    ivlf1 = f2/ivlf3;
    clog << f2 << " / " << ivlf3 << " = " << ivlf1 << endl;
    clog << endl;
    

    ivlq1 = ivlq2+ivlq3;
    clog << ivlq2 << " + " << ivlq3 << " = " << ivlq1 << endl;
    ivlq1 = ivlq2-ivlq3;
    clog << ivlq2 << " - " << ivlq3 << " = " << ivlq1 << endl;
    ivlq1 = ivlq2*ivlq3;
    clog << ivlq2 << " * " << ivlq3 << " = " << ivlq1 << endl;
    ivlq1 = ivlq2/ivlq3; 
    clog << ivlq2 << " / " << ivlq3 << " =" << ivlq1 << endl;
    clog << endl;
    //ivlr1 = sin(ivlr2);

  }
  catch(exception& e) {
    cout << "EXCEPTION " << e.what() << "\n";
    return 1;
  }
  
  return 0;
}
