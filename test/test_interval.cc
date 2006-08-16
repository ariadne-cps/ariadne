#include "ariadne.h"
#include "real_typedef.h"
#include "numeric/numerical_types.h"
#include "numeric/interval.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "test.h"

using namespace Ariadne;
using std::ofstream; using std::stringstream; using std::string; 
using std::cout; using std::endl; using std::flush;
using std::exception;

template class Interval<double>;
template class Interval<MPFloat>;
template class Interval<Dyadic>;
template class Interval<Rational>;

int main() {
  cout << "test_interval: " << flush;
  ofstream clog("test_interval.log");
  
  Interval<double> ivld1(1.125,2.25);
  Interval<double> ivld2;
  Interval<double> ivld3(2.1,3.2);
  
  Interval<Rational> ivlq1(1.1,2.2);
  Interval<Rational> ivlq2(Rational(11,10),Rational(22,10));
  Interval<Rational> ivlq3(Rational(21,10),Rational(32,10));
  
  MPFloat f0=-2.25;
  MPFloat f1=1.5;
  MPFloat f2=2.25;
  MPFloat f3=3.125;
  MPFloat f4=4.0625;
  Interval<MPFloat> ivlf1;
  Interval<MPFloat> ivlf2(f1,f2);
  Interval<MPFloat> ivlf3(f3,f4);
  
  Dyadic r1=1.5;
  Dyadic r2=2.25;
  Dyadic r3=3.125;
  Dyadic r4=4.0625;
  Interval<Dyadic> ivlr1;
  Interval<Dyadic> ivlr2(r1,r2);
  //Interval<Dyadic> ivlr3(r3,r4);
  Interval<Dyadic> ivlr3(r1,r4);

  clog << "ivlr1=" << ivlr1 << ", ivlr2=" << ivlr2 << ", ivlr3=" << ivlr3 << endl;
  
  Interval<MPFloat> ivlf0(f2,f1);
  clog << "ivlf0=" << ivlf0 << endl;
  clog << "ivlf0.empty()=" << ivlf0.empty() << ", ivlf1.empty()=" << ivlf1.empty() << endl;
  
  ivlf1=Interval<MPFloat>(f1,f3);
  ivlf2=Interval<MPFloat>(f0,f1);
  ivlf3=Interval<MPFloat>(f2,f4);
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
    
    Interval<MPFloat>& ivlf1ref=ivlf1;
    ivlf1ref=Interval<MPFloat>(5.25,7.375);
    clog << "ivlf1ref=" << ivlf1ref << endl;
    test_assert(ivlf1ref.lower()==MPFloat(5.25),"copy assignment");
    
    clog << "double arithmetic" << endl;
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
    
    clog << "mpfloat arithmetic" << endl;
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
    

    clog << "rational arithmetic" << endl;
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
  
  clog.close();
  cout << "INCOMPLETE - Transcendental functions not fully implemented\n";
  
  return 0;
}
