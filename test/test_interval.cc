#include "ariadne.h"
#include "numerical_type.h"
#include "interval.h"

#include <iostream>
#include <string>
#include <sstream>

#include "test.h"

using namespace Ariadne;
using namespace std;

template class Interval<double>;
template class Interval<Dyadic>;
template class Interval<Rational>;

int main() {
  cout << "test_interval: " << flush;
  
  Interval<double> ivld1(1.1,2.2);
  Interval<double> ivld2;
  Interval<double> ivld3(2.1,3.2);
  
  Interval<Rational> ivlr1(1.1,2.2);
  Interval<Rational> ivlr2(Rational(11,10),Rational(22,10));
  Interval<Rational> ivlr3(Rational(21,10),Rational(32,10));
  
  Dyadic f1=1.5;
  Dyadic f2=2.25;
  Dyadic f3=3.125;
  Dyadic f4=4.0625;
  Interval<Dyadic> ivlf1;
  Interval<Dyadic> ivlf2(f1,f2);
  Interval<Dyadic> ivlf3(f3,f4);
  
  try {
    string input("[1.1,2.2] ");
    stringstream iss(input);
    
    iss >> ivld2;
    test_assert(ivld1.lower()==ivld2.lower() && ivld1.upper()==ivld2.upper(),
		"construction from stream");
    
    try {
      test_assert(ivld1==ivld2, "operator==");
    }
    catch(...) { }
    
    ivld1 = ivld2+ivld3;
    ivld1 = ivld2-ivld3;
    ivld1 = ivld2*ivld3;
    ivld1 = ivld2/ivld3;
    //ivld1 = cos(ivld2);
    
    ivlr1 = ivlr2+ivlr3;
    ivlr1 = ivlr2-ivlr3;
    ivlr1 = ivlr2*ivlr3;
    ivlr1 = ivlr2/ivlr3;
    //  ivlr1 = sin(ivlr2);
    
    ivlf1 = ivlf2+ivlf3;
    ivlf1 = ivlf2-ivlf3;
    ivlf1 = ivlf2*ivlf3;
    ivlf1 = ivlf2/ivlf3;
    //     cerr << ivlf1  << " = " << ivlf2  << " / " << ivlf3 << "\n";
    //  ivlf1 = sin(ivlf2);
   
    cout << "INCOMPLETE\n";
  }
  catch(exception& e) {
    cout << "EXCEPTION " << e.what() << "\n";
    return 1;
  }
  
  return 0;
}
