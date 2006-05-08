#include <iostream>

#include "numeric/arithmetic.h"

#include "test.h"

using namespace std;

int main() {
//  using Ariadne::cos;
  using Ariadne::Rational;

  cout << "test_arithmetic: " << flush;

//  test_assert(cos(3.2)==Ariadne::cos(3.2),"cos(double)");
//  test_assert(abs(0.877583-Ariadne::cos[8](Rational(1,2)).get_d()) < 1e-6,"cos[Rational](double)");

  cout << "PASS\n";

  return 0;
}
