/***************************************************************************
 *            test_first_derivative.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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

#include "test_float.h"
#include "numeric/rational.h"
#include "linear_algebra/vector.h"
#include "function/first_derivative.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Function;
using namespace std;

template class FirstDerivative< Flt, Vector<Flt> >;
template class FirstDerivative< Rational, Vector<Rational> >;

template<class X, class V>
class TestFirstDerivative {
 private:
  FirstDerivative<X,V> x1,x2,x3,x4,x5;
 public:
  TestFirstDerivative() {
    x1=FirstDerivative<X,V>(2.0,2,0);
    x2=FirstDerivative<X,V>(2.0,2,1);
    x3=FirstDerivative<X,V>(2.0,"[3,2]");
    x4=FirstDerivative<X,V>(-1.0,"[1,4]");
    x5=FirstDerivative<X,V>(2,"[1,-1]");

    test_add();
    test_sub();
    test_mul();
    test_div();
    test_pow();
  }

  void test_add() {
    cout << x1 << "+" << x2 << " = " << x1+x2 << std::endl;
    assert(((x1+x2)==FirstDerivative<X,V>("4","[1,1]")));
  }

  void test_sub() {
    cout << x1 << "-" << x2 << " = " << x1-x2 << std::endl;
    assert(((x1-x2)==FirstDerivative<X,V>("0","[1,-1]")));
  }

  void test_mul() {
    cout << x3 << "*" << x4 << " = " << x3*x4 << std::endl;
    assert(((x3*x4)==FirstDerivative<X,V>("-2","[-1,6]")));
  }

  void test_div() {
    cout << x3 << "/" << x4 << " = " << x3/x4 << std::endl;
    assert(((x3/x4)==FirstDerivative<X,V>("-2","[-5,-10]")));
  }

  void test_pow() {
    cout << x5 << "^5 = " << pow(x5,5) << std::endl;
    assert((pow(x5,5)==FirstDerivative<X,V>("32","[80,-80]")));
  }

};



int main() {
  TestFirstDerivative<Rational, Vector<Rational> > t1;
  
  return 0;
}
