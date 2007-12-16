/***************************************************************************
 *            test_taylor_variable.cc
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
#include "function/taylor_series.h"
#include "function/taylor_variable.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Function;
using namespace std;

template<class X>
class TestTaylorVariable {
 private:
  TaylorVariable<X> x1,x2,x3;
 public:
  TestTaylorVariable() {
    double a1[15]={ 2.0, 1.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double a2[15]={ 3.0, 1.0, 0.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double a3[15]={ 2.0, 1.0, 0.0, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    x1=TaylorVariable<X>(2,4,a1);
    x2=TaylorVariable<X>(2,4,a2);
    x3=TaylorVariable<X>(1,4,a3);
    
    ARIADNE_TEST_CALL(test_degree());
    ARIADNE_TEST_CALL(test_add());
    ARIADNE_TEST_CALL(test_sub());
    ARIADNE_TEST_CALL(test_mul());
    ARIADNE_TEST_CALL(test_div());
    ARIADNE_TEST_CALL(test_rec());
    ARIADNE_TEST_CALL(test_pow());
    ARIADNE_TEST_CALL(test_compose());
  }

  void test_degree() {
    ARIADNE_TEST_ASSERT(x1.degree()==4);
  }

  void test_add() {
    cout << x1 << "+" << x2 << " = " << x1+x2 << std::endl;
    //assert((x1+x2)==TaylorVariable<X>("[3,2,0,0]"));
  }

  void test_sub() {
    cout << x1 << "-" << x2 << " = " << x1-x2 << std::endl;
    //assert((x1-x2)==TaylorVariable<X>("[-1,0,0,0]"));
  }

  void test_mul() {
    cout << x1 << "*" << x2 << " = " << x1*x2 << std::endl;
    //assert((x1*x2)==TaylorVariable<X>("[2,3,2,0]"));
  }

  void test_div() {
    /*
    TaylorVariable<X> x3("[2,3,4]");
    TaylorVariable<X> x4("[1,0,0]");
    cout << x3 << "/" << x4 << " = " << x3/x4 << std::endl;
    assert((x3/x4)==x3);
    cout << x4 << "/" << x3 << " = " << x4/x3 << std::endl;
    assert((x4/x3)==TaylorVariable<X>("[0.5,-0.75,1.25]"));
    cout << 1 << "/" << x2 << " = " << 1/x2 << std::endl;
    assert((1/x2)==TaylorVariable<X>("[0.5,-0.25,0.25,-0.375]"));
    cout << x1 << "/" << x2 << " = " << x1/x2 << std::endl;
    assert((x1/x2)==TaylorVariable<X>("[0.5,0.25,-0.25,0.375]"));
    */
  }

  void test_rec() {
    TaylorVariable<X> x=TaylorVariable<X>::variable(1,5,1.0,0);
    assert(rec(rec(x))==x);
  }

  void test_pow() {
    cout << x2 << "^5 = " << pow(x2,5) << std::endl;
    //    assert(pow(x2,5)==TaylorVariable<X>("[32,80,160,240]"));
  }

  void test_compose() {
    //double ax[10] = { 3.0, 1.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double ax[10] = { 3.0, 1.0, 0.0, 0.0, 0.125, 0.25, 0.0, 0.0, 0.0, 0.0 };
    double ay[4] = { 1.0, -1.0, 0.5, -0.25 };
    double aid[4] = { ax[0], 1.0, 0.0, 0.0 };
    TaylorVariable<X> x(2,3,ax);
    TaylorSeries<X> y(3,ay);
    TaylorSeries<X> id(3,aid);
    cout << "x=" << x << endl;
    cout << "y=" << y << endl;
    cout << "compose(y,x)=" << compose(y,x) << endl;
    cout << "compose(id,x)=" << compose(id,x) << endl;
    cout << "compose(id,x)-x=" << compose(id,x)-x << endl;
  }


};


int main() {
  TestTaylorVariable<Rational> t1;
  cout << "INCOMPLETE " << flush;
  return 0;
}
