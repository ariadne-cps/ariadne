/***************************************************************************
 *            test_linear_program.cc
 *
 *  Copyright  2006  Pieter Collins, Alberto Casagrande
 *  Email Pieter.Collins@cwi.nl, casagrande@dimi.uniud.it
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
#include <iomanip>
#include <fstream>
#include <string>

#include "output/logging.h"
#include "test/test.h"
#include "test/test_float.h"
#include "numeric/rational.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "linear_programming/linear_program.h"

using namespace Ariadne;
using namespace std;

template<class R>
class TestLinearProgram 
{
  Matrix<R> A;
  Vector<R> b;
  Vector<R> c;
 public:
  TestLinearProgram() 
    : A("[2, 1; 1, 4; 5, 6]"),
      b("[30, 64, 110]"),
      c("[10,10]")
  { }
  
  void test_tableau() const {
    Matrix<R> T1("[1,0,-1,0,1; 0,1,0,-1,1; 1,1,0,0,3; -1,-1,1,1,-2]");
    Matrix<R> T2("[1,0,-1,0,1; 0,1,0,-1,1; 1,1,0,0,3/2; -1,-1,1,1,-2]");
    LinearProgram<R> LP;
    
    LP=LinearProgram<R>(T1);
    cout << LP << endl;
    cout << "feasible = " << LP.is_feasible() << endl;
    cout << LP << endl;
    cout << "feasible_point = " << LP.feasible_point() << endl;
    cout << "optimal_value = " << LP.optimal_value() <<endl;
    
    LP=LinearProgram<R>(T2);
    cout << LP << endl;
    cout << "feasible = " << LP.is_feasible() << endl;
    cout << LP << endl;
  }

  void
  test_feasibility() {
    cout << feasible(A,b);
  }

  void
  test_dual_feasibility() {
    cout << feasible(A,b);
  }

  void
  test_solve() {
    Matrix<R> A("[2, 1; 1, 4; 5, 6]");
    Vector<R> b("[30, 64, 110]");
    Vector<R> c1("[10, 10]");
    Vector<R> c2("[10, 20]");
    LinearProgram<R> LP;
    cout << A << endl;
    cout << b << endl;
    cout << c1 << endl;
    cout << c2 << endl;
    
    cout << feasible(A,b);
    cout << solve(A,b,c1);
    
    LP=LinearProgram<R>(A, b, c1);
    cout << LP << endl;
    cout << LP.is_feasible() << endl;
    cout << LP.feasible_point() << endl;
    cout << LP.optimal_value() << endl;
    cout << LP << endl;
    
    cout << solve(A,b,c1);
  } 
  
  void test() {
    ARIADNE_TEST_CALL(test_feasibility());
    ARIADNE_TEST_CALL(test_dual_feasibility());
    ARIADNE_TEST_CALL(test_solve());
    ARIADNE_TEST_CALL(test_tableau());
  }

};

int main() {
  set_linear_algebra_verbosity(0);
  
  TestLinearProgram<Rational>().test();
  cerr << "INCOMPLETE ";
  return 0;
}


