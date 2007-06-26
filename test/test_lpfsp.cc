/***************************************************************************
 *            test_lpfsp.cc
 *
 *  31 Jan 2006
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
#include "test/test_float.h"
#include "numeric/rational.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "linear_algebra/lp.h"
#include "linear_algebra/permutation.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace std;

template<class R> void test_lpfsp();

int main() {
  
  set_linear_algebra_verbosity(1);
  
#ifdef FLOAT_TEST
  //  test_lpfsp<double>();
  test_lpfsp<Float64>();
#else
  test_lpfsp<Rational>();
#endif
  
  cerr << "COMPLETE ";
  return 0;
}

template<class R> bool doit(Matrix<R> &A, Vector<R> &b, Vector<R> &c) {
  bool ok = true;
  Vector<R> x;
  Vector<R> y;
  Permutation p(b.size() + c.size());
  bool ans;
  
  // call lpfsp
    ans = lpfsp(A, b, p, x, y);

    // verify if the output vector x has positive coefficients
  for (int i = 0; i < x.size(); i++) {
    if (x[i] < 0) {
      ok = false;
      cout << "Test_lpfsp -- error: negative coefficient x[" << i << "] = " << x[i] << endl;
    }
  }
  
  // verify if the output vector y has positive coefficients
  for (int i = 0; i < y.size(); i++) {
    if (y[i] < 0) {
      ok = false;
      cout << "Test_lpfsp -- error: negative coefficient y[" << i << "] = " << y[i] << endl;
    }
  }
  
  // verify primal inequalities
  for (int i = 0; i < y.size(); i++) {
    R sum = 0;
    for (int j = 0; j < x.size(); j++)
      sum += A(i, j)*x(j);
    if (sum > b(i)) {
      ok = false;
      cout << "Test_lpfsp -- error in inequalities: sum=" << sum << " > b(" << i << ")=" << b(i) << endl;
    }
  }
  
//  if (!ok)
    cout << "Test_lpfsp -- A = " << A << "\nb = " << b <<
    ", lpfsp = " << ans << ",  x = " << x << ", y = " << y << endl << endl;
  return ok;
}

template<class R> void test_lpfsp() {
  
  Matrix<R> mat;
  Vector<R> b;
  Vector<R> c;

  // cycling possible
  mat = Matrix<R>("[0.5, -5.5, -2.5, 9; 0.5, -1.5, -0.5, 1; 1, 0, 0, 0]");
  b = Vector<R>("[0, 0, 1]");
  c = Vector<R>("[10, -57, -9, -24]");
  doit(mat, b, c);
  
  // infeasible origin
  mat = Matrix<R>("[2, -1, 2; 2, -3, 1; -1, 1, -2]");
  b = Vector<R>("[4, -5, -1]");
  c = Vector<R>("[1, -1, 1]");
  doit(mat, b, c);
  
  mat = Matrix<R>("[1, -1; -1, -1; 2, 1]");
  b = Vector<R>("[-1, -3, 4]");
  c = Vector<R>("[3, 1]");
  doit(mat, b, c);
  
  // unbounded
  mat = Matrix<R>("[-2, 1; -1, -2]");
  b = Vector<R>("[-1, -2]");
  c = Vector<R>("[1, -1]");
  doit(mat, b, c);
  
  mat = Matrix<R>("[1, -1; -1, -1; 2, -1]");
  b = Vector<R>("[-1, -3, 2]");
  c = Vector<R>("[3, 1]");
  doit(mat, b, c);
  
  // infeasible
  mat = Matrix<R>("[1, 1; -2, -2]");
  b = Vector<R>("[2, -10]");
  c = Vector<R>("[3, -1]");
  doit(mat, b, c);
  
  mat = Matrix<R>("[1, -1; -1, -1; 2, 1]");
  b = Vector<R>("[-1, -3, 2]");
  c = Vector<R>("[3, 1]");
  doit(mat, b, c);
  
}

