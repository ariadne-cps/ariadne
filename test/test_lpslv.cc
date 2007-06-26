/***************************************************************************
 *            test_lpslv.cc
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

//#define FLOAT_TEST

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace std;

template<class R> void test_lpslv();

int main() {
  
  set_linear_algebra_verbosity(2);
  
#ifdef FLOAT_TEST
  //  test_lpslv<double>();
  test_lpslv<Float64>();
#else
  test_lpslv<Rational>();
#endif
  
  cerr << "COMPLETE ";
  return 0;
}

template<class R> void dumperror(Matrix<R> &A, Vector<R> &b, Vector<R> &c, Permutation p, Vector<R>x, Vector<R>y, R optval, char* str) {
    cout << str << endl;
    cout << "Test_lpslv -- A = " << A << " b = " << b << ", c = " << c <<
    ",\n\tlpslv = " << optval << ", p = " << p << ", x = " << x << ", y = " << y << endl << endl;
}

template<class R> bool doit(Matrix<R> &A, Vector<R> &b, Vector<R> &c) {
  bool ok = true;
  Vector<R> x;
  Vector<R> y;
  Permutation p(b.size() + c.size());
  R optval;
  
  // call lpslv
  try {
    optval = lpslv(A, b, c, p, x, y);
  } catch (InfeasibleSolution ) {
    dumperror(A, b, c, p, x, y, optval, "Caught Infeasible Solution exception");
    return false;
  } catch (UnboundProblem ) {
    dumperror(A, b, c, p, x, y, optval, "Caught UnboundProblem exception");
    return false;
  } catch (InfeasibleOrigin ) {
    dumperror(A, b, c, p, x, y, optval, "Caught InfeasibleOrigin exception");
    return false;
  } catch ( ... ) {
    dumperror(A, b, c, p, x, y, optval, "Caught default exception");
    return false;
  }
  
  // verify if the output vector x has positive coefficients and produces the returned value
  R check = 0;
  for (int i = 0; i < x.size(); i++) {
    if (x[i] < 0) {
      ok = false;
      cout << "Test_lpslv -- error: negative coefficient x[" << i << "] = " << x[i] << endl;
    }
    check += x[i]*c[i];
  }
  if (check != optval) {
    ok = false;
    cout << "Test_lpslv -- error: optimal value lpslv (=" << optval << ") != x[i]*c[i] (=" << check << ")\n";
  }
  
  // verify if the output vector y has positive coefficients and produces the returned value
  check = 0;
  for (int i = 0; i < y.size(); i++) {
    if (y[i] < 0) {
      ok = false;
      cout << "Test_lpslv -- error: negative coefficient y[" << i << "] = " << y[i] << endl;
    }
    check += y[i]*b[i];
  }
  if (check != optval) {
    ok = false;
    cout << "Test_lpslv -- error: optimal value lpslv (=" << optval << ") != y[i]*b[i] (=" << check << ")\n";
  }
  
  // verify primal inequalities
  for (int i = 0; i < y.size(); i++) {
    R sum = 0;
    for (int j = 0; j < x.size(); j++)
      sum += A(i, j)*x(j);
    if (sum > b(i)) {
      ok = false;
      cout << "Test_lpslv -- error in inequalities: sum=" << sum << " > b(" << i << ")=" << b(i) << endl;
    }
  }
  
  // verify dual inequalities
  for (int j = 0; j < x.size(); j++) {
    R sum = 0;
    for (int i = 0; i < y.size(); i++)
      sum += A(i, j)*y(i);
    if (sum < c(j)) {
      ok = false;
      cout << "Test_lpslv -- error in inequalities: sum=" << sum << " < c(" << j << ")=" << c(j) << endl;
    }
  }
  
  if (!ok)
    cout << "Test_lpslv -- A = " << A << " b = " << b << ", c = " << c <<
    ",\n\tlpslv = " << optval << ", p = " << p << ", x = " << x << ", y = " << y << endl << endl;
  return ok;
}

template<class R> void test_lpslv() {
  
  Matrix<R> mat;
  Vector<R> b;
  Vector<R> c;
  
  mat = Matrix<R>("[2, 1; 1, 4; 5, 6]");
  b = Vector<R>("[30, 64, 110]");
  c = Vector<R>("[10, 10]");
  doit(mat, b, c);
  
  c = Vector<R>("[10, 20]");
  doit(mat, b, c);
  
  mat = Matrix<R>("[-1, 2; 1, 6; 3, 5; 5, 3; 6, 1; 2, -1]");
  b = Vector<R>("[36, 132, 136, 136, 132, 36]");
  c = Vector<R>("[10, 10]");
  doit(mat, b, c);
  
  mat = Matrix<R>("[2, 1; -3, 1; 1, -2]");
  b = Vector<R>("[1, 1, 1]");
  c = Vector<R>("[1, 1]");
  doit(mat, b, c);
  
  mat = Matrix<R>("[1, 3, 1; -1, 0, 3; 2, -1, 2; 2, 3, -1]");
  b = Vector<R>("[3, 2, 4, 2]");
  c = Vector<R>("[5, 5, 3]");
  doit(mat, b, c);
  
  mat = Matrix<R>("[2, 3, 1; 4, 1, 2; 3, 4, 2]");
  b = Vector<R>("[5, 11, 8]");
  c = Vector<R>("[5, 4, 3]");
  doit(mat, b, c);
  
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

