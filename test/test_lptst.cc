/***************************************************************************
 *            test_lptst.cc
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
#include "test/test_float.h"
#include "numeric/rational.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "linear_algebra/permutation.h"
#include "linear_programming/lptst.h"


using namespace Ariadne;
using namespace std;

template<class R> void test_lptst();

int main() {
  
  set_linear_algebra_verbosity(1);
  
  //test_lptst<Flt>();
  test_lptst<Rational>();
  
  cerr << "COMPLETE ";
  return 0;
}

template<class R> void test_lptst() {
  
  Matrix<R> A("[2, 1; 1, 4; 5, 6]");
  Vector<R> b("[30, 64, 110]");
  
  Matrix<R> B = A;
  Vector<R> x1("[5, 10]");    // in
  Vector<R> x2("[2, 16]");    // out
  Vector<R> x3("[2, 15.5]"); // edge
  Vector<R> x4("[13, 5]");    // out
  Vector<R> x5("[12.5, 5]");  // edge
  Vector<R> x6("[7, 13]");     // out
  Vector<R> x7("[7, 12.5]");  // edge
  Vector<R> x8("[10, 10]");  // point
  Interval<R> i = Interval<R>(9, 11);
  Vector< Interval <R> > x9 = Vector<Interval < R > >(2, i);  // point
  Vector< Interval <R> > x10 = x9 +Vector<R>("[-2, 4]");   // out
  
  cout << "lptst result x1 = " << lptst<R, R>(A, b, x1) << " (inside)" << endl;
  cout << "lptst result x2 = " << lptst<R, R>(A, b, x2) << " (outside)" << endl;
  cout << "lptst result x3 = " << lptst<R, R>(A, b, x3) << " (edge)" << endl;
  cout << "lptst result x4 = " << lptst<R, R>(A, b, x4) << " (outside)" << endl;
  cout << "lptst result x5 = " << lptst<R, R>(A, b, x5) << " (edge)" << endl;
  cout << "lptst result x6 = " << lptst<R, R>(A, b, x6) << " (outside)" << endl;
  cout << "lptst result x7 = " << lptst<R, R>(A, b, x7) << " (edge)" << endl;
  cout << "lptst result x8 = " << lptst<R, R>(A, b, x8) << " (point)" << endl;
  cout << "lptst result x9 = " << lptst<R, Interval <R> >(A, b, x9) << " (~point)" << endl;
  cout << "lptst result x10 = " << lptst<R, Interval <R> >(A, b, x10) << " (outside)" << endl;
  
}
