/***************************************************************************
 *            test_lpstp.cc
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
#include "linear_algebra/permutation.h"
#include "linear_programming/lp.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace std;

template<class R> void test_lpstp();

int main() {
  
  set_linear_algebra_verbosity(0);
   
  test_lpstp<Rational>();
  cerr << "INCOMPLETE ";
  return 0;
}

template<class R> void test_lpstp() {
  
  Matrix<R> A("[2, 1; 1, 4; 5, 6]");
  Vector<R> b("[30, 64, 110]");
  Vector<R> c1("[10, 10]");
  Vector<R> c2("[10, 20]");
  Permutation p(b.size()+c1.size());

  Matrix<R> B = A;
//  B.inverse();
  Vector<R> x;
  Vector<R> y;
  Vector<R> z;
  
  cout << p << endl;
  
  cout << "lpstp result c1 = " << lpstp(A, b, c1, p, B, x, y, z) << endl;
  cout << "lpstp result c2 = " << lpstp(A, b, c2, p, B, x, y, z) << endl;
   
}
