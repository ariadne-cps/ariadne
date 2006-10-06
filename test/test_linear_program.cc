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

#include "real_typedef.h"
#include "numeric/rational.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "linear_algebra/linear_program.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace std;

template<typename R>
void test_linear_program();


int main() {
  test_linear_program<Rational>();
  cout << "Incomplete";
  return 0;
}

template<typename R>
void 
test_linear_program()
{
  Matrix<R> A("[1,1,1,1; 1,1,1,1; 1,1,1,1]");
  Vector<R> b("[1,1,1]");
  Vector<R> c("[1,1,1,1]");
  
  LinearProgram<R> LP(A,b,c);
  cout << LP << endl;

  LP.solve();
  cout << LP << endl;

  
}
