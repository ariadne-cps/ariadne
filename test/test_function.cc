/***************************************************************************
 *            test_function.cc
 *
 *  Copyright  2007  Alberto Casagrande,  Pieter Collins
 *  Email  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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
#include <fstream>
#include <string>

#include "ariadne.h"
#include "test_float.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "function/interpreted_function.h"
#include "output/logging.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::Function;
using namespace Ariadne::Output;
using namespace std;

template<class R> int test_function();
  
int main() {
  return test_function<Float>();
}

template<class R>
int test_function() 
{
  set_input_verbosity(0);

  typedef typename Numeric::traits<R>::arithmetic_type A;
  ifstream fis("test_function.dat");
  
  do {
    cout << "Inputting function" << endl;
    InterpretedFunction<R> f;
    fis >> f;
    cout << "Done read" << endl;
    cout << f << endl;
    Vector<A> x;
    fis >> x;
    cout << "x=" << x << endl;
    cout << "f.evaluate(x)=" << flush; cout << f.evaluate(x) << endl;
    cout << "f(x)=" << flush; cout << f(x) << endl;
    cout << "f.jacobian(x)=" << flush; cout << f.jacobian(x) << endl;
    cout << endl;

    // Look for end-of-file
    char c; fis>>c; fis.putback(c);
  } while(!fis.eof());
  
  return 0;
}
