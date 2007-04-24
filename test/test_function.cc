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
#include "system/function.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace std;

template<class R> int test_function();
  
int main() {
  return test_function<Float>();
}

template<class R>
int test_function() 
{
  typedef typename Numeric::traits<R>::arithmetic_type A;
  std::string f1str="y[0] = ( x[1] + x[2] ) * x[0];";

  GeneralFunction<R> f1(f1str);
  cout << "Done read" << endl;
  cout << f1 << endl;

  Vector<R> x1("[2,3,5]");
  Vector<A> ix1(x1);
  ARIADNE_EVALUATE(f1.image(ix1));
  
  std::string f2str="y[0] = x[0] - sin(x[1]) + x[2];";
  GeneralFunction<R> f2(f2str);
  cout << "Done read" << endl;
  cout << f2 << endl;

  Vector<R> x2("[2,3,5]");
  Vector<A> ix2(x2);
  ARIADNE_EVALUATE(f2.image(ix2));
  
  std::string f3str=
    "y[0] = x[0] - x[1] + x[2]; "
    "y[1] = (x[0] - x[1]) * x[2] + x[1]; ;";
  GeneralFunction<Rational> f3(f3str);
  cout << "Done read" << endl;
  cout << f3 << endl;

  Vector<Rational> qx3("[2,3,5]");
  cout << f3str << endl << f3 << endl << qx3 << endl;
  ARIADNE_EVALUATE(f3.image(qx3));
  
  return 0;
}
