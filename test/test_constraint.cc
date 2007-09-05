/***************************************************************************
 *            test_constraint.cc
 *
 *  Copyright  2007  Pieter Collins
 *  Email Pieter.Collins@cwi.nl
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
#include <fstream>
#include <string>

#include <vector>

#include "test.h"
#include "test_float.h"

#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/grid_set.h"
#include "geometry/constraint.h"
#include "function/function_interface.h"
#include "function/interpreted_function.h"
#include "output/epsstream.h"
#include "output/logging.h"


using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Function;
using namespace Ariadne::Geometry;
using namespace Ariadne::Output;
using namespace std;

template<class R> int test_constraint();

int main() {
  
  cout << boolalpha;
  
  test_constraint<Float>();
}



template<class R>
int 
test_constraint() 
{
  cout << "test_constraint<" << name<R>() << ">" << endl;

  InterpretedFunction<R> f("function disc output Real y; input Real[2] x; algorithm y=1-(x[0]^2+x[1]^2); end disc;");

  Constraint<R> c(f);

  Rectangle<R> r;
  Zonotope<Interval<R>,R> z;
  
  
  // satisfies
  r=Rectangle<R>("[0.4,0.5]x[0.45,0.65]");
  ARIADNE_TEST_ASSERT(satisfies(r,c));
  z=Zonotope<Interval<R>,R>(r);
  ARIADNE_TEST_ASSERT(satisfies(z,c));

  // does not satisfy
  r=Rectangle<R>("[0.8,0.9]x[0.95,0.95]");
  ARIADNE_TEST_ASSERT(!satisfies(r,c));
  z=Zonotope<Interval<R>,R>(r);
  ARIADNE_TEST_ASSERT(!satisfies(z,c));

  
  LinearConstraint<R> lc(Vector<R>("[1,1]"),Geometry::less,R(2));
  return 0;
}


