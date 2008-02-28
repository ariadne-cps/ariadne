/***************************************************************************
 *            test_constraint_set.cc
 *
 *  Copyright  2007  Alberto Casagrande,  Pieter Collins
 *  Email casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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
#include "geometry/constraint_set.h"
#include "function/function_interface.h"
#include "function/build_function.h"
#include "output/epsstream.h"
#include "output/logging.h"


using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Function;
using namespace Ariadne::Geometry;
using namespace Ariadne::Output;
using namespace std;

template<class R> int test_constraint_set();

template<class R, class A, class P> 
void radius_function(R& r, const A& a, const P& p) {
  r[0] = 1 - ( pow(a[0],2u) + pow(a[1],2u) );
}

ARIADNE_BUILD_FUNCTION(Radius,radius_function,1,2,0,255);

int main() {
  
  cout << boolalpha;
  
  test_constraint_set<Flt>();
}



template<class R>
int 
test_constraint_set() 
{
  cout << "test_constraint_set<" << Numeric::name<R>() << ">" << endl;

  Radius<R> f=Radius<R>(Vector<R>());
  Box<R> po=Box<R>::positive_orthant(1);

  ConstraintSet<R> s(f,po);

  Box<R> r;
  
  // superset
  r=Box<R>("[0.4,0.5]x[0.45,0.65]");
  ARIADNE_TEST_ASSERT(s.superset(r));

  // disjoint
  r=Box<R>("[0.8,0.9]x[0.95,0.95]");
  ARIADNE_TEST_ASSERT(s.disjoint(r));


  // approximations
  Grid<R> g(Vector<R>("[0.1875,0.125]"));
  Box<R> bb("[-1.01,1.01]x[-1.01,1.01]");
  FiniteGrid<R> fg(g,bb);
  GridMaskSet<R> gmsia(fg); 
  GridMaskSet<R> gmsoa(fg); 
  gmsia.adjoin_inner_approximation(static_cast<const SetInterface<R>&>(s));
  gmsoa.adjoin_outer_approximation(static_cast<const SetInterface<R>&>(s));

  // graphical output
  epsfstream eps;
  eps.open("test_constraint_set-1.eps",bb.neighbourhood(0.25));
  eps << fill_colour(red) << gmsoa;
  eps << fill_colour(green);
  eps << line_style(false) << static_cast<const SetInterface<R>&>(s);
  eps << line_style(true) << fill_colour(blue) << gmsia;
  eps.close();
  
  cout << endl;
  
  return 0;
}


