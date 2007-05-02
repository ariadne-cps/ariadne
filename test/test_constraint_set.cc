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

#include "test_float.h"
#include "test_macros.h"

#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/grid_set.h"
#include "geometry/constraint_set.h"
#include "system/function_interface.h"
#include "system/function.h"
#include "output/epsfstream.h"
#include "output/logging.h"


using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Output;
using namespace std;

template<class R> int test_constraint_set();

int main() {
  
  cout << boolalpha;
  
  test_constraint_set<Float>();
}



template<class R>
int 
test_constraint_set() 
{
  cout << "test_constraint_set<" << name<R>() << ">" << endl;

  Function<R> f("function disc output Real y; input Real[2] x; algorithm y=1-(x[0]^2+x[1]^2); end disc;");

  ConstraintSet<R> s(f);

  Rectangle<R> r;
  
  // superset
  r=Rectangle<R>("[0.4,0.5]x[0.45,0.65]");
  ARIADNE_ASSERT(s.superset(r));

  // disjoint
  r=Rectangle<R>("[0.8,0.9]x[0.95,0.95]");
  ARIADNE_ASSERT(s.disjoint(r));


  // approximations
  Grid<R> g(Vector<R>("[0.1875,0.125]"));
  Rectangle<R> bb("[-1.01,1.01]x[-1.01,1.01]");
  FiniteGrid<R> fg(g,bb);
  GridMaskSet<R> gmsia(fg); 
  GridMaskSet<R> gmsoa(fg); 
  gmsia.adjoin_inner_approximation(static_cast<const SetInterface<R>&>(s));
  gmsoa.adjoin_outer_approximation(static_cast<const SetInterface<R>&>(s));

  // graphical output
  epsfstream eps;
  eps.open("test_constraint_set-1.eps",bb.neighbourhood(0.25));
  eps.set_fill_colour("red");
  eps << gmsoa;
  eps.set_fill_colour("green");
  eps.set_line_style(false);
  eps << static_cast<const SetInterface<R>&>(s);
  eps.set_line_style(true);
  eps.set_fill_colour("blue");
  eps << gmsia;
  eps.close();
  
  cout << endl;
  
  return 0;
}


