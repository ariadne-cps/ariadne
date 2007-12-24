/***************************************************************************
 *            test_epsstream.cc
 *
 *  Copyright  2005-7  Pieter Collins
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
#include <fstream>
#include <sstream>
#include <string>

#include "ariadne.h"
#include "test_float.h"

#include "base/utility.h"
#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/zonotope.h"
#include "geometry/list_set.h"
#include "output/epsstream.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::Output;
using namespace std;

int main() {

  Box<Flt> bbox(2);

  Point<Flt> pt("(0.0,0.0)");

  Box<Flt> r1,r2,r3,r4;
  string input("[-0.125,1.125]x[-0.25, 3.25] "
               "[ 0.0125,1.0]x[0.0,2.0] "
               "[ 0.5,1.0]x[1.0,3.0] "
               "[ 0,0.3333333]x[2.3333,3] "
               "[ 0.06125,0.125]x[0.5,2.75] "
               );
  stringstream iss(input);

  iss >> bbox >> r1 >> r2 >> r3 >> r4;
  Zonotope<Flt> z3(r3);
  Polytope<Flt> p4(r4);
  
  cout << "bbox=" << bbox << "\n";

  
  // Test output of basic sets
  epsfstream eps;
  eps.open("test_epsstream-1.eps",bbox);
  cout << "r1=" << r1 << endl;
  eps << r1;
  eps << fill_colour(blue);
  cout << "r2=" << r2 << endl;
  eps << r2;
  eps << fill_colour(red);
  cout << "z3=" << z3 << endl;
  eps << z3;
  cout << "p4=" << p4 << endl;
  eps << p4;
  eps.close();
  eps << pt;

  cout << endl;

  // Test output of grid mask set
  Box<Flt> bb("[0,1]x[0,1]x[0,1]");
  Grid<Flt> g(Vector<Flt>("[0.25,0.25,0.25]"));
  GridMaskSet<Flt> gms(g,bb);
  Box<Flt> r("[0.33,0.66]x[0.125,0.375]x[0.25,0.75]");
  gms.adjoin_outer_approximation(r);
  cout << "gms.size()=" << gms.size() << endl;
  eps.open("test_epsstream-2.eps",bb.neighbourhood(0.125));
  eps << gms;
  eps.close();
  
  return 0;
}
