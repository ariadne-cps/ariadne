/***************************************************************************
 *            test_simplex.cc
 *
 *  Copyright  2006  Pieter Collins
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
#include <string>

#include "ariadne.h"
#include "real_typedef.h"
#include "base/exception.h"
#include "base/utility.h"
#include "numeric/numerical_types.h"
#include "geometry/point.h"
#include "geometry/simplex.h"
#include "geometry/polyhedron.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;

int main() {
  cout << "test_simplex: " << flush;
  ofstream clog("test_simplex.log");
  
  Point<Real> v0("(0.125,-0.25)");
  Point<Real> v1("(1.0,-0.5)");
  Point<Real> v2("(-0.25,0.75)");

  std::vector< Point<Real> > pl;
  pl.push_back(v0);
  pl.push_back(v1);
  pl.push_back(v2);

  Simplex<Real> s(pl);
  clog << s << endl;

  Point<Real> pt1,pt2,pt3;

  pt1=Point<Real>("(17.0,15.0)");
  pt2=Point<Real>("(0.125,-0.25)");
  pt3=Point<Real>("(0.0,0.0)");

  assert(!s.contains(pt1));
  assert(s.contains(pt2));
  assert(!s.interior_contains(pt2));
  assert(s.interior_contains(pt3));

  Polyhedron<Real> ply(s);
  clog << ply << endl;
  clog.close();
  cout << "INCOMPLETE\n";

  return 0;
}
