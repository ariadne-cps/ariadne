/***************************************************************************
 *            test_zonotope.cc
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
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;

int main() {
  cout << "test_zonotope: " << flush;
  ofstream clog("test_zonotope.log");
  
  Point<Real> c("(0.125,-0.25,0.5)");
  LinearAlgebra::Vector<Real> v("[0.0,1.0,-1.5]");
  LinearAlgebra::Matrix<Real> a("[2.0,1.0,-1.5; 1.0,1.0,0.5; 0.0,0.0,0.375]");
  
  Rectangle<Real> r1=Rectangle<Real>("[9,11]x[5,11]x[0,0]");
  Rectangle<Real> r2=Rectangle<Real>("[5,6]x[3,4]x[-0.5,0.5]");

  Parallelotope<Real> p2=Parallelotope<Real>(r2);

  Zonotope<Real> z1=Zonotope<Real>(r1);
  Zonotope<Real> z2=Zonotope<Real>(p2);

  Zonotope<Real> z3=minkowski_sum(z1,z2);

  clog << z1 << z2 << z3 << endl;

  Point<Real> pt1,pt2,pt3,pt4,pt5,pt6;

  pt1=Point<Real>("(17.0,15.0,-0.5)");
  pt2=Point<Real>("(17.0,8.0,-0.5)");
  pt3=Point<Real>("(14.0,15.0,-0.5)");
  pt4=Point<Real>("(14.0,8.0,-0.5)");
  pt5=Point<Real>("(14.0,5.0,-0.5)");
  pt6=Point<Real>("(15.5,11.5,0)");

  assert(z3.contains(pt1));
  assert(z3.contains(pt2));
  assert(z3.contains(pt3));
  assert(z3.contains(pt4));
  assert(!z3.contains(pt5));
  assert(z3.contains(pt6));

  clog << pt1 << " " << pt2 << " " << pt3 << " " << pt4 << " " << pt5 << " " << pt6 << endl;
  assert(z3.interior_contains(pt1));
  assert(z3.interior_contains(pt2));
  assert(z3.interior_contains(pt3));
  assert(z3.interior_contains(pt4));
  assert(!z3.interior_contains(pt5));
  assert(z3.interior_contains(pt6));

  clog.close();
  cout << "INCOMPLETE\n";

  return 0;
}
