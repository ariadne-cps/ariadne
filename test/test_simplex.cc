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

#include "test_float.h"

#include "ariadne.h"
#include "geometry/point.h"
#include "geometry/simplex.h"
#include "geometry/polyhedron.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Geometry;
using namespace std;

template<class R> int test_simplex();

int main() {
  test_simplex<Flt>();

  tribool x=indeterminate;
  cout << x << endl;
  cout << boolalpha;
  cout << x << endl;
  cout << x << endl;
  cout << indeterminate(x) << endl;
  cout << endl;
  cerr << "INCOMPLETE ";
}
  
template<class R> 
int 
test_simplex() 
{
  Point<R> v0("(-0.25,-0.125)");
  Point<R> v1("(1.0,-0.5)");
  Point<R> v2("(-0.375,0.75)");

  cout << v0 << " " << v1 << " " << v2 << endl;
  
  PointList<R> pl;
  pl.push_back(v0);
  pl.push_back(v1);
  pl.push_back(v2);

  cout << "pl.dimension()=" << pl.dimension() << endl;
  cout << "pl.size()=" << pl.size() << endl;
  cout << "pl.capacity()=" << pl.capacity() << endl;
  cout << "pl=" << pl << std::endl;
  
  Simplex<R> s(pl);
  cout << "s=" << flush << s << endl;

  cout << "s.generators()=" << s.generators() << endl;
  cout << "s.vertices()=" << s.vertices() << endl;
  
  Point<R> pt1,pt2,pt3;

  pt1=Point<R>("(17.0,15.0)");
  pt2=Point<R>("(-0.25,-0.125)");
  pt3=Point<R>("(0.0,0.0)");

  cout << "s.coordinates(pt1)=" << flush;
  cout << s.coordinates(pt1) << endl;
  cout << "s.coordinates(pt2)=" << s.coordinates(pt2) << endl;
  cout << "s.coordinates(pt3)=" << s.coordinates(pt3) << endl;
  
  assert((bool)(!s.contains(pt1)));
  assert((bool)(s.contains(pt3)));
  assert(indeterminate(s.contains(pt2)));
  
  return 0;
}
