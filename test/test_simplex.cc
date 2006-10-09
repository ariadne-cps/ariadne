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

#include "real_typedef.h"

#include "ariadne.h"
#include "base/exception.h"
#include "base/utility.h"
#include "geometry/point.h"
#include "geometry/simplex.h"
#include "geometry/polyhedron.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;

template<class R> int test_simplex();

int main() {
  test_simplex<Real>();
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

  Point<R> pt1,pt2,pt3;

  pt1=Point<R>("(17.0,15.0)");
  pt2=Point<R>("(-0.25,-0.125)");
  pt3=Point<R>("(0.0,0.0)");

  assert(!s.contains(pt1));
  assert(s.contains(pt2));
  assert(!s.interior_contains(pt2));
  assert(s.interior_contains(pt3));

  
  cerr << "INCOMPLETE ";

  return 0;
}
