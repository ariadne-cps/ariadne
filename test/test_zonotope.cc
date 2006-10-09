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

#include "real_typedef.h"

#include "ariadne.h"
#include "base/exception.h"
#include "base/utility.h"
#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"

#include "output/epsfstream.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;

template<typename R> int test_zonotope();
 
int main() {
  test_zonotope<Real>();
  
  cerr << "INCOMPLETE ";
  return 0;
}


template<typename R>
int 
test_zonotope()
{
  Point<R> c("(0.125,-0.25,0.5)");
  LinearAlgebra::Vector<R> v("[0.0,1.0,-1.5]");
  LinearAlgebra::Matrix<R> a("[2.0,1.0,-1.5; 1.0,1.0,0.5; 0.0,0.0,0.375]");
  
  cout << "c=" << c << "\nv=" << v << "\na=" << a << endl;
  
  Rectangle<R> r1=Rectangle<R>("[9,11]x[5,11]x[0,0]");
  Rectangle<R> r2=Rectangle<R>("[5,6]x[3,4]x[-0.5,0.5]");
  cout << "r1=" << r1 << "\nr2=" << r2 << endl;

  Parallelotope<R> p2=Parallelotope<R>(r2);
  cout << "p2=" << p2 << endl;

  Zonotope<R> z1=Zonotope<R>(r1);
  cout << "z1=" << z1 << endl;
  Zonotope<R> z2=Zonotope<R>(p2);
  cout << "z2=" << z2 << endl;

  Zonotope<R> z3=minkowski_sum(z1,z2);
  cout << "z3=" << z3 << endl;

  Point<R> pt1,pt2,pt3,pt4,pt5,pt6;

  pt1=Point<R>("(17.0,15.0,-0.5)");
  pt2=Point<R>("(17.0,8.0,-0.5)");
  pt3=Point<R>("(14.0,15.0,-0.5)");
  pt4=Point<R>("(14.0,8.0,-0.5)");
  pt5=Point<R>("(14.0,5.0,-0.5)");
  pt6=Point<R>("(15.5,11.5,0)");

  Postscript::epsfstream eps("test_zonotope.eps",z3.bounding_box(),0,1);
  eps << z3 << pt1 << pt2 << pt3 << pt4 << pt5 << pt6 << endl;
  eps.close();
  
  assert(z3.contains(pt1));
  assert(z3.contains(pt2));
  assert(z3.contains(pt3));
  assert(z3.contains(pt4));
  assert(!z3.contains(pt5));
  assert(z3.contains(pt6));

  assert(z3.contains(z3.centre()));

  cout << pt1 << " " << pt2 << " " << pt3 << " " << pt4 << " " << pt5 << " " << pt6 << endl;
  assert(!z3.interior_contains(pt1));
  assert(!z3.interior_contains(pt2));
  assert(!z3.interior_contains(pt3));
  assert(!z3.interior_contains(pt4));
  assert(!z3.interior_contains(pt5));
  assert(z3.interior_contains(pt6));

  assert(z3.interior_contains(z3.centre()));
  
  assert(!z3.empty());
  assert(!z3.empty_interior());
  
  return 0;
}
