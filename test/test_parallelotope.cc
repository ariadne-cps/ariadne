/***************************************************************************
 *            test_parallelotope.cc
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

#include <vector>

#include "real_typedef.h"

#include "ariadne.h"
#include "base/exception.h"
#include "base/utility.h"
#include "numeric/numerical_types.h"
#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;

template<typename R> int test_parallelotope();

int main() {
  test_parallelotope<Real>();
  cerr << "INCOMPLETE ";
  return 0;
}

template<typename R>
int
test_parallelotope()
{ 
  Rectangle<R> r1=Rectangle<R>("[9,11]x[5,11]x[-1,1]");
  Rectangle<R> r2=Rectangle<R>("[5,6]x[3,4]x[-0.5,0.5]");

  Parallelotope<R> p1=Parallelotope<R>(r1);
  Parallelotope<R> p2=Parallelotope<R>(r2);

  cout << p1 << p2 << endl;

  Point<R> pt1,pt2,pt3,pt4,pt5,pt6;
  pt1=Point<R>("(17.0,15.0,-0.5) ");
  pt2=Point<R>("(17.0,8.0,-0.5) ");
  pt3=Point<R>("(14.0,15.0,-0.5) ");
  pt4=Point<R>("(14.0,15.0,-0.5) ");
  pt5=Point<R>("(14.0,15.0,-0.5) ");
  pt6=Point<R>("(15.5,11.5,0) ");


  cout << pt1 << pt2 << pt3 << pt4 << pt5 << pt6 << endl;
  
  assert(!p1.contains(pt1));
  assert(!p1.contains(pt2));
  assert(!p1.contains(pt3));
  assert(!p1.contains(pt4));
  assert(!p1.contains(pt5));
  assert(!p1.contains(pt6));

  assert(!p1.interior_contains(pt1));
  assert(!p1.interior_contains(pt2));
  assert(!p1.interior_contains(pt3));
  assert(!p1.interior_contains(pt4));
  assert(!p1.interior_contains(pt5));
  assert(!p1.interior_contains(pt6));

  return 0;
}
