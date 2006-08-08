/***************************************************************************
 *            test_polyhedron.cc
 *
 *  2 May 2005
 *  Copyright  2005  Pieter Collins
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

#include "ariadne.h"
#include "real_typedef.h"
#include "base/exception.h"
#include "base/utility.h"
#include "numeric/numerical_types.h"
#include "geometry/point.h"
#include "geometry/polyhedron.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;

int main() {
  cout << "test_polyhedron: " << flush;
  ofstream clog("test_polyhedron.log");
  
  Point<Real> s1(2);
  s1[0]=Real(1);
  s1[1]=Real(1);
  Point<Real> s2(2);
  s2[0]=Real(4,3);
  Point<Real> s3(2);
  s3[1]=Real(4,3);
  
  std::vector< Point<Real> > ptl;
  ptl.push_back(s1);
  ptl.push_back(s2);
  ptl.push_back(s3);
  
  Polyhedron<Real> p(ptl);
  clog << p << endl;
  
  Point<Rational> qs1(s1);
  Point<Rational> qs2(s2);
  Point<Rational> qs3(s3);
  std::vector< Point<Rational> > qptl;
  qptl.push_back(qs1);
  qptl.push_back(qs2);
  qptl.push_back(qs3);
  
  Polyhedron<Rational> qp(qptl);
  clog << qp << endl;
  
  clog.close();
  cout << "INCOMPLETE\n";

  return 0;
}
