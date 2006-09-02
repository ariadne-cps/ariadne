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

#include <ppl.hh>

#include "ariadne.h"
#include "real_typedef.h"
#include "base/exception.h"
#include "base/utility.h"
#include "numeric/numerical_types.h"
#include "geometry/point.h"
#include "geometry/point_list.h"
#include "geometry/polyhedron.h"

#include "geometry/polyhedron.tpl"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;

using namespace Parma_Polyhedra_Library::IO_Operators;

Parma_Polyhedra_Library::C_Polyhedron ppl_polyhedron(const Point<Real>& pt) {
  return ppl_polyhedron(LinearAlgebra::Vector<Rational>(pt.position_vector())); }
  
int main() {
  cout << "test_polyhedron: " << flush;
  ofstream clog("test_polyhedron.log");
  
  Point<Real> s1("(1.0,0.875)");
  Point<Real> s2("(1.375,0.5)");
  Point<Real> s3("(1.50,1.25)");
  Point<Real> s4("(1.25,1.03125)");
  Point<Real> s5("(0.75,1.25)");
  
  PointList<Real> ptl;
  ptl.push_back(s1);
  ptl.push_back(s2);
  ptl.push_back(s3);
  clog << ptl << endl;
 
  Polyhedron<Real> p(ptl);
  clog << p << endl;
  
  Parma_Polyhedra_Library::C_Polyhedron ppl_p(p);
  clog << ppl_p.generators() << endl;
  clog << PointList<Real>(generators(ppl_p)) << endl;
  
  clog << ppl_p.constraints() << endl;
  clog << constraints_matrix(ppl_p) << " " << constraints_vector(ppl_p) << endl;
  
  Parma_Polyhedra_Library::NNC_Polyhedron ppl_ip=interior(ppl_p);
  clog << ppl_ip.constraints() << endl;
  
  Parma_Polyhedra_Library::C_Polyhedron ppl_s=ppl_polyhedron(s4);
  clog << ppl_s.generators() << endl;
  
  assert(p.contains(s1));
  assert(!p.interior_contains(s1));
  assert(p.contains(s4));
  assert(p.interior_contains(s4));
  assert(!p.contains(s5));
  assert(!p.interior_contains(s5));
  
  Point<Rational> qs1(s1);
  Point<Rational> qs2(s2);
  Point<Rational> qs3(s3);
  Point<Rational> qs4(s4);
  Point<Rational> qs5(s5);
  PointList<Rational> qptl;
  qptl.push_back(qs1);
  qptl.push_back(qs2);
  qptl.push_back(qs3);
  
  Polyhedron<Rational> qp(qptl);
  clog << qp << endl;
  
  assert(qp.contains(qs1));
  assert(!qp.interior_contains(qs1));
  assert(qp.contains(qs4));
  assert(qp.interior_contains(qs4));
  assert(!qp.contains(qs5));
  assert(!qp.interior_contains(qs5));

  clog.close();
  cout << "INCOMPLETE\n";

  return 0;
}
