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
#include "base/exception.h"
#include "base/utility.h"
#include "numeric/float64.h"
#include "numeric/mpfloat.h"
#include "numeric/rational.h"
#include "geometry/point.h"
#include "geometry/point_list.h"
#include "geometry/polyhedron.h"

#include "geometry/polyhedron.tpl"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;

using namespace Parma_Polyhedra_Library::IO_Operators;


int
test_ppl_polyhedron(ostream& clog)
{
  return 0;
}  

template<typename R>
int 
test_polyhedron(ostream& clog) 
{
  clog << "test_polyhedron<" << name<R>() << ">" << endl;
  Point<R> s1("(1.0,0.875)");
  Point<R> s2("(1.375,0.5)");
  Point<R> s3("(1.50,1.25)");
  Point<R> s4("(1.25,1.03125)");
  Point<R> s5("(0.75,1.25)");
  
  PointList<R> ptl;
  ptl.push_back(s1);
  ptl.push_back(s2);
  ptl.push_back(s3);
  clog << ptl << endl;
 
  Polyhedron<R> p(ptl);
  clog << p << endl;
  
  Parma_Polyhedra_Library::C_Polyhedron ppl_p(p);
  clog << ppl_p.generators() << endl;
  
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
  
  return 0;
}



template<typename R>
int 
test_polytope(ostream& clog) 
{
  clog << "test_polytope<" << name<R>() << ">" << endl;
  LinearAlgebra::Matrix<R> A("(1.0,0.875;-1,1.125;0.125,-2.25)");
  LinearAlgebra::Vector<R> b("(1.375,0.5,0.25)");
  Polytope<R> ptp;
  
  ptp=Polytope<R>(A,b);
  
  return 0;
}



int main() {
  cout << "test_polyhedron: " << flush;
  ofstream clog("test_polyhedron.log");

  test_ppl_polyhedron(clog);
  
  test_polyhedron<Float64>(clog);
  test_polyhedron<MPFloat>(clog);
  test_polyhedron<Rational>(clog);
  
  test_polytope<Float64>(clog);

  clog.close();
  cout << "INCOMPLETE\n";
}
