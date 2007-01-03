/***************************************************************************
 *            test_polytope.cc
 *
 *  Copyright  2005-6  Alberto Casagrande,  Pieter Collins
 *  Email casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>

#include <vector>

#include <ppl.hh>

#include "numeric/float64.h"
#include "numeric/mpfloat.h"
#include "numeric/rational.h"
#include "geometry/point.h"
#include "geometry/point_list.h"
#include "geometry/polytope.h"
#include "geometry/ppl_polyhedron.h"
#include "output/epsfstream.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::Output;
using namespace std;

using namespace Parma_Polyhedra_Library::IO_Operators;

template<class R> int test_polytope();
template<> int test_polytope<Rational>();

int main() {
  //test_polytope<Float64>();
  test_polytope<MPFloat>();
  test_polytope<Rational>();
  
  cerr << "INCOMPLETE ";
}


template<class R>
int 
test_polytope() 
{
  cout << "test_polytope<" << name<R>() << ">" << endl;
  Point<R> pt1("(1.0,0.875)");
  Point<R> pt2("(1.375,0.5)");
  Point<R> pt3("(1.50,1.25)");
  Point<R> pt4("(1.25,1.03125)");
  Point<R> pt5("(0.75,1.25)");
  
  PointList<R> ptl;
  ptl.push_back(pt1);
  ptl.push_back(pt2);
  ptl.push_back(pt3);
  cout << ptl << endl;
 
  Polytope<R> p(ptl);
  cout << "p=" << p << endl;
  
  cout << "p.generators()=" << p.generators() << endl;
  cout << "p.vertices()=" << p.vertices() << endl;
  
  for(size_type i=0; i!=p.number_of_vertices(); ++i) {
    cout << "p.vertex(" << i << ")=" << p.vertex(i) << endl;
  }
  
  cout << "p.vertices_begin() to p.vertices_end() = \n  " << flush;
  for(typename Polytope<R>::vertices_const_iterator v=p.vertices_begin();
      v!=p.vertices_end(); ++v) 
  {
    cout << *v << ", " << flush;
  }
  cout << endl;
  
  
  //assert(p.contains(pt4));
  //assert(!p.contains(pt5));
  //assert(indeterminate(p.contains(pt1)));
  
  cout << endl;
  
  Rectangle<R> bbox=p.bounding_box();
  cout << "p.bounding_box()=" << bbox << endl;
  
  epsfstream eps("test_polytope-1.eps",bbox,0,1);
  eps << p << pt1 << pt2 << pt3 << pt4 << pt5 << endl;
  eps.close();
  
  cout << endl;
  
  Rectangle<R> r("[-0.125,0.25]x[0.125,0.75]");
  cout << "r=" << r << endl;
  p=Polytope<R>(r);
  cout << "Polytope(r)=" << p << endl;
  
  cout << endl;


  Polytope<R> pltp2=Polytope<R>(Matrix<R>("[3.125,1.125,-2.875,-0.875;1.75,-0.25,-2.25,-0.25]"));
  cout << "pltp2=" << pltp2 << endl << "pltp2.bounding_box()=" << pltp2.bounding_box() << endl;
  RegularGrid<R> gr2(2,0.125);
  GridCellListSet<R> oap2=over_approximation(pltp2,gr2);
  GridCellListSet<R> uap2=under_approximation(pltp2,gr2);
  cout << "oap2.size()=" << oap2.size() << endl;
  cout << "uap2.size()=" << uap2.size() << endl;
  Rectangle<R> bbox2=pltp2.bounding_box().expand_by(0.25);
  eps.open("test_polytope-2.eps",bbox2);
  eps.set_fill_colour("white");
  eps << pltp2.bounding_box();
  eps.set_fill_colour("red");
  eps << oap2;
  eps.set_fill_colour("green");
  eps << pltp2;
  eps.set_fill_colour("blue");
  eps << uap2;
  eps.close();

  cout << endl;

  return 0;
}

template<>
int 
test_polytope<Rational>() 
{
  typedef Rational R;
  
  cout << "test_polytope<" << name<R>() << ">" << endl;
  Point<R> s1("(1,7/8)");
  Point<R> s2("(11/8,1/2)");
  Point<R> s3("(3/2,5/4)");
  Point<R> s4("(5/4,33/32)");
  Point<R> s5("(3/4,5/4)");
  
  PointList<R> ptl;
  ptl.push_back(s1);
  ptl.push_back(s2);
  ptl.push_back(s3);
  cout << ptl << endl;
 
  Polytope<R> p(ptl);
  cout << "p=" << p << endl;
  
  
  Parma_Polyhedra_Library::C_Polyhedron ppl_p=ppl_polyhedron(p);
  cout << ppl_p.generators() << endl;
  
  cout << ppl_p.constraints() << endl;
  cout << constraints_matrix(ppl_p) << " " << constraints_vector(ppl_p) << endl;
  
  Parma_Polyhedra_Library::NNC_Polyhedron ppl_ip=interior(ppl_p);
  cout << ppl_ip.constraints() << endl;
  
  Parma_Polyhedra_Library::C_Polyhedron ppl_s=ppl_polyhedron(s4.position_vector());
  cout << ppl_s.generators() << endl;

  //assert(p.contains(s4));
  //assert(!p.contains(s5));
  //assert(indeterminate(p.contains(s1)));
 
  cout << endl;
  return 0;
}
