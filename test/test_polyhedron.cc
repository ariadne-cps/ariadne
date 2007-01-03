/***************************************************************************
 *            test_polyhedron.cc
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

#include "real_typedef.h"

#include "ariadne.h"
#include "numeric/float64.h"
#include "numeric/mpfloat.h"
#include "numeric/rational.h"
#include "geometry/point.h"
#include "geometry/point_list.h"
#include "geometry/rectangle.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "geometry/polyhedron.h"
#include "geometry/grid.h"
#include "geometry/grid_set.h"
#include "output/epsfstream.h"

#include "geometry/ppl_polyhedron.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::Output;
using namespace std;

template<class R> int test_polyhedron();
template<> int test_polyhedron<Rational>();

int main() {
  
  cout << boolalpha;
  
  test_polyhedron<MPFloat>();
  test_polyhedron<Rational>();
   
  cerr << "INCOMPLETE ";
}



template<class R>
int 
test_polyhedron() 
{
  cout << "test_polyhedron<" << name<R>() << ">" << endl;
  Matrix<R> A("[1.0,0.875;-1,1.125;0.125,-2.25]");
  Vector<R> b("[1.375,0.5,0.25]");
  Polyhedron<R> phd1;
  
  phd1=Polyhedron<R>(A,b);
  cout << "phd1=" << phd1 << endl;
  cout << "phd1.dimension()=" << phd1.dimension() << ", " 
       << "phd1.number_of_constraints()=" << phd1.number_of_constraints() << endl;
  //assert(!phd1.empty());

  Point<R> pt1("(0.25,0.375)");
  cout << "pt1=" << pt1 << endl;

  Constraint<R> c=Constraint<R> (phd1.dimension(),phd1.A().begin(),phd1.b().begin()[0]);
  cout << c << flush;
  cout << "  " << c.satisfied_by(pt1) << endl;
  
  typename Polyhedron<R>::constraints_const_iterator iter=phd1.constraints_begin();
  for(typename Polyhedron<R>::constraints_const_iterator c_iter=phd1.constraints_begin();
      c_iter!=phd1.constraints_end(); ++c_iter)
  {
    const Constraint<R>& c=*c_iter;
    cout << c << flush;
    cout << "  " << c.satisfied_by(pt1) << endl;
  }
  
  cout << "phd1.contains(pt1)=" << flush; cout << phd1.contains(pt1) << endl;
  assert(phd1.contains(pt1));
  
  Rectangle<R> r1("[-0.06125,0.25]x[0.125,0.375]");
  cout << "r1=" << r1 << endl;
  for(class Rectangle<R>::vertices_const_iterator v_iter=r1.vertices_begin();
      v_iter!=r1.vertices_end(); ++v_iter)
  {
    for(typename Polyhedron<R>::constraints_const_iterator c_iter=phd1.constraints_begin();
        c_iter!=phd1.constraints_end(); ++c_iter)
    {
      cout << *c_iter << ".satisfied_by" << *v_iter << "=" << c_iter->satisfied_by(*v_iter) << endl;
    }
  }
  cout << "subset(r1,phd1)=" << subset(r1,phd1) << endl;
  assert(subset(r1,phd1));
  
  Rectangle<R> r2("[-0.125,0.25]x[0.125,0.75]");
  cout << "r2=" << r1 << endl;
  cout << "subset(r2,phd1)=" << subset(r2,phd1) << endl;
  assert(!subset(r2,phd1));
  
  Zonotope<R> z1(r1);
  cout << "z1=" << z1 << endl;
  cout << "subset(z1,phd1)=" << subset(z1,phd1) << endl;
  assert(subset(z1,phd1));
  
  Polytope<R> p2(r1);
  cout << "p2=" << flush; cout << p2 << endl;
  cout << "subset(p2,phd1)=" << subset(p2,phd1) << endl;
  assert(subset(p2,phd1));
  
  cout << endl;


  Polyhedron<R> phd2(Matrix<R>("[1,-1;1,-2;-1,1;-1,2]"),Vector<R>("[1.375,1.625,0.625,0.375]"));
  cout << "phd2=" << phd2 << endl << "phd2.vertices()=" << phd2.vertices() << endl;
  cout << "phd2.bounding_box()=" << phd2.bounding_box() << endl;

  Polytope<Rational> qpltp2=Polytope<Rational>(Polyhedron<Rational>(phd2));
  cout << "qpltp2=" << qpltp2 << endl << "qpltp2.bounding_box()=" << qpltp2.bounding_box() << endl;
  Matrix<R>  G=Matrix<R>("[0.889,0.615,-1.600,-2.667;-0.889,-1.231,1.600,5.333]");
  G=Matrix<R>("[3.125,1.125,-2.875,-0.875;1.75,-0.25,-2.25,-0.25]");
  cout << "G=" << G << endl;
  Polytope<R> pltp2=Polytope<R>(G);
  cout << "pltp2=" << pltp2 << endl << "pltp2.bounding_box()=" << pltp2.bounding_box() << endl;

  Rectangle<R> bbox2=phd2.bounding_box().expand_by(0.25);
  bbox2=Rectangle<R>("[-4,4]x[-4,4]");
  RegularGrid<R> gr2(2,0.125);
  GridCellListSet<R> uap2=under_approximation(phd2,gr2);
  cout << "uap2.size()=" << uap2.size() << endl;
  epsfstream eps;
  eps.open("test_polyhedron.eps",bbox2);
  eps.set_fill_colour("white");
  eps << bbox2;
  eps.set_fill_colour("green");
  eps << phd2;
  eps.set_fill_colour("blue");
  eps << uap2;
  eps.close();
  
  cout << endl;
  
  return 0;
}


template<>
int 
test_polyhedron<Rational>() 
{
  typedef Rational R;
  cout << "test_polyhedron<" << name<R>() << ">" << endl;
  LinearAlgebra::Matrix<R> A("[1,7/8;-1,9/8;1/8,-9/4]");
  LinearAlgebra::Vector<R> b("[11/8,1/2,1/4]");
  Polyhedron<R> phd(A,b);
  cout << "phd=" << phd << endl;
    
  Polytope<R> pltp(phd);
  cout << "Polytope(phd)=" << pltp << endl;
    
  Polyhedron<R> phd2(pltp);
  cout << "Polyhedron(Polytope(phd))=" << phd2 << endl;
  
  return 0;
}
