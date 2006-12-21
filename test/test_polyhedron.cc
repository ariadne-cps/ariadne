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

#include <ppl.hh>

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

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;

using namespace Parma_Polyhedra_Library::IO_Operators;

int test_ppl_polyhedron();
template<class R> int test_polyhedron();
template<> int test_polyhedron<Rational>();

int main() {
  
  cout << boolalpha;
  
  test_ppl_polyhedron();
  test_polyhedron<MPFloat>();
  test_polyhedron<Rational>();
   
  cerr << "INCOMPLETE ";
}


int
test_ppl_polyhedron()
{
  return 0;
}  


template<class R>
int 
test_polyhedron() 
{
  cout << "test_polyhedron<" << name<R>() << ">" << endl;
  LinearAlgebra::Matrix<R> A("[1.0,0.875;-1,1.125;0.125,-2.25]");
  LinearAlgebra::Vector<R> b("[1.375,0.5,0.25]");
  Polyhedron<R> phd;
  
  phd=Polyhedron<R>(A,b);
  cout << "phd=" << phd << endl;
  cout << "phd.dimension()=" << phd.dimension() << ", " 
       << "phd.number_of_constraints()=" << phd.number_of_constraints() << endl;
  //assert(!phd.empty());

  Point<R> pt1("(0.25,0.375)");
  cout << "pt1=" << pt1 << endl;

  Constraint<R> c=Constraint<R> (phd.dimension(),phd.A().begin(),phd.b().begin()[0]);
  cout << c << flush;
  cout << "  " << c.satisfied_by(pt1) << endl;
  
  typename Polyhedron<R>::constraints_const_iterator iter=phd.constraints_begin();
  for(typename Polyhedron<R>::constraints_const_iterator c_iter=phd.constraints_begin();
      c_iter!=phd.constraints_end(); ++c_iter)
  {
    const Constraint<R>& c=*c_iter;
    cout << c << flush;
    cout << "  " << c.satisfied_by(pt1) << endl;
  }
  
  cout << "phd.contains(pt1)=" << flush; cout << phd.contains(pt1) << endl;
  assert(phd.contains(pt1));
  
  Rectangle<R> r1("[-0.06125,0.25]x[0.125,0.375]");
  cout << "r1=" << r1 << endl;
  for(class Rectangle<R>::vertices_iterator v_iter=r1.vertices_begin();
      v_iter!=r1.vertices_end(); ++v_iter)
  {
    for(typename Polyhedron<R>::constraints_const_iterator c_iter=phd.constraints_begin();
        c_iter!=phd.constraints_end(); ++c_iter)
    {
      cout << *c_iter << ".satisfied_by" << *v_iter << "=" << c_iter->satisfied_by(*v_iter) << endl;
    }
  }
  cout << "subset(r1,phd)=" << subset(r1,phd) << endl;
  assert(subset(r1,phd));
  
  Rectangle<R> r2("[-0.125,0.25]x[0.125,0.75]");
  cout << "r2=" << r1 << endl;
  cout << "subset(r2,phd)=" << subset(r2,phd) << endl;
  assert(!subset(r2,phd));
  
  Zonotope<R> z1(r1);
  cout << "z1=" << z1 << endl;
  cout << "subset(z1,phd)=" << subset(z1,phd) << endl;
  assert(subset(z1,phd));
  
  Polytope<R> p2(r1);
  cout << "p2=" << flush; cout << p2 << endl;
  cout << "subset(p2,phd)=" << subset(p2,phd) << endl;
  assert(subset(p2,phd));
 
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
