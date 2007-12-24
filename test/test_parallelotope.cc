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

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>

#include <vector>

#include "test_float.h"

#include "base/utility.h"
#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "output/epsstream.h"
#include "output/logging.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::Output;
using namespace std;

template<class R> int test_parallelotope();

int main() {
  cout << boolalpha;
  test_parallelotope<Flt>();
  cerr << "INCOMPLETE ";
  return 0;
}

template<class R>
int
test_parallelotope()
{ 
  typedef Numeric::Interval<R> I;

  Box<R> r1=Box<R>("[9,11]x[5,11]x[-1,1]");
  Box<R> r2=Box<R>("[4.875,5.125]x[2.875,3.125]x[1.75,2.25]");
  Box<R> r3=Box<R>("[-1,9]x[-1,6]x[0,4]");
  Box<R> r4=Box<R>("[-1,9]x[-1,6]x[1,4]");

  Point<R> c2("(5,3,2)");
  Matrix<R> G2("[2,1,0;1,1,0;0,0,1]");

  Point<R> c3("(0.125,-0.25)");
  Matrix<R> G3("[2,1;1,1]");

  Parallelotope<R> p1=Parallelotope<R>(r1);
  Parallelotope<R> p2=Parallelotope<R>(c2,G2);
  Parallelotope<R> p3=Parallelotope<R>(c3,G3);

  cout << "r1=" << r1 << "\nr2=" << r2<< "\nr3=" << r3 << endl;
  cout << "p1=" << p1 << "\np2=" << p2 << "\np3=" << p3 << endl;

  cout << "p1.empty()=" << p1.empty() << endl;
//  cout << "disjoint(p1,p2)=" << disjoint(p1,p2) << endl;
//  cout << "subset(p1,p2)=" << subset(p1,p2) << endl;
  for(typename Box<R>::vertices_const_iterator v_iter=r3.vertices_begin();
      v_iter!=r3.vertices_end(); ++v_iter)
  {
    cout << "p2.contains(" << *v_iter << ")=" << p2.contains(*v_iter) << endl;
  }
  
  cout << "p2.volume()=" << p2.volume() << endl;
  assert(encloses(Interval<R>(0.99,1.01),p2.volume()));
  
  Zonotope<R>& z2=p2;
  cout << "disjoint(r1,z2)=" << flush; cout << disjoint(r1,z2) << endl;
  //cout << "subset(r1,z2)=" << flush; cout << subset(r1,z2) << endl;
  //cout << "subset(r2,z2)=" << flush; cout << subset(r2,z2) << endl;
  cout << "subset(z2,r3)=" << flush; cout << subset(z2,r3) << endl;
  cout << "subset(z2,r4)=" << flush; cout << subset(z2,r4) << endl;
  cout << "disjoint(r1,p2)=" << disjoint(r1,p2) << endl;
  cout << "subset(r1,p2)=" << subset(r1,p2) << endl; 
  cout << "subset(r2,p2)=" << subset(r2,p2) << endl; 
  cout << "subset(p2,r3)=" << subset(p2,r3) << endl; 
  cout << "subset(p2,r4)=" << subset(p2,r4) << endl; 

  Point<R> pts[6];
  pts[0]=Point<R>("(17.0,15.0,-0.5) ");
  pts[1]=Point<R>("(17.0,8.0,-0.5) ");
  pts[2]=Point<R>("(14.0,15.0,-0.5) ");
  pts[3]=Point<R>("(14.0,15.0,-0.5) ");
  pts[4]=Point<R>("(14.0,15.0,-0.5) ");
  pts[5]=Point<R>("(15.5,11.5,0) ");

  
  cout << pts[0] << pts[1] << pts[2] << pts[3]<< pts[4] << pts[5]<< endl;
  
  assert((bool)(!p1.contains(pts[0])));
  assert((bool)(!p1.contains(pts[1])));
  assert((bool)(!p1.contains(pts[2])));
  assert((bool)(!p1.contains(pts[3])));
  assert((bool)(!p1.contains(pts[4])));
  assert((bool)(!p1.contains(pts[5])));

  // Test grid approximation routines
  Parallelotope<I,I> ip(p1);
  Parallelotope<I,R> ep=over_approximation(ip);
  Parallelotope<R,R> oap=over_approximation(ep);
  Parallelotope<I,I> qoap=orthogonal_over_approximation(ip);

  Box<R> r5("[1.125,1.25]x[1.625,1.75]");
  Box<R> r6("[0.875,1.00]x[0.125,0.25]");
  assert((bool)(!disjoint(p3,r6)));
  assert((bool)(disjoint(p3,r5)));
 
  Box<R> bbox3=p3.bounding_box().neighbourhood(0.25);
  Grid<R> gr3(2,0.125);
  
  GridCellListSet<R> oap3=outer_approximation(p3,gr3);
  GridCellListSet<R> uap3=inner_approximation(p3,gr3);
  epsfstream eps;
  eps.open("test_parallelotope.eps",bbox3,0,1);
  eps << fill_colour(white) << bbox3;
  eps << fill_colour(red) << oap3;
  eps << fill_colour(green) << p3;
  eps << fill_colour(blue) << uap3;
  eps << fill_colour(yellow) << r5 << r6;
  eps.close();
  
  return 0;
}
