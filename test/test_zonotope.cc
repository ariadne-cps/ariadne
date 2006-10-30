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

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>

#include "real_typedef.h"

#include "ariadne.h"
#include "base/exceptions.h"
#include "base/utility.h"
#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"

#include "output/epsfstream.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Output;
using namespace std;

template<class R> int test_zonotope();
 
int main() {
  test_zonotope<Real>();
  
  cerr << "INCOMPLETE ";
  return 0;
}


template<class R>
int 
test_zonotope()
{
  typedef typename Numeric::traits<R>::arithmetic_type F;
  
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
  
  typename Zonotope<R>::vertices_const_iterator vend=z2.vertices_end();
  typename Zonotope<R>::vertices_const_iterator vi=z2.vertices_begin();
  while(vi!=vend) {
    cout << *vi << endl;
    cout << vi << endl;
    ++vi;
  }
  cout << "z2.vertices()=" << z2.vertices() << endl;
  
  Zonotope<F> z3=minkowski_sum(z1,z2);
  cout << "z3=" << z3 << endl;

  Point<R> pts[6];
  
  pts[0]=Point<R>("(17.0,15.0,-0.5)");
  pts[1]=Point<R>("(17.0,8.0,-0.5)");
  pts[2]=Point<R>("(14.0,15.0,-0.5)");
  pts[3]=Point<R>("(14.0,8.0,-0.5)");
  pts[4]=Point<R>("(14.0,5.0,-0.5)");
  pts[5]=Point<R>("(15.5,11.5,0)");

  Interval<R> unit=Interval<R>(-1,1);
  
  cout << "Writing to eps stream" << endl;
  
  Rectangle<R> bbox=z2.bounding_box().expand_by(R(0.5));
  epsfstream eps("test_zonotope.eps",bbox,0,1);
  eps << z2;
  for(uint i=0; i!=6; ++i) {
    if(z2.contains(pts[i])) {
      eps.set_fill_colour("black");
    } else {
      eps.set_fill_colour("red");
    }
  }
  eps.close();
  
  
  
  try {
    //Rectangle<R> r("[1,17/16]x[19/16,5/4]");
    //Zonotope<R> z(Point<R>("(1/2, 1/10)"),Matrix<R>("[1,1/2;1/2,3/5]"));
    Rectangle<R> r1("[1,1.0625]x[1.1875,1.25]");
    Zonotope<R> z1(r1);
    Zonotope<R> z2(Point<R>("(0.5,0.1)"),Matrix<R>("[1.0,0.5;0.5,0.6]"));
    cout << "r1=" << r1 << "\nz1=" << z1 << "\nz2=" << z2 << endl;
    cout << "disjoint(r1,z2)=" << disjoint(r1,z2) << endl;
    cout << "subset(r1,z2)=" << subset(r1,z2) << endl;
    cout << "subset(z2,r1)=" << subset(z2,r1) << endl;
    cout << "disjoint(z1,z2)=" << disjoint(z1,z2) << endl;
    cout << "subset(z1,z2)=" << subset(z1,z2) << endl;
    assert(disjoint(r1,z2));
  }
  catch(NotImplemented e) {
    cerr << "WARNING: " << e.what();
  }
  
  /*
  assert(z3.contains(pts[0]));
  assert(z3.contains(pts[1]));
  assert(z3.contains(pts[2]));
  assert(z3.contains(pts[3]));
  assert(!z3.contains(pts[4]));
  assert(z3.contains(pts[5]));

  assert(z3.contains(z3.centre()));
  
  assert(!z3.empty());
  */

  return 0;
}
