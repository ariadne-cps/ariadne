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

#include "test_float.h"

#include "ariadne.h"

#include "base/utility.h"
#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/polyhedron.h"

#include "output/logging.h"
#include "output/epsstream.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Geometry;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Output;
using namespace std;

template<class R> int test_zonotope();
 
int main() {
  test_zonotope<Flt>();
  
  cerr << "INCOMPLETE ";
  return 0;
}


template<class R>
int 
test_zonotope()
{
  typedef typename Numeric::traits<R>::arithmetic_type F;
  typedef Interval<R> I;

  Point<R> c("(0.125,-0.25,0.5)");
  LinearAlgebra::Vector<R> v("[0.0,1.0,-1.5]");
  LinearAlgebra::Matrix<R> a("[2.0,1.0,-1.5; 1.0,1.0,0.5; 0.0,0.0,0.375]");
  
  cout << "c=" << c << "\nv=" << v << "\na=" << a << endl;
  
  Rectangle<R> r1=Rectangle<R>("[9,11]x[5,11]x[0,0]");
  Rectangle<R> r2=Rectangle<R>("[5,6]x[3,4]x[-0.5,0.5]");
  cout << "r1=" << r1 << "\nr2=" << r2 << endl;

  Parallelotope<R> p2=Parallelotope<R>(r2);
  cout << "p2=" << p2 << endl;
  
  Point<R> pt;
  cout << pt << endl;
  Rectangle<R> r;
  cout << r << endl;
  Zonotope<R> z;
  cout << z << endl;
  Zonotope<F> iz;
  cout << iz << endl;
  
  Zonotope<R> z1=Zonotope<R>(r1);
  cout << "z1=" << z1 << endl;
  Zonotope<R> z2=Zonotope<R>(r1);
  cout << "z2=" << z2 << endl;
  Zonotope<R> z3=Zonotope<R>(Point<R>("(0.125,-0.25)"),Matrix<R>("[2,1;1,1]"));
  cout << "z3=" << z3 << endl;
  Zonotope<R> z4=Zonotope<R>(Point<R>("(0,0)"),Matrix<R>("[1,0,1;0,1,1]"));
  cout << "z4=" << z4 << endl;
  cout << endl;

  // construct from a string literal
  Zonotope<R> zs("{(0,1),[1,2,3;4,5,6]}");
  cout << "zs=" << zs << endl;
  
  cout << "z2.dimension()=" << z2.dimension() << endl;
  cout << "z2.number_of_generators()=" << z2.number_of_generators() << endl;
  cout << "z2.centre()=" << z2.centre() << endl;
  for(size_type i=0; i!=z2.number_of_generators(); ++i) {
    cout << "z2.generator(" << i << ")=" << z2.generator(i) << endl;
  }
  cout << "z2.vertices()=" << z2.vertices() << endl;
  cout << endl;

  typename Zonotope<R>::vertices_const_iterator vend=z2.vertices_end();
  typename Zonotope<R>::vertices_const_iterator vi=z2.vertices_begin();
  while(vi!=vend) {
    cout << vi << flush;
    cout << "  " << *vi << endl;
    ++vi;
  }
  cout << endl;
  
  cout << "z3.vertices()=" << z3.vertices() << endl;
  cout << "z4.vertices()=" << z4.vertices() << endl;

  // Minkowski sum
  cout << "z1=" << z1 << "\nz2=" << z2 << endl;
  iz=minkowski_sum(z1,z2);
  cout << "minkowski_sum(z1,z2)=" << iz << endl;
  iz=minkowski_difference(z1,z2);
  cout << "minkowski_difference(z1,z2)=" << iz << endl;
  cout << endl;
  
  ListSet< Zonotope<R> > zls;
  cout << "z2=" << z2 << std::endl;
  
  zls.clear();
  zls=z2.subdivide();
  cout << "z2.subdivide()=" << zls << std::endl;
  
  zls.clear();
  zls=z2.divide();
  cout << "z2.divide()=" << zls << std::endl;
  
  
  Point<R> pts[6];
  
  pts[0]=Point<R>("(17.0,15.0,-0.5)");
  pts[1]=Point<R>("(17.0,8.0,-0.5)");
  pts[2]=Point<R>("(14.0,15.0,-0.5)");
  pts[3]=Point<R>("(14.0,8.0,-0.5)");
  pts[4]=Point<R>("(14.0,5.0,-0.5)");
  pts[5]=Point<R>("(15.5,11.5,0)");

  I unit=I(-1,1);
  
  cout << "Writing to eps stream" << endl;
  
  Rectangle<R> bbox=z2.bounding_box().expand_by(R(0.5));
  epsfstream eps;
  eps.open("test_zonotope-1.eps",bbox,0,1);
  eps << z2;
  for(uint i=0; i!=6; ++i) {
    if(z2.contains(pts[i])) {
      eps << fill_colour(black);
    } else {
      eps << fill_colour(red);
    }
  }
  eps.close();
  cout << "   done" << endl;
  

  Rectangle<R> bbox3=z3.bounding_box().expand_by(0.25);
  cout << "z3=" << z3 << std::endl;
  Grid<R> gr3(2,0.125);
  GridCellListSet<R> oaz3=outer_approximation(z3,gr3);
  cout << "outer_approximation(z,gr)="<<oaz3<<std::endl;
  assert(oaz3.size()==372);
  GridCellListSet<R> iaz3=inner_approximation(z3,gr3);
  cout << "inner_approximation(z,gr)="<<iaz3<<std::endl;
  cout << "iaz3.size()="<<iaz3.size() << endl;

  bbox=z3.bounding_box().expand_by(R(0.5));
  eps.open("test_zonotope-2.eps",bbox);
  eps << fill_colour(white) << bbox3;
  //eps << fill_colour(red) << oaz3;
  eps << fill_colour(green) << z3;
  //eps << fill_colour(yellow) << polyhedron(Zonotope<Rational>(z3));
  eps << fill_colour(blue) << iaz3;
  eps.close();
  assert(iaz3.size()==210);
  
  bbox=z4.bounding_box().expand_by(R(0.5));
  eps.open("test_zonotope-3.eps",bbox);
  eps << z4;
  eps.close();

  // Approximations of real zonotope
  {
    cout << "Inner and outer approximations of Zonotope<R,R>" << endl;
    //Zonotope<R,R> z(Point<R>("(0.25,-0.125)"),Matrix<R>("[2,1,1;1,-1,1]"));
    Zonotope<R,R> z(Point<R>("(0.25,-0.125)"),Matrix<R>("[2,1;1,-1]"));
    
    Grid<R> gr(2,0.5);
    GridCellListSet<R> foaz=fuzzy_outer_approximation(z,gr);
    GridCellListSet<R> oaz=outer_approximation(z,gr);
    GridCellListSet<R> iaz=inner_approximation(z,gr);
    ARIADNE_TEST_ASSERT(subset(iaz,oaz));
    ARIADNE_TEST_ASSERT(subset(oaz,foaz));
    r=Rectangle<R>("[-2.5,-2.0]x[0.0,0.5]");
    pt=Point<R>("(-2.5,0.5)");
    bbox=z.bounding_box().expand_by(R(0.5));
    eps.open("test_zonotope-real.eps",bbox);
    eps << fill_colour(white) << z.bounding_box();
    eps << fill_colour(red) << foaz;
    eps << fill_colour(yellow) << oaz;
    eps << fill_colour(blue) << iaz;
    eps << fill_colour(transparant) << z;
    eps.close();
  }    

  // over approximations
  {
    Zonotope<I> fz(Point<I>("([1.99,2.01],2)"),Matrix<I>("[[2,2.25],1;1,1]"));
    Zonotope<R> z;
    over_approximate_(z,fz);
    cout << "fz="<<fz<<endl;
    cout << "over_approximation(fz)=" << z << endl;
  }


  // Uniform error zonotope
  {
    cout << "Inner and outer approximations of Zonotope<I,R>" << endl;
    Zonotope<R,R> z(Point<R>("(0.25,-0.125)"),Matrix<R>("[2,1,1,0.125,0;1,-1,1,0,0.125]"));
    Zonotope<I,R> ez(Point<I>("([0.125,0.375],[-0.25,0.0])"),Matrix<R>("[2,1,1;1,-1,1]"));
    Zonotope<R,R> aez=over_approximation(ez);
    
    Grid<R> gr(2,0.5);
    GridCellListSet<R> foaez=fuzzy_outer_approximation(ez,gr);
    GridCellListSet<R> oaez=outer_approximation(ez,gr);
    GridCellListSet<R> iaez=inner_approximation(ez,gr);
    ARIADNE_TEST_ASSERT(subset(iaez,oaez));
    ARIADNE_TEST_ASSERT(subset(oaez,foaez));
    r=Rectangle<R>("[-2.5,-2.0]x[0.0,0.5]");
    pt=Point<R>("(-2.5,0.5)");
    bbox=z.bounding_box().expand_by(R(0.5));
    eps.open("test_zonotope-uniform.eps",bbox);
    eps << fill_colour(white) << ez.bounding_box();
    eps << fill_colour(red) << foaez;
    eps << fill_colour(yellow) << oaez;
    eps << fill_colour(blue) << iaez;
    eps << fill_colour(transparant) << ez;
    eps.close();
  }    


  {
    cout << "Orthogonal over approximations of Zonotope<I,R>" << endl;
    // Interval zonotope
    //Zonotope<I,R> ez(Point<I>("([0.125,0.375],[-0.25,0.0])"),Matrix<R>("[2,1,1;1,-1,1]"));
    //Zonotope<I,R> ez(Point<I>("([0.125,0.375],[-0.25,0.0])"),Matrix<R>("[2,1,1;1,-1,4]"));
    Zonotope<I,R> ez(Point<I>("([0.125,0.375],[-0.25,0.0])"),Matrix<R>("[4,1,1,1,1;1,1,3,3,3]"));
    Zonotope<R,R> oaez(over_approximation(ez));
    Zonotope<I,R> ooaez(orthogonal_over_approximation(ez));
    
    eps.open("test_zonotope-orthogonal.eps",ooaez.bounding_box());
    eps << fill_colour(white) << ooaez.bounding_box();
    eps << fill_colour(red) << ooaez;
    eps << fill_colour(yellow) << oaez;
    eps << fill_colour(green) << ez;
    eps.close();
  }
 

  {
    cout << "Outer approximations of Zonotope<I,I>" << endl;
    // Interval zonotope
    Zonotope<R,R> z(Point<R>("(0.25,-0.125)"),Matrix<R>("[2,1,1,0.125,0;1,-1,1,0,0.125]"));
    Zonotope<I,I> iz(Point<I>("([0.125,0.375],[-0.25,0.0])"),Matrix<R>("[2,1,1;1,-1,1]"));
    Zonotope<I,R> ez=over_approximation(iz);
    
    Grid<R> gr(2,0.5);
    GridCellListSet<R> oaz=outer_approximation(iz,gr);
    r=Rectangle<R>("[-2.5,-2.0]x[0.0,0.5]");
    pt=Point<R>("(-2.5,0.5)");
    bbox=z.bounding_box().expand_by(R(0.5));
    eps.open("test_zonotope-interval.eps",bbox);
    eps << fill_colour(white) << z.bounding_box();
    eps << fill_colour(red) << oaz;
    eps << fill_colour(yellow) << z;
    eps << fill_colour(green) << approximation(iz);
    eps.close();
    
    
    cout << "z=" << z << endl;
    cout << "r=" << r << endl;
    cout << "pt=" << pt << endl;
    cout << "z.bounding_box()=" << z.bounding_box() << endl;
    cout << "disjoint(z,r)=" << disjoint(r,z) << endl;
    cout << "subset(r,z)=" << subset(r,z) << endl;
    cout << "contains(z,pt)=" << z.contains(pt) << endl;
    pt=Point<R>("(-2.0,0.0)");
    cout << "contains(z,"<<pt<<")=" << z.contains(pt) << endl;
    pt=Point<R>("(-2.5,0.0)");
    cout << "contains(z,"<<pt<<")=" << z.contains(pt) << endl;
    pt=Point<R>("(-2.5,0.5)");
    cout << "contains(z,"<<pt<<")=" << z.contains(pt) << endl;
    pt=Point<R>("(-2.0,0.5)");
    cout << "contains(z,"<<pt<<")=" << z.contains(pt) << endl;
  }


  try {
    //Rectangle<R> r("[1,17/16]x[19/16,5/4]");
    //Zonotope<R> z(Point<R>("(1/2, 1/10)"),Matrix<R>("[1,1/2;1/2,3/5]"));
    Rectangle<R> r1("[1,1.0625]x[1.1875,1.25]");
    Zonotope<R> z1(r1);
    Zonotope<R> z2(Point<R>("(0.5,0.125)"),Matrix<R>("[1.0,0.5;0.5,0.625]"));
    cout << "r1=" << r1 << "\nz1=" << z1 << "\nz2=" << z2 << endl;
    cout << "disjoint(r1,z2)=" << disjoint(r1,z2) << endl;
    cout << "subset(r1,z2)=" << subset(r1,z2) << endl;
    cout << "subset(z2,r1)=" << subset(z2,r1) << endl;
    cout << "disjoint(z1,z2)=" << disjoint(z1,z2) << endl;
    cout << "subset(z1,z2)=" << subset(z1,z2) << endl;
    assert((bool)(disjoint(r1,z2)));
  }
  catch(NotImplemented e) {
    cerr << "Warning: " << e.what() << " not implemented\n";
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


