/***************************************************************************
 *            test_orthotope.cc
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
#include "geometry/orthotope.h"
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

template<class R> int test_orthotope();
 
int main() {
  test_orthotope<Flt>();
  
  cerr << "INCOMPLETE ";
  return 0;
}


template<class R>
int 
test_orthotope()
{
  typedef typename Numeric::traits<R>::arithmetic_type F;
  typedef Interval<R> I;

  Point<R> c("(0.125,-0.25,0.5)");
  LinearAlgebra::Vector<R> v("[0.0,1.0,-1.5]");
  LinearAlgebra::Matrix<R> a("[2.0,1.0,-1.5; 1.0,1.0,0.5; 0.0,0.0,0.375]");
  
  cout << "c=" << c << "\nv=" << v << "\na=" << a << endl;
  
  Box<R> r1=Box<R>("[9,11]x[5,11]x[0,0]");
  Box<R> r2=Box<R>("[5,6]x[3,4]x[-0.5,0.5]");
  cout << "r1=" << r1 << "\nr2=" << r2 << endl;

  Point<R> pt;
  cout << pt << endl;
  Box<R> r;
  cout << r << endl;
  
  Orthotope<F,R> o1=Orthotope<F,R>(r1);
  cout << "o1=" << o1 << endl;
  Orthotope<F,R> o2=Orthotope<F,R>(r1);
  cout << "o2=" << o2 << endl;
  cout << endl;


  ListSet< Orthotope<F,R> > ols;
  
  ols.clear();
  ols=o2.subdivide();
  cout << "o2.subdivide()=" << ols << std::endl;
  
  ols.clear();
  ols=o2.divide();
  cout << "o2.divide()=" << ols << std::endl;
  
  return 0;
}


