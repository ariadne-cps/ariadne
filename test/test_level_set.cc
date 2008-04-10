/***************************************************************************
 *            test_level_set.cc
 *
 *  Copyright  2007  Alberto Casagrande,  Pieter Collins
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

#include "test.h"
#include "test_float.h"

#include "numeric/traits.h"
#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/grid_set.h"
#include "geometry/level_set.h"
#include "function/function_interface.h"
#include "function/build_function.h"
#include "output/epsstream.h"
#include "output/logging.h"


using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Function;
using namespace Ariadne::Geometry;
using namespace Ariadne::Output;
using namespace std;

template<class R> int test_level_set();

int main() {
  
  cout << boolalpha;
  
  test_level_set<Flt>();
}


template<class R, class A, class P> 
void radius_function(R& r, const A& a, const P& p) {
  r[0] =1 - ( pow(a[0],2u) + pow(a[1],2u) );
}

ARIADNE_BUILD_FUNCTION(Radius,radius_function,1,2,0,255);

template<class R>
int 
test_level_set() 
{
  cout << "test_level_set<" << Numeric::name<R>() << ">" << endl;
  typedef typename Numeric::traits<R>::arithmetic_type A;

  Radius<R> f=Radius<R>(Vector<R>());

  LevelSet<R> s(f);

  Box<R> r;
  Point<A> pt1,pt2;
  
  // separates
  pt1=Point<A>("(0.77,0.30)");
  pt2=Point<A>("(-0.79,0.59)");
  ARIADNE_TEST_ASSERT(!s.separates(pt1,pt2));
  pt2=Point<A>("(0.81,0.63)");
  ARIADNE_TEST_ASSERT(s.separates(pt1,pt2));

  // disjoint
  r=Box<R>("[0.8,0.9]x[0.95,0.95]");
  ARIADNE_TEST_ASSERT(s.disjoint(r));
  r=Box<R>("[0.59,0.61]x[-0.81,-0.79]");
  ARIADNE_TEST_ASSERT(!bool(s.disjoint(r)));


  // approximations
  Grid<R> g(Vector<R>("[0.1875,0.125]"));
  Box<R> bb("[-1.01,1.01]x[-1.01,1.01]");
  FiniteGrid<R> fg(g,bb);
  GridMaskSet<R> gmsia(fg); 
  GridMaskSet<R> gmsoa(fg); 
  gmsia.adjoin_inner_approximation(static_cast<const SetInterface< Box<R> >&>(s));
  gmsoa.adjoin_outer_approximation(static_cast<const SetInterface< Box<R> >&>(s));

  // graphical output
  epsfstream eps;
  eps.open("test_level_set-1.eps",bb.neighbourhood(0.25));
  eps << fill_colour(red) << gmsoa;
  eps << line_style(false) << fill_colour(green) << static_cast<const SetInterface< Box<R> >&>(s);
  eps << line_style(true) << fill_colour(blue) << gmsia;
  eps.close();
  
  cout << endl;
  
  return 0;
}


