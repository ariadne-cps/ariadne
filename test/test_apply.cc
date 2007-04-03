/***************************************************************************
 *            test_apply.cc
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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

#include "test_float.h"

#include "base/pointer.h"
#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "geometry/rectangular_set.h"
#include "evaluation/applicator.h"
#include "output/epsfstream.h"

#include "models/henon.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Output;
using namespace std;

template<class R> int test_apply();

int main() {
  return test_apply<Float>();
}

template<class R> 
int 
test_apply()
{
  Point<R> params=Point<R>("(1.5,0.875)");
  R a=params[0];
  R b=params[1];

  HenonMap<R> h=HenonMap<R>(a,b);
  Rectangle<R> gbb=Rectangle<R>("[-11,5]x[-8,8]") ;
  FiniteGrid<R> fg=FiniteGrid<R>(gbb,128); // grid
  const Grid<R>& g=fg.grid(); // grid
  Rectangle<R> cb=Rectangle<R>("[-4,4]x[-4,4]"); // cutoff box
  Rectangle<R> epsbb=Rectangle<R>("[-4.1,4.1]x[-4.1,4.1]"); // eps bounding box
  
  Rectangle<R> ir=Rectangle<R>("[1.499,1.501]x[0.499,0.501]"); // initial state
  Parallelotope<R> ip=Parallelotope<R>(Rectangle<R>("[1.499,1.501]x[0.499,0.501]")); // initial state
  Zonotope<R> iz=Zonotope<R>(Rectangle<R>("[1.499,1.501]x[0.499,0.501]")); // initial state
  Polytope<R> ipl=Polytope<R>(Rectangle<R>("[1.499,1.501]x[0.499,0.501]")); // initial state
  
  cb=Rectangle<R>(gbb); // cutoff box
  epsbb=Rectangle<R>(gbb); // eps bounding box
  
  Applicator<R> apply;
  
  Rectangle<R> fr=apply.image(h,ir);
  Parallelotope<R> fp=apply.image(h,ip);
  //Zonotope<R> fz=apply(h,iz);
  //Polytope<R> fpl=apply(h,ipl);
  
  RectangularSet<R> bs(cb);
  RectangularSet<R> is(ir);
  
/*
  shared_ptr< SetInterface<R> > ims(apply.image(h,is));
  shared_ptr< SetInterface<R> > prims=apply.preimage(h,is);
*/
  shared_ptr< SetInterface<R> > reach_set(apply.reach(h,is));
  shared_ptr< SetInterface<R> > chain_reach_set(apply.chainreach(h,is,bs));
  
  cout << "rs.dimension()=" << reach_set->dimension() << endl;
  const ListSet< Parallelotope<R> >& lrs=dynamic_cast<const ListSet< Parallelotope<R> >&>(*reach_set);
  cout << "lrs.dimension()=" << lrs.dimension() << endl;
  cout << "lrs.size()=" << lrs.size() << std::endl;
  
  cout << "crs.dimension()=" << chain_reach_set->dimension() << endl;
  const GridMaskSet<R>& gmcrs=dynamic_cast<const GridMaskSet<R>&>(*chain_reach_set);
  cout << "gmcrs.dimension()=" << gmcrs.dimension() << endl;
  cout << "gmcrs.grid()=" << gmcrs.grid() << endl;
  cout << "gmcrs.block()=" << gmcrs.block() << endl;
  cout << "gmcrs.size()=" << gmcrs.size() << " of " << gmcrs.capacity() << std::endl;
  
  epsfstream eps;
  eps.open("test_apply.eps",epsbb);
  eps.set_line_style(false);
  eps.set_fill_colour("green");
  eps << gmcrs;
  eps.set_fill_colour("blue");
  eps << lrs;
  eps.set_fill_colour("yellow");
  eps << ir;
  eps.close();
  
  
  return 0;
}
