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

#include "real_typedef.h"

#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "system/henon_map.h"
#include "evaluation/apply.h"
#include "output/epsfstream.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Output;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace std;

template<typename R> int test_apply();

int main() {
  return test_apply<Real>();
}

template<typename R> 
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
  
  
  Rectangle<R> fr=apply(h,ir);
  Parallelotope<R> fp=apply(h,ip);
  //Zonotope<R> fz=apply(h,iz);
  //Polytope<R> fpl=apply(h,ipl);
  
  return 0;
}
