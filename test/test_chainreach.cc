/***************************************************************************
 *            test_chainreach.cc
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

#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/grid.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"
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

template<class R> int test_chainreach();

int main() {
  return test_chainreach<Float>();
}

template<class R> 
int 
test_chainreach()
{
  Point<R> params=Point<R>("(1.5,0.875)");
  R a=params[0];
  R b=params[1];

  HenonMap<R> h=HenonMap<R>(a,b);
  Rectangle<R> gbb=Rectangle<R>("[-11.00000,5.00000]x[-8.00000,8.00000]") ;
  cout << "gbb=" << gbb << endl;
  FiniteGrid<R> fg=FiniteGrid<R>(gbb,128); // grid
  cout << "fg=" << fg << endl;
  cout << fg.subdivision_coordinate(0,0)<< std::endl;
  cout << fg.subdivision_coordinate(0,-88) << std::endl;
  cout << fg.subdivision_coordinate(0,-88)+11.000 << std::endl;
  const Grid<R>& g=fg.grid(); // grid
  Rectangle<R> ir=Rectangle<R>("[1.499,1.501]x[0.499,0.501]"); // initial state
  Rectangle<R> cb=Rectangle<R>("[-4,4]x[-4,4]"); // cutoff box
  Rectangle<R> epsbb=Rectangle<R>("[-4.1,4.1]x[-4.1,4.1]"); // eps bounding box
  
  cb=Rectangle<R>(gbb); // cutoff box
  epsbb=Rectangle<R>(gbb); // eps bounding box
  
  GridMaskSet<R> in=GridMaskSet<R>(fg);
  GridMaskSet<R> bd=GridMaskSet<R>(fg);
  in.adjoin(over_approximation(ir,g));
  bd.adjoin(over_approximation(gbb,g));

  Applicator<R> apply;

  GridMaskSet<R> cr=apply.chainreach(h,in,bd);
  PartitionTreeSet<R> ptcr=PartitionTreeSet<R>(cr);

  cout << cr <<endl;
  cout << ptcr <<endl;
  
  cout << "cr.size()=" << cr.size() << endl;
  cout << "ptcr.size()=" << ptcr.size() << "  " << flush;
  cout << "ptcr.capacity()=" << ptcr.capacity() << endl;

  epsfstream eps;
  eps.open("test_chainreach-1.eps",epsbb);
  eps.set_pen_colour("black");
  eps.set_fill_colour("white");
  eps << cb;
  eps.set_line_style(false);
  eps.set_fill_colour("green");
  eps << cr;
  eps.set_line_style(true);
  eps.set_fill_style(false);
  eps << ptcr.partition_tree();
  eps.set_line_style(true);
  eps.set_fill_colour("blue");
  eps << ir;
  eps.close();

  GridMaskSet<R> gmcr=GridMaskSet<R>(cr);
  eps.open("test_chainreach-2.eps",epsbb);
  eps.set_line_style(false);
  eps.set_fill_colour("red");
  eps << difference(gmcr.neighbourhood(),gmcr);
  eps.set_fill_colour("blue");
  eps << gmcr.adjoining();
  eps.set_fill_colour("green");
  eps << gmcr;
  eps.close();

  return 0;
}
