/***************************************************************************
 *            henon_chainreach.cc
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

#include "numeric/float.h"
#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/grid.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"
#include "geometry/rectangular_set.h"
#include "evaluation/map_evolver.h"
#include "output/epsstream.h"
#include "output/logging.h"
#include "models/henon.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Output;
using namespace std;

template<class R> int henon_chainreach();

int main() {
  return henon_chainreach<Float64>();
}

template<class R> 
int 
henon_chainreach()
{
  set_applicator_verbosity(0);
  
  typedef Numeric::Interval<R> I;
  typedef Geometry::Zonotope<I,R> BS;

  double maximum_basic_set_radius=0.25;
  double grid_length=0.125;

  int subdivisions=128;

  Point<R> params=Point<R>("(1.5,0.875)");
  R a=params[0];
  R b=params[1];


  HenonMap<R> h=HenonMap<R>(params);
  Point<R> pt("(3,5)"); smoothness_type s=3;
  cout << "pt="<<pt<<" s="<<s<<endl;
  cout << "h(pt)="<<h(pt)<<endl;
  cout << "h.jacobian(pt)="<<h.jacobian(pt) << endl; 
  cout << "h.derivative(pt,s)="<<h.derivative(pt,s) << endl; 

  Rectangle<R> gbb=Rectangle<R>("[-11.0,5.0]x[-8.0,8.0]") ;
  cout << "gbb=" << gbb << endl;
  FiniteGrid<R> fg=FiniteGrid<R>(gbb,subdivisions); // grid
  cout << "fg=" << fg << endl;
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

  EvolutionParameters<R> parameters;
  parameters.set_maximum_basic_set_radius(maximum_basic_set_radius);
  parameters.set_grid_length(grid_length);
  
  MapEvolver<R> evolver(parameters);
  evolver.parameters().set_grid_length(16.0/subdivisions);

  GridMaskSet<R> gmcr=evolver.chainreach(h,in,bd);
  PartitionTreeSet<R> ptcr=PartitionTreeSet<R>(gmcr);

  cout << gmcr <<endl;
  cout << ptcr <<endl;
  
  cout << "gmcr.size()=" << gmcr.size() << endl;
  cout << "ptcr.size()=" << ptcr.size() << "  " << flush;
  cout << "ptcr.capacity()=" << ptcr.capacity() << endl;

  RectangularSet<R> ins(ir);
  RectangularSet<R> cbs(cb);
  SetInterface<R>* crs=evolver.chainreach(h,ins,cbs);
  cout << "*crs=" << *crs <<endl;
  GridMaskSet<R>& gmcrs=*dynamic_cast<GridMaskSet<R>*>(crs);
  cout << "gmcrs=" << gmcrs <<endl;
  cout << "gmcrs.size()=" << gmcrs.size() << " of " << gmcrs.capacity() << endl;

  epsfstream eps;
  eps.open("henon_chainreach-1.eps",epsbb);
  eps << fill_colour(white) << cb;
  eps << line_style(false);
  eps << fill_colour(green) << gmcr;
  eps << fill_colour(blue) << ir;
  eps << line_style(true);
  eps << ptcr.partition_tree();
  eps.close();

  eps.open("henon_chainreach-2.eps",epsbb);
  eps << line_style(false);
  eps << fill_colour(red) << difference(gmcr.neighbourhood(),gmcr);
  eps << fill_colour(blue) << gmcr.adjoining();
  eps << fill_colour(green) << gmcr;
  eps.close();

  epsbb=Rectangle<R>("[-4.1,4.1]x[-4.1,4.1]"); // eps bounding box
  eps.open("henon_chainreach-3.eps",epsbb);
  eps << line_style(false);
  eps << fill_colour(green) << gmcr;
  eps.close();

  return 0;
}
