/***************************************************************************
 *            henon_chainreach.cc
 *
 *  Copyright  2005-8  Pieter Collins
 *
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
#include "numeric/approximate_float.h"
#include "function/function_interface.h"
#include "geometry/point.h"
#include "geometry/box.h"
#include "geometry/grid.h"
#include "geometry/grid_set.h"
#include "geometry/image_set.h"
#include "geometry/constraint_set.h"
#include "geometry/partition_tree_set.h"
#include "geometry/zonotope.h"
#include "evaluation/evolution_parameters.h"
#include "evaluation/map_evolver.h"
#include "evaluation/default_evolver.h"
#include "evaluation/default_approximator.h"
#include "evaluation/reachability_analyser.h"
#include "output/epsstream.h"
#include "output/logging.h"
#include "models/henon.h"

using namespace Ariadne;
using namespace Ariadne::Models;
using namespace std;

template<class R> int henon_chainreach();

int main() {
  return henon_chainreach<Float64>();
}

template<class R> 
int 
henon_chainreach()
{
  typedef Zonotope<R> ES;

  uint steps = 5;
  
  double maximum_enclosure_radius=0.25;
  double grid_length=0.0625;
  EvolutionParameters<R> evolution_parameters; evolution_parameters.set_maximum_enclosure_radius(maximum_enclosure_radius);
  evolution_parameters.set_grid_length(grid_length);

  EvolverInterface<Map<R>,ES>& evolver=*default_evolver<Map<R>,ES>(evolution_parameters);
  ApproximatorInterface<GridApproximationScheme<R>,ES>& approximator=*default_approximator<Box<R>,ES>();
  ReachabilityAnalyser< Map<R>, GridApproximationScheme<R> > analyser(evolution_parameters,evolver,approximator);

  Point<R> params=Point<R>("(1.5,-0.375)");
  R a=params[0];
  R b=params[1];

  HenonMap<R> henon=HenonMap<R>(params);
  Point<R> pt("(3,5)"); 
  smoothness_type sm=3;
  cout << "h="<<henon<<"\npt="<<pt<<"\nsm="<<sm<<endl;
  cout << "h(pt)="<<henon(pt)<<endl;
  cout << "h.jacobian(pt)="<<henon.jacobian(pt) << endl; 
  cout << "h.derivative(pt,sm)="<<henon.derivative(pt,sm) << endl; 

  Box<R> initial_box("[1.49,1.51]x[0.49,0.51]");
  Box<R> bounding_box("[-4,4]x[-3,3]");
  Box<R> graphics_bounding_box("[-4.1,4.1]x[-3.1,3.1]");

  Zonotope<R> initial_zonotope(initial_box);
  ListSet< Zonotope<R> > reach=evolver.reach(henon,initial_zonotope,steps,upper_semantics);
  ListSet< Zonotope<R> > evolve=evolver.evolve(henon,initial_zonotope,steps,upper_semantics);


  ImageSet<R> initial_set(initial_box);
  ConstraintSet<R> bounding_set(bounding_box);

  SetInterface< Box<R> >* chain_reach=analyser.chain_reach(henon,initial_set);
  GridMaskSet<R> grid_chain_reach=*dynamic_cast<GridMaskSet<R>*>(chain_reach);
  PartitionTreeSet<R> tree_chain_reach=PartitionTreeSet<R>(grid_chain_reach);

  cout << grid_chain_reach <<endl;
  cout << tree_chain_reach <<endl;
  
  cout << "evolve(" << steps << ").size()=" << evolve.size() << endl;
  cout << "reach(" << steps << ").size()=" << reach.size() << endl;

  cout << "grid_chain_reach.size()=" << grid_chain_reach.size() << endl;
  cout << "tree_chain_reach.size()=" << tree_chain_reach.size() << "  " << flush;
  cout << "tree_chain_reach.capacity()=" << tree_chain_reach.capacity() << endl;

  epsfstream eps;
  eps.open("henon_map-reach.eps",graphics_bounding_box);
  eps << fill_colour(white) << bounding_box;
  eps << line_style(true);
  eps << fill_colour(green) << reach;
  eps << fill_colour(yellow) << evolve;
  eps << fill_colour(blue) << initial_zonotope;
  eps.close();

  eps.open("henon_map_chainreach.eps",graphics_bounding_box);
  eps << fill_colour(white) << bounding_box;
  eps << line_style(false);
  eps << fill_colour(green) << grid_chain_reach;
  eps << fill_colour(blue) << initial_box;
  eps << line_style(true);
  eps << tree_chain_reach.partition_tree();
  eps.close();

  return 0;
}
